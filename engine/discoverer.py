import pysam
import collections
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

# Install htslib warning filter early so it's active before any BAM opens
try:
    from utils.hts_filter import install_hts_warning_filter
    install_hts_warning_filter()
except ImportError:
    pass  # If utils not available, continue without filter


def _read_to_dict(read):
    """Serialize a pysam AlignedSegment to a dict (picklable for process-based discovery). reference_end = exclusive (pysam style)."""
    end = read.reference_end
    if end is None and read.cigartuples:
        end = read.reference_start
        for op, length in read.cigartuples:
            if op in (0, 2, 3, 7, 8):  # M, D, N, =, X
                end += length
    return {
        'reference_name': read.reference_name,
        'reference_start': read.reference_start,
        'reference_end': end,
        'query_sequence': read.query_sequence or '',
        'cigartuples': tuple(read.cigartuples) if read.cigartuples else (),
        'next_reference_name': read.next_reference_name or '',
        'next_reference_start': read.next_reference_start,
        'query_name': read.query_name,
        'tags': dict(read.get_tags()) if read.get_tags() else {},
    }


class ReadLike:
    """Minimal read-like object for the assembler (works with process-based discovery)."""
    __slots__ = ('reference_name', 'reference_start', 'reference_end', 'query_sequence', 'cigartuples',
                  'next_reference_name', 'next_reference_start', 'query_name', '_tags')

    def __init__(self, d):
        self.reference_name = d.get('reference_name', '')
        self.reference_start = int(d.get('reference_start', 0))
        self.reference_end = d.get('reference_end')  # pysam-style exclusive end (assembler does -1 for last base)
        self.query_sequence = d.get('query_sequence') or ''
        self.cigartuples = d.get('cigartuples') or ()
        self.next_reference_name = d.get('next_reference_name') or ''
        self.next_reference_start = int(d.get('next_reference_start', 0))
        self.query_name = d.get('query_name') or ''
        self._tags = d.get('tags') or {}

    def has_tag(self, tag):
        return tag in self._tags

    def get_tag(self, tag):
        return self._tags[tag]


def _scan_chunk_process(args):
    """Top-level function for process workers: scan BAM chunk, return (gene_pair, [read_dict, ...])."""
    from utils.hts_filter import install_hts_warning_filter
    install_hts_warning_filter()
    bam_path, ref_list, idx, min_mapq = args
    local = collections.defaultdict(list)
    sam = pysam.AlignmentFile(bam_path, "rb")
    try:
        for ref in ref_list:
            try:
                for read in sam.fetch(ref):
                    if read.is_unmapped or read.is_duplicate or read.mapping_quality < min_mapq:
                        continue
                    if read.has_tag("SA"):
                        genes_a = idx.get_gene_at(read.reference_name, read.reference_start)
                        if not genes_a:
                            continue
                        for part in read.get_tag("SA").split(';'):
                            if not part:
                                continue
                            sa_chrom, sa_pos = part.split(',')[0], int(part.split(',')[1])
                            genes_b = idx.get_gene_at(sa_chrom, sa_pos)
                            if genes_b:
                                for g_a in genes_a:
                                    for g_b in genes_b:
                                        if g_a != g_b:
                                            key = tuple(sorted([str(g_a), str(g_b)]))
                                            local[key].append(_read_to_dict(read))
                    elif not read.is_proper_pair and not read.mate_is_unmapped:
                        genes_a = idx.get_gene_at(read.reference_name, read.reference_start)
                        if not genes_a:
                            continue
                        genes_b = idx.get_gene_at(read.next_reference_name, read.next_reference_start)
                        if genes_b:
                            for g_a in genes_a:
                                for g_b in genes_b:
                                    if g_a != g_b:
                                        key = tuple(sorted([str(g_a), str(g_b)]))
                                        local[key].append(_read_to_dict(read))
            except Exception:
                continue
    finally:
        sam.close()
    return dict(local)


class FusionDiscoverer:
    def __init__(self, bam_path, idx, user_targets=None, n_workers=1, min_mapq=20, use_process_discovery=False):
        self.bam_path = bam_path
        self.idx = idx
        self.targets = user_targets if user_targets else set()
        self.candidates = collections.defaultdict(list)
        self.n_workers = max(1, int(n_workers))
        self.min_mapq = max(0, int(min_mapq))
        self.use_process_discovery = bool(use_process_discovery)

    def _process_read(self, read):
        """Process one read; returns list of (gene_pair_key, read) or empty."""
        if read.is_unmapped or read.is_duplicate or read.mapping_quality < self.min_mapq:
            return []
        out = []
        if read.has_tag("SA"):
            genes_a = self.idx.get_gene_at(read.reference_name, read.reference_start)
            if not genes_a:
                return []
            for part in read.get_tag("SA").split(';'):
                if not part:
                    continue
                sa_chrom, sa_pos = part.split(',')[0], int(part.split(',')[1])
                genes_b = self.idx.get_gene_at(sa_chrom, sa_pos)
                if genes_b:
                    for g_a in genes_a:
                        for g_b in genes_b:
                            if g_a != g_b:
                                key = tuple(sorted([str(g_a), str(g_b)]))
                                out.append((key, read))
            return out
        if not read.is_proper_pair and not read.mate_is_unmapped:
            genes_a = self.idx.get_gene_at(read.reference_name, read.reference_start)
            if not genes_a:
                return []
            genes_b = self.idx.get_gene_at(read.next_reference_name, read.next_reference_start)
            if genes_b:
                for g_a in genes_a:
                    for g_b in genes_b:
                        if g_a != g_b:
                            key = tuple(sorted([str(g_a), str(g_b)]))
                            out.append((key, read))
        return out

    def _collect_refs(self, refs):
        """Scan given reference names; returns dict gene_pair -> [reads]. Each thread uses its own BAM handle."""
        local = collections.defaultdict(list)
        sam = pysam.AlignmentFile(self.bam_path, "rb")
        try:
            for ref in refs:
                try:
                    for read in sam.fetch(ref):
                        for key, r in self._process_read(read):
                            local[key].append(r)
                except Exception:
                    continue
        finally:
            sam.close()
        return dict(local)

    def collect_seeds(self):
        if self.n_workers <= 1:
            return self._collect_refs(self._get_all_refs())

        refs = self._get_all_refs()
        if not refs:
            return self.candidates

        n = len(refs)
        chunk_size = max(1, (n + self.n_workers - 1) // self.n_workers)
        chunks = [refs[i:i + chunk_size] for i in range(0, n, chunk_size)]
        n_workers_use = min(len(chunks), self.n_workers)

        if self.use_process_discovery:
            # Processes: --cores N really uses N cores (threads are GIL-bound and often show as 1 core)
            print(f"[*] Discovery: using {n_workers_use} processes for BAM scan ({n} refs)", flush=True)
            args_list = [(self.bam_path, chunk, self.idx, self.min_mapq) for chunk in chunks]
            with ProcessPoolExecutor(max_workers=n_workers_use) as executor:
                for result in executor.map(_scan_chunk_process, args_list):
                    for key, read_dicts in result.items():
                        self.candidates[key].extend([ReadLike(d) for d in read_dicts])
        else:
            # Threads: when multiple BAMs, one process per BAM already uses multiple cores
            print(f"[*] Discovery: using {n_workers_use} threads for BAM scan ({n} refs)", flush=True)
            with ThreadPoolExecutor(max_workers=n_workers_use) as executor:
                futures = [executor.submit(self._collect_refs, chunk) for chunk in chunks]
                for fut in as_completed(futures):
                    for key, reads in fut.result().items():
                        self.candidates[key].extend(reads)
        return self.candidates

    def _get_all_refs(self):
        """Return list of reference names in the BAM."""
        with pysam.AlignmentFile(self.bam_path, "rb") as sam:
            return list(sam.references)