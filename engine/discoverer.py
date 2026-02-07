import pysam
import collections
from concurrent.futures import ThreadPoolExecutor, as_completed


class FusionDiscoverer:
    def __init__(self, bam_path, idx, user_targets=None, n_workers=1, min_mapq=20):
        self.bam_path = bam_path
        self.idx = idx
        self.targets = user_targets if user_targets else set()
        self.candidates = collections.defaultdict(list)
        self.n_workers = max(1, int(n_workers))
        self.min_mapq = max(0, int(min_mapq))

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

        # Split refs across workers (roughly equal by count)
        n = len(refs)
        chunk_size = max(1, (n + self.n_workers - 1) // self.n_workers)
        chunks = [refs[i:i + chunk_size] for i in range(0, n, chunk_size)]
        n_threads = min(len(chunks), self.n_workers)
        print(f"[*] Discovery: using {n_threads} threads for BAM scan ({n} refs)", flush=True)
        # ThreadPoolExecutor: pysam AlignedSegment is not picklable, so we cannot use processes.

        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            futures = [executor.submit(self._collect_refs, chunk) for chunk in chunks]
            for fut in as_completed(futures):
                for key, reads in fut.result().items():
                    self.candidates[key].extend(reads)
        return self.candidates

    def _get_all_refs(self):
        """Return list of reference names in the BAM."""
        with pysam.AlignmentFile(self.bam_path, "rb") as sam:
            return list(sam.references)