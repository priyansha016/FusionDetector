import collections
from collections import Counter
import statistics
import pysam

class FusionAssembler:
    def __init__(self, bam_path=None):
        self.bam_path = bam_path

    def find_breakpoint(self, reads, bam=None):
        """Pass open bam handle to reuse (one open per sample); otherwise opens/closes per call."""
        if not reads:
            return None, 0

        evidence_a = []
        evidence_a_sa = []
        evidence_b = []
        evidence_b_sa = []
        mate_end_cache = {}
        open_here = bam is None
        if bam is None and self.bam_path:
            bam = pysam.AlignmentFile(self.bam_path, "rb")
        try:
            for r in reads:
                if r.has_tag("SA"):
                    sa = r.get_tag("SA").split(';')[0].split(',')
                    c_b, p_b = sa[0], int(sa[1])
                    evidence_b.append((c_b, p_b))
                    evidence_b_sa.append((c_b, p_b))
                    evidence_a.append((r.reference_name, self._primary_segment_end(r)))
                    evidence_a_sa.append((r.reference_name, self._primary_segment_end(r)))
                else:
                    evidence_a.append((r.reference_name, self._primary_segment_end(r)))
                    mate_end = self._get_mate_end(r, mate_end_cache, bam)
                    evidence_b.append((r.next_reference_name, mate_end))
        finally:
            if open_here and bam:
                bam.close()

        if not evidence_a or not evidence_b:
            return None, 0

        # Chromosome consensus
        use_a = evidence_a_sa if evidence_a_sa else evidence_a
        chrom_a = Counter(c for c, _ in use_a).most_common(1)[0][0]
        evidence_b_other_chr = [(c, p) for (c, p) in evidence_b if c != chrom_a]
        if not evidence_b_other_chr:
            evidence_b_other_chr = evidence_b
        other_chr_sa = [(c, p) for (c, p) in evidence_b_sa if c != chrom_a] if chrom_a else evidence_b_sa
        use_b = other_chr_sa if other_chr_sa else evidence_b_other_chr
        chrom_b = Counter(c for c, _ in use_b).most_common(1)[0][0]

        # Same-chromosome: two clusters or one cluster + noise (e.g. ALK-EML4 real at 29Mb, noise at 42Mb)
        MIN_CLUSTER_GAP = 1_000_000   # 1 Mb
        NOISE_CLUSTER_GAP = 5_000_000  # gap > 5 Mb: smaller cluster likely noise, use dominant only
        if chrom_a == chrom_b:
            all_a = [p for c, p in (evidence_a_sa if evidence_a_sa else evidence_a) if c == chrom_a]
            all_b = [p for c, p in evidence_b_sa if c == chrom_b] or [p for c, p in evidence_b if c == chrom_b]
            combined = sorted(all_a + all_b)
            if len(combined) >= 2:
                gaps = [(combined[i + 1] - combined[i], i + 1) for i in range(len(combined) - 1)]
                gap_size, split_at = max(gaps, key=lambda x: x[0])
                if gap_size >= MIN_CLUSTER_GAP:
                    left, right = combined[:split_at], combined[split_at:]
                    # Very large gap + one cluster much smaller => treat smaller as noise (e.g. ALK-EML4 29Mb vs 42Mb)
                    if gap_size > NOISE_CLUSTER_GAP and (len(left) >= 3 * len(right) or len(right) >= 3 * len(left)):
                        use = left if len(left) >= len(right) else right
                        med_use = self._consensus_pos(use)
                        within = 2_000_000
                        a_near = [p for p in all_a if abs(p - med_use) <= within] or all_a
                        b_near = [p for p in all_b if abs(p - med_use) <= within] or all_b
                        bp_a = (chrom_a, self._consensus_pos(a_near))
                        bp_b = (chrom_b, self._consensus_pos(b_near))
                    else:
                        med_left = self._consensus_pos(left)
                        med_right = self._consensus_pos(right)
                        n_a_left = sum(1 for p in all_a if abs(p - med_left) <= abs(p - med_right))
                        if n_a_left >= len(all_a) / 2:
                            bp_a, bp_b = (chrom_a, med_left), (chrom_b, med_right)
                        else:
                            bp_a, bp_b = (chrom_a, med_right), (chrom_b, med_left)
                    return (bp_a, bp_b), len(reads)

        # Default: primary vs other (and for same-chr with no clear gap, use original logic)
        pos_a_sa = [p for c, p in evidence_a_sa if c == chrom_a]
        pos_a = self._consensus_pos(pos_a_sa if pos_a_sa else [p for c, p in use_a if c == chrom_a])
        bp_a = (chrom_a, pos_a)
        pos_b_sa = [p for c, p in evidence_b_sa if c == chrom_b]
        pos_b = self._consensus_pos(pos_b_sa if pos_b_sa else [p for c, p in use_b if c == chrom_b])
        bp_b = (chrom_b, pos_b)
        return (bp_a, bp_b), len(reads)

    @staticmethod
    def _consensus_pos(positions):
        """Median position (robust); fallback to mode if only one value."""
        if not positions:
            return 0
        if len(positions) == 1:
            return positions[0]
        return int(statistics.median(positions))

    def _get_mate_end(self, read, cache, bam=None):
        """Mate's last ref base (0-based). Uses mate from BAM when available; else estimate."""
        qname = read.query_name
        if qname in cache:
            return cache[qname]
        mate_end = read.next_reference_start + max(0, len(read.query_sequence) - 1)  # fallback
        if bam:
            try:
                for mate in bam.fetch(
                    read.next_reference_name,
                    max(0, read.next_reference_start - 50),
                    read.next_reference_start + 600,
                ):
                    if mate.query_name == qname and mate.reference_start != read.reference_start:
                        if mate.reference_end is not None:
                            mate_end = mate.reference_end - 1
                        break
            except Exception:
                pass
        cache[qname] = mate_end
        return mate_end

    @staticmethod
    def _primary_segment_end(read):
        """Last reference base (0-based) of the primary aligned segment."""
        if read.reference_end is not None:
            return read.reference_end - 1
        end = read.reference_start
        for op, length in (read.cigartuples or []):
            if op in (0, 2, 3, 7, 8):  # M, D, N, =, X
                end += length
        return end - 1