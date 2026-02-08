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

        # Prefer split-read (SA) evidence for breakpoint positions; discordant mates give fragment ends, not junction (Factera/JuLI use clip boundaries)
        use_a_for_bp = evidence_a_sa if evidence_a_sa else evidence_a
        use_b_for_bp = evidence_b_sa if evidence_b_sa else evidence_b

        # Same-chromosome: two distinct breakpoint clusters (e.g. ALK ~29Mb, EML4 ~42Mb on chr2)
        MIN_CLUSTER_GAP = 1_000_000   # 1 Mb
        if chrom_a == chrom_b:
            all_a = [p for c, p in use_a_for_bp if c == chrom_a]
            all_b = [p for c, p in use_b_for_bp if c == chrom_b]
            combined = sorted(all_a + all_b)
            if len(combined) >= 2:
                gaps = [(combined[i + 1] - combined[i], i + 1) for i in range(len(combined) - 1)]
                gap_size, split_at = max(gaps, key=lambda x: x[0])
                if gap_size >= MIN_CLUSTER_GAP:
                    left, right = combined[:split_at], combined[split_at:]
                    med_left = self._consensus_pos(left)
                    med_right = self._consensus_pos(right)
                    # Assign bp_a to cluster with more evidence_a, bp_b to the other (so we get two distinct positions)
                    n_a_left = sum(1 for p in all_a if abs(p - med_left) <= abs(p - med_right))
                    if n_a_left >= len(all_a) / 2:
                        bp_a, bp_b = (chrom_a, med_left), (chrom_b, med_right)
                    else:
                        bp_a, bp_b = (chrom_a, med_right), (chrom_b, med_left)
                    return (bp_a, bp_b), len(reads)

        # Default: use split-read positions only when available (so breakpoints are junction boundaries, not mate fragment ends)
        pos_a_candidates = [p for c, p in use_a_for_bp if c == chrom_a]
        pos_b_candidates = [p for c, p in use_b_for_bp if c == chrom_b]
        pos_a = self._consensus_pos(pos_a_candidates)
        pos_b = self._consensus_pos(pos_b_candidates)
        bp_a = (chrom_a, pos_a)
        bp_b = (chrom_b, pos_b)
        return (bp_a, bp_b), len(reads)

    @staticmethod
    def _consensus_pos(positions, window_bp=50, trim_frac=0.15):
        """Best estimate of breakpoint from a list of supporting positions.
        - Use exact mode if one position has >= 20% of evidence (true junction often has many reads at same base).
        - Else use densest window: find the window of size window_bp with the most positions, return median inside it (robust to outliers).
        - Fallback: trimmed median (drop trim_frac from each tail) then median.
        """
        if not positions:
            return 0
        if len(positions) == 1:
            return positions[0]
        cnt = Counter(positions)
        mode_pos, mode_count = cnt.most_common(1)[0]
        if mode_count >= max(2, len(positions) * 0.20):
            return mode_pos
        sorted_pos = sorted(positions)
        n = len(sorted_pos)
        # Densest window: slide a window of size window_bp, keep the one with max count
        best_start = 0
        best_count = 0
        left = 0
        for right in range(n):
            while sorted_pos[right] - sorted_pos[left] > window_bp:
                left += 1
            if right - left + 1 > best_count:
                best_count = right - left + 1
                best_start = left
        if best_count >= max(2, n * 0.15):
            cluster = sorted_pos[best_start : best_start + best_count]
            return int(statistics.median(cluster))
        # Trimmed median: drop trim_frac from each tail to reduce outlier effect
        drop = max(0, int(n * trim_frac))
        if drop > 0 and n > 2 * drop:
            trimmed = sorted_pos[drop : n - drop]
            if trimmed:
                return int(statistics.median(trimmed))
        return int(statistics.median(sorted_pos))

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