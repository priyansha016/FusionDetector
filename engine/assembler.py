import collections
from collections import Counter
import statistics
import pysam

# CIGAR ops: M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8
_CIGAR_REF = (0, 2, 3, 7, 8)   # consume reference
_CIGAR_QUERY = (0, 1, 4, 7, 8)  # consume query (M,I,S,=,X)


class FusionAssembler:
    def __init__(self, bam_path=None, fasta_path=None):
        self.bam_path = bam_path
        self.fasta_path = fasta_path
        self._fasta_handle = None  # Cache FASTA handle to avoid repeated open/close
    
    def _get_fasta(self):
        """Get cached FASTA handle, opening it if needed."""
        if self._fasta_handle is None and self.fasta_path:
            try:
                self._fasta_handle = pysam.FastaFile(self.fasta_path)
            except Exception:
                return None
        return self._fasta_handle
    
    def __del__(self):
        """Clean up FASTA handle on deletion."""
        if self._fasta_handle is not None:
            try:
                self._fasta_handle.close()
            except Exception:
                pass

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
                    if self.fasta_path and evidence_a_sa:
                        bp_a, bp_b = self._refine_with_clips(reads, chrom_a, chrom_b, bp_a, bp_b)
                        bp_a, bp_b = self._bp_correction_factera(reads, chrom_a, bp_a[1], chrom_b, bp_b[1])
                    return (bp_a, bp_b), len(reads)

        # Default: use split-read positions only when available (so breakpoints are junction boundaries, not mate fragment ends)
        pos_a_candidates = [p for c, p in use_a_for_bp if c == chrom_a]
        pos_b_candidates = [p for c, p in use_b_for_bp if c == chrom_b]
        # Factera/JuLI: use exact inner boundary from CIGAR – min(primary ends) and max(SA starts) – whenever we have SA evidence
        if pos_a_candidates and pos_b_candidates and use_a_for_bp is evidence_a_sa and use_b_for_bp is evidence_b_sa:
            pos_a = min(pos_a_candidates)
            pos_b = max(pos_b_candidates)
        else:
            pos_a = self._consensus_pos(pos_a_candidates)
            pos_b = self._consensus_pos(pos_b_candidates)
        bp_a = (chrom_a, pos_a)
        bp_b = (chrom_b, pos_b)
        # Reference-based refinement (breakpoints only; fusion detection unchanged)
        if self.fasta_path and evidence_a_sa:
            bp_a, bp_b = self._refine_with_clips(reads, chrom_a, chrom_b, bp_a, bp_b)
            bp_a, bp_b = self._bp_correction_factera(reads, chrom_a, bp_a[1], chrom_b, bp_b[1])
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

    @staticmethod
    def _get_clip_sequences(read):
        """Return (five_prime_clip, three_prime_clip) as (start_offset, length), (start_offset, length) or (None, None)."""
        cig = read.cigartuples or []
        if not cig:
            return (None, None), (None, None)
        query_len = sum(l for op, l in cig if op in _CIGAR_QUERY)
        seq = read.query_sequence or ""
        if not seq:
            return (None, None), (None, None)
        # 5' clip: first op is S
        five = (None, None)
        if cig[0][0] == 4:  # S
            five = (0, cig[0][1])
        # 3' clip: last op is S
        three = (None, None)
        if cig[-1][0] == 4:
            three = (query_len - cig[-1][1], cig[-1][1])
        return five, three

    def _refine_with_clips(self, reads, chrom_a, chrom_b, bp_a, bp_b, window=150, min_clip_len=15):
        """Refine breakpoints by aligning soft-clipped sequences to the reference. Returns (bp_a, bp_b) refined."""
        fasta = self._get_fasta()
        if fasta is None:
            return bp_a, bp_b
        pos_a, pos_b = bp_a[1], bp_b[1]
        # Normalize chrom for FASTA (some have chr, some don't)
        def ref_chrom(c):
            c = str(c)
            if c in fasta.references:
                return c
            if c.startswith("chr") and c[3:] in fasta.references:
                return c[3:]
            if ("chr" + c) in fasta.references:
                return "chr" + c
            return c
        chrom_a_ref = ref_chrom(chrom_a)
        chrom_b_ref = ref_chrom(chrom_b)
        if chrom_a_ref not in fasta.references or chrom_b_ref not in fasta.references:
            return bp_a, bp_b
        clips_for_a = []  # (approx_pos, sequence) for clips that map to chrom_a
        clips_for_b = []  # (approx_pos, sequence) for clips that map to chrom_b
        for r in reads:
            if not r.has_tag("SA"):
                continue
            seq = r.query_sequence
            if not seq:
                continue
            five, three = self._get_clip_sequences(r)
            primary_chrom = r.reference_name
            sa = r.get_tag("SA").split(";")[0].split(",")
            sa_chrom, sa_pos = sa[0], int(sa[1])
            is_primary_a = ref_chrom(primary_chrom) == chrom_a_ref
            is_sa_b = ref_chrom(sa_chrom) == chrom_b_ref
            # 3' clip maps to the SA chrom
            if three[0] is not None and three[1] >= min_clip_len:
                clip_seq = seq[three[0] : three[0] + three[1]].upper()
                if is_primary_a and is_sa_b:
                    clips_for_b.append((sa_pos, clip_seq))
                elif not is_primary_a and ref_chrom(sa_chrom) == chrom_a_ref:
                    clips_for_a.append((sa_pos, clip_seq))
            # 5' clip maps to the SA chrom (segment before primary)
            if five[0] is not None and five[1] >= min_clip_len:
                clip_seq = seq[five[0] : five[0] + five[1]].upper()
                if is_primary_a and is_sa_b:
                    clips_for_b.append((sa_pos, clip_seq))
                elif not is_primary_a and ref_chrom(sa_chrom) == chrom_a_ref:
                    clips_for_a.append((sa_pos, clip_seq))
        refined_a, refined_b = pos_a, pos_b
        # Refinement: among high-scoring clip alignments, use inner boundary (min for 5' end, max for 3' start)
        if clips_for_a:
            scored = self._align_clips_to_ref_scored(fasta, chrom_a_ref, pos_a, clips_for_a, window, use_end=True)
            if scored:
                best_score_a = max(s for _, s in scored)
                top_a = [p for p, s in scored if s >= 0.85 * best_score_a]
                if top_a:
                    refined_a = min(top_a)
        if clips_for_b:
            scored = self._align_clips_to_ref_scored(fasta, chrom_b_ref, pos_b, clips_for_b, window, use_end=False)
            if scored:
                best_score_b = max(s for _, s in scored)
                top_b = [p for p, s in scored if s >= 0.85 * best_score_b]
                if top_b:
                    refined_b = max(top_b)
        # Don't close fasta - it's cached for reuse
        return (chrom_a, refined_a), (chrom_b, refined_b)

    def _bp_correction_factera(self, reads, chrom_a, pos_a, chrom_b, pos_b, window=200, min_clip=10):
        """Factera-style: compare read sequence to reference at the other breakpoint; if clip matches ref at a different offset, use that position (breakpoints only)."""
        fasta = self._get_fasta()
        if fasta is None:
            return (chrom_a, pos_a), (chrom_b, pos_b)
        def ref_chrom(c):
            c = str(c)
            if c in fasta.references:
                return c
            if c.startswith("chr") and c[3:] in fasta.references:
                return c[3:]
            if ("chr" + c) in fasta.references:
                return "chr" + c
            return c
        ca, cb = ref_chrom(chrom_a), ref_chrom(chrom_b)
        if ca not in fasta.references or cb not in fasta.references:
            return (chrom_a, pos_a), (chrom_b, pos_b)
        # Collect clips that map to each chrom (same as refinement)
        clips_b, clips_a = [], []
        for r in reads:
            if not r.has_tag("SA") or not r.query_sequence:
                continue
            seq = r.query_sequence.upper()
            five, three = self._get_clip_sequences(r)
            pc, sa = r.reference_name, r.get_tag("SA").split(";")[0].split(",")
            sa_chrom = sa[0]
            if ref_chrom(pc) == ca and ref_chrom(sa_chrom) == cb and three[0] is not None and three[1] >= min_clip:
                clips_b.append(seq[three[0] : three[0] + three[1]])
            if ref_chrom(pc) == cb and ref_chrom(sa_chrom) == ca and three[0] is not None and three[1] >= min_clip:
                clips_a.append(seq[three[0] : three[0] + three[1]])
            if five[0] is not None and five[1] >= min_clip:
                if ref_chrom(pc) == ca and ref_chrom(sa_chrom) == cb:
                    clips_b.append(seq[five[0] : five[0] + five[1]])
                elif ref_chrom(pc) == cb and ref_chrom(sa_chrom) == ca:
                    clips_a.append(seq[five[0] : five[0] + five[1]])
        rc_table = str.maketrans("ACGT", "TGCA")
        def best_pos(chrom_ref, center, clips, use_end):
            ref_len = fasta.get_reference_length(chrom_ref)
            start = max(0, center - window)
            end = min(ref_len, center + window)
            ref_seq = fasta.fetch(chrom_ref, start, end).upper()
            if len(ref_seq) < min_clip:
                return center
            best_pos_val = center
            best_score = -1
            for clip in clips:
                if len(clip) > len(ref_seq):
                    continue
                for is_rc, seq in enumerate((clip, clip[::-1].translate(rc_table))):
                    for offset in range(0, len(ref_seq) - len(seq) + 1):
                        score = sum(1 for a, b in zip(seq, ref_seq[offset : offset + len(seq)]) if a == b)
                        if score > best_score and score >= max(5, int(len(seq) * 0.5)):
                            best_score = score
                            if use_end:
                                best_pos_val = (start + offset) if is_rc else (start + offset + len(seq) - 1)
                            else:
                                best_pos_val = start + offset
            return best_pos_val if best_score >= 5 else center
        new_a = best_pos(ca, pos_a, clips_a, use_end=True)
        new_b = best_pos(cb, pos_b, clips_b, use_end=False)
        # Don't close fasta - it's cached for reuse
        return (chrom_a, new_a), (chrom_b, new_b)

    def _align_clips_to_ref_scored(self, fasta, chrom, center_pos, clips, window, use_end=False, min_identity=0.5):
        """Align each clip to reference; return [(pos, score), ...]. Caller picks best (e.g. max score)."""
        ref_len = fasta.get_reference_length(chrom)
        start = max(0, center_pos - window)
        end = min(ref_len, center_pos + window)
        ref_seq = fasta.fetch(chrom, start, end).upper()
        if len(ref_seq) < 10:
            return []
        rc_table = str.maketrans("ACGT", "TGCA")
        scored = []
        for _approx, clip_seq in clips:
            if len(clip_seq) > len(ref_seq):
                continue
            best_offset = None
            best_score = -1
            best_was_rc = False
            for is_rc, seq in enumerate((clip_seq, clip_seq[::-1].translate(rc_table))):
                for offset in range(0, len(ref_seq) - len(seq) + 1):
                    score = sum(1 for a, b in zip(seq, ref_seq[offset : offset + len(seq)]) if a == b)
                    if score > best_score:
                        best_score = score
                        best_offset = offset
                        best_was_rc = bool(is_rc)
            min_score = max(5, int(len(clip_seq) * min_identity))
            if best_offset is not None and best_score >= min_score:
                if use_end:
                    pos = (start + best_offset) if best_was_rc else (start + best_offset + len(clip_seq) - 1)
                else:
                    pos = start + best_offset
                scored.append((pos, best_score))
        return scored
