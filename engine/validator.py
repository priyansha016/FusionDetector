import pysam

class FusionValidator:
    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.fasta = None
        self.blacklist_prefixes = ("NBPF", "PCDH", "LOC", "HLA", "OR", "MUC", "KRT", "GOLG", "LINC", "RPS")

    def _ensure_fasta(self):
        """Opens the FASTA handle locally in the worker process."""
        if self.fasta is None:
            self.fasta = pysam.FastaFile(self.fasta_path)

    def get_sequence_context(self, chrom, pos, window=50):
        """Fetches genomic sequence around the breakpoint."""
        self._ensure_fasta()
        try:
            start = max(0, pos - window)
            end = pos + window
            return self.fasta.fetch(chrom, start, end).upper()
        except Exception:
            return ""

    def is_likely_fp(self, gene_a, gene_b, chrom_a, chrom_b, pos_a, pos_b, support=None, sa_fraction=None):
        """Factera-style filter for technical artifacts. Enhanced with support and SA tag checks.
        Returns: (is_fp: bool, reason: str)
        """
        def clean(g):
            return str(g).split('(')[0].replace('[','').replace(']','').replace("'","").strip().upper()
        
        g1, g2 = clean(gene_a), clean(gene_b)

        # 1. Same Gene or Paralog family

        if g1 == g2 or g1[:4] == g2[:4]:
            return (True, f"FILTERED: Same gene or paralog family ({g1} == {g2})")
        # Same gene family (ZNF, GYP) — often aligner confusion
        if len(g1) >= 3 and len(g2) >= 3 and g1[:3] == g2[:3]:
            if g1.startswith("ZNF") or g1.startswith("GYP"):
                return (True, f"FILTERED: Same gene family ({g1[:3]})")

        # 2. Intra-chromosomal Neighbor Filter (tiered filtering by distance)
        # Very close breakpoints (<100KB) are almost always artifacts, even with high support
        # Moderate distance (100KB-1MB) might be real inversions if high quality
        if chrom_a == chrom_b:
            try:
                dist = abs(int(pos_a) - int(pos_b))
                
                # Tier 1: Very close (<100KB) - filter almost all, very strict exceptions
                if dist < 100_000:  # 100KB
                    # Exception: Extremely high support (>100) AND high SA fraction (>0.5) might be real
                    # But most <100KB fusions are read-through transcripts or mapping artifacts
                    if support is None or support <= 100 or sa_fraction is None or sa_fraction <= 0.5:
                        return (True, f"FILTERED: Intra-chromosomal distance {dist:,}bp < 100KB (likely read-through or mapping artifact)")
                
                # Tier 2: Moderate distance (100KB-1MB) - filter unless high quality
                elif dist < 1_000_000:  # 1MB
                    # Exception: High support (>50) AND high SA fraction (>0.3) might be real inversion
                    if support is None or support <= 50 or sa_fraction is None or sa_fraction <= 0.3:
                        return (True, f"FILTERED: Intra-chromosomal distance {dist:,}bp < 1MB (likely mapping artifact)")
                
                # Do NOT require SA tags for same-chromosome fusions (e.g. EML4-ALK on chr2);
                # legacy BAMs may have low SA fraction but true fusions.
                # Same-chromosome fusions >1MB apart are handled by other filters (support thresholds, etc.)
            except (ValueError, TypeError):
                pass
        
        # 3. Blacklist Check
        # Filter ALL fusions involving blacklisted genes (LINC, LOC, etc.) regardless of support
        # These are typically repetitive elements or non-coding RNAs that cause mapping artifacts
        # No exceptions - if any gene is blacklisted, filter the fusion
        g1_blacklisted = any(g1.startswith(p) for p in self.blacklist_prefixes)
        g2_blacklisted = any(g2.startswith(p) for p in self.blacklist_prefixes)
        if g1_blacklisted or g2_blacklisted:
            matched_prefix = next((p for p in self.blacklist_prefixes if g1.startswith(p) or g2.startswith(p)), "unknown")
            return (True, f"FILTERED: Blacklisted prefix ({matched_prefix}) - these genes cause mapping artifacts")
        
        # 4. Adaptive support threshold: stricter when combined with suspicious signals
        # Base threshold: 20 reads (true fusions can have 20-25 reads)
        # Most FPs (61%) have support < 25, so we filter more aggressively
        if support is not None:
            if support < 20:  # Absolute minimum
                return (True, f"FILTERED: Support={support} < 20 (absolute minimum)")
            # Low support + low SA tags = suspicious (discordant-only artifacts)
            if support < 25 and sa_fraction is not None and sa_fraction < 0.12:
                return (True, f"FILTERED: Support={support} < 25 and SA fraction={sa_fraction:.2f} < 0.12")
            # Very low support + inter-chromosomal = likely artifact
            if support < 22 and chrom_a != chrom_b:
                return (True, f"FILTERED: Support={support} < 22 (inter-chromosomal)")
        
        # 5. Stricter SA requirements for inter-chromosomal fusions (most FPs are inter-chr)
        # Most inter-chromosomal FPs have support 22-24, so we need stricter thresholds
        if chrom_a != chrom_b:
            if sa_fraction is not None:
                # Support 22-24: require SA >= 0.20 (stricter than before)
                if support is not None and 22 <= support <= 24 and sa_fraction < 0.20:
                    return (True, f"FILTERED: Inter-chromosomal, support={support} <= 24 and SA fraction={sa_fraction:.2f} < 0.20")
                # Support 25-28: require SA >= 0.15
                if support is not None and 25 <= support <= 28 and sa_fraction < 0.15:
                    return (True, f"FILTERED: Inter-chromosomal, support={support} <= 28 and SA fraction={sa_fraction:.2f} < 0.15")
                # Support < 26: require SA >= 0.10 (catch low SA cases)
                if support is not None and support < 26 and sa_fraction < 0.10:
                    return (True, f"FILTERED: Inter-chromosomal, support={support} < 26 and SA fraction={sa_fraction:.2f} < 0.10")
        
        # 6. Minimum absolute SA read count for inter-chromosomal (reduces discordant-only artifacts)
        # Re-enabled with stricter thresholds
        # Exception: High support (>= 35) fusions are preserved (e.g., ROS1-SLC34A2 with support 42)
        if chrom_a != chrom_b and support is not None and sa_fraction is not None:
            # High support fusions are preserved regardless of SA count
            if support >= 35:
                pass  # Skip SA count check for high-support fusions
            else:
                estimated_sa_count = support * sa_fraction
                # Filter if support <= 28 and SA count < 2
                if support <= 28 and estimated_sa_count < 2:
                    return (True, f"FILTERED: Inter-chromosomal, support={support} <= 28 and estimated SA count={estimated_sa_count:.1f} < 2")
                # Filter if support <= 35 and SA count < 3
                if 28 < support <= 35 and estimated_sa_count < 3:
                    return (True, f"FILTERED: Inter-chromosomal, support={support} <= 35 and estimated SA count={estimated_sa_count:.1f} < 3")
        
        # 7. Same-chromosome with very low support and no SA evidence
        if chrom_a == chrom_b and support is not None and 20 <= support <= 24:
            if sa_fraction is not None and sa_fraction < 0.15:
                return (True, f"FILTERED: Same-chromosome, support={support} <= 24 and SA fraction={sa_fraction:.2f} < 0.15")
        
        # 8. Moderate support but very low SA fraction (inter-chromosomal)
        if chrom_a != chrom_b and support is not None and 26 <= support <= 32:
            if sa_fraction is not None and sa_fraction < 0.08:
                return (True, f"FILTERED: Inter-chromosomal, support={support} <= 32 and SA fraction={sa_fraction:.2f} < 0.08")
        
        # 9. Sequence Complexity Filter (conservative: filter fusions with low complexity sequences)
        # Only applies if FASTA is available
        # Apply to moderate-support fusions (support < 50) to catch artifacts in repeat regions
        if self.fasta_path and support is not None and support < 50:
            try:
                seq_a = self.get_sequence_context(chrom_a, pos_a, window=50)
                seq_b = self.get_sequence_context(chrom_b, pos_b, window=50)
                # Check complexity at both breakpoints
                if seq_a and seq_b:
                    complex_a = self.is_complex(seq_a)
                    complex_b = self.is_complex(seq_b)
                    # Very conservative: only filter if BOTH breakpoints have low complexity AND support is moderate/low
                    if not complex_a and not complex_b:
                        # Filter if support is low-moderate (more aggressive for very low support)
                        if support < 30:
                            return (True, f"FILTERED: Low sequence complexity at both breakpoints (support={support} < 30)")
                        elif support < 40:
                            return (True, f"FILTERED: Low sequence complexity at both breakpoints (support={support} < 40)")
                    # Also filter if one breakpoint is low complexity AND support is very low
                    elif support < 25 and (not complex_a or not complex_b):
                        return (True, f"FILTERED: Low sequence complexity at breakpoint (support={support} < 25)")
            except Exception as e:
                # If sequence extraction fails, skip complexity check (don't filter)
                # Silently continue - complexity check is optional
                pass
        
        # Passed all filters
        chrom_type = "inter-chromosomal" if chrom_a != chrom_b else "intra-chromosomal"
        sa_str = f"{sa_fraction:.2f}" if sa_fraction is not None else "N/A"
        reason = f"PASSED: Support={support}, SA={sa_str}, {chrom_type}"
        return (False, reason)

    def is_complex(self, sequence):
        """Checks if sequence has sufficient complexity (not a simple repeat).
        Returns True if sequence is complex enough, False if too simple.
        Conservative thresholds to avoid filtering true positives in repeat regions.
        """
        if not sequence or len(sequence) < 20: 
            return False  # Too short to assess
        
        # Check for extreme AT or GC bias (common in mapping artifacts)
        # Conservative: only flag if extremely biased (>90% or <10%)
        at_content = (sequence.count('A') + sequence.count('T')) / len(sequence)
        if at_content > 0.90 or at_content < 0.10:
            return False  # Too simple/extreme bias
            
        # Ensure at least 3 different nucleotides are present
        unique_nucleotides = len(set(sequence))
        if unique_nucleotides < 3:
            return False  # Too few unique nucleotides
        
        # Additional check: ensure no single nucleotide dominates (>80% of sequence)
        for nucleotide in ['A', 'T', 'G', 'C']:
            if sequence.count(nucleotide) / len(sequence) > 0.80:
                return False  # Single nucleotide dominates
        
        return True  # Sequence has sufficient complexity