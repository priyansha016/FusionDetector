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
        """Factera-style filter for technical artifacts. Enhanced with support and SA tag checks."""
        def clean(g):
            return str(g).split('(')[0].replace('[','').replace(']','').replace("'","").strip().upper()
        
        g1, g2 = clean(gene_a), clean(gene_b)

        # 1. Same Gene or Paralog family
        if g1 == g2 or g1[:4] == g2[:4]:
            return True
        # Same gene family (ZNF, GYP) — often aligner confusion
        if len(g1) >= 3 and len(g2) >= 3 and g1[:3] == g2[:3]:
            if g1.startswith("ZNF") or g1.startswith("GYP"):
                return True

        # 2. Intra-chromosomal Neighbor Filter (500KB: likely read-through)
        if chrom_a == chrom_b:
            try:
                dist = abs(int(pos_a) - int(pos_b))
                if dist < 500000:  # 500KB: likely neighbor/read-through (e.g. GYPA-GYPE)
                    return True
                # Do NOT require SA tags for same-chromosome fusions (e.g. EML4-ALK on chr2);
                # legacy BAMs may have low SA fraction but true fusions.
            except (ValueError, TypeError):
                pass
        
        # 3. Blacklist Check
        if any(g1.startswith(p) or g2.startswith(p) for p in self.blacklist_prefixes):
            return True
        
        # 4. Adaptive support threshold: stricter when combined with suspicious signals
        # Base threshold: 20 reads (true fusions can have 20-25 reads)
        # Stricter thresholds apply when other suspicious signals are present
        if support is not None:
            if support < 20:  # Absolute minimum
                return True
            # Low support + low SA tags = suspicious (discordant-only artifacts)
            if support < 25 and sa_fraction is not None and sa_fraction < 0.15:
                return True
            # Very low support + inter-chromosomal = likely artifact
            if support < 22 and chrom_a != chrom_b:
                return True
        
        # 5. Require minimum SA tag fraction for inter-chromosomal fusions (if available)
        if chrom_a != chrom_b and sa_fraction is not None and sa_fraction < 0.15:
            if support is not None and support < 28:
                return True
        if chrom_a != chrom_b and support is not None and 20 <= support <= 21 and sa_fraction is not None and sa_fraction < 0.25:
            return True  # Very low support + low SA → filter
        
        # 6. Minimum absolute SA read count for inter-chromosomal (reduces discordant-only artifacts)
        if chrom_a != chrom_b and support is not None and sa_fraction is not None:
            estimated_sa_count = support * sa_fraction
            if support <= 30 and estimated_sa_count < 3:
                return True
            if 30 < support <= 45 and estimated_sa_count < 4:
                return True  # Keep 24-read fusions with ≥3 SA (e.g. CD74-ROS1)
        
        # 7. Same-chromosome with very low support and no SA evidence
        if chrom_a == chrom_b and support is not None and 20 <= support <= 23:
            if sa_fraction is not None and sa_fraction < 0.15:
                return True  # Discordant-only same-chromosome at 20–23 reads
        
        # 8. Moderate support but very low SA fraction (inter-chromosomal)
        if chrom_a != chrom_b and support is not None and 24 <= support <= 32:
            if sa_fraction is not None and sa_fraction < 0.1:
                return True  # Require at least ~10% SA when support is moderate
                
        return False

    def is_complex(self, sequence):
        """Filters out simple repeats (low AT/GC diversity)."""
        if not sequence or len(sequence) < 20: 
            return False
        
        # Check for extreme AT or GC bias (common in mapping artifacts)
        at_content = (sequence.count('A') + sequence.count('T')) / len(sequence)
        if at_content > 0.85 or at_content < 0.15:
            return False # Too simple
            
        # Ensure at least 3 different nucleotides are present
        return len(set(sequence)) >= 3