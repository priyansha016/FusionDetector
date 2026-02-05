import pysam
import collections

class BamSniffer:
    """
    Inspects the BAM to decide the best detection strategy.
    """
    def __init__(self, bam_path):
        self.bam_path = bam_path

    def sniff_strategy(self, sample_size=10000):
        """
        Scans the beginning of the BAM to detect the 'Alignment Signature'.
        Returns: 'LEGACY' or 'MODERN'
        """
        sam = pysam.AlignmentFile(self.bam_path, "rb")
        
        sa_tag_count = 0
        discordant_count = 0
        total_mapped = 0
        
        # Look at a sample of reads to identify the aligner behavior
        for i, read in enumerate(sam):
            if i >= sample_size:
                break
            
            if read.is_unmapped or read.is_secondary:
                continue
                
            total_mapped += 1
            
            # Evidence of Modern Aligners (BWA-MEM, etc.)
            if read.has_tag("SA"):
                sa_tag_count += 1
                
            # Evidence of Discordant Pairs (Factera's bread and butter)
            if not read.is_proper_pair and not read.mate_is_unmapped:
                # Check if chromosomes are different
                if read.reference_id != read.next_reference_id:
                    discordant_count += 1
        
        sam.close()
        
        if total_mapped == 0:
            return "UNKNOWN"

        # Logic: If SA tags are rare (<1%) but discordants exist, it's Legacy.
        sa_ratio = sa_tag_count / total_mapped
        
        if sa_ratio < 0.01:
            print(f"[*] Sniffer: Low SA tags ({sa_ratio:.2%}). Activating LEGACY mode (Factera-style).")
            return "LEGACY"
        else:
            print(f"[*] Sniffer: High SA tags ({sa_ratio:.2%}). Activating MODERN mode (JuLI/Genefuse-style).")
            return "MODERN"

# --- Testing Trail ---
# sniffer = BamSniffer("your_sample.bam")
# mode = sniffer.sniff_strategy()