import pysam
import collections


class FusionDiscoverer:
    def __init__(self, bam_path, idx, mode, user_targets=None):
        self.bam_path = bam_path
        self.idx = idx
        self.mode = mode
        # Store targets as a set for fast lookup
        self.targets = user_targets if user_targets else set()
        self.candidates = collections.defaultdict(list)

    def collect_seeds(self):
        sam = pysam.AlignmentFile(self.bam_path, "rb")
        
        for read in sam.fetch():
            if read.is_unmapped or read.is_duplicate or read.mapping_quality < 20:
                continue

            # 1. ALWAYS process Split-Reads (SA Tags)
            if read.has_tag("SA"):
                genes_a = self.idx.get_gene_at(read.reference_name, read.reference_start)
                if not genes_a: continue
                
                for part in read.get_tag("SA").split(';'):
                    if not part: continue
                    sa_chrom, sa_pos = part.split(',')[0], int(part.split(',')[1])
                    genes_b = self.idx.get_gene_at(sa_chrom, sa_pos)
                    if genes_b:
                        for g_a in genes_a:
                            for g_b in genes_b:
                                if g_a != g_b:
                                    key = tuple(sorted([str(g_a), str(g_b)]))
                                    self.candidates[key].append(read)

            # 2. TARGETED Discordant Reads (The Speed Hack)
            elif not read.is_proper_pair and not read.mate_is_unmapped:
                genes_a = self.idx.get_gene_at(read.reference_name, read.reference_start)
                if not genes_a: continue
                
                # Only check Side B if Side A is in your genes.txt
                if any(str(g).upper() in self.targets for g in genes_a):
                    genes_b = self.idx.get_gene_at(read.next_reference_name, read.next_reference_start)
                    if genes_b:
                        for g_a in genes_a:
                            for g_b in genes_b:
                                if g_a != g_b:
                                    key = tuple(sorted([str(g_a), str(g_b)]))
                                    self.candidates[key].append(read)
        sam.close()
        return self.candidates