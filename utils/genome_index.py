import collections
import gzip
from intervaltree import IntervalTree

class GenomeIndex:
    def __init__(self):
        # We'll store everything internally WITHOUT the 'chr' prefix 
        # to stay consistent regardless of the source.
        self.trees = collections.defaultdict(IntervalTree)

    def load_refflat(self, refflat_path):
        opener = gzip.open(refflat_path, 'rt') if refflat_path.endswith('.gz') else open(refflat_path, 'r')
        count = 0
        with opener as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 6: continue
                
                gene_name = parts[0]
                chrom = parts[2].replace('chr', '')
                start = int(parts[4])
                end = int(parts[5])
                
                if start < end:
                    self.trees[chrom].addi(start, end, gene_name)
                    count += 1
        
        # DO NOT call merge_overlaps() here.
        print(f"[*] Loaded {count} gene intervals into memory.")

    def get_gene_at(self, chrom, pos):
        clean_chrom = str(chrom).replace('chr', '')
        matches = self.trees[clean_chrom][pos]
        
        # Extract unique gene names from the data field of the matches
        return list(set(m.data for m in matches if m.data))

    def is_known_chrom(self, chrom):
        clean_chrom = str(chrom).replace('chr', '')
        return clean_chrom in self.trees