import collections
import gzip
import re
from intervaltree import IntervalTree

def _normalize_chrom(chrom):
    """Strip leading chr/Chr so 'chr2', 'Chr2', '2' all become '2'."""
    s = str(chrom).strip()
    return re.sub(r'^[Cc]hr', '', s) if s else s

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
                chrom = _normalize_chrom(parts[2])
                start = int(parts[4])
                end = int(parts[5])
                
                if start < end:
                    self.trees[chrom].addi(start, end, gene_name)
                    count += 1
        
        # DO NOT call merge_overlaps() here.
        print(f"[*] Loaded {count} gene intervals into memory.")

    def get_gene_at(self, chrom, pos):
        clean_chrom = _normalize_chrom(chrom)
        matches = self.trees[clean_chrom][pos]
        
        # Extract unique gene names from the data field of the matches
        return list(set(m.data for m in matches if m.data))

    def is_known_chrom(self, chrom):
        return _normalize_chrom(chrom) in self.trees


def load_exon_boundaries(refflat_path):
    """
    Load exon boundaries from refFlat file for post-processing filter.
    Returns: dict mapping chrom -> set of exon boundary positions
    """
    exon_boundaries = collections.defaultdict(set)
    opener = gzip.open(refflat_path, 'rt') if refflat_path.endswith('.gz') else open(refflat_path, 'r')
    count = 0
    with opener as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 11: continue  # Need at least 11 columns for exon info
            
            chrom = _normalize_chrom(parts[2])
            
            # Parse exon boundaries (columns 9 and 10: exonStarts, exonEnds)
            try:
                exon_starts = [int(x) for x in parts[9].rstrip(',').split(',') if x]
                exon_ends = [int(x) for x in parts[10].rstrip(',').split(',') if x]
                
                # Add all exon boundaries (both starts and ends)
                for pos in exon_starts + exon_ends:
                    exon_boundaries[chrom].add(pos)
                    count += 1
            except (ValueError, IndexError):
                # Skip if exon parsing fails
                pass
    
    print(f"[*] Loaded {count} exon boundaries for post-processing filter.")
    return exon_boundaries


def get_nearest_exon_boundary(exon_boundaries, chrom, pos, window=20):
    """
    Find the distance to the nearest exon boundary.
    Returns: (distance, is_near_boundary) where is_near_boundary is True if within window bp
    """
    clean_chrom = _normalize_chrom(chrom)
    if clean_chrom not in exon_boundaries:
        return (None, False)
    
    boundaries = exon_boundaries[clean_chrom]
    if not boundaries:
        return (None, False)
    
    # Find nearest boundary
    min_dist = None
    for boundary_pos in boundaries:
        dist = abs(pos - boundary_pos)
        if min_dist is None or dist < min_dist:
            min_dist = dist
    
    is_near = min_dist is not None and min_dist <= window
    return (min_dist, is_near)