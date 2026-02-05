import collections
import gzip
from intervaltree import IntervalTree

refflat_path = 'Test/refFlat.txt.gz'
target_chrom = 'chr4'
target_pos = 25666628

# Normalize target chrom
clean_target = target_chrom.replace('chr', '')

tree = IntervalTree()
with gzip.open(refflat_path, 'rt') as f:
    for line in f:
        parts = line.strip().split('\t')
        if parts[0] == 'SLC34A2':
            chrom = parts[2].replace('chr', '')
            start, end = int(parts[4]), int(parts[5])
            print(f"DEBUG: Found SLC34A2 in file: {chrom}:{start}-{end}")
            tree.addi(start, end, parts[0])

matches = tree[target_pos]
print(f"DEBUG: Search Result for {target_pos}: {matches}")