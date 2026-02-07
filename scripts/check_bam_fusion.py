#!/usr/bin/env python3
"""
Quick check: does this BAM have evidence of EML4-ALK (or any fusion)?
Usage: python scripts/check_bam_fusion.py Test/H3122/H3122.bam [--ref Test/refFlat.txt.gz]
"""
import argparse
import sys
import os

def main():
    ap = argparse.ArgumentParser(description="Check BAM for fusion evidence (e.g. EML4-ALK)")
    ap.add_argument("bam", help="Path to BAM")
    ap.add_argument("--ref", default="Test/refFlat.txt.gz", help="RefFlat (for gene lookup)")
    args = ap.parse_args()
    if not os.path.exists(args.bam):
        print(f"BAM not found: {args.bam}", file=sys.stderr)
        sys.exit(1)

    try:
        import pysam
    except ImportError:
        print("Install pysam: pip install pysam", file=sys.stderr)
        sys.exit(1)

    # Load gene index (same as pipeline)
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from utils.genome_index import GenomeIndex
    idx = GenomeIndex()
    idx.load_refflat(args.ref)

    sa_eml4_alk = 0
    discordant_eml4_alk = 0
    sa_any = 0
    disc_any = 0
    low_mapq_sa = 0
    low_mapq_disc = 0

    with pysam.AlignmentFile(args.bam, "rb") as sam:
        refs = list(sam.references)
        for ref in refs:
            try:
                for read in sam.fetch(ref):
                    if read.is_duplicate:
                        continue
                    if read.has_tag("SA"):
                        sa_any += 1
                        if read.mapping_quality < 20:
                            low_mapq_sa += 1
                        genes_a = idx.get_gene_at(read.reference_name, read.reference_start)
                        if not genes_a:
                            continue
                        for part in read.get_tag("SA").split(';'):
                            if not part:
                                continue
                            sa_chrom, sa_pos = part.split(',')[0], int(part.split(',')[1])
                            genes_b = idx.get_gene_at(sa_chrom, sa_pos)
                            if genes_b:
                                names_a = {str(g).split('(')[0].upper() for g in genes_a}
                                names_b = {str(g).split('(')[0].upper() for g in genes_b}
                                if ('EML4' in names_a and 'ALK' in names_b) or ('ALK' in names_a and 'EML4' in names_b):
                                    sa_eml4_alk += 1
                    elif not read.is_proper_pair and not read.mate_is_unmapped:
                        disc_any += 1
                        if read.mapping_quality < 20:
                            low_mapq_disc += 1
                        genes_a = idx.get_gene_at(read.reference_name, read.reference_start)
                        genes_b = idx.get_gene_at(read.next_reference_name, read.next_reference_start)
                        if genes_a and genes_b:
                            names_a = {str(g).split('(')[0].upper() for g in genes_a}
                            names_b = {str(g).split('(')[0].upper() for g in genes_b}
                            if ('EML4' in names_a and 'ALK' in names_b) or ('ALK' in names_a and 'EML4' in names_b):
                                discordant_eml4_alk += 1
            except Exception as e:
                continue

    print(f"BAM: {args.bam}")
    print(f"  Split reads (SA) total: {sa_any} (of those, MAPQ<20: {low_mapq_sa})")
    print(f"  Discordant pairs total: {disc_any} (of those, MAPQ<20: {low_mapq_disc})")
    print(f"  EML4-ALK split reads:   {sa_eml4_alk}")
    print(f"  EML4-ALK discordant:   {discordant_eml4_alk}")
    total_eml4_alk = sa_eml4_alk + discordant_eml4_alk
    if total_eml4_alk > 0:
        print(f"  => EML4-ALK evidence: YES (total {total_eml4_alk})")
    else:
        print(f"  => EML4-ALK evidence: NONE in this BAM with current refFlat.")
        if low_mapq_sa or low_mapq_disc:
            print(f"  Tip: try running the pipeline with --min-mapq 10 to include MAPQ<20 reads.")

if __name__ == "__main__":
    main()
