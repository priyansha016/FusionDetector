#!/usr/bin/env python3
"""
Compare FusionDetector output to expected breakpoints (truth).
Usage:
  python scripts/check_breakpoints.py Test/results/H3122_fusions.tsv --truth scripts/truth_H3122.txt
  python scripts/check_breakpoints.py Test/results/H3122_fusions.tsv --expected "ALK,EML4,chr2,29446607,chr2,42526889"
Truth file format (one line per fusion, tab-separated):
  sample  gene_a  gene_b  chrom_a  pos_a_expected  chrom_b  pos_b_expected
Positions in truth are 1-based (standard). TSV may be 0- or 1-based; we normalize.
"""
import argparse
import csv
import sys


def parse_truth_line(line):
    toks = [t.strip() for t in line.split("\t")]
    if len(toks) < 7:
        return None
    return {
        "sample": toks[0],
        "gene_a": toks[1].upper(),
        "gene_b": toks[2].upper(),
        "chrom_a": toks[3],
        "pos_a": int(toks[4]),
        "chrom_b": toks[5],
        "pos_b": int(toks[6]),
    }


def load_truth_from_file(path):
    truth = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            t = parse_truth_line(line)
            if t:
                truth.append(t)
    return truth


def load_truth_from_args(sample, gene_a, gene_b, chrom_a, pos_a, chrom_b, pos_b):
    return [{
        "sample": sample,
        "gene_a": gene_a.upper(),
        "gene_b": gene_b.upper(),
        "chrom_a": chrom_a,
        "pos_a": int(pos_a),
        "chrom_b": chrom_b,
        "pos_b": int(pos_b),
    }]


def norm_chrom(c):
    c = str(c).strip()
    if c.upper().startswith("CHR"):
        return c
    return f"chr{c}"


def match_fusion(row, truth_entry):
    ga = str(row.get("gene_a", "")).upper()
    gb = str(row.get("gene_b", "")).upper()
    ta = truth_entry["gene_a"]
    tb = truth_entry["gene_b"]
    return (ga == ta and gb == tb) or (ga == tb and gb == ta)


def main():
    ap = argparse.ArgumentParser(description="Compare fusion TSV to expected breakpoints")
    ap.add_argument("tsv", help="Path to sample_fusions.tsv")
    ap.add_argument("--truth", help="Path to truth file (tab: sample gene_a gene_b chrom_a pos_a chrom_b pos_b)")
    ap.add_argument("--expected", help="Inline: gene_a,gene_b,chrom_a,pos_a,chrom_b,pos_b (1-based)")
    ap.add_argument("--tolerance", type=int, default=500, help="Max bp distance to consider 'pass' (default 500)")
    ap.add_argument("--sample", default=None, help="Sample name if not from TSV")
    args = ap.parse_args()

    # Load truth
    if args.truth:
        truth_list = load_truth_from_file(args.truth)
    elif args.expected:
        parts = [p.strip() for p in args.expected.split(",")]
        if len(parts) < 6:
            print("--expected must be gene_a,gene_b,chrom_a,pos_a,chrom_b,pos_b", file=sys.stderr)
            sys.exit(1)
        gene_a, gene_b, chrom_a, pos_a, chrom_b, pos_b = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
        sample = args.sample or "H3122"
        truth_list = load_truth_from_args(sample, gene_a, gene_b, chrom_a, pos_a, chrom_b, pos_b)
    else:
        print("Provide --truth <file> or --expected <gene_a,gene_b,chr,pos,chr,pos>", file=sys.stderr)
        sys.exit(1)

    # Load TSV (assume 0-based if values look smaller than truth; we'll compare in 1-based)
    rows = []
    with open(args.tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)

    tolerance = args.tolerance
    found_any = False
    for te in truth_list:
        # Expect 1-based truth
        t_pos_a = te["pos_a"]
        t_pos_b = te["pos_b"]
        t_chr_a = norm_chrom(te["chrom_a"])
        t_chr_b = norm_chrom(te["chrom_b"])

        for row in rows:
            if not match_fusion(row, te):
                continue
            found_any = True
            # FusionDetector reports 1-based in TSV (genomics standard). Truth is 1-based.
            pos_a = int(row["pos_a"])
            pos_b = int(row["pos_b"])
            chrom_a = norm_chrom(row["chrom_a"])
            chrom_b = norm_chrom(row["chrom_b"])
            # If TSV was from an old run (0-based), convert for comparison
            if pos_a < t_pos_a - 500:
                pos_a += 1
            if pos_b < t_pos_b - 500:
                pos_b += 1
            dist_a = abs(pos_a - t_pos_a)
            dist_b = abs(pos_b - t_pos_b)
            pass_a = dist_a <= tolerance
            pass_b = dist_b <= tolerance
            status = "PASS" if (pass_a and pass_b) else "FAIL"
            print(f"Fusion: {row['gene_a']}-{row['gene_b']}  [{status}]")
            print(f"  Break1: {chrom_a}  reported={pos_a}  expected={t_pos_a}  dist={dist_a} bp  {'ok' if pass_a else 'OFF'}")
            print(f"  Break2: {chrom_b}  reported={pos_b}  expected={t_pos_b}  dist={dist_b} bp  {'ok' if pass_b else 'OFF'}")
            print()

    if not found_any:
        print("No matching fusion found in TSV for the given truth.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
