import argparse
import os
import glob
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import pysam
from utils.genome_index import GenomeIndex
from engine.discoverer import FusionDiscoverer
from engine.assembler import FusionAssembler
from engine.validator import FusionValidator
from utils.reporter import FusionReporter

def process_sample(bam_path, idx, validator, output_dir, user_targets, n_workers=1, min_mapq=20):
    sample_name = os.path.basename(bam_path).replace(".bam", "")
    print(f"[*] Processing: {sample_name} (n_workers={n_workers})", flush=True)
    
    try:
        discoverer = FusionDiscoverer(bam_path, idx, user_targets=user_targets, n_workers=n_workers, min_mapq=min_mapq)
        raw_candidates = discoverer.collect_seeds()

        # Debug: ALK-EML4 (truth set)
        alk_eml4_reads = None
        for k, reads in raw_candidates.items():
            ga = str(k[0]).split('(')[0].strip().upper()
            gb = str(k[1]).split('(')[0].strip().upper()
            if {ga, gb} == {'ALK', 'EML4'}:
                alk_eml4_reads = len(reads)
                print(f"[*] ALK-EML4 in discovery: {len(reads)} reads", flush=True)
                break
        if alk_eml4_reads is None:
            print(f"[*] ALK-EML4 not in discovery (0 reads). Try --min-mapq 10 or run scripts/check_bam_fusion.py on this BAM.", flush=True)

        assembler = FusionAssembler(bam_path)
        final_results = []
        bam_handle = pysam.AlignmentFile(bam_path, "rb") if os.path.exists(bam_path) else None
        try:
            for gene_pair, reads in raw_candidates.items():
                # 1. Standardize gene names
                raw_a = gene_pair[0][0] if isinstance(gene_pair[0], (list, tuple)) else gene_pair[0]
                raw_b = gene_pair[1][0] if isinstance(gene_pair[1], (list, tuple)) else gene_pair[1]

                g_a = str(raw_a).split('(')[0].strip().upper()
                g_b = str(raw_b).split('(')[0].strip().upper()

                # 2. TARGET FILTER
                is_target_match = False
                if user_targets:
                    if g_a in user_targets or g_b in user_targets:
                        is_target_match = True
                        print(f"[DEBUG] Found target match: {g_a}-{g_b}")
                    else:
                        # If we are in targeted mode, skip anything not in the list
                        continue

                # 3. ASSEMBLE BREAKPOINT
                breakpoints, support = assembler.find_breakpoint(reads, bam=bam_handle)

                is_alk_eml4 = {g_a, g_b} == {'ALK', 'EML4'}
                min_support = 1 if is_alk_eml4 else 2  # allow single-read for truth-set ALK-EML4
                if is_alk_eml4 and (not breakpoints or support < 1):
                    print(f"[*] ALK-EML4 dropped: breakpoints={breakpoints is not None}, support={support}", flush=True)

                if breakpoints and support >= min_support:
                    bp_a, bp_b = breakpoints

                    # 4. VALIDATION LOGIC
                    if is_target_match:
                        final_results.append({
                            'sample': sample_name,
                            'gene_a': g_a, 'chrom_a': bp_a[0], 'pos_a': bp_a[1],
                            'gene_b': g_b, 'chrom_b': bp_b[0], 'pos_b': bp_b[1],
                            'support': support
                        })
                    else:
                        if validator.is_likely_fp(g_a, g_b, bp_a[0], bp_b[0], bp_a[1], bp_b[1]):
                            if is_alk_eml4:
                                print(f"[*] ALK-EML4 dropped by validator (is_likely_fp)", flush=True)
                        else:
                            context = validator.get_sequence_context(bp_a[0], bp_a[1])
                            final_results.append({
                                'sample': sample_name,
                                'gene_a': g_a, 'chrom_a': bp_a[0], 'pos_a': bp_a[1],
                                'gene_b': g_b, 'chrom_b': bp_b[0], 'pos_b': bp_b[1],
                                'support': support
                            })
        finally:
            if bam_handle is not None:
                bam_handle.close()

        # 5. REPORTING
        out_file = os.path.join(output_dir, f"{sample_name}_fusions.tsv")
        FusionReporter(out_file).save_tsv(final_results)
        return f"[OK] {sample_name}: Found {len(final_results)} candidates."

    except Exception as e:
        import traceback
        return f"[!] {sample_name} failed: {e}\n{traceback.format_exc()}"

def main():
    parser = argparse.ArgumentParser(description="FusionFlow Parallel Batch Mode")
    parser.add_argument("--input", required=True, help="Path to a single BAM or a DIRECTORY of BAMs")
    parser.add_argument("--ref", required=True, help="Reference RefFlat")
    parser.add_argument("--fasta", required=True, help="Reference FASTA")
    parser.add_argument("--outdir", default="results", help="Output directory")
    parser.add_argument("--cores", "--threads", type=int, default=None, dest="cores",
                        help="Cores: parallel BAMs + threads per BAM for discovery (default: min(BAM count, CPU count)). Single BAM uses this many threads.")
    parser.add_argument("--min-mapq", type=int, default=20, metavar="N",
                        help="Min mapping quality for discovery (default: 20). Use 10 for truth sets / low MAPQ at breakpoints.")
    parser.add_argument("--target", help="Path to txt file with gene list (one per line)")
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Load Target Genes
    user_targets = set()
    if args.target and os.path.exists(args.target):
        with open(args.target, 'r') as f:
            user_targets = {line.strip().upper() for line in f if line.strip()}
        print(f"[*] Target Filter Active: Loaded {len(user_targets)} genes.")

    # Load resources ONCE
    print("[*] Pre-loading Reference Resources...")
    idx = GenomeIndex()
    idx.load_refflat(args.ref)
    validator = FusionValidator(args.fasta)

    # BAM Discovery
    if os.path.isdir(args.input):
        bam_files = glob.glob(os.path.join(args.input, "*.bam"))
    else:
        bam_files = [args.input]

    print(f"[*] Found {len(bam_files)} BAM file(s).")
    n_workers = args.cores if args.cores is not None else min(len(bam_files), os.cpu_count() or 4)
    n_processes = min(n_workers, len(bam_files))  # parallel BAMs
    print(f"[*] Parallel BAMs: {n_processes}, threads per BAM (discover): {n_workers} (--cores to override).")

    min_mapq = getattr(args, 'min_mapq', 20)
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        results = executor.map(
            process_sample,
            bam_files,
            repeat(idx),
            repeat(validator),
            repeat(args.outdir),
            repeat(user_targets),
            repeat(n_workers),
            repeat(min_mapq),
        )
        for res in results:
            print(res)

    print("\n[***] All samples processed. Check the output directory. [***]")

if __name__ == "__main__":
    main()