import argparse
import os
import glob
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from utils.genome_index import GenomeIndex
from engine.sniffer import BamSniffer
from engine.discoverer import FusionDiscoverer
from engine.assembler import FusionAssembler
from engine.validator import FusionValidator
from utils.reporter import FusionReporter

def process_sample(bam_path, idx, validator, output_dir, user_targets):
    sample_name = os.path.basename(bam_path).replace(".bam", "")
    print(f"[*] Processing: {sample_name}")
    
    try:
        sniffer = BamSniffer(bam_path)
        # Pass user_targets to the discoverer for the "Speed Hack"
        discoverer = FusionDiscoverer(bam_path, idx, sniffer.sniff_strategy(), user_targets=user_targets)
        raw_candidates = discoverer.collect_seeds()

        assembler = FusionAssembler()
        final_results = []
        
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
            breakpoints, support = assembler.find_breakpoint(reads)
            
            if breakpoints and support >= 3:
                bp_a, bp_b = breakpoints
                
                # 4. VALIDATION LOGIC
                # If it's a target (like CD74-ROS1), we trust it and bypass filters
                if is_target_match:
                    final_results.append({
                        'sample': sample_name,
                        'gene_a': g_a, 'chrom_a': bp_a[0], 'pos_a': bp_a[1],
                        'gene_b': g_b, 'chrom_b': bp_b[0], 'pos_b': bp_b[1],
                        'support': support
                    })
                else:
                    # For non-targets (Discovery mode), use the full validator
                    # This kills neighbors, paralogs, and simple repeats
                    if not validator.is_likely_fp(g_a, g_b, bp_a[0], bp_b[0], bp_a[1], bp_b[1]):
                        context = validator.get_sequence_context(bp_a[0], bp_a[1])
                        # if validator.is_complex(context):
                        final_results.append({
                            'sample': sample_name,
                            'gene_a': g_a, 'chrom_a': bp_a[0], 'pos_a': bp_a[1],
                            'gene_b': g_b, 'chrom_b': bp_b[0], 'pos_b': bp_b[1],
                            'support': support
                        })

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
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel samples to process")
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
    print(f"[*] Launching Parallel Pool with {args.threads} workers...")

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(
            process_sample, 
            bam_files, 
            repeat(idx), 
            repeat(validator), 
            repeat(args.outdir),
            repeat(user_targets)
        )
        for res in results:
            print(res)

    print("\n[***] All samples processed. Check the output directory. [***]")

if __name__ == "__main__":
    main()