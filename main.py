import argparse
import os
import glob
import multiprocessing
import time
import pysam
from utils.genome_index import GenomeIndex
from engine.discoverer import FusionDiscoverer, ReadLike, _read_to_dict
from engine.assembler import FusionAssembler
from engine.validator import FusionValidator
from utils.reporter import FusionReporter

from utils.hts_filter import install_hts_warning_filter

# Install filter early so it catches all warnings from the start
install_hts_warning_filter()

MIN_SUPPORT = 20  # require support >= 20 reads for all fusions


def _serialize_read(read):
    """Convert a read (ReadLike or pysam AlignedSegment) to a picklable dict."""
    if hasattr(read, "_tags"):
        # ReadLike
        return {
            "reference_name": read.reference_name,
            "reference_start": read.reference_start,
            "reference_end": read.reference_end,
            "query_sequence": read.query_sequence or "",
            "cigartuples": read.cigartuples or (),
            "next_reference_name": read.next_reference_name or "",
            "next_reference_start": read.next_reference_start,
            "query_name": read.query_name or "",
            "tags": read._tags,
        }
    return _read_to_dict(read)


def _process_one_item(item):
    """Unpack tuple for imap_unordered (worker receives single argument)."""
    return _process_one_gene_pair(*item)


def _normalize_gene_name(gene):
    """Normalize gene name: remove suffixes like 'OS', handle common aliases."""
    g = str(gene).split("(")[0].strip().upper()
    # Handle common aliases (e.g., MYCNOS -> MYCN)
    aliases = {
        "MYCNOS": "MYCN",
        "ZNF765-ZNF761": "ZNF765",  # Handle hyphenated names
    }
    return aliases.get(g, g)


def _filter_breakpoint_clusters(results):
    """
    Filter breakpoint clusters: when multiple fusions share the same breakpoint region
    (within 10KB), it's likely a mapping artifact. This detects cases like multiple
    genes all "fusing" to chr2:33141xxx (false positive cluster).
    """
    if not results or len(results) < 2:
        return results
    
    # Cluster breakpoints within 10KB window
    CLUSTER_WINDOW = 10000
    bp_clusters = {}  # (chrom, cluster_pos) -> set of fusion indices
    fusion_list = list(results)  # Keep reference to original list
    
    # First pass: build clusters
    for idx, r in enumerate(fusion_list):
        bp_a = (r.get("chrom_a", ""), r.get("pos_a", 0))
        bp_b = (r.get("chrom_b", ""), r.get("pos_b", 0))
        
        # Check both breakpoints for clustering
        assigned_cluster = None
        for chrom, pos in [bp_a, bp_b]:
            # Find existing cluster
            for (c, cluster_pos) in list(bp_clusters.keys()):
                if c == chrom and abs(pos - cluster_pos) <= CLUSTER_WINDOW:
                    assigned_cluster = (c, cluster_pos)
                    break
            
            if assigned_cluster:
                break
        
        if assigned_cluster:
            # Add to existing cluster
            bp_clusters[assigned_cluster].add(idx)
        else:
            # Create new cluster using bp_a as anchor
            bp_clusters[bp_a] = {idx}
    
    # Filter clusters with >=3 fusions (likely artifacts)
    suspicious_clusters = {}
    for cluster_key, fusion_indices in bp_clusters.items():
        if len(fusion_indices) >= 3:
            suspicious_clusters[cluster_key] = fusion_indices
    
    if not suspicious_clusters:
        return results
    
    # Filter fusions in suspicious clusters, but keep the highest support one per cluster
    filtered = []
    cluster_best = {}  # cluster_key -> best fusion
    
    for idx, r in enumerate(fusion_list):
        bp_a = (r.get("chrom_a", ""), r.get("pos_a", 0))
        bp_b = (r.get("chrom_b", ""), r.get("pos_b", 0))
        support = r.get("support", 0)
        
        # Check if this fusion is in any suspicious cluster
        in_cluster = None
        for cluster_key, fusion_indices in suspicious_clusters.items():
            if idx in fusion_indices:
                in_cluster = cluster_key
                break
        
        if in_cluster:
            # Keep only the highest support fusion per cluster
            if in_cluster not in cluster_best or support > cluster_best[in_cluster].get("support", 0):
                # Remove previous best if exists
                if in_cluster in cluster_best:
                    if cluster_best[in_cluster] in filtered:
                        filtered.remove(cluster_best[in_cluster])
                cluster_best[in_cluster] = r
                filtered.append(r)
        else:
            # Not in a suspicious cluster, keep it
            filtered.append(r)
    
    if len(filtered) < len(results):
        print(f"[*] Filtered {len(results) - len(filtered)} breakpoint cluster fusions", flush=True)
    return filtered


def _filter_duplicate_breakpoints(results):
    """
    Filter duplicate breakpoints: if the same chromosome:position appears in multiple
    fusions, keep only the one with highest support. This handles cases like GLI4-LRMDA
    and EP400-GLI4 sharing the same breakpoint (chr8:144351031).
    """
    if not results:
        return results
    
    # Track which breakpoints we've seen and their best fusion
    bp_to_result = {}
    for r in results:
        bp_a = (r.get("chrom_a", ""), r.get("pos_a", 0))
        bp_b = (r.get("chrom_b", ""), r.get("pos_b", 0))
        fusion_key = (bp_a, bp_b)
        
        # Check if either breakpoint matches an existing fusion's breakpoints
        found_dup = False
        for existing_key, existing_r in list(bp_to_result.items()):
            existing_bp_a, existing_bp_b = existing_key
            # Check if any breakpoint overlaps
            if (bp_a == existing_bp_a or bp_a == existing_bp_b or 
                bp_b == existing_bp_a or bp_b == existing_bp_b):
                # Same breakpoint found - keep the one with higher support
                if r.get("support", 0) > existing_r.get("support", 0):
                    bp_to_result[existing_key] = r
                found_dup = True
                break
        
        if not found_dup:
            bp_to_result[fusion_key] = r
    
    deduped = list(bp_to_result.values())
    if len(deduped) < len(results):
        print(f"[*] Filtered {len(results) - len(deduped)} duplicate breakpoint fusions", flush=True)
    return deduped


def _deduplicate_gene_aliases(results):
    """Remove duplicates where same breakpoint has different gene names (e.g., MYCN vs MYCNOS)."""
    if not results:
        return results
    
    # Group by breakpoint, keep the one with highest support
    bp_to_result = {}
    for r in results:
        bp_key = (
            r.get("chrom_a", ""),
            r.get("pos_a", 0),
            r.get("chrom_b", ""),
            r.get("pos_b", 0),
        )
        # Normalize gene names for comparison
        g_a_norm = _normalize_gene_name(r.get("gene_a", ""))
        g_b_norm = _normalize_gene_name(r.get("gene_b", ""))
        norm_key = (bp_key, g_a_norm, g_b_norm)
        
        if norm_key not in bp_to_result:
            bp_to_result[norm_key] = r
        else:
            # Keep the one with higher support
            if r.get("support", 0) > bp_to_result[norm_key].get("support", 0):
                bp_to_result[norm_key] = r
    
    deduped = list(bp_to_result.values())
    if len(deduped) < len(results):
        print(f"[*] Deduplicated {len(results) - len(deduped)} gene alias duplicates", flush=True)
    return deduped


def _filter_hotspots(results, user_targets):
    """
    Filter fusions where both genes are hotspots in this sample's results.
    Hotspots are computed per sample: a gene is a hotspot only if it appears
    in >3 fusions in this sample. So SYPL1 is not a global hotspot—only in
    samples where it shows up in many fusions. Preserves target fusions and
    high-support calls.
    """
    if not results:
        return results
    
    # Count gene appearances in this sample only (using normalized names)
    gene_counts = {}
    for r in results:
        g_a = _normalize_gene_name(r.get("gene_a", ""))
        g_b = _normalize_gene_name(r.get("gene_b", ""))
        gene_counts[g_a] = gene_counts.get(g_a, 0) + 1
        gene_counts[g_b] = gene_counts.get(g_b, 0) + 1
    
    # Per-sample hotspots: genes appearing in >=3 fusions in this sample only
    HOTSPOT_THRESHOLD = 3
    hotspots = {g for g, count in gene_counts.items() if count >= HOTSPOT_THRESHOLD}
    
    # Filter duplicate breakpoints (same chr:pos appearing in multiple fusions)
    filtered = _filter_duplicate_breakpoints(results)
    
    if not hotspots:
        return filtered
    
    # Filter: remove fusions where BOTH genes are hotspots (unless one is a target or high support)
    user_targets = set(user_targets) if user_targets else set()
    final_filtered = []
    for r in filtered:
        g_a = _normalize_gene_name(r.get("gene_a", ""))
        g_b = _normalize_gene_name(r.get("gene_b", ""))
        is_hotspot_fusion = (g_a in hotspots and g_b in hotspots)
        is_target_fusion = (g_a in user_targets or g_b in user_targets)
        support = r.get("support", 0)
        
        # Keep if: not a hotspot fusion, OR it involves a target gene, OR high support
        if not is_hotspot_fusion or is_target_fusion or support >= 40:
            final_filtered.append(r)
        # Filter if ONE gene is a hotspot AND support is very low (< 25) AND low SA tags
        elif (g_a in hotspots or g_b in hotspots):
            sa_frac = r.get("sa_fraction", 1.0)  # Default to 1.0 if not present (conservative)
            # Only filter if low support AND low SA tags (likely artifact)
            if support < 25 and sa_frac < 0.2:
                continue  # Filter out low-support, low-SA-tag fusions involving hotspots
            else:
                final_filtered.append(r)  # Keep if decent SA tags even with hotspot
        else:
            final_filtered.append(r)
    
    if len(final_filtered) < len(filtered):
        print(f"[*] Filtered {len(filtered) - len(final_filtered)} hotspot fusions", flush=True)
    
    return final_filtered


def _process_one_gene_pair(bam_path, fasta_path, sample_name, gene_pair, reads_dicts, user_targets, min_support):
    """
    Worker for one gene pair: open BAM/FASTA, run assembler + validator, return result dict or None.
    Must be module-level and picklable; only receives picklable args.
    """
    install_hts_warning_filter()
    if not reads_dicts:
        return None
    reads = [ReadLike(d) for d in reads_dicts]
    user_targets = set(user_targets) if user_targets else set()

    raw_a = gene_pair[0][0] if isinstance(gene_pair[0], (list, tuple)) else gene_pair[0]
    raw_b = gene_pair[1][0] if isinstance(gene_pair[1], (list, tuple)) else gene_pair[1]
    g_a = str(raw_a).split("(")[0].strip().upper()
    g_b = str(raw_b).split("(")[0].strip().upper()

    if user_targets and g_a not in user_targets and g_b not in user_targets:
        return None

    is_target_match = bool(user_targets and (g_a in user_targets or g_b in user_targets))
    bam_handle = None
    try:
        if not os.path.exists(bam_path):
            return None
        bam_handle = pysam.AlignmentFile(bam_path, "rb")
        assembler = FusionAssembler(bam_path, fasta_path=fasta_path)
        validator = FusionValidator(fasta_path)
        breakpoints, support = assembler.find_breakpoint(reads, bam=bam_handle)
        if not breakpoints or support < min_support:
            return None
        bp_a, bp_b = breakpoints
        
        # Calculate SA tag fraction (split-read evidence)
        sa_count = sum(1 for r in reads if r.has_tag("SA"))
        sa_fraction = sa_count / len(reads) if reads else 0.0
        
        if is_target_match:
            return {
                "sample": sample_name,
                "gene_a": g_a,
                "chrom_a": bp_a[0],
                "pos_a": bp_a[1],
                "gene_b": g_b,
                "chrom_b": bp_b[0],
                "pos_b": bp_b[1],
                "support": support,
                "sa_fraction": sa_fraction,
            }
        # Non-target fusions: apply stricter filters
        if validator.is_likely_fp(g_a, g_b, bp_a[0], bp_b[0], bp_a[1], bp_b[1], support=support, sa_fraction=sa_fraction):
            return None
        return {
            "sample": sample_name,
            "gene_a": g_a,
            "chrom_a": bp_a[0],
            "pos_a": bp_a[1],
            "gene_b": g_b,
            "chrom_b": bp_b[0],
            "pos_b": bp_b[1],
            "support": support,
            "sa_fraction": sa_fraction,
        }
    except Exception:
        return None
    finally:
        if bam_handle is not None:
            bam_handle.close()


def process_sample(bam_path, ref_path, fasta_path, output_dir, user_targets, n_workers=1, min_mapq=20, use_process_discovery=False):
    sample_name = os.path.basename(bam_path).replace(".bam", "")
    print(f"[*] Processing: {sample_name} (n_workers={n_workers}, discovery={'processes' if use_process_discovery else 'threads'})", flush=True)

    try:
        idx = GenomeIndex()
        idx.load_refflat(ref_path)
        validator = FusionValidator(fasta_path)
        discoverer = FusionDiscoverer(bam_path, idx, user_targets=user_targets, n_workers=n_workers, min_mapq=min_mapq, use_process_discovery=use_process_discovery)
        t0 = time.perf_counter()
        raw_candidates = discoverer.collect_seeds()
        print(f"[*] Discovery done in {time.perf_counter() - t0:.1f}s ({len(raw_candidates)} gene pairs)", flush=True)

        # Build list of (gene_pair, reads_dicts) for pool workers; serialize reads to dicts for pickling
        items = []
        for gene_pair, reads in raw_candidates.items():
            if not reads:
                continue
            reads_dicts = [_serialize_read(r) for r in reads]
            items.append((bam_path, fasta_path, sample_name, gene_pair, reads_dicts, tuple(user_targets) if user_targets else (), MIN_SUPPORT))

        if not items:
            final_results = []
        else:
            total = len(items)
            t0 = time.perf_counter()
            if n_workers > 1:
                # Check if we're running inside a multiprocessing worker (nested pool)
                current_process = multiprocessing.current_process()
                is_worker_process = current_process.name != "MainProcess"
                if is_worker_process:
                    final_results = []
                    for i, args in enumerate(items):
                        r = _process_one_gene_pair(*args)
                        if r is not None:
                            final_results.append(r)
                        if (i + 1) % max(1, total // 10) == 0 or i + 1 == total:
                            print(f"[*] Assembly progress: {i + 1}/{total} gene pairs", flush=True)
                else:
                    n_workers_use = min(n_workers, len(items))
                    with multiprocessing.Pool(n_workers_use) as pool:
                        results = []
                        for done, r in enumerate(pool.imap_unordered(_process_one_item, items, chunksize=1), 1):
                            results.append(r)
                            if done % max(1, total // 10) == 0 or done == total:
                                print(f"[*] Assembly progress: {done}/{total} gene pairs", flush=True)
                    final_results = [r for r in results if r is not None]
            else:
                final_results = []
                for i, args in enumerate(items):
                    r = _process_one_gene_pair(*args)
                    if r is not None:
                        final_results.append(r)
                    if (i + 1) % max(1, total // 10) == 0 or i + 1 == total:
                        print(f"[*] Assembly progress: {i + 1}/{total} gene pairs", flush=True)
            
            # Post-process: deduplicate gene aliases, filter breakpoint clusters, filter duplicates, then filter hotspots
            final_results = _deduplicate_gene_aliases(final_results)
            final_results = _filter_breakpoint_clusters(final_results)
            final_results = _filter_duplicate_breakpoints(final_results)
            final_results = _filter_hotspots(final_results, user_targets)
            print(f"[*] Assembly done in {time.perf_counter() - t0:.1f}s ({len(final_results)} fusions after filtering)", flush=True)

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
        with open(args.target, "r") as f:
            user_targets = {line.strip().upper() for line in f if line.strip()}
        print(f"[*] Target Filter Active: Loaded {len(user_targets)} genes.")

    if os.path.isdir(args.input):
        bam_files = glob.glob(os.path.join(args.input, "*.bam"))
    else:
        bam_files = [args.input]

    print(f"[*] Found {len(bam_files)} BAM file(s).")
    if not bam_files:
        print("[*] No BAM files found. Exiting.")
        return
    # Filter is already installed at module level, but ensure it's active
    install_hts_warning_filter()
    n_workers = args.cores if args.cores is not None else min(len(bam_files), os.cpu_count() or 4)
    n_processes = min(n_workers, len(bam_files))
    use_process_discovery = len(bam_files) == 1
    print(f"[*] Parallel BAMs: {n_processes}, discovery: {'processes' if use_process_discovery else 'threads'} (--cores={n_workers}).")

    min_mapq = getattr(args, "min_mapq", 20)
    if n_processes == 1:
        # Single BAM: call directly to avoid nested pools (daemonic processes can't spawn children)
        # The per-gene-pair pool inside process_sample will still create N worker processes
        result = process_sample(bam_files[0], args.ref, args.fasta, args.outdir, user_targets, n_workers, min_mapq, use_process_discovery)
        print(result)
    else:
        # Multiple BAMs: use Pool for parallel BAM processing
        # Discovery will use threads (use_process_discovery=False) to avoid nested pools
        with multiprocessing.Pool(n_processes) as pool:
            results = pool.starmap(
                process_sample,
                [(f, args.ref, args.fasta, args.outdir, user_targets, n_workers, min_mapq, use_process_discovery) for f in bam_files],
            )
        for res in results:
            print(res)

    print("\n[***] All samples processed. Check the output directory. [***]")

if __name__ == "__main__":
    main()