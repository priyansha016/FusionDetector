import argparse
import os
import sys
import glob
import multiprocessing
import time
from collections import Counter
import pysam
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from utils.genome_index import GenomeIndex, load_exon_boundaries, get_nearest_exon_boundary
from engine.discoverer import FusionDiscoverer, ReadLike, _read_to_dict
from engine.assembler import FusionAssembler
from engine.validator import FusionValidator
from utils.reporter import FusionReporter

from utils.hts_filter import install_hts_warning_filter

# Install filter early so it catches all warnings from the start
install_hts_warning_filter()

MIN_SUPPORT = 20  # require support >= 20 reads for all fusions


class TraceLogger:
    """Thread-safe logger for writing trace information to a log file."""
    _initialized_files = set()
    _init_lock = threading.Lock()
    _file_handles = {}  # Cache open file handles per process
    _handle_lock = threading.Lock()
    
    def __init__(self, log_file_path):
        self.log_file_path = log_file_path
        self.lock = threading.Lock()
        # Only initialize file once (in main process)
        with TraceLogger._init_lock:
            if log_file_path not in TraceLogger._initialized_files:
                with open(self.log_file_path, 'w') as f:
                    f.write(f"# Fusion Tracer Log\n")
                    f.write(f"# Log file: {log_file_path}\n\n")
                TraceLogger._initialized_files.add(log_file_path)
    
    def _get_file_handle(self):
        """Get or create a buffered file handle for this logger."""
        # Use process ID to avoid sharing handles across processes
        import os
        pid = os.getpid()
        key = (pid, self.log_file_path)
        
        with TraceLogger._handle_lock:
            if key not in TraceLogger._file_handles:
                # Open in append mode with buffering
                TraceLogger._file_handles[key] = open(self.log_file_path, 'a', buffering=8192)
            return TraceLogger._file_handles[key]
    
    def log(self, message):
        """Write a trace message to the log file (thread-safe, buffered)."""
        with self.lock:
            f = self._get_file_handle()
            f.write(f"{message}\n")
            # Flush periodically but not every line (buffering helps performance)
            if f.tell() % 65536 < len(message):  # Flush every ~64KB
                f.flush()
    
    def log_filter_result(self, fusion_name, filter_name, passed, reason=None):
        """Log a filter result for a traced fusion."""
        status = "PASSED" if passed else "FILTERED"
        msg = f"[TRACE] {fusion_name}: Post-filter: {filter_name} - {status}"
        if reason:
            msg += f" ({reason})"
        self.log(msg)


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
    # Current format: 10 args (bam_path, fasta_path, ref_path, sample_name, gene_pair, reads_dicts, user_targets, min_support, traced_pairs, trace_log_path)
    if len(item) == 10:
        return _process_one_gene_pair(*item)
    # Legacy compatibility (should not be needed in current codebase)
    elif len(item) == 7:
        return _process_one_gene_pair(*item[:2], None, *item[2:], None, None)
    elif len(item) == 9:
        return _process_one_gene_pair(item[0], item[1], None, *item[2:])
    else:
        raise ValueError(f"Unexpected item length: {len(item)} (expected 10)")


def _normalize_gene_name(gene):
    """Normalize gene name: remove suffixes like 'OS', handle common aliases."""
    g = str(gene).split("(")[0].strip().upper()
    # Handle common aliases (e.g., MYCNOS -> MYCN)
    aliases = {
        "MYCNOS": "MYCN",
        "ZNF765-ZNF761": "ZNF765",  # Handle hyphenated names
    }
    return aliases.get(g, g)


def _filter_sink_breakpoints(results, user_targets=None, traced_pairs=None, trace_logger=None, filter_log_file=None):
    """
    Sink filter: a breakpoint region that appears in many different fusions is likely
    a mapping artifact (e.g. chr2:33141xxx "fusing" with many genes). 
    
    Smart filtering: Only filters low-quality fusions in sink regions. Preserves:
    - Target genes (user_targets)
    - High-support fusions (support >= 30)
    - High SA-fraction fusions (SA >= 0.25)
    
    This prevents filtering true fusions like ROS1-SLC34A2 that have high support
    but happen to be in regions with many false positives.
    """
    if not results or len(results) < 2:
        return results, []
    BIN_KB = 10  # 10kb bin (same order as cluster window)
    SINK_THRESHOLD = 3  # bin in >3 fusions = sink
    MIN_SUPPORT_EXEMPT = 30  # Fusions with support >= this are exempt from sink filtering
    MIN_SA_EXEMPT = 0.25  # Fusions with SA fraction >= this are exempt from sink filtering
    
    bin_counts = Counter()
    # Pre-compute bins and normalized names
    fusion_data = []
    user_targets = set(user_targets) if user_targets else set()
    
    for r in results:
        chrom_a, pos_a = r.get("chrom_a", ""), r.get("pos_a", 0)
        chrom_b, pos_b = r.get("chrom_b", ""), r.get("pos_b", 0)
        bin_a = (chrom_a, pos_a // (BIN_KB * 1000)) if chrom_a and pos_a is not None else None
        bin_b = (chrom_b, pos_b // (BIN_KB * 1000)) if chrom_b and pos_b is not None else None
        g_a_norm = _normalize_gene_name(r.get("gene_a", ""))
        g_b_norm = _normalize_gene_name(r.get("gene_b", ""))
        fusion_name = f"{g_a_norm}-{g_b_norm}"
        support = r.get("support", 0)
        sa_frac = r.get("sa_fraction", 0.0)
        is_target = g_a_norm in user_targets or g_b_norm in user_targets
        fusion_data.append((r, bin_a, bin_b, fusion_name, g_a_norm, g_b_norm, support, sa_frac, is_target))
        
        if bin_a:
            bin_counts[bin_a] += 1
        if bin_b and bin_b != bin_a:
            bin_counts[bin_b] += 1
    
    sink_bins = {b for b, c in bin_counts.items() if c > SINK_THRESHOLD}
    if not sink_bins:
        return results, []
    
    # Identify severe sinks (appear in many fusions) - these filter ALL fusions regardless of support
    SEVERE_SINK_THRESHOLD = 5  # If a bin appears in >= 5 fusions, it's a severe sink
    severe_sink_bins = {b for b, c in bin_counts.items() if c >= SEVERE_SINK_THRESHOLD}
    
    # Chromosome-specific sink detection: chr1, chr2, chr3 are large chromosomes prone to mapping artifacts
    # If a breakpoint is on these chromosomes AND in a sink bin, apply stricter filtering
    LARGE_CHROMOSOMES = {"chr1", "chr2", "chr3"}
    
    filtered = []
    filtered_names = []
    for r, bin_a, bin_b, fusion_name, g_a_norm, g_b_norm, support, sa_frac, is_target in fusion_data:
        is_traced = traced_pairs and ((g_a_norm, g_b_norm) in traced_pairs or (g_b_norm, g_a_norm) in traced_pairs)
        
        # Check if breakpoints are in sink bins
        bp_a_in_sink = bin_a in sink_bins if bin_a else False
        bp_b_in_sink = bin_b in sink_bins if bin_b else False
        bp_a_in_severe_sink = bin_a in severe_sink_bins if bin_a else False
        bp_b_in_severe_sink = bin_b in severe_sink_bins if bin_b else False
        
        # Exemptions: don't filter if...
        # 1. Target gene (user wants to see these) - always exempt
        if is_target:
            if is_traced and trace_logger:
                trace_logger.log_filter_result(fusion_name, "sink_breakpoints", True, "Target gene (bypasses filter)")
            filtered.append(r)
            continue
        
        # 2. Severe sink: Filter most fusions, but keep very high-quality ones
        # (e.g., chr2:33141xxx appearing in 6+ fusions is likely artifact, but keep best ones)
        # Stricter: Require very high support OR both moderate-high support AND high SA fraction
        # This filters more mapping artifacts while still preserving true very high-quality fusions
        if bp_a_in_severe_sink or bp_b_in_severe_sink:
            sink_bin = bin_a if bp_a_in_severe_sink else bin_b
            sink_count = bin_counts[sink_bin]
            
            # In severe sinks, keep only very high-quality fusions:
            # - Very high support (>= 80) alone, OR
            # - Moderate-high support (>= 40) AND high SA fraction (>= 0.30) together
            # This filters mapping artifacts like chr2:33141xxx while preserving true fusions like ROS1-SLC34A2 (support 42, SA 0.30)
            # Note: Default min_support is 20, so 40-80 is 2-4x the minimum, balanced for severe sinks
            if support >= 80 or (support >= 40 and sa_frac >= 0.30):
                if is_traced and trace_logger:
                    if support >= 80:
                        quality_reason = f"very high support={support} >= 80"
                    else:
                        quality_reason = f"high support={support} >= 50 AND SA fraction={sa_frac:.2f} >= 0.30"
                    trace_logger.log_filter_result(fusion_name, "sink_breakpoints", True, f"Severe sink bin {sink_bin[0]}:{sink_bin[1]*BIN_KB}KB ({sink_count} fusions) but kept due to {quality_reason}")
                filtered.append(r)
                continue
            else:
                # Filter low-quality fusions in severe sinks
                filtered_names.append(fusion_name)
                if is_traced and trace_logger:
                    trace_logger.log_filter_result(fusion_name, "sink_breakpoints", False, f"Breakpoint in severe sink bin {sink_bin[0]}:{sink_bin[1]*BIN_KB}KB (appears in {sink_count} fusions, support={support} < 80 and (support < 40 or SA < 0.30))")
                continue
        
        # 3. Regular sink: Apply exemptions for high-quality fusions
        if bp_a_in_sink or bp_b_in_sink:
            sink_bin = bin_a if bp_a_in_sink else bin_b
            sink_count = bin_counts[sink_bin]
            
            # Chromosome-specific stricter filtering: chr1/chr2/chr3 are prone to mapping artifacts
            # If breakpoint is on large chromosome AND in sink, require higher support/SA
            is_large_chr_sink = sink_bin[0] in LARGE_CHROMOSOMES
            
            if is_large_chr_sink:
                # Large chromosome sink: require support >= 30 OR SA >= 0.30 (stricter than regular sinks)
                if support >= 30 or sa_frac >= 0.30:
                    if is_traced and trace_logger:
                        quality_reason = f"support={support} >= 30" if support >= 30 else f"SA fraction={sa_frac:.2f} >= 0.30"
                        sink_info = f" (breakpoint on {sink_bin[0]} in sink bin {sink_bin[0]}:{sink_bin[1]*BIN_KB}KB but {quality_reason})"
                        trace_logger.log_filter_result(fusion_name, "sink_breakpoints", True, f"Large chromosome sink exempt{sink_info}")
                    filtered.append(r)
                    continue
                else:
                    # Filter low-quality fusions in large chromosome sinks
                    filtered_names.append(fusion_name)
                    if is_traced and trace_logger:
                        trace_logger.log_filter_result(fusion_name, "sink_breakpoints", False, f"Breakpoint on {sink_bin[0]} in sink bin {sink_bin[0]}:{sink_bin[1]*BIN_KB}KB (appears in {sink_count} fusions, support={support} < 30 and SA={sa_frac:.2f} < 0.30)")
                    continue
            
            # Regular sink (not on large chromosome): standard exemptions
            # High support exemption (for moderate sinks with 3-4 fusions)
            # Lowered threshold to 25 to catch ROS1-SLC34A2 (support 42) and similar high-quality fusions
            if support >= 25:
                if is_traced and trace_logger:
                    sink_info = f" (breakpoint in sink bin {sink_bin[0]}:{sink_bin[1]*BIN_KB}KB but support={support} >= 25)"
                    trace_logger.log_filter_result(fusion_name, "sink_breakpoints", True, f"High support exempt{sink_info}")
                filtered.append(r)
                continue
            
            # High SA fraction exemption
            if sa_frac >= MIN_SA_EXEMPT:
                if is_traced and trace_logger:
                    sink_info = f" (breakpoint in sink bin {sink_bin[0]}:{sink_bin[1]*BIN_KB}KB but SA fraction={sa_frac:.2f} >= {MIN_SA_EXEMPT})"
                    trace_logger.log_filter_result(fusion_name, "sink_breakpoints", True, f"High SA fraction exempt{sink_info}")
                filtered.append(r)
                continue
            
            # Filter low-quality fusions in regular sinks
            filtered_names.append(fusion_name)
            if is_traced and trace_logger:
                trace_logger.log_filter_result(fusion_name, "sink_breakpoints", False, f"Breakpoint in sink bin {sink_bin[0]}:{sink_bin[1]*BIN_KB}KB (appears in {sink_count} fusions, support={support} < 25)")
            continue
        
        # Not filtered - keep it
        if is_traced and trace_logger:
            trace_logger.log_filter_result(fusion_name, "sink_breakpoints", True, "Breakpoints not in sink bins or exempted")
        filtered.append(r)
    
    return filtered, filtered_names


def _filter_breakpoint_clusters(results, user_targets=None, traced_pairs=None, trace_logger=None, filter_log_file=None):
    """
    Filter breakpoint clusters: when multiple fusions share the same breakpoint region
    (within 10KB), it's likely a mapping artifact. This detects cases like multiple
    genes all "fusing" to chr2:33141xxx (false positive cluster).
    """
    if not results or len(results) < 2:
        return results, []
    
    # Cluster breakpoints within 10KB window
    CLUSTER_WINDOW = 10000
    bp_clusters = {}  # (chrom, cluster_pos) -> set of fusion indices
    fusion_list = list(results)  # Keep reference to original list
    
    # Pre-normalize gene names once
    normalized_genes = {}
    for idx, r in enumerate(fusion_list):
        g_a = r.get("gene_a", "")
        g_b = r.get("gene_b", "")
        normalized_genes[idx] = (_normalize_gene_name(g_a), _normalize_gene_name(g_b))
    
    # First pass: build clusters (optimized - check existing clusters more efficiently)
    for idx, r in enumerate(fusion_list):
        bp_a = (r.get("chrom_a", ""), r.get("pos_a", 0))
        bp_b = (r.get("chrom_b", ""), r.get("pos_b", 0))
        
        # Check both breakpoints for clustering
        assigned_cluster = None
        for chrom, pos in [bp_a, bp_b]:
            # Find existing cluster - check only clusters for this chromosome
            for cluster_key in bp_clusters.keys():
                c, cluster_pos = cluster_key
                if c == chrom and abs(pos - cluster_pos) <= CLUSTER_WINDOW:
                    assigned_cluster = cluster_key
                    break
            
            if assigned_cluster:
                break
        
        if assigned_cluster:
            # Add to existing cluster
            bp_clusters[assigned_cluster].add(idx)
        else:
            # Create new cluster using bp_a as anchor
            bp_clusters[bp_a] = {idx}
    
    # Filter clusters with >=2 fusions (stricter: clusters are strong FP signals)
    suspicious_clusters = {k: v for k, v in bp_clusters.items() if len(v) >= 2}
    
    if not suspicious_clusters:
        return results, []
    
    # Build index: fusion_idx -> cluster_key for O(1) lookup
    fusion_to_cluster = {}
    for cluster_key, fusion_indices in suspicious_clusters.items():
        for idx in fusion_indices:
            fusion_to_cluster[idx] = cluster_key
    
    # Filter fusions in suspicious clusters, but keep high-quality ones
    # Strategy: Keep top N fusions per cluster (by support), or all fusions with support >= threshold
    # This preserves true fusions like ROS1-SLC34A2 even if they're in clusters
    # Stricter: Require higher support OR both moderate support AND good SA evidence
    CLUSTER_KEEP_SUPPORT_THRESHOLD = 50  # Keep all fusions with support >= this, even in clusters (stricter)
    CLUSTER_KEEP_SA_THRESHOLD = 0.30  # Or keep if support >= 40 AND SA >= this
    CLUSTER_KEEP_TOP_N = 1  # Keep only top 1 fusion per cluster (stricter - most clusters are artifacts)
    
    # First pass: find best fusions per cluster
    cluster_fusions = {}  # cluster_key -> list of (idx, support, sa_frac)
    for idx, r in enumerate(fusion_list):
        in_cluster = fusion_to_cluster.get(idx)
        if in_cluster:
            support = r.get("support", 0)
            sa_frac = r.get("sa_fraction", 0.0)
            if in_cluster not in cluster_fusions:
                cluster_fusions[in_cluster] = []
            cluster_fusions[in_cluster].append((idx, support, sa_frac))
    
    # Determine which fusions to keep per cluster
    cluster_keep_indices = set()
    for cluster_key, fusions in cluster_fusions.items():
        # Sort by support (descending), then SA fraction
        sorted_fusions = sorted(fusions, key=lambda x: (x[1], x[2]), reverse=True)
        
        # Keep all fusions with support >= threshold
        for idx, support, sa_frac in sorted_fusions:
            if support >= CLUSTER_KEEP_SUPPORT_THRESHOLD:
                cluster_keep_indices.add(idx)
            # Or keep if moderate support AND good SA evidence
            elif support >= 40 and sa_frac >= CLUSTER_KEEP_SA_THRESHOLD:
                cluster_keep_indices.add(idx)
        
        # Also keep top N fusions per cluster (even if below threshold) - but only if they meet quality criteria
        kept_top = 0
        for idx, support, sa_frac in sorted_fusions:
            if idx not in cluster_keep_indices and kept_top < CLUSTER_KEEP_TOP_N:
                # Only keep top N if they have reasonable quality
                if support >= 35 or (support >= 30 and sa_frac >= 0.25):
                    cluster_keep_indices.add(idx)
                    kept_top += 1
    
    # Second pass: build filtered list and track removals
    filtered = []
    removed_fusions = []
    
    for idx, r in enumerate(fusion_list):
        g_a_norm, g_b_norm = normalized_genes[idx]
        fusion_name = f"{g_a_norm}-{g_b_norm}"
        is_traced = traced_pairs and ((g_a_norm, g_b_norm) in traced_pairs or (g_b_norm, g_a_norm) in traced_pairs)
        
        in_cluster = fusion_to_cluster.get(idx)
        if in_cluster:
            if idx in cluster_keep_indices:
                # This fusion is kept (high support or top N in cluster)
                filtered.append(r)
                if is_traced and trace_logger:
                    support = r.get("support", 0)
                    reason = f"In cluster but kept (support={support} >= {CLUSTER_KEEP_SUPPORT_THRESHOLD} or in top {CLUSTER_KEEP_TOP_N})"
                    trace_logger.log_filter_result(fusion_name, "breakpoint_clusters", True, reason)
            else:
                # This fusion is in a cluster but not kept, remove it
                removed_fusions.append((fusion_name, in_cluster, r.get("support", 0)))
        else:
            # Not in a suspicious cluster, keep it
            if is_traced and trace_logger:
                trace_logger.log_filter_result(fusion_name, "breakpoint_clusters", True, "Not in a suspicious cluster")
            filtered.append(r)
    
    # Log trace for removed fusions
    if traced_pairs and trace_logger:
        for fusion_name, cluster_key, support_val in removed_fusions:
            g_a_norm, g_b_norm = fusion_name.split('-', 1)
            if ((g_a_norm, g_b_norm) in traced_pairs or (g_b_norm, g_a_norm) in traced_pairs):
                trace_logger.log_filter_result(fusion_name, "breakpoint_clusters", False, f"In cluster {cluster_key[0]}:{cluster_key[1]}, support={support_val} not in top {CLUSTER_KEEP_TOP_N} or < {CLUSTER_KEEP_SUPPORT_THRESHOLD}")
    
    filtered_names = [name for name, _, _ in removed_fusions]
    return filtered, filtered_names


def _filter_duplicate_breakpoints(results, user_targets=None, traced_pairs=None, trace_logger=None, filter_log_file=None):
    """
    Filter duplicate breakpoints: if the same chromosome:position appears in multiple
    fusions, keep only the one with highest support. This handles cases like GLI4-LRMDA
    and EP400-GLI4 sharing the same breakpoint (chr8:144351031).
    Optimized: uses breakpoint -> result mapping for O(1) lookup instead of O(n) iteration.
    """
    if not results:
        return results, []
    
    # Track which breakpoints we've seen and their best fusion
    bp_to_result = {}  # (bp_a, bp_b) -> result dict
    bp_set = set()  # Set of all breakpoints seen (for fast lookup)
    removed_fusions = []
    
    # Pre-normalize gene names once
    normalized_names = {}
    for idx, r in enumerate(results):
        g_a_norm = _normalize_gene_name(r.get("gene_a", ""))
        g_b_norm = _normalize_gene_name(r.get("gene_b", ""))
        normalized_names[idx] = (g_a_norm, g_b_norm, f"{g_a_norm}-{g_b_norm}")
    
    for idx, r in enumerate(results):
        bp_a = (r.get("chrom_a", ""), r.get("pos_a", 0))
        bp_b = (r.get("chrom_b", ""), r.get("pos_b", 0))
        fusion_key = (bp_a, bp_b)
        g_a_norm, g_b_norm, fusion_name = normalized_names[idx]
        
        # Check if any breakpoint matches an existing fusion's breakpoints (O(1) lookup)
        found_dup = False
        for bp in [bp_a, bp_b]:
            if bp in bp_set:
                # Find the fusion(s) with this breakpoint
                for existing_key, existing_r in bp_to_result.items():
                    existing_bp_a, existing_bp_b = existing_key
                    if bp == existing_bp_a or bp == existing_bp_b:
                        # Same breakpoint found - keep the one with higher support
                        existing_support = existing_r.get("support", 0)
                        current_support = r.get("support", 0)
                        if current_support > existing_support:
                            # Remove old one, keep new one
                            existing_idx = existing_r.get("_idx")
                            if existing_idx is not None:
                                existing_g_a, existing_g_b, existing_fusion_name = normalized_names[existing_idx]
                                removed_fusions.append((existing_fusion_name, existing_support))
                            # Update mapping
                            bp_to_result[existing_key] = r
                            r["_idx"] = idx
                        else:
                            removed_fusions.append((fusion_name, current_support))
                        found_dup = True
                        break
                if found_dup:
                    break
        
        if not found_dup:
            bp_to_result[fusion_key] = r
            r["_idx"] = idx
            bp_set.add(bp_a)
            bp_set.add(bp_b)
    
    # Log trace for removed fusions
    if traced_pairs and trace_logger:
        for fusion_name, support_val in removed_fusions:
            g_a_norm, g_b_norm = fusion_name.split('-', 1)
            if ((g_a_norm, g_b_norm) in traced_pairs or (g_b_norm, g_a_norm) in traced_pairs):
                trace_logger.log_filter_result(fusion_name, "duplicate_breakpoints", False, f"Exact duplicate breakpoint with support={support_val} (not highest)")
    
    filtered_names = [name for name, _ in removed_fusions]
    if filtered_names and filter_log_file:
        filter_log_file.write(f"  Filtered fusions: {', '.join(filtered_names)}\n")
    deduped = list(bp_to_result.values())
    return deduped, filtered_names


def _filter_gene_name_artifacts(results, user_targets=None, traced_pairs=None, trace_logger=None, filter_log_file=None):
    """
    Filter gene naming artifacts (e.g., BORCS8-MEF2B-MEF2B where gene name contains hyphen).
    These are typically annotation issues, not real fusions.
    """
    if not results:
        return results, []
    
    filtered = []
    filtered_names = []
    
    for r in results:
        g_a = str(r.get("gene_a", ""))
        g_b = str(r.get("gene_b", ""))
        fusion_name = f"{g_a}-{g_b}"
        
        # Check for problematic patterns:
        # 1. Gene name contains hyphen (e.g., BORCS8-MEF2B-MEF2B)
        # 2. Gene name is very long (likely annotation artifact)
        # 3. Same gene appears twice (e.g., MEF2B-MEF2B)
        
        g_a_has_hyphen = '-' in g_a and g_a.count('-') > 0
        g_b_has_hyphen = '-' in g_b and g_b.count('-') > 0
        g_a_too_long = len(g_a) > 30  # Very long gene names are often artifacts
        g_b_too_long = len(g_b) > 30
        same_gene = g_a.upper() == g_b.upper()
        
        if same_gene:
            filtered_names.append(fusion_name)
            is_traced = traced_pairs and ((g_a, g_b) in traced_pairs or (g_b, g_a) in traced_pairs)
            if is_traced and trace_logger:
                trace_logger.log_filter_result(fusion_name, "gene_name_artifacts", False, "Same gene on both sides")
            continue
        
        # Filter if gene name has hyphens (annotation artifact) unless it's a known pattern
        # Some real genes have hyphens (e.g., C8ORF37-AS1), so be careful
        # Filter if BOTH genes have hyphens OR if one has multiple hyphens
        if (g_a_has_hyphen and g_b_has_hyphen) or g_a.count('-') > 1 or g_b.count('-') > 1:
            filtered_names.append(fusion_name)
            is_traced = traced_pairs and ((g_a, g_b) in traced_pairs or (g_b, g_a) in traced_pairs)
            if is_traced and trace_logger:
                trace_logger.log_filter_result(fusion_name, "gene_name_artifacts", False, f"Gene name artifact (hyphens in both or multiple hyphens)")
            continue
        
        # Filter if gene name is suspiciously long (likely annotation artifact)
        if g_a_too_long or g_b_too_long:
            filtered_names.append(fusion_name)
            is_traced = traced_pairs and ((g_a, g_b) in traced_pairs or (g_b, g_a) in traced_pairs)
            if is_traced and trace_logger:
                trace_logger.log_filter_result(fusion_name, "gene_name_artifacts", False, f"Gene name too long (likely annotation artifact)")
            continue
        
        filtered.append(r)
    
    if filtered_names and filter_log_file:
        filter_log_file.write(f"  Filtered fusions: {', '.join(filtered_names)}\n")
    
    return filtered, filtered_names


def _deduplicate_gene_aliases(results, user_targets=None, traced_pairs=None, trace_logger=None, filter_log_file=None):
    """Remove duplicates where same breakpoint has different gene names (e.g., MYCN vs MYCNOS)."""
    if not results:
        return results, []
    
    # Pre-normalize gene names once
    normalized_names = {}
    for idx, r in enumerate(results):
        g_a_norm = _normalize_gene_name(r.get("gene_a", ""))
        g_b_norm = _normalize_gene_name(r.get("gene_b", ""))
        normalized_names[idx] = (g_a_norm, g_b_norm, f"{g_a_norm}-{g_b_norm}")
    
    # Group by breakpoint, keep the one with highest support
    bp_to_result = {}  # (bp_key, g_a_norm, g_b_norm) -> result dict
    removed_fusions = set()
    for idx, r in enumerate(results):
        bp_key = (
            r.get("chrom_a", ""),
            r.get("pos_a", 0),
            r.get("chrom_b", ""),
            r.get("pos_b", 0),
        )
        g_a_norm, g_b_norm, fusion_name = normalized_names[idx]
        norm_key = (bp_key, g_a_norm, g_b_norm)
        
        if norm_key not in bp_to_result:
            bp_to_result[norm_key] = r
        else:
            # Keep the one with higher support
            existing_support = bp_to_result[norm_key].get("support", 0)
            current_support = r.get("support", 0)
            if current_support > existing_support:
                # Remove old one, keep new one
                existing_idx = bp_to_result[norm_key].get("_idx")
                if existing_idx is not None:
                    _, _, existing_fusion_name = normalized_names[existing_idx]
                    removed_fusions.add(existing_fusion_name)
                bp_to_result[norm_key] = r
                r["_idx"] = idx
            else:
                removed_fusions.add(fusion_name)
    
    # Log trace for removed fusions
    if traced_pairs and trace_logger:
        for fusion_name in removed_fusions:
            g_a_norm, g_b_norm = fusion_name.split('-', 1)
            if ((g_a_norm, g_b_norm) in traced_pairs or (g_b_norm, g_a_norm) in traced_pairs):
                trace_logger.log_filter_result(fusion_name, "deduplicate_gene_aliases", False, "Duplicate breakpoint with lower support")
    
    filtered_names = sorted(removed_fusions) if removed_fusions else []
    deduped = list(bp_to_result.values())
    return deduped, filtered_names


def _filter_low_confidence(results, user_targets, traced_pairs=None, trace_logger=None, filter_log_file=None):
    """
    Drop low-confidence calls with stricter thresholds for inter-chromosomal fusions.
    Most FPs (61%) have support < 25, so we filter more aggressively.
    """
    if not results:
        return results, []
    user_targets = set(user_targets) if user_targets else set()
    filtered = []
    filtered_names = []
    for r in results:
        support = r.get("support", 0)
        sa_frac = r.get("sa_fraction", 1.0)
        chrom_a = r.get("chrom_a", "")
        chrom_b = r.get("chrom_b", "")
        g_a = _normalize_gene_name(r.get("gene_a", ""))
        g_b = _normalize_gene_name(r.get("gene_b", ""))
        fusion_name = f"{g_a}-{g_b}"
        is_target = g_a in user_targets or g_b in user_targets
        is_traced = traced_pairs and ((g_a, g_b) in traced_pairs or (g_b, g_a) in traced_pairs)
        is_inter_chr = chrom_a != chrom_b
        
        if is_target:
            if is_traced and trace_logger:
                trace_logger.log_filter_result(fusion_name, "low_confidence", True, "Target gene (bypasses filter)")
            filtered.append(r)
            continue
        
        # Stricter filtering for inter-chromosomal fusions (most FPs are inter-chr)
        # But preserve high-support fusions like ROS1-SLC34A2 (support 42)
        if is_inter_chr:
            # High support (>= 35) - always keep (true fusions often have high support)
            if support >= 35:
                if is_traced and trace_logger:
                    trace_logger.log_filter_result(fusion_name, "low_confidence", True, f"Inter-chromosomal: high support={support} >= 35 (preserved)")
                filtered.append(r)
                continue
            
            # Inter-chromosomal: require higher support OR good SA evidence
            # Support 24: slightly more lenient (SA >= 0.15) to preserve borderline true fusions
            # Support < 24: stricter (SA >= 0.20) to filter low-quality artifacts
            if support < 24:
                # Very low support (< 24) - filter unless high SA
                if sa_frac < 0.20:
                    filtered_names.append(fusion_name)
                    if is_traced and trace_logger:
                        trace_logger.log_filter_result(fusion_name, "low_confidence", False, f"Inter-chromosomal: support={support} < 24 and SA fraction={sa_frac:.2f} < 0.20")
                    continue
            elif support == 24:
                # Support exactly 24: slightly more lenient SA requirement (0.15) to preserve borderline true fusions
                if sa_frac < 0.15:
                    filtered_names.append(fusion_name)
                    if is_traced and trace_logger:
                        trace_logger.log_filter_result(fusion_name, "low_confidence", False, f"Inter-chromosomal: support={support} == 24 and SA fraction={sa_frac:.2f} < 0.15")
                    continue
            elif support < 35:
                # Moderate support (25-34) - require decent SA evidence
                if sa_frac < 0.15:
                    filtered_names.append(fusion_name)
                    if is_traced and trace_logger:
                        trace_logger.log_filter_result(fusion_name, "low_confidence", False, f"Inter-chromosomal: support={support} < 35 and SA fraction={sa_frac:.2f} < 0.15")
                    continue
        else:
            # Same-chromosome: less strict (true inversions can have lower SA)
            if support < 22:
                filtered_names.append(fusion_name)
                if is_traced and trace_logger:
                    trace_logger.log_filter_result(fusion_name, "low_confidence", False, f"Same-chromosome: support={support} < 22")
                continue
        
        if is_traced and trace_logger:
            trace_logger.log_filter_result(fusion_name, "low_confidence", True, f"Support={support} and SA fraction={sa_frac:.2f} meet thresholds")
        filtered.append(r)
    
    return filtered, filtered_names


# Global threshold: same for all samples (genes appearing this many times = repeating)
# Lowered from 3 to 2 to catch genes like LINC00486 that appear in many fusions
REPEATING_GENE_THRESHOLD = 2


def _filter_exon_boundaries(results, exon_boundaries, user_targets, traced_pairs=None, trace_logger=None, filter_log_file=None):
    """
    Post-processing filter: Check if breakpoints are near exon boundaries.
    Real fusions typically occur at exon boundaries. Filter breakpoints far from exon boundaries
    if they also have low support/SA evidence (conservative: only filter if multiple conditions met).
    
    NOTE: Currently disabled by default - may be too aggressive. Enable by uncommenting in process_sample.
    """
    if not results or not exon_boundaries:
        return results, []
    
    EXON_BOUNDARY_WINDOW = 20  # Consider "near" if within 20bp of exon boundary
    MIN_DIST_FROM_BOUNDARY = 150  # Less strict: filter if >150bp from boundary (was 80bp - too aggressive)
    MIN_SUPPORT_THRESHOLD = 22  # Less strict: require support >= 22 if far from boundaries (was 28 - too aggressive)
    
    filtered = []
    filtered_names = []
    for r in results:
        chrom_a = r.get("chrom_a", "")
        pos_a = r.get("pos_a", 0)
        chrom_b = r.get("chrom_b", "")
        pos_b = r.get("pos_b", 0)
        support = r.get("support", 0)
        sa_frac = r.get("sa_fraction", 1.0)
        g_a = _normalize_gene_name(r.get("gene_a", ""))
        g_b = _normalize_gene_name(r.get("gene_b", ""))
        fusion_name = f"{g_a}-{g_b}"
        is_target = g_a in user_targets or g_b in user_targets
        is_traced = traced_pairs and ((g_a, g_b) in traced_pairs or (g_b, g_a) in traced_pairs)
        
        # Skip target genes
        if is_target:
            if is_traced and trace_logger:
                trace_logger.log_filter_result(fusion_name, "exon_boundaries", True, "Target gene (bypasses filter)")
            filtered.append(r)
            continue
        
        # Check both breakpoints
        dist_a, near_a = get_nearest_exon_boundary(exon_boundaries, chrom_a, pos_a, EXON_BOUNDARY_WINDOW)
        dist_b, near_b = get_nearest_exon_boundary(exon_boundaries, chrom_b, pos_b, EXON_BOUNDARY_WINDOW)
        
        # More lenient filter: only filter if BOTH breakpoints are very far from boundaries AND support is very low
        if dist_a is not None and dist_b is not None:
            both_far_from_boundary = (dist_a > MIN_DIST_FROM_BOUNDARY and dist_b > MIN_DIST_FROM_BOUNDARY)
            if both_far_from_boundary:
                # Only filter if support is very low (less strict approach)
                if support < MIN_SUPPORT_THRESHOLD:
                    reason = f"Both breakpoints far from boundaries (dist_a={dist_a}bp, dist_b={dist_b}bp) and support={support} < {MIN_SUPPORT_THRESHOLD}"
                    filtered_names.append(f"{fusion_name} (support={support}, dist_a={dist_a}bp, dist_b={dist_b}bp)")
                    if is_traced and trace_logger:
                        trace_logger.log_filter_result(fusion_name, "exon_boundaries", False, reason)
                    continue
                # Also filter if moderate support but very low SA fraction
                if support < 30 and sa_frac < 0.10:  # Less strict: support < 30 (was 38) and SA < 0.10 (was 0.12)
                    reason = f"Both breakpoints far from boundaries (dist_a={dist_a}bp, dist_b={dist_b}bp), support={support} < 30 and SA fraction={sa_frac:.2f} < 0.10"
                    filtered_names.append(f"{fusion_name} (support={support}, SA={sa_frac:.2f}, dist_a={dist_a}bp, dist_b={dist_b}bp)")
                    if is_traced and trace_logger:
                        trace_logger.log_filter_result(fusion_name, "exon_boundaries", False, reason)
                    continue
        
        # Passed exon boundary check
        if is_traced and trace_logger:
            if dist_a is not None and dist_b is not None:
                trace_logger.log_filter_result(fusion_name, "exon_boundaries", True, f"Near exon boundaries (dist_a={dist_a}bp, dist_b={dist_b}bp)")
            else:
                trace_logger.log_filter_result(fusion_name, "exon_boundaries", True, "Exon boundary data not available")
        filtered.append(r)
    
    return filtered, filtered_names


def _filter_target_gene_quality(results, user_targets, traced_pairs=None, trace_logger=None, filter_log_file=None):
    """
    Filter target gene fusions to keep only high-quality ones.
    When a target gene (e.g., ROS1) appears in many fusions, only keep the best ones.
    
    Strategy:
    - Keep top N fusions per target gene (by support)
    - Apply quality thresholds (support, SA fraction)
    - Keep unique breakpoints (if same breakpoint, keep highest support)
    """
    if not results or not user_targets:
        return results, []
    
    user_targets = set(user_targets)
    filtered = []
    filtered_names = []
    
    # Group target gene fusions by target gene
    target_fusions = {}  # target_gene -> list of (fusion_dict, fusion_name, support, sa_frac)
    non_target_fusions = []
    
    for r in results:
        g_a = _normalize_gene_name(r.get("gene_a", ""))
        g_b = _normalize_gene_name(r.get("gene_b", ""))
        fusion_name = f"{g_a}-{g_b}"
        support = r.get("support", 0)
        sa_frac = r.get("sa_fraction", 0.0)
        
        is_target_a = g_a in user_targets
        is_target_b = g_b in user_targets
        
        if is_target_a or is_target_b:
            # This fusion involves a target gene
            target_gene = g_a if is_target_a else g_b
            if target_gene not in target_fusions:
                target_fusions[target_gene] = []
            target_fusions[target_gene].append((r, fusion_name, support, sa_frac))
        else:
            # Non-target fusion - keep as-is
            non_target_fusions.append(r)
    
    # Process each target gene's fusions
    MIN_TARGET_SUPPORT = 25  # Minimum support for target gene fusions
    MIN_TARGET_SA = 0.15  # Minimum SA fraction for target gene fusions
    MAX_TARGET_FUSIONS_PER_GENE = 5  # Keep top 5 fusions per target gene
    
    for target_gene, fusions in target_fusions.items():
        # Sort by support (descending), then by SA fraction
        sorted_fusions = sorted(fusions, key=lambda x: (x[2], x[3]), reverse=True)
        
        # Filter by quality thresholds
        quality_fusions = [
            (r, name, sup, sa) for r, name, sup, sa in sorted_fusions
            if sup >= MIN_TARGET_SUPPORT and sa >= MIN_TARGET_SA
        ]
        
        # If no fusions meet quality threshold, keep top 2 anyway (might be real but low quality)
        if not quality_fusions and sorted_fusions:
            quality_fusions = sorted_fusions[:2]
        
        # Keep top N fusions per target gene
        kept_fusions = quality_fusions[:MAX_TARGET_FUSIONS_PER_GENE]
        kept_names = {name for _, name, _, _ in kept_fusions}
        
        # Track filtered ones
        for r, name, sup, sa in sorted_fusions:
            if name in kept_names:
                filtered.append(r)
            else:
                filtered_names.append(name)
                is_traced = traced_pairs and ((name.split('-')[0], name.split('-')[1]) in traced_pairs or (name.split('-')[1], name.split('-')[0]) in traced_pairs)
                if is_traced and trace_logger:
                    reason = f"Target gene {target_gene}: kept top {len(kept_fusions)} fusions (support={sup}, SA={sa:.2f} not in top {MAX_TARGET_FUSIONS_PER_GENE})"
                    trace_logger.log_filter_result(name, "target_gene_quality", False, reason)
    
    # Add non-target fusions
    filtered.extend(non_target_fusions)
    
    if filtered_names and filter_log_file:
        filter_log_file.write(f"  Filtered target gene fusions: {', '.join(filtered_names)}\n")
    
    return filtered, filtered_names


def _filter_repeating_genes(results, user_targets, traced_pairs=None, trace_logger=None, filter_log_file=None):
    """
    Global filter: genes appearing >= REPEATING_GENE_THRESHOLD times (in either column)
    are treated as repeating. Drop fusions where both genes are repeating, or one is
    repeating with low support/SA. Same threshold for every sample.
    """
    if not results:
        return results, []
    
    # Pre-normalize gene names once
    normalized_names = {}
    for idx, r in enumerate(results):
        g_a = _normalize_gene_name(r.get("gene_a", ""))
        g_b = _normalize_gene_name(r.get("gene_b", ""))
        normalized_names[idx] = (g_a, g_b, f"{g_a}-{g_b}")
    
    # Count each gene in both columns (global rule: use both gene_a and gene_b)
    gene_counts = {}
    for idx, r in enumerate(results):
        g_a, g_b, _ = normalized_names[idx]
        gene_counts[g_a] = gene_counts.get(g_a, 0) + 1
        gene_counts[g_b] = gene_counts.get(g_b, 0) + 1
    
    repeating = {g for g, c in gene_counts.items() if c >= REPEATING_GENE_THRESHOLD}
    if not repeating:
        return results, []
    
    user_targets = set(user_targets) if user_targets else set()
    filtered = []
    filtered_names = []
    
    for idx, r in enumerate(results):
        g_a, g_b, fusion_name = normalized_names[idx]
        support = r.get("support", 0)
        sa_frac = r.get("sa_fraction", 1.0)
        is_target = g_a in user_targets or g_b in user_targets
        is_traced = traced_pairs and ((g_a, g_b) in traced_pairs or (g_b, g_a) in traced_pairs)
        both_repeating = (g_a in repeating and g_b in repeating)
        one_repeating = (g_a in repeating or g_b in repeating)
        
        if is_target:
            if is_traced and trace_logger:
                trace_logger.log_filter_result(fusion_name, "repeating_genes", True, "Target gene (bypasses filter)")
            filtered.append(r)
            continue
        # More aggressive filtering for repeating genes (catches LINC00486, ALK clusters, etc.)
        if both_repeating:
            # Both genes repeating - filter unless very high support
            if support < 50:
                filtered_names.append(fusion_name)
                if is_traced and trace_logger:
                    trace_logger.log_filter_result(fusion_name, "repeating_genes", False, f"Both genes repeating ({g_a} count={gene_counts[g_a]}, {g_b} count={gene_counts[g_b]}), support={support} < 50")
                continue
        elif one_repeating:
            # One gene repeating - filter unless moderate-high support AND good SA evidence
            repeating_gene = g_a if g_a in repeating else g_b
            if support < 40 or (support < 50 and sa_frac < 0.30):
                filtered_names.append(fusion_name)
                if is_traced and trace_logger:
                    trace_logger.log_filter_result(fusion_name, "repeating_genes", False, f"{repeating_gene} repeating (count={gene_counts[repeating_gene]}), support={support} < 40 or (support < 50 and SA fraction={sa_frac:.2f} < 0.30)")
                continue
        if is_traced and trace_logger:
            reason_parts = []
            if both_repeating:
                reason_parts.append(f"Both repeating but support={support} >= 40")
            elif one_repeating:
                repeating_gene = g_a if g_a in repeating else g_b
                reason_parts.append(f"{repeating_gene} repeating but support={support} >= 35 or SA={sa_frac:.2f} >= 0.25")
            else:
                reason_parts.append(f"Neither gene repeating ({g_a} count={gene_counts.get(g_a, 0)}, {g_b} count={gene_counts.get(g_b, 0)})")
            trace_logger.log_filter_result(fusion_name, "repeating_genes", True, ", ".join(reason_parts))
        filtered.append(r)
    
    return filtered, filtered_names


def _process_one_gene_pair(bam_path, fasta_path, ref_path, sample_name, gene_pair, reads_dicts, user_targets, min_support, traced_pairs=None, trace_log_path=None):
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
    
    # Normalize for tracing
    g_a_norm = _normalize_gene_name(g_a)
    g_b_norm = _normalize_gene_name(g_b)
    fusion_name = f"{g_a_norm}-{g_b_norm}"
    is_traced = traced_pairs and ((g_a_norm, g_b_norm) in traced_pairs or (g_b_norm, g_a_norm) in traced_pairs)
    
    # Create trace logger if needed (each worker opens file in append mode)
    trace_logger = None
    if is_traced and trace_log_path:
        trace_logger = TraceLogger(trace_log_path)  # Will append to existing file
    
    if is_traced and trace_logger:
        trace_logger.log(f"[TRACE] {fusion_name}: Discovery - {len(reads)} reads found")

    # Target filter: if user_targets is provided and non-empty, only keep fusions involving target genes
    if user_targets and len(user_targets) > 0:
        if g_a not in user_targets and g_b not in user_targets:
            if is_traced and trace_logger:
                trace_logger.log(f"[TRACE] {fusion_name}: ✗ DROPPED: Not in target gene list")
            return None

    is_target_match = bool(user_targets and (g_a in user_targets or g_b in user_targets))
    bam_handle = None
    assembler = None
    try:
        if not os.path.exists(bam_path):
            if is_traced and trace_logger:
                trace_logger.log(f"[TRACE] {fusion_name}: ✗ DROPPED: BAM file not found")
            return None
        bam_handle = pysam.AlignmentFile(bam_path, "rb")
        assembler = FusionAssembler(bam_path, fasta_path=fasta_path)
        validator = FusionValidator(fasta_path)
        breakpoints, support = assembler.find_breakpoint(reads, bam=bam_handle)
        if not breakpoints or support < min_support:
            if is_traced and trace_logger:
                reason = "No breakpoints found" if not breakpoints else f"Support={support} < min_support={min_support}"
                trace_logger.log(f"[TRACE] {fusion_name}: ✗ DROPPED: {reason}")
            return None
        bp_a, bp_b = breakpoints
        
        # Calculate SA tag fraction (split-read evidence)
        sa_count = sum(1 for r in reads if r.has_tag("SA"))
        sa_fraction = sa_count / len(reads) if reads else 0.0
        
        # JuLI-style filtering: For inter-chromosomal fusions, require BOTH discordant pairs AND split-reads
        # This significantly reduces FPs while preserving true fusions like ROS1-SLC34A2
        # (True fusions typically have both types of evidence)
        chrom_a = bp_a[0]
        chrom_b = bp_b[0]
        is_inter_chr = chrom_a != chrom_b
        
        if is_inter_chr:
            # Count discordant pairs (reads without SA tags)
            discordant_count = len(reads) - sa_count
            # JuLI-style filtering: Require BOTH discordant pairs AND split-reads for inter-chromosomal fusions
            # But allow exceptions for high-support fusions (support >= 40) with good SA fraction (>= 0.15)
            # This preserves true fusions like ROS1-SLC34A2 while filtering most FPs
            MIN_DISCORDANT = 3
            MIN_SPLIT_READS = 2
            
            # High-support exception: if support >= 40 AND SA >= 0.15, allow even if split-reads < 2
            # (True fusions can have moderate SA fraction but high overall support)
            high_support_exception = support >= 40 and sa_fraction >= 0.15
            
            if not high_support_exception:
                # Standard requirement: BOTH discordant pairs >= 3 AND split-reads >= 2
                if discordant_count < MIN_DISCORDANT or sa_count < MIN_SPLIT_READS:
                    if is_traced and trace_logger:
                        reason = f"Inter-chromosomal: requires BOTH discordant pairs >= {MIN_DISCORDANT} AND split-reads >= {MIN_SPLIT_READS}, but got discordant={discordant_count}, split-reads={sa_count}"
                        trace_logger.log(f"[TRACE] {fusion_name}: ✗ DROPPED: {reason}")
                    return None
            else:
                # High-support exception: still require discordant pairs >= 3
                if discordant_count < MIN_DISCORDANT:
                    if is_traced and trace_logger:
                        reason = f"Inter-chromosomal: high-support exception (support={support} >= 40, SA={sa_fraction:.2f} >= 0.15) but discordant pairs {discordant_count} < {MIN_DISCORDANT}"
                        trace_logger.log(f"[TRACE] {fusion_name}: ✗ DROPPED: {reason}")
                    return None
                if is_traced and trace_logger:
                    trace_logger.log(f"[TRACE] {fusion_name}: Dual evidence - high-support exception (support={support} >= 40, SA={sa_fraction:.2f} >= 0.15)")
        
        if is_traced and trace_logger:
            trace_logger.log(f"[TRACE] {fusion_name}: Assembly - support={support}, {bp_a[0]}:{bp_a[1]} + {bp_b[0]}:{bp_b[1]}, sa_fraction={sa_fraction:.2f}")
        
        if is_target_match:
            if is_traced and trace_logger:
                trace_logger.log(f"[TRACE] {fusion_name}: Validator - PASSED: Target gene match (bypasses validator)")
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
        try:
            is_fp, validator_reason = validator.is_likely_fp(g_a, g_b, bp_a[0], bp_b[0], bp_a[1], bp_b[1], support=support, sa_fraction=sa_fraction)
            
        except Exception as validator_error:
            # If validator throws an exception, log it and treat as FP
            import traceback
            validator_reason = f"VALIDATOR EXCEPTION: {str(validator_error)}\n{traceback.format_exc()}"
            is_fp = True
            if is_traced and trace_logger:
                trace_logger.log(f"[TRACE] {fusion_name}: Validator exception - {validator_reason}")
            # Also print to stderr for debugging
            print(f"[!] Validator exception for {fusion_name}: {validator_error}", file=sys.stderr, flush=True)
            print(f"[!] Full traceback:\n{traceback.format_exc()}", file=sys.stderr, flush=True)
        
        if is_traced and trace_logger:
            trace_logger.log(f"[TRACE] {fusion_name}: Validator - {validator_reason}")
        if is_fp:
            if is_traced and trace_logger:
                trace_logger.log(f"[TRACE] {fusion_name}: ✗ DROPPED by: validator - {validator_reason}")
            return None
        if is_traced and trace_logger:
            trace_logger.log(f"[TRACE] {fusion_name}: ✓ PASSED validator, proceeding to post-processing filters")
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
    except Exception as e:
        if is_traced and trace_logger:
            trace_logger.log(f"[TRACE] {fusion_name}: ✗ DROPPED: Exception during processing - {str(e)}")
        return None
    finally:
        if bam_handle is not None:
            bam_handle.close()
        # Clean up assembler's cached FASTA handle
        if assembler and assembler._fasta_handle is not None:
            try:
                assembler._fasta_handle.close()
                assembler._fasta_handle = None
            except Exception:
                pass


def process_sample(bam_path, ref_path, fasta_path, output_dir, user_targets, n_workers=1, min_mapq=20, use_process_discovery=False, min_support=None, traced_pairs=None):
    if min_support is None:
        min_support = MIN_SUPPORT
    sample_name = os.path.basename(bam_path).replace(".bam", "")
    print(f"[*] Processing: {sample_name} (n_workers={n_workers}, min_support={min_support}, discovery={'processes' if use_process_discovery else 'threads'})", flush=True)

    # Create trace logger if tracing is enabled
    trace_logger = None
    trace_log_path = None
    if traced_pairs:
        trace_log_path = os.path.join(output_dir, f"{sample_name}_trace.log")
        trace_logger = TraceLogger(trace_log_path)
        print(f"[*] Trace log: {trace_log_path}", flush=True)
    
    # Create general filter log file (shows what filters are removing)
    filter_log_path = os.path.join(output_dir, f"{sample_name}_filter_log.txt")
    filter_log_file = None

    try:
        idx = GenomeIndex()
        idx.load_refflat(ref_path)
        # Load exon boundaries once for post-processing filter (only for final fusions)
        exon_boundaries = load_exon_boundaries(ref_path)
        validator = FusionValidator(fasta_path)
        discoverer = FusionDiscoverer(bam_path, idx, user_targets=user_targets, n_workers=n_workers, min_mapq=min_mapq, use_process_discovery=use_process_discovery)
        t0 = time.perf_counter()
        raw_candidates = discoverer.collect_seeds()
        print(f"[*] [{sample_name}] Discovery done in {time.perf_counter() - t0:.1f}s ({len(raw_candidates)} gene pairs)", flush=True)

        # Build list of (gene_pair, reads_dicts) for pool workers; serialize reads to dicts for pickling
        items = []
        for gene_pair, reads in raw_candidates.items():
            if not reads:
                continue
            reads_dicts = [_serialize_read(r) for r in reads]
            items.append((bam_path, fasta_path, ref_path, sample_name, gene_pair, reads_dicts, tuple(user_targets) if user_targets else (), min_support, traced_pairs, trace_log_path))

        if not items:
            final_results = []
            # Still create filter log even if no items
            filter_log_file = open(filter_log_path, 'w')
            filter_log_file.write(f"# Filter Log for {sample_name}\n")
            filter_log_file.write(f"# Shows how many fusions pass/fail each filter\n\n")
            filter_log_file.write(f"Total gene pairs from discovery: {len(raw_candidates)}\n")
            filter_log_file.write(f"Total gene pairs processed (with reads): 0\n")
            filter_log_file.write(f"Fusions passing validator (before post-processing): 0\n\n")
            filter_log_file.write(f"Before post-processing filters: 0 fusions\n")
            filter_log_file.write(f"Final result: 0 fusions reported\n")
            filter_log_file.close()
            print(f"[*] No gene pairs with reads found. Filter log saved: {filter_log_path}", flush=True)
        else:
            total = len(items)
            print(f"[*] [{sample_name}] Processing {total} gene pairs (assembly + validation)...", flush=True)
            t0 = time.perf_counter()
            if n_workers > 1:
                # Check if we're running inside a multiprocessing worker (nested pool)
                current_process = multiprocessing.current_process()
                is_worker_process = current_process.name != "MainProcess"
                if is_worker_process:
                    # In a worker process: use threads for gene-pair processing (can't create nested pools)
                    # This allows cores to be utilized even when processing multiple BAMs
                    n_workers_use = min(n_workers, len(items))
                    validator_passed = 0
                    validator_failed = 0
                    with ThreadPoolExecutor(max_workers=n_workers_use) as executor:
                        futures = [executor.submit(_process_one_gene_pair, *args) for args in items]
                        results = []
                        # Process results as they complete (like imap_unordered) to show progress
                        for done, fut in enumerate(as_completed(futures), 1):
                            r = fut.result()
                            results.append(r)
                            if r is not None:
                                validator_passed += 1
                            else:
                                validator_failed += 1
                            # Print progress every 10% or every 1000 items, whichever is more frequent
                            if done % max(1, min(1000, total // 10)) == 0 or done == total:
                                print(f"[*] [{sample_name}] Processed {done}/{total} gene pairs ({validator_passed} passed, {validator_failed} failed)...", flush=True)
                    final_results = [r for r in results if r is not None]
                    # Progress logged to filter log file only (no console spam)
                    # Validator summary goes to log file only
                    print(f"[*] [{sample_name}] Assembly/validation done in {time.perf_counter() - t0:.1f}s ({validator_passed} passed, {validator_failed} failed)", flush=True)
                    if validator_failed > 0 and validator_passed == 0:
                        print(f"[!] [{sample_name}] WARNING: All fusions failed validator! Use --trace to debug specific pairs.", flush=True)
                else:
                    # Not in a worker process: use processes for gene-pair processing (better for CPU-bound work)
                    n_workers_use = min(n_workers, len(items))
                    validator_passed = 0
                    validator_failed = 0
                    with multiprocessing.Pool(n_workers_use) as pool:
                        results = []
                        for done, r in enumerate(pool.imap_unordered(_process_one_item, items, chunksize=1), 1):
                            results.append(r)
                            if r is not None:
                                validator_passed += 1
                            else:
                                validator_failed += 1
                            # Progress logged to filter log file only (no console spam)
                    final_results = [r for r in results if r is not None]
                    # Validator summary goes to log file only
                    print(f"[*] [{sample_name}] Assembly/validation done in {time.perf_counter() - t0:.1f}s ({validator_passed} passed, {validator_failed} failed)", flush=True)
                    if validator_failed > 0 and validator_passed == 0:
                        print(f"[!] WARNING: All fusions failed validator! Use --trace to debug specific pairs.", flush=True)
            else:
                final_results = []
                validator_passed = 0
                validator_failed = 0
                failure_reasons = {}  # Track why fusions are failing
                sample_failures = []  # Store first 10 failure reasons for debugging
                for i, args in enumerate(items):
                    r = _process_one_gene_pair(*args)
                    if r is not None:
                        final_results.append(r)
                        validator_passed += 1
                    else:
                        validator_failed += 1
                        # Sample first 10 failures to understand why they're failing
                        if len(sample_failures) < 10:
                            gene_pair = args[4]  # gene_pair is at index 4
                            raw_a = gene_pair[0][0] if isinstance(gene_pair[0], (list, tuple)) else gene_pair[0]
                            raw_b = gene_pair[1][0] if isinstance(gene_pair[1], (list, tuple)) else gene_pair[1]
                            g_a = str(raw_a).split("(")[0].strip().upper()
                            g_b = str(raw_b).split("(")[0].strip().upper()
                            sample_failures.append(f"{g_a}-{g_b}")
                    # Progress logged to filter log file only (no console spam)
                # Validator summary goes to log file only
                print(f"[*] [{sample_name}] Assembly/validation done in {time.perf_counter() - t0:.1f}s ({validator_passed} passed, {validator_failed} failed)", flush=True)
                if sample_failures and validator_passed == 0:
                    print(f"[!] [{sample_name}] Sample failed gene pairs (first 10): {', '.join(sample_failures)}", flush=True)
                    print(f"[!] [{sample_name}] Use --trace to debug specific pairs, e.g., --trace {sample_failures[0]}", flush=True)
            
            # Post-process: global filters only (same thresholds for all samples)
            print(f"[*] [{sample_name}] Applying post-processing filters ({len(final_results)} fusions)...", flush=True)
            # Create filter log file to track what's being filtered
            filter_log_file = open(filter_log_path, 'w')
            filter_log_file.write(f"# Filter Log for {sample_name}\n")
            filter_log_file.write(f"# Shows how many fusions pass/fail each filter\n\n")
            filter_log_file.write(f"Validator summary: {validator_passed} passed, {validator_failed} failed (out of {total} gene pairs)\n\n")
            filter_log_file.write(f"Total gene pairs from discovery: {len(raw_candidates)}\n")
            filter_log_file.write(f"Total gene pairs processed (with reads): {total}\n")
            filter_log_file.write(f"Fusions passing validator (before post-processing): {len(final_results)}\n\n")
            
            filter_log_file.write(f"Before post-processing filters: {len(final_results)} fusions\n")
            before_count = len(final_results)
            final_results, filtered_names = _filter_gene_name_artifacts(final_results, user_targets, traced_pairs, trace_logger, filter_log_file)
            removed_count = before_count - len(final_results)
            filter_log_file.write(f"After gene_name_artifacts: {len(final_results)} fusions (removed {removed_count})\n")
            if filtered_names:
                filter_log_file.write(f"  Filtered fusions: {', '.join(filtered_names)}\n")
            before_count = len(final_results)
            final_results, filtered_names = _deduplicate_gene_aliases(final_results, user_targets, traced_pairs, trace_logger, filter_log_file)
            removed_count = before_count - len(final_results)
            filter_log_file.write(f"After deduplicate_gene_aliases: {len(final_results)} fusions (removed {removed_count})\n")
            if filtered_names:
                filter_log_file.write(f"  Filtered fusions: {', '.join(filtered_names)}\n")
            before_count = len(final_results)
            final_results, filtered_names = _filter_sink_breakpoints(final_results, user_targets, traced_pairs, trace_logger, filter_log_file)
            removed_count = before_count - len(final_results)
            filter_log_file.write(f"After sink_breakpoints: {len(final_results)} fusions (removed {removed_count})\n")
            if filtered_names:
                filter_log_file.write(f"  Filtered fusions: {', '.join(filtered_names)}\n")
            before_count = len(final_results)
            final_results, filtered_names = _filter_breakpoint_clusters(final_results, user_targets, traced_pairs, trace_logger, filter_log_file)
            removed_count = before_count - len(final_results)
            filter_log_file.write(f"After breakpoint_clusters: {len(final_results)} fusions (removed {removed_count})\n")
            if filtered_names:
                filter_log_file.write(f"  Filtered fusions: {', '.join(filtered_names)}\n")
            before_count = len(final_results)
            final_results, filtered_names = _filter_duplicate_breakpoints(final_results, user_targets, traced_pairs, trace_logger, filter_log_file)
            removed_count = before_count - len(final_results)
            filter_log_file.write(f"After duplicate_breakpoints: {len(final_results)} fusions (removed {removed_count})\n")
            if filtered_names:
                filter_log_file.write(f"  Filtered fusions: {', '.join(filtered_names)}\n")
            before_count = len(final_results)
            final_results, filtered_names = _filter_low_confidence(final_results, user_targets, traced_pairs, trace_logger, filter_log_file)
            removed_count = before_count - len(final_results)
            filter_log_file.write(f"After low_confidence: {len(final_results)} fusions (removed {removed_count})\n")
            if filtered_names:
                filter_log_file.write(f"  Filtered fusions: {', '.join(filtered_names)}\n")
            before_count = len(final_results)
            final_results, filtered_names = _filter_repeating_genes(final_results, user_targets, traced_pairs, trace_logger, filter_log_file)
            removed_count = before_count - len(final_results)
            filter_log_file.write(f"After repeating_genes: {len(final_results)} fusions (removed {removed_count})\n")
            if filtered_names:
                filter_log_file.write(f"  Filtered fusions: {', '.join(filtered_names)}\n")
            # Target gene quality filter: keep only high-quality target gene fusions
            before_count = len(final_results)
            final_results, filtered_names = _filter_target_gene_quality(final_results, user_targets, traced_pairs, trace_logger, filter_log_file)
            removed_count = before_count - len(final_results)
            filter_log_file.write(f"After target_gene_quality: {len(final_results)} fusions (removed {removed_count})\n")
            if filtered_names:
                filter_log_file.write(f"  Filtered fusions: {', '.join(filtered_names)}\n")
            # Exon boundary check: only check final fusions (post-processing, faster)
            before_count = len(final_results)
            final_results, filtered_names = _filter_exon_boundaries(final_results, exon_boundaries, user_targets, traced_pairs, trace_logger, filter_log_file)
            removed_count = before_count - len(final_results)
            filter_log_file.write(f"After exon_boundaries: {len(final_results)} fusions (removed {removed_count})\n")
            if filtered_names:
                # Write fusion names with details (support, distances) if available
                filter_log_file.write(f"  Filtered fusions:\n")
                for fusion_info in filtered_names:
                    filter_log_file.write(f"    {fusion_info}\n")
            filter_log_file.write(f"\nFinal result: {len(final_results)} fusions reported\n")
            filter_log_file.write(f"\n# Summary:\n")
            filter_log_file.write(f"# - Discovery found {len(raw_candidates)} gene pairs\n")
            filter_log_file.write(f"# - {total} gene pairs had reads and were processed\n")
            filter_log_file.write(f"# - Validator passed: {len(final_results)} fusions (before post-processing: {len(final_results) + sum([len(raw_candidates) - total, 0])})\n")
            filter_log_file.write(f"# - {len(final_results)} fusions passed validator and all post-processing filters\n")
            filter_log_file.close()
            print(f"[*] [{sample_name}] Post-processing done: {len(final_results)} fusions remaining", flush=True)
            
            # Log final results for traced fusions
            if traced_pairs and trace_logger:
                for r in final_results:
                    g_a_norm = _normalize_gene_name(r.get("gene_a", ""))
                    g_b_norm = _normalize_gene_name(r.get("gene_b", ""))
                    fusion_name = f"{g_a_norm}-{g_b_norm}"
                    if ((g_a_norm, g_b_norm) in traced_pairs or (g_b_norm, g_a_norm) in traced_pairs):
                        trace_logger.log(f"[TRACE] {fusion_name}: ✓ REPORTED in final output")

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
    parser.add_argument("--min-support", type=int, default=MIN_SUPPORT, metavar="N",
                        help=f"Minimum read support for a fusion to be reported (default: {MIN_SUPPORT}).")
    parser.add_argument("--target", help="Path to txt file with gene list (one per line)")
    parser.add_argument("--trace", help="Comma-separated gene pairs to trace (e.g., ROS1-SLC34A2,CCDC6-RET)")
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Load Target Genes
    user_targets = set()
    if args.target and os.path.exists(args.target):
        with open(args.target, "r") as f:
            user_targets = {line.strip().upper() for line in f if line.strip()}
        print(f"[*] Target Filter Active: Loaded {len(user_targets)} genes.")
    
    # Parse traced gene pairs
    traced_pairs = set()
    if args.trace:
        for pair_str in args.trace.split(','):
            pair_str = pair_str.strip()
            if '-' in pair_str:
                parts = pair_str.split('-', 1)
                if len(parts) == 2:
                    g_a = _normalize_gene_name(parts[0].strip())
                    g_b = _normalize_gene_name(parts[1].strip())
                    traced_pairs.add((g_a, g_b))
        if traced_pairs:
            print(f"[*] Trace Mode Active: Tracing {len(traced_pairs)} gene pair(s).")

    if os.path.isdir(args.input):
        bam_files = glob.glob(os.path.join(args.input, "*.bam"))
    else:
        bam_files = [args.input]

    print(f"[*] Found {len(bam_files)} BAM file(s).")
    if not bam_files:
        print("[*] No BAM files found. Exiting.")
        return
    
    # Check for BAM index files (.bai) - required for efficient BAM access
    missing_indices = []
    for bam_file in bam_files:
        bai_file = bam_file + ".bai"
        if not os.path.exists(bai_file):
            missing_indices.append(bam_file)
    
    if missing_indices:
        print("\n[!] ERROR: Missing BAM index files (.bai). Index files are required for fusion detection.")
        print("[!] Missing indices for:")
        for bam in missing_indices:
            print(f"[!]   - {bam}")
        print(f"\n[!] To create index files, run:")
        print(f"[!]   samtools index {missing_indices[0]}")
        print(f"[!] (Repeat for each BAM file)\n")
        sys.exit(1)
    
    # Filter is already installed at module level, but ensure it's active
    install_hts_warning_filter()
    n_workers = args.cores if args.cores is not None else os.cpu_count() or 4
    use_process_discovery = len(bam_files) == 1
    
    if len(bam_files) == 1:
        # Single BAM: use all cores for gene-pair processing
        n_processes = 1
        print(f"[*] Single BAM: using {n_workers} cores for gene-pair processing, discovery: processes (--cores={n_workers}).")
    else:
        # Multiple BAMs: distribute cores between BAM-level and gene-pair-level parallelism
        # Use all available cores for BAM processing to keep cores busy
        n_processes = min(n_workers, len(bam_files))
        # Each BAM process will use threads for discovery (to avoid nested pools)
        # and will use threads for gene-pair processing (since it's in a worker process)
        print(f"[*] Multiple BAMs: {n_processes} parallel BAM processes, discovery: threads, gene-pair processing: threads (--cores={n_workers}).")

    min_mapq = getattr(args, "min_mapq", 20)
    min_support = getattr(args, "min_support", MIN_SUPPORT)
    if n_processes == 1:
        # Single BAM: call directly to avoid nested pools (daemonic processes can't spawn children)
        # The per-gene-pair pool inside process_sample will still create N worker processes
        result = process_sample(bam_files[0], args.ref, args.fasta, args.outdir, user_targets, n_workers, min_mapq, use_process_discovery, min_support, traced_pairs)
        print(result)
    else:
        # Multiple BAMs: use Pool for parallel BAM processing
        # Discovery will use threads (use_process_discovery=False) to avoid nested pools
        # The pool automatically keeps cores busy by starting new tasks as workers become available
        with multiprocessing.Pool(n_processes) as pool:
            # Create argument tuples for each BAM file
            args_list = [(f, args.ref, args.fasta, args.outdir, user_targets, n_workers, min_mapq, use_process_discovery, min_support, traced_pairs) for f in bam_files]
            # Use starmap_async to start all tasks immediately
            # The pool will automatically start new tasks as workers become available
            async_result = pool.starmap_async(process_sample, args_list)
            # Get all results (pool keeps cores busy by starting new tasks as workers finish)
            results = async_result.get()
            for res in results:
                print(res, flush=True)

    print("\n[***] All samples processed. Check the output directory. [***]")

if __name__ == "__main__":
    main()