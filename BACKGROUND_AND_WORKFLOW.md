# FusionDetector: Background, How It Works, and Workflow

This document describes the background of the project, how the pipeline works, how it differs from other fusion-detection tools, and the role of each file in the workflow from input to output.

---

## 1. Background

**FusionDetector** (FusionFlow Chimeric Read Detector) was built to detect gene fusions from BAM files in a way that works across different alignment types. Existing tools (e.g. Factera, GeneFuse, JuLI) often work well on one kind of BAM but fail or underperform on others: for example, some find fusions on modern alignments but miss them on legacy BAMs, and vice versa.

The goal was a **single, BAM-type-aware pipeline** that:

- Uses **both** discordant read pairs and split-read (SA tag) evidence
- Adapts to the BAM at hand (modern vs legacy) without requiring different tools
- Reduces false positives with filters tuned for real fusion vs artifact patterns
- Produces precise breakpoints when a reference FASTA is provided

The design was informed by how other tools work (discordant pairs, soft-clips, SA tags, breakpoint refinement), but FusionDetector is **not** a port of any one tool; it implements its own strategy and filters.

---

## 2. How It Works (High-Level)

The pipeline has four main stages:

1. **Discovery**  
   Scan the BAM and assign reads to **gene pairs**. A read contributes to a fusion candidate if:
   - It has an **SA tag** (split-read): primary alignment in one gene, secondary (SA) in another, or  
   - It is **discordant**: mates map to different chromosomes/genes (improper pair, mate mapped).

2. **Assembly**  
   For each gene-pair candidate with enough supporting reads (≥20 by default), compute **consensus breakpoints**:
   - Use primary and SA (or mate) positions to get two breakpoint regions (chromosome + position).
   - For same-chromosome fusions, separate two clusters (e.g. EML4 vs ALK on chr2) by a minimum gap.
   - If a reference FASTA is given and there are SA tags, **refine** breakpoints by aligning soft-clipped sequences to the reference (Factera-style refinement).

3. **Validation**  
   Apply **false-positive filters** (e.g. same-gene/paralog, intra-chromosomal neighbors &lt;200 kb, blacklisted gene families, adaptive support and SA-fraction thresholds). Only candidates that pass are kept.

4. **Post-processing and reporting**  
   Deduplicate gene aliases, filter **breakpoint clusters** (many fusions sharing one breakpoint region = likely artifact), filter **hotspots** (genes appearing in many fusions in one sample), then write the final list to a TSV file.

So: **BAM + RefFlat + optional FASTA** → **Discovery** → **Assembly** → **Validation** → **Post-processing** → **TSV output**.

---

## 3. How FusionDetector Differs from Other Tools

| Aspect | FusionDetector | Typical alternatives |
|--------|----------------|----------------------|
| **BAM type** | Single pipeline for both modern (SA-rich) and legacy (discordant-heavy) BAMs. | Many tools assume one BAM/aligner type and fail on others. |
| **Evidence** | Combines discordant pairs and split-reads (SA) in one flow; no need to run separate tools. | Some tools are mainly discordant-based or mainly split-read-based. |
| **False positives** | Per-sample hotspot detection, breakpoint-cluster filtering, adaptive support/SA thresholds, no global “hotspot list” that might drop true fusions. | Often rely on fixed thresholds or global blacklists. |
| **Breakpoints** | When FASTA is provided and SA tags exist, refines breakpoints using reference and soft-clips. | Not all tools do reference-based refinement. |
| **Input** | BAM + RefFlat + FASTA (optional but recommended). Optional target gene list. | Some tools need FASTQ; some need BAM + different annotation formats. |
| **Scalability** | Multiprocessing for many BAMs and for per-sample gene-pair parallelism. | Varies by tool. |

---

## 4. Role of Each File (From Input to Output)

### Entry and orchestration

| File | Role |
|------|------|
| **`main.py`** | CLI entry point. Parses arguments (`--input`, `--ref`, `--fasta`, `--outdir`, `--cores`, `--target`, etc.), resolves BAM list (single file or directory), loads target genes if provided. For each sample, calls `process_sample()`. Handles parallelism: single BAM → one process with internal pool for gene pairs; multiple BAMs → one pool of processes (one per BAM). Builds per-gene-pair work items (with serialized reads), runs assembly/validation in workers, then runs post-processing (dedup, breakpoint-cluster filter, duplicate-breakpoint filter, hotspot filter) and reporting. |

### Engine (core logic)

| File | Role |
|------|------|
| **`engine/discoverer.py`** | **Discovery.** Loads BAM and, using `GenomeIndex`, assigns reads to gene pairs. For each read: if it has an SA tag, gets genes at primary and SA positions; if it’s discordant (improper pair, mate mapped), gets genes at read and mate positions. Output: dictionary `gene_pair → list of reads`. Supports multithreaded or multiprocess BAM scan. Defines `ReadLike` and `_read_to_dict` for serializing reads for multiprocessing. |
| **`engine/assembler.py`** | **Breakpoint assembly.** For a list of reads for one gene pair, computes consensus breakpoints (chromosome + position for both sides). Handles inter- vs intra-chromosomal cases; for same chromosome, enforces a minimum gap between the two breakpoint clusters. If `fasta_path` is set and there are SA tags, refines breakpoints via `_refine_with_clips()` and `_bp_correction_factera()` (align soft-clips to reference). Returns `(bp_a, bp_b), support`. |
| **`engine/validator.py`** | **Validation / FP filtering.** Implements `is_likely_fp()`: same gene or 4-letter prefix, intra-chromosomal &lt;200 kb, blacklisted gene prefixes, adaptive support and SA-fraction thresholds. Used after assembly to drop low-confidence or artifact-like candidates. Also provides `get_sequence_context()` and `is_complex()` (e.g. for potential future sequence-based filters). |
| **`engine/sniffer.py`** | **BAM-type detection (optional).** Can sample the BAM to classify it as e.g. “LEGACY” (few SA tags) vs “MODERN” (many SA tags). Currently not used to change pipeline behavior; discovery always uses both SA and discordant evidence. Kept for possible future tuning. |

### Utilities

| File | Role |
|------|------|
| **`utils/genome_index.py`** | **Gene annotation index.** Loads RefFlat (or similar) into interval trees per chromosome. Provides `get_gene_at(chrom, pos)` to map a position to gene name(s). Used by the discoverer to assign reads to genes. |
| **`utils/reporter.py`** | **Output.** Writes the final fusion table to a TSV file (sample, gene_a, gene_b, chrom_a, pos_a, chrom_b, pos_b, support). Converts 0-based breakpoints to 1-based for the report. |
| **`utils/hts_filter.py`** | **Stderr filtering.** Suppresses repeated htslib BAM-index warnings (e.g. “index older than data file”) so they appear at most once. Used so logs stay readable when many BAMs or workers open files. |

### Configuration and scripts

| File | Role |
|------|------|
| **`requirements.txt`** | Python dependencies (e.g. pysam, intervaltree). |
| **`fusiondetect.yml`** | Conda environment definition. |
| **`Makefile`** | Convenience targets to run the pipeline on test BAMs with fixed ref/FASTA/outdir. |
| **`scripts/check_breakpoints.py`** | Compares pipeline output to a truth file (expected breakpoints). |
| **`scripts/check_bam_fusion.py`** | Utility to check BAMs for fusion-related evidence. |
| **`design-log/01-project-origin-and-context.md`** | Design log: origin, constraints, architecture, and reference-tool notes. |

---

## 5. Workflow From Input to Output

```
INPUTS
------
  • BAM file(s)          (--input: file or directory)
  • RefFlat              (--ref: gene annotations)
  • Reference FASTA      (--fasta: for breakpoint refinement)
  • Optional: target genes (--target: one gene per line)
  • Optional: --cores, --min-mapq

STEP 1 – SETUP (main.py)
-------------------------
  • Resolve list of BAM paths.
  • Load target genes (if --target).
  • Create output directory.

STEP 2 – PER-SAMPLE (process_sample in main.py)
-----------------------------------------------
  For each BAM:

  2a. Load genome index (utils/genome_index.py)
      • Load RefFlat into interval trees.
      • Used to map (chrom, pos) → gene name(s).

  2b. Discovery (engine/discoverer.py)
      • Open BAM; scan by reference (optionally in parallel threads/processes).
      • For each read: if SA tag → genes at primary + SA; if discordant → genes at read + mate.
      • Build candidates: gene_pair → [list of reads].
      • Output: raw_candidates (dict gene_pair → reads).

  2c. Prepare work items (main.py)
      • For each gene_pair with enough reads, serialize reads to picklable dicts.
      • Build list of (bam_path, fasta_path, sample_name, gene_pair, reads_dicts, user_targets, MIN_SUPPORT).

  2d. Assembly + validation (main.py workers → engine/assembler.py, engine/validator.py)
      • Each worker: open BAM/FASTA, run FusionAssembler.find_breakpoint(reads) → breakpoints, support.
      • If support &lt; MIN_SUPPORT or target filter excludes pair → skip.
      • If target fusion → keep. Else run FusionValidator.is_likely_fp(...) → drop if True.
      • Return fusion dict (sample, gene_a, gene_b, chrom_a, pos_a, chrom_b, pos_b, support, sa_fraction).

  2e. Post-processing (main.py)
      • Deduplicate gene aliases (e.g. MYCN/MYCNOS same breakpoint).
      • Filter breakpoint clusters (≥3 fusions in same small region → keep best).
      • Filter duplicate breakpoints (same chr:pos in multiple fusions → keep best).
      • Filter hotspots (per-sample: genes in &gt;3 fusions, drop low-support hotspot–hotspot pairs with exceptions).

  2f. Report (utils/reporter.py)
      • Write final list to <outdir>/<sample_name>_fusions.tsv (1-based positions).

STEP 3 – DONE
-------------
  • One TSV per sample in --outdir.
  • Columns: sample, gene_a, gene_b, chrom_a, pos_a, chrom_b, pos_b, support.
```

---

## 6. Data Flow Summary

- **RefFlat** → **GenomeIndex** → used in **Discoverer** to assign reads to gene pairs.
- **BAM** → **Discoverer** → gene_pair → reads; then **Assembler** → breakpoints + support; **Validator** → keep/drop; **main.py** post-processing → final list; **Reporter** → TSV.
- **FASTA** → **Assembler** (and optionally **Validator**) for breakpoint refinement and sequence context when SA tags exist.
- **Target genes** → applied in **main.py** (only report fusions involving at least one target gene when --target is used).

This is the full path from your inputs to the reported fusions.
