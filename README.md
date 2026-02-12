# FusionDetector

**FusionFlow Chimeric Read Detector** – A robust, BAM-type-aware pipeline for detecting gene fusions from next-generation sequencing (NGS) alignments.

## Overview

FusionDetector identifies gene fusions from BAM files by leveraging both discordant read pairs and split-read evidence. It automatically adapts to different BAM types (modern aligners with SA tags vs. legacy discordant-only BAMs) and provides high-precision fusion detection with advanced false-positive filtering.

## Key Features

- **BAM-Type Adaptive**: Automatically detects and adapts to modern alignments (with SA tags) or legacy BAMs (discordant pairs only)
- **Dual Evidence Support**: Uses both discordant read pairs and split-read (SA tag) evidence for comprehensive fusion detection
- **High-Precision Filtering**: Advanced false-positive reduction through global filters (same thresholds for all samples):
  - **Sink breakpoint filter**: Removes fusions where a breakpoint region appears in many different fusions (mapping artifacts)
  - **Repeating genes filter**: Filters genes appearing ≥3 times across fusions (both columns checked)
  - **Breakpoint clustering**: Groups nearby breakpoints and keeps only the best per cluster
  - **Duplicate breakpoint removal**: Removes exact duplicate breakpoints
  - **Low-confidence filter**: Filters very low support (20-21 reads) with minimal SA evidence
  - **Adaptive support thresholds**: Stricter requirements when combined with low SA tag fraction
  - **Intra-chromosomal distance filtering**: Filters same-chromosome fusions <500KB apart
  - **Gene family filtering**: Filters same-family genes (ZNF, GYP) and blacklisted prefixes (RPS, HLA, etc.)
  - **Exon boundary check**: Filters breakpoints far from exon boundaries (>50bp) with low support/SA evidence (inspired by Factera/GeneFuse)
- **Breakpoint Refinement**: Uses reference FASTA to refine breakpoint positions by aligning soft-clipped sequences
- **Parallel Processing**: Efficient multiprocessing support for both single and batch BAM processing
- **Target Gene Filtering**: Optional focus on specific genes of interest
- **Comprehensive Output**: TSV format with fusion details including breakpoint positions and read support

## Installation

```bash
# Clone the repository
git clone https://github.com/priyansha016/FusionDetector.git
cd FusionDetector

# Create virtual environment (optional but recommended)
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
# Or use conda:
# conda env create -f fusiondetect.yml
# conda activate fusiondetect
```

## Usage

### Basic Usage

```bash
python main.py --input <BAM_FILE> --ref <REFFLAT> --fasta <FASTA> --outdir <OUTPUT_DIR>
```

### Required Arguments

- `--input`: Path to a single BAM file or directory containing BAM files
- `--ref`: Path to reference RefFlat file (gene annotations)
- `--fasta`: Path to reference FASTA file (for breakpoint refinement)
- `--outdir`: Output directory for results

### Required Index Files

FusionDetector requires index files for efficient BAM and FASTA access:
- **BAM index (`.bai`)**: Must be present for each BAM file (e.g., `sample.bam.bai` or `sample.bai`)
- **FASTA index (`.fai`)**: Must be present for the reference FASTA file (e.g., `hg19.fa.fai`)

**Creating index files:**
```bash
# Create BAM index
samtools index your_file.bam

# Create FASTA index
samtools faidx your_reference.fa
```

The pipeline will exit with an error if index files are missing.

### Optional Arguments

- `--cores`: Number of CPU cores to use for parallel processing (default: auto-detect)
- `--min-mapq`: Minimum mapping quality threshold (default: 20)
- `--min-support`: Minimum read support for a fusion to be reported (default: 20)
- `--target`: Path to text file with target genes (one per line) for focused analysis

### Examples

```bash
# Single BAM file
python main.py --input sample.bam --ref refFlat.txt.gz --fasta hg19.fa --outdir results/

# Multiple BAM files (batch processing)
python main.py --input bam_directory/ --ref refFlat.txt.gz --fasta hg19.fa --outdir results/ --cores 8

# With target gene filtering
python main.py --input sample.bam --ref refFlat.txt.gz --fasta hg19.fa --outdir results/ --target target_genes.txt

# Custom minimum read support
python main.py --input sample.bam --ref refFlat.txt.gz --fasta hg19.fa --outdir results/ --min-support 25
```

## Output Format

Results are written to `<outdir>/<sample_name>_fusions.tsv` with the following columns:

- `sample`: Sample name
- `gene_a`, `gene_b`: Fusion partner genes
- `chrom_a`, `chrom_b`: Chromosomes
- `pos_a`, `pos_b`: Breakpoint positions (1-based)
- `support`: Number of supporting reads

## How FusionDetector Differs from Other Tools

1. **BAM-Type Awareness**: Unlike tools that assume a specific BAM format, FusionDetector automatically adapts to both modern aligners (BWA-MEM, STAR) and legacy alignments (bwa-sampe), ensuring broad compatibility.

2. **Unified Evidence Integration**: Combines discordant pairs and split-read evidence in a single pipeline, rather than requiring separate tools for different evidence types.

3. **Advanced Artifact Filtering**: Implements comprehensive false-positive reduction with **global filters only** (same thresholds for all samples):
   - **Sink breakpoint detection**: Identifies breakpoint regions that appear in many fusions (e.g. chr2:33141xxx "fusing" with multiple genes) and filters all fusions involving such sinks
   - **Repeating genes**: Filters genes appearing ≥3 times across fusions (checks both gene_a and gene_b columns)
   - **Breakpoint clustering**: Groups breakpoints within 10KB and keeps only the highest-support fusion per cluster
   - **Duplicate breakpoint removal**: Removes exact duplicate breakpoints (same chr:pos in multiple fusions)
   - **Low-confidence filtering**: Filters fusions with support 20-21 and low SA tag fraction (<0.2)
   - **Adaptive support/SA thresholds**: Requires higher support or SA fraction when other suspicious signals are present
   - **Intra-chromosomal filtering**: Filters same-chromosome fusions <500KB apart (read-through artifacts)
   - **Gene family/blacklist**: Filters same-family genes (ZNF, GYP) and blacklisted prefixes (RPS, HLA, LOC, etc.)
   - **Exon boundary check**: Filters breakpoints far from exon boundaries (>50bp) when combined with low support/SA evidence (conservative filter inspired by Factera/GeneFuse)

4. **Breakpoint Precision**: Uses reference FASTA for breakpoint refinement, improving accuracy especially for split-read evidence.

5. **Scalable Processing**: Efficient multiprocessing architecture handles both single samples and batch processing with optimal resource utilization.

## Dependencies

- Python 3.8–3.12 (tested; 3.x should work)
- `pysam` (for BAM/FASTA file handling)
- `intervaltree` (for gene annotation indexing)

See `requirements.txt` or `fusiondetect.yml` for complete dependency lists.

## Notes

- **Index files required**: BAM (`.bai`) and FASTA (`.fai`) index files are mandatory. The pipeline will check for missing indices and exit with an error message if they are not found.
- **BAM index warning**: If you see `[W::hts_idx_load3] The index file is older than the data file`, the pipeline suppresses it to at most once. To avoid it entirely, refresh BAM indexes: `samtools index your.bam`
- **Cloning on another system**: Use Python 3.8+ and reinstall dependencies (`pip install -r requirements.txt`). Ensure index files are present for all BAM and FASTA files.
- **All filters are global**: Same thresholds and rules apply to every sample—no sample-specific tuning.

## Project Structure

```
FusionDetector/
├── main.py              # CLI entry point
├── engine/
│   ├── discoverer.py   # Fusion candidate discovery
│   ├── assembler.py     # Breakpoint assembly and refinement
│   └── validator.py     # False-positive filtering
├── utils/
│   ├── genome_index.py # Gene annotation indexing
│   └── reporter.py      # Output formatting
└── design-log/          # Design documentation
```
