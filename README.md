# FusionDetector

**FusionFlow Chimeric Read Detector** – A robust, BAM-type-aware pipeline for detecting gene fusions from next-generation sequencing (NGS) alignments.

## Overview

FusionDetector identifies gene fusions from BAM files by leveraging both discordant read pairs and split-read evidence. It automatically adapts to different BAM types (modern aligners with SA tags vs. legacy discordant-only BAMs) and provides high-precision fusion detection with advanced false-positive filtering.

## Key Features

- **BAM-Type Adaptive**: Automatically detects and adapts to modern alignments (with SA tags) or legacy BAMs (discordant pairs only)
- **Dual Evidence Support**: Uses both discordant read pairs and split-read (SA tag) evidence for comprehensive fusion detection
- **High-Precision Filtering**: Advanced false-positive reduction through:
  - Breakpoint clustering detection (filters mapping artifacts)
  - Hotspot gene filtering (per-sample artifact detection)
  - Adaptive support thresholds based on evidence quality
  - Intra-chromosomal distance filtering
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

### Optional Arguments

- `--cores`: Number of CPU cores to use for parallel processing (default: auto-detect)
- `--min-mapq`: Minimum mapping quality threshold (default: 20)
- `--target`: Path to text file with target genes (one per line) for focused analysis

### Examples

```bash
# Single BAM file
python main.py --input sample.bam --ref refFlat.txt.gz --fasta hg19.fa --outdir results/

# Multiple BAM files (batch processing)
python main.py --input bam_directory/ --ref refFlat.txt.gz --fasta hg19.fa --outdir results/ --cores 8

# With target gene filtering
python main.py --input sample.bam --ref refFlat.txt.gz --fasta hg19.fa --outdir results/ --target target_genes.txt
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

3. **Advanced Artifact Filtering**: Implements sophisticated false-positive reduction including:
   - Per-sample hotspot detection (not global blacklists)
   - Breakpoint clustering to identify mapping artifacts
   - Adaptive thresholds based on evidence quality (SA tag fraction)

4. **Breakpoint Precision**: Uses reference FASTA for breakpoint refinement, improving accuracy especially for split-read evidence.

5. **Scalable Processing**: Efficient multiprocessing architecture handles both single samples and batch processing with optimal resource utilization.

## Dependencies

- Python 3.x
- `pysam` (for BAM/FASTA file handling)
- `intervaltree` (for gene annotation indexing)

See `requirements.txt` or `fusiondetect.yml` for complete dependency lists.

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

## License

[Add your license information here]

## Citation

[Add citation information if applicable]
