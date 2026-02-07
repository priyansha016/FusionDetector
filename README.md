# FusionDetector

**FusionFlow Chimeric Read Detector** – BAM-type-aware pipeline to find gene fusions from NGS alignments.

## Why this project

Three existing tools (Factera, GeneFuse, JuLI) were tried on two sample sets: VA/F1cdx BAMs (bwa-sampe) and literature benchmarking data. **All fail on different BAM types:** JuLI and GeneFuse work on external samples but find no fusions in VA BAMs; Factera fails on most external samples but works ~50% on VA BAMs. Their logic was studied as **reference only**; this project implements a **different, BAM-type–aware** pipeline that:

1. **Detects BAM type** (e.g. high SA-tag/split-read density vs legacy discordant-only).
2. **Finds fusions** from BAM using a strategy that fits that type.

For **new-type BAMs** (modern alignments with good chimeric/SA support), FusionDetector uses **both** discordant pairs and split reads; for legacy/VA-style BAMs it can emphasize discordant logic. See `design-log/01-project-origin-and-context.md` for reference-tool summaries and design details.

## Modes

- **Split-reads (SA tags):** Chimeric reads with the junction in the read; high precision.
- **Discordant pairs (Legacy):** Mates map to different genes/chromosomes; junction not sequenced; lower precision.
- **New-type BAMs:** Both discordant and split-read evidence are used.

## Layout

- `main.py` – CLI, parallel batch (`ProcessPoolExecutor`).
- `engine/` – Sniffer (BAM type) → Discoverer (candidates) → Assembler (breakpoints) → Validator (filters).
- `utils/` – Genome index (e.g. refFlat), reporter (e.g. CSV).
- `design-log/01-project-origin-and-context.md` – Initial prompt, architecture, milestones, and next steps.

## Run

```bash
# Virtual env
source .venv/bin/activate   # or: conda activate fusion

# Example
python main.py --bam path/to/sample.bam --refFlat refFlat.txt.gz [--target genes.txt]
```

## Dependencies

- Python 3.x, `pysam`, `intervaltree` (see `requirements.txt` or `fusiondetect.yml`).
