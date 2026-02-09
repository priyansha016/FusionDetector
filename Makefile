# FusionDetector — run pipeline for BAM files in Test/
# Usage: make all              — run for every BAM
#        make SRR948904        — run one sample
#        make all CORES=8      — use 8 workers
#        make index            — refresh BAM indexes

REF    := Test/refFlat.txt.gz
FASTA  := Test/hg19.fa
OUTDIR := Test/results
CORES  ?= 6
THREADS := $(CORES)

BAMS := Test/SRR948904.bam Test/SRR25536655.bam Test/SRR25536660.bam Test/HorizonDx_HD753.bam Test/H3122.bam

.PHONY: all index $(BAMS)

# Refresh BAM indexes so they are newer than the BAM (avoids htslib warnings when using make)
index:
	samtools index Test/SRR948904.bam
	samtools index Test/SRR25536655.bam
	samtools index Test/SRR25536660.bam
	samtools index Test/HorizonDx_HD753.bam
	samtools index Test/H3122.bam

all: $(BAMS)

# One target per BAM; pass CORES=N to limit workers (e.g. make all CORES=8)
Test/SRR948904.bam:
	python main.py --input Test/SRR948904.bam --ref $(REF) --fasta $(FASTA) --outdir $(OUTDIR) $(if $(THREADS),--cores $(THREADS),)

Test/SRR25536655.bam:
	python main.py --input Test/SRR25536655.bam --ref $(REF) --fasta $(FASTA) --outdir $(OUTDIR) $(if $(THREADS),--cores $(THREADS),)

Test/SRR25536660.bam:
	python main.py --input Test/SRR25536660.bam --ref $(REF) --fasta $(FASTA) --outdir $(OUTDIR) $(if $(THREADS),--cores $(THREADS),)

Test/HorizonDx_HD753.bam:
	python main.py --input Test/HorizonDx_HD753.bam --ref $(REF) --fasta $(FASTA) --outdir $(OUTDIR) $(if $(THREADS),--cores $(THREADS),)

Test/H3122.bam:
	python main.py --input Test/H3122.bam --ref $(REF) --fasta $(FASTA) --outdir $(OUTDIR) $(if $(THREADS),--cores $(THREADS),)

# Run H3122 and check ALK-EML4 breakpoints vs truth (chr2:29446607, chr2:42526889)
check-H3122: Test/H3122.bam
	python scripts/check_breakpoints.py $(OUTDIR)/H3122_fusions.tsv --truth scripts/truth_H3122.txt --tolerance 100

# Run SRR25536660 and check CD74-ROS1 breakpoints vs truth (chr6:117647499, chr5:149783394); use tolerance 2 for 1–2 bp accuracy
check-SRR25536660: Test/SRR25536660.bam
	python scripts/check_breakpoints.py $(OUTDIR)/SRR25536660_fusions.tsv --truth scripts/truth_SRR25536660.txt --tolerance 2

# Convenience targets by sample name (run e.g. make SRR948904)
SRR948904:       Test/SRR948904.bam
SRR25536655:     Test/SRR25536655.bam
SRR25536660:     Test/SRR25536660.bam
HorizonDx_HD753: Test/HorizonDx_HD753.bam
H3122: Test/H3122.bam
