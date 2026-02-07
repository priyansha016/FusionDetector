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

BAMS := Test/SRR948904.bam Test/SRR25536655.bam Test/SRR25536660.bam Test/HorizonDx_HD753.bam

.PHONY: all index $(BAMS)

# Refresh BAM indexes so they are newer than the BAM (avoids htslib warnings when using make)
index:
	samtools index Test/SRR948904.bam
	samtools index Test/SRR25536655.bam
	samtools index Test/SRR25536660.bam
	samtools index Test/HorizonDx_HD753.bam

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

# Convenience targets by sample name (run e.g. make SRR948904)
SRR948904:       Test/SRR948904.bam
SRR25536655:     Test/SRR25536655.bam
SRR25536660:     Test/SRR25536660.bam
HorizonDx_HD753: Test/HorizonDx_HD753.bam
