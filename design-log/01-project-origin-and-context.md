# Design Log #01: Project Origin & Context

## Background

FusionDetector (FusionFlow Chimeric Read Detector) was generated to solve a real-world fusion-calling gap across different BAM types.

**Reference tools (all failing on different BAM types – used only for logic/knowledge, not as “shortlisted” runners):** Factera, GeneFuse, JuLI.

**Sample sets:**
- **F1cdx / VA BAMs** – older alignments (bwa-sampe).
- **Benchmarking data** – from literature.

**Observed behavior:**
- JuLI and GeneFuse: work on external samples; **no fusions** on VA BAM samples.
- Factera: fails on most external samples; **~50%** on VA BAMs.

**Goal:** Use the three tools’ logic as reference to build a **single, BAM-type-aware** pipeline that (1) **detects what kind of BAM** it is and (2) **finds fusions** from BAM with an approach suited to that type.

---

## FusionDetector’s Approach (Not a Port of Any One Tool)

- **BAM-type detection first** (sniffer): e.g. high SA-tag / split-read density vs legacy (discordant-only) BAMs.
- **For “new type” BAMs** (modern alignments with good chimeric/SA support): use **both** discordant pairs and split reads.
- **For legacy/VA-style BAMs:** pipeline can emphasize discordant-pair logic where split/SA evidence is weak or absent.
- Logic of Factera, GeneFuse, and JuLI was studied for **how** fusion tools can work (discordant pairs, soft-clips, SA tags, breakpoint refinement, validation); FusionDetector implements a **different, type-based strategy** rather than reimplementing any one of them.

---

## Input Data Types (Design Constraint)

| Mode | Description | Precision |
|------|-------------|-----------|
| **Split-reads (SA tags)** | Chimeric reads with junction in the read; SA tag gives second location | High |
| **Discordant pairs (PE / Legacy)** | Mates map to different genes/chromosomes; junction not sequenced (Factera-style) | Lower |

The pipeline must support both and choose behavior based on BAM characteristics.

---

## Current Architecture

```
main.py          → CLI, parallel batch (ProcessPoolExecutor)
engine/discoverer→ Gene-pair candidates (SA + discordant) using GenomeIndex
engine/assembler → Consensus breakpoint (bpa, bpb), read support
engine/validator → Technical filters (neighbors, paralogs, low complexity)
utils/genome_index → Gene intervals (e.g. refFlat)
utils/reporter   → Output (e.g. CSV)
```

**Discovery:** We always collect **both** split-read (SA) and discordant pair candidates for every BAM. Mode (LEGACY/MODERN) was removed; pipeline runs both SA and discordant regardless of BAM type. If this captures too much noise in practice, we can reintroduce the sniffer and use mode to tune discovery (e.g. restrict discordant in MODERN or emphasize one signal). Sniffer module left in repo for that possibility.

**Target filter:** `--target genes.txt` is applied **on top** of discovery: only fusions involving at least one listed gene are reported (in main.py). It does not change which reads are collected; it only filters which candidates are written to output.

---

## Key Fixes Already Applied

1. **“Same chromosome” bug (critical)** – **Fixed in assembler**  
   CD74–ROS1 (Chr5/Chr6) was mis-assigned as intra-chromosomal (chr6, 1 bp apart) and then removed by the neighbor filter.  
   **Cause:** Consensus in `FusionAssembler.find_breakpoint` took `most_common(1)` of evidence_a and evidence_b independently, so when more reads mapped to one side (e.g. ROS1 on chr6), both bp_a and bp_b could be (chr6, …).  
   **Fix (engine/assembler.py):** After setting bp_a = most_common(evidence_a), bp_b is chosen from evidence_b **only on chromosomes ≠ bp_a[0]**. So inter-chromosomal fusions get correct chr/pos for both sides; same-chromosome fusions unchanged (fallback to most_common(evidence_b)).

2. **Filtering balance**  
   - `is_complex`: disabled (was killing fusions in intronic/repeat regions).  
   - `is_likely_fp`: active (intra-chromosomal neighbor filter <500KB, paralog-prefix filter). With correct chr assignment, the neighbor filter no longer removes CD74–ROS1.

---

## Open Technical Challenges

- **False positives:** With `is_complex` off, noise increases; balance “VIP gene list” vs general filters.
- **Optional:** Local mini-assembly for 1 bp-off discordant reads to refine junction.

---

## Recommended Next Steps

1. ~~Assembler: fix inter-chromosomal reporting~~ **Done.**
2. Keep the intra-chromosomal distance filter in validator; it now correctly retains inter-chromosomal fusions.
3. **Local assembly (optional):** Consider mini-assembly for discordant-only support to get exact junction.

---

## Reference Tool Logic (For Knowledge Only)

Summaries of how each reference tool works – **not** what FusionDetector runs; used to inform our BAM-type–based design.

### Factera (Perl)

- **Input:** BAM (indexed), exons.bed, genome.2bit; optional targets.bed.
- **Method:** (1) Find improperly paired reads (samtools `-F 2`, FLAG filter), assign to gene pairs via exon coords. (2) Find soft-clipped reads in fusion-target regions (`-L` fusetargets), min clip length (e.g. 16 bp). (3) Breakpoint from soft-clips: k-mer comparison of clipped vs non-clipped segment, consensus, orientation (NC/CN). (4) Build fusion template (twoBitToFa), BLAST reads (improper, soft-clipped, unmapped) against it; require min breakpoint-spanning reads and similarity.
- **Dependencies:** SAMtools, twoBitToFa, blastn, makeblastdb, Statistics::Descriptive.
- **Note:** Heavy reliance on discordant pairs to define candidates, then soft-clips + BLAST to confirm; works better on some VA BAMs, fails on many external BAMs.

### GeneFuse

- **Input:** FASTQ (not BAM); refFlat for gene/transcript structure.
- **Shown snippet:** Utility to generate a “fusion gene” file from a gene list (gene ± transcript) using refFlat (longest transcript if not specified). Main fusion caller logic not in the snippet; tool is FASTQ-based.
- **Note:** Works on external samples, finds no fusions in VA BAMs; different input and alignment assumptions than BAM-first pipelines.

### JuLI (R)

- **Input:** Case BAM(s), optional control BAM, refgene, Gap file, reference FASTA; TargetBed optional; AnalysisType e.g. `'DP'` (Discordant + Proper/split).
- **Method:** (1) BAM stats (read length, insert size, split-read count). (2) Candidate breaks from CIGAR (S/H): start/end break positions, overlap with Gap filtered out. (3) Optional control/panel filter. (4) **Discordant (D):** improper pairs + SA tag; count supporting reads per break and mate region; require discordant + split support; build consensus contigs; compare contigs between breaks (pairwise alignment); nucleotide diversity filter. (5) **Proper/split (P):** properly paired reads with soft-clips; count support; consensus; compare contigs within chromosome. (6) Output per sample with chr, pos, ori, spl/dis counts, counter (partner break).
- **Note:** Uses both split-reads and discordant pairs; SA tag and CIGAR; works on external, no fusions in VA BAMs – likely tuned to different BAM/aligner behavior.

---

## Reference

- Initial prompt and project summary derived from Gemini conversation (project generation).
- Focus: **finding fusions from BAM** with BAM-type-aware logic; for new-type BAMs, use **both** discordant and split reads.
