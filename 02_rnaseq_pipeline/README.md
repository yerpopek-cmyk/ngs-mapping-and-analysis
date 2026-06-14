# 🧪 Project 2: RNA-seq Alignment Suite

> Splice-aware alignment (HISAT2) + complete RSeQC quality control battery.

## Scripts

| Script | Purpose |
|--------|---------|
| `rnaseq_hisat2.sh` | Full alignment pipeline: fastp → HISAT2 → Picard → Qualimap → MultiQC |
| `rseqc_suite.sh` | Standalone comprehensive RSeQC QC battery (8 modules) |

## Quick start

### 1. Build HISAT2 index (one-time, per genome)

```bash
hisat2-build -p 8 reference/genome.fa reference/genome_index
```

### 2. Run the alignment pipeline

```bash
bash rnaseq_hisat2.sh \
  -s SAMPLE001 \
  -1 SAMPLE001_R1.fastq.gz \
  -2 SAMPLE001_R2.fastq.gz \
  -i reference/genome_index \
  -a annotation.gtf \
  -o results/ \
  -S RF \
  -t 16
```

`-S` strandedness options:
- `RF` — dUTP / Illumina stranded (most common, default)
- `FR` — standard stranded forward
- `unstranded` — non-strand-specific libraries

### 3. Run the full RSeQC battery

```bash
# First, convert GTF → BED12 (required by RSeQC)
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons annotation.gtf annotation.genePred
genePredToBed annotation.genePred annotation.bed

# Then run the suite on your aligned BAM
bash rseqc_suite.sh \
  -b results/aligned/SAMPLE001.dedup.bam \
  -r annotation.bed \
  -o results/rseqc/ \
  -s SAMPLE001
```

## RSeQC modules included

| Module | What it tells you |
|--------|-------------------|
| `read_distribution.py` | % reads in CDS/UTR/intron/intergenic |
| `geneBody_coverage.py` | 5'→3' coverage uniformity (degradation check) |
| `junction_saturation.py` | Known vs novel splice junction discovery |
| `inner_distance.py` | Fragment size between paired reads |
| `infer_experiment.py` | Confirms/detects library strandedness |
| `clipping_profile.py` | Adapter/quality clipping rates |
| `bam_stat.py` | General BAM-level statistics |
| `tin.py` | Transcript Integrity Number (RNA degradation) |

## Pipeline output structure

```
results/
├── aligned/
│   ├── SAMPLE001.markdup.bam   (kept for quantification)
│   └── SAMPLE001.dedup.bam     (for RSeQC, dup-removed)
├── qc/
│   ├── SAMPLE001_hisat2_summary.txt
│   ├── SAMPLE001.markdup_metrics.txt
│   ├── SAMPLE001_fastp.html
│   └── qualimap_SAMPLE001/
├── logs/
└── multiqc/
```

## Expected QC ranges (RNA-seq)

| Metric | Good range |
|--------|-----------|
| Overall alignment rate | >70% |
| Duplication rate | <30% |
| rRNA contamination | <5% |
| Intronic mapping (read_distribution) | <30% |
| Gene body coverage | flat profile (5'/3' bias <2x) |
