# 🧬 Project 1: WGS/WES Alignment Pipeline

> GATK Best Practices alignment: BWA-MEM2 → Picard MarkDuplicates → mosdepth → MultiQC

Available in two flavours:
- **`wgs_align.sh`** — single-sample Bash pipeline, ideal for local/WSL runs
- **`nextflow/`** — DSL2 pipeline with local, SLURM, and AWS Batch profiles

## Bash usage

```bash
bash wgs_align.sh \
  -s SAMPLE001 \
  -1 SAMPLE001_R1.fastq.gz \
  -2 SAMPLE001_R2.fastq.gz \
  -r reference/hg38.fa \
  -o results/ \
  -t 16

# WES mode (with target BED)
bash wgs_align.sh \
  -s SAMPLE001 -1 R1.fastq.gz -2 R2.fastq.gz \
  -r reference/hg38.fa -o results/ \
  -b targets.bed
```

### Pipeline steps

1. **Reference validation** — auto-builds BWA-MEM2 index if missing
2. **Pre-alignment QC** — fastp quality report
3. **Alignment** — BWA-MEM2 with read groups, sorted output
4. **Duplicate marking** — Picard MarkDuplicates (optical pixel distance 2500)
5. **Alignment statistics** — samtools flagstat + stats
6. **Coverage** — mosdepth (genome-wide or WES target-restricted)
7. **MultiQC** — aggregate HTML report

### Output structure

```
results/
├── aligned/
│   ├── SAMPLE001.markdup.bam
│   └── SAMPLE001.markdup.bai
├── qc/
│   ├── SAMPLE001.flagstat.txt
│   ├── SAMPLE001.markdup_metrics.txt
│   ├── SAMPLE001.mosdepth.summary.txt
│   └── SAMPLE001_fastp.html
├── logs/
└── multiqc/
    └── SAMPLE001_multiqc_report.html
```

## Nextflow usage

```bash
cd nextflow/

# Local run (Docker)
nextflow run main.nf -profile local \
  --reads "data/*_{R1,R2}*.fastq.gz" \
  --reference reference/hg38.fa \
  --outdir results

# SLURM HPC
nextflow run main.nf -profile slurm \
  --reads "data/*_{R1,R2}*.fastq.gz" \
  --reference /shared/refs/hg38.fa \
  --outdir /scratch/results

# AWS Batch
nextflow run main.nf -profile aws \
  --reads "s3://bucket/fastq/*_{R1,R2}*.fastq.gz" \
  --reference s3://bucket/refs/hg38.fa \
  --outdir s3://bucket/results
```

### Available parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--reads` | FASTQ glob pattern (paired) | `data/*_{R1,R2}*.fastq.gz` |
| `--reference` | Reference FASTA | `reference/genome.fa` |
| `--outdir` | Output directory | `results` |
| `--threads` | CPUs per process | 8 |
| `--remove_dups` | Remove (vs mark) duplicates | false |
| `--bed` | WES target BED for mosdepth | none |

## Requirements

All tools come from `environment.yml` at the repo root, or use the
Nextflow container profiles (Docker/Singularity) which pull pinned images
automatically.

## Expected QC ranges (GATK / ENCODE)

| Metric | WGS | WES |
|--------|-----|-----|
| Mapping rate | >90% | >90% |
| Duplication rate | <20% | <40% |
| Mean coverage | 30-40X | 100-150X |
