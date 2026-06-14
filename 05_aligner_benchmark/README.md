# 📊 Project 5: Aligner Benchmark

> Run BWA-MEM2, minimap2, and HISAT2 on the same data, get a composite score, pick a winner.

## What it does

`run_benchmark.sh` subsamples your FASTQ to a manageable size, then runs each
available aligner with `/usr/bin/time` wrapping to capture runtime and peak
memory (MaxRSS). `compare_aligners.py` turns the raw TSV into bar charts and
an HTML report with a weighted composite score.

## Quick start

```bash
pip install -r requirements.txt

# DNA-seq only (BWA-MEM2 vs minimap2)
bash run_benchmark.sh \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -r reference/hg38.fa \
  -o benchmark_results \
  -n 2000000 -t 16

# Include HISAT2 (RNA-seq capable comparison)
bash run_benchmark.sh \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -r reference/hg38.fa \
  -x reference/hg38_hisat2_index \
  -o benchmark_results \
  -n 2000000 -t 16
```

## Options

| Flag | Description | Default |
|------|-------------|---------|
| `-1 / -2` | Paired FASTQ files | required |
| `-r` | Reference FASTA (BWA-MEM2 + minimap2) | required |
| `-x` | HISAT2 index prefix (optional) | none |
| `-o` | Output directory | `benchmark_results` |
| `-t` | Threads | 8 |
| `-n` | Reads to subsample (for fast comparison) | 2,000,000 |

## Composite score formula

```
score = 0.40 × (alignment_rate / max_alignment_rate)
      + 0.25 × (1/runtime / max(1/runtime))
      + 0.20 × (properly_paired_pct / 100)
      + 0.15 × (1/memory / max(1/memory))
```

Higher is better. The aligner with the top score is flagged with ⭐ in the report.

## Output

```
benchmark_results/
├── bam/                     # sorted, indexed BAM per aligner
├── logs/                    # per-aligner stderr + timing
├── metrics/
│   └── benchmark_raw.tsv    # raw runtime/memory/alignment metrics
├── subsample/               # subsampled FASTQs (reused across aligners)
├── benchmark_plots.png      # 4-panel comparison chart
└── benchmark_report.html    # full HTML report with recommendation
```

## Interpreting results

- **BWA-MEM2**: gold standard for DNA-seq variant calling, GATK-compatible
- **minimap2**: very fast, excellent for long reads, slightly different
  scoring for short-read SNV calling
- **HISAT2**: splice-aware — only meaningful comparison for RNA-seq data;
  will show lower "properly paired" on genomic DNA due to splice handling

The "best" aligner depends on your downstream application — the composite
score is a general-purpose heuristic, not a substitute for benchmarking
against your specific variant-calling or quantification pipeline.

## Standalone comparison script

You can re-run just the comparison/report step on an existing TSV:

```bash
python compare_aligners.py \
  --metrics benchmark_results/metrics/benchmark_raw.tsv \
  --output benchmark_results/benchmark_report.html \
  --plot benchmark_results/benchmark_plots.png \
  --reads 2000000
```
