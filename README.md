# NGS Mapping WGS, RNA-seq, Panorama Suite

A production-grade bioinformatics workspace that implements standard genomic alignment workflows, quality control batteries, coverage visualization, and multi-aligner benchmarking. Designed for whole-genome sequencing (WGS), whole-exome sequencing (WES), and RNA-seq datasets, it provides a structured, reproducible framework for processing high-throughput sequencing reads from raw FASTQ to final reports.

To facilitate local testing and development without terabyte-scale storage, the suite ships with an automated data preparation script that downloads human **Chromosome 22** from the UCSC Genome Browser and simulates 50,000 paired-end reads, allowing the entire pipeline to run end-to-end on a standard laptop or WSL environment.

---

## Why This Suite Exists

Aligning high-throughput sequencing data is the critical first step in nearly every genomic analysis — variant calling, gene expression quantification, structural variant detection, and copy number profiling all depend on accurate, well-characterized alignments. However, the standard bioinformatics workflow involves running multiple independent tools (aligners, duplicate markers, coverage calculators, QC aggregators), each with its own input formats, parameters, and failure modes. This suite wraps those tools into automated, reproducible pipelines and adds custom downstream utilities to visualize results, enforce quality standards, and compare tool performance objectively.

---

## Repository Structure

```
ngs-mapping-and-analysis/
├── 01_wgs_pipeline/          # WGS/WES DNA alignment pipeline
│   ├── wgs_align.sh          # Single-sample Bash pipeline
│   └── nextflow/             # DSL2 Nextflow pipeline (local, SLURM, AWS)
├── 02_rnaseq_pipeline/       # RNA-seq splice-aware alignment suite
│   ├── rnaseq_hisat2.sh      # HISAT2 alignment + Qualimap + MultiQC
│   └── rseqc_suite.sh        # Comprehensive RSeQC QC battery (8 modules)
├── 03_coverage_panorama/     # Genome-wide coverage visualization
│   ├── coverage_panorama.py  # Karyotype-style coverage renderer
│   └── karyotype_painter.py  # Synthetic data generator for demos
├── 04_qc_sentinel/           # Unified QC dashboard with GATK/ENCODE thresholds
│   ├── qc_sentinel.py        # Cohort QC scanner and reporter
│   └── parsers/              # Modular parsers for flagstat, Picard, mosdepth
├── 05_aligner_benchmark/     # Multi-aligner benchmarking harness
│   ├── run_benchmark.sh      # Subsample, align, measure runtime/memory
│   └── compare_aligners.py   # Composite scoring and HTML report generator
├── download_and_prep.py      # Chromosome 22 reference + read simulator
├── simulate_reads.py         # Standalone paired-end read simulator
├── generate_mock_results.py  # Mock QC files for testing without alignment
├── run_pipelines.sh          # Master orchestrator: WGS + RNA-seq pipelines
├── run_downstream.sh         # Master orchestrator: Projects 03, 04, 05
├── environment.yml           # Conda environment specification
└── LICENSE
```

---

## Projects

### Project 1 — WGS/WES Alignment Pipeline

**Directory**: `01_wgs_pipeline/`

This pipeline implements the GATK Best Practices alignment workflow for DNA sequencing data. Starting from raw paired-end FASTQ files and a reference genome, it produces a duplicate-marked, indexed BAM file accompanied by alignment statistics and a coverage profile.

**What the pipeline does, step by step:**

1. **Reference validation** — checks for (and auto-builds) the BWA-MEM2 index if it does not already exist.
2. **Pre-alignment QC** — runs `fastp` to generate a read quality report (adapter content, quality distribution, GC bias).
3. **Alignment** — maps reads to the reference using BWA-MEM2, the SIMD-accelerated successor to BWA-MEM. Read group tags (`@RG`) are injected during alignment so downstream tools can identify samples and libraries.
4. **Duplicate marking** — Picard MarkDuplicates identifies PCR and optical duplicates. PCR duplicates arise when the same DNA template molecule is amplified multiple times during library preparation; optical duplicates are spurious clusters on the flow cell surface. These duplicates inflate variant allele frequencies if not flagged.
5. **Alignment statistics** — `samtools flagstat` and `samtools stats` calculate mapping rate, properly paired percentage, singleton rate, and other core metrics.
6. **Coverage calculation** — `mosdepth` computes per-region and genome-wide mean coverage depth. In WES mode, coverage is restricted to a user-supplied target BED file.
7. **Report aggregation** — MultiQC combines all individual tool reports into a single interactive HTML dashboard.

The pipeline is available in two execution modes: a single-sample Bash script (`wgs_align.sh`) for local or WSL runs, and a Nextflow DSL2 pipeline with profiles for local Docker, SLURM HPC, and AWS Batch execution.

---

### Project 2 — RNA-seq Alignment Suite

**Directory**: `02_rnaseq_pipeline/`

RNA-seq alignment differs fundamentally from DNA-seq alignment because RNA reads span exon-exon junctions — the aligner must be aware of introns and splice sites to map reads correctly. This suite handles the full RNA-seq alignment workflow and includes a comprehensive transcriptomic QC battery.

**Alignment pipeline** (`rnaseq_hisat2.sh`):

1. **Pre-alignment QC** — `fastp` quality and adapter report.
2. **Splice-aware alignment** — HISAT2 uses a hierarchical indexing scheme (global FM index + local indexes for splice sites) to efficiently map reads across exon-intron boundaries. The user can specify library strandedness (`RF` for dUTP/Illumina stranded, `FR` for standard forward, or `unstranded`).
3. **Duplicate marking** — Picard MarkDuplicates, as in the WGS pipeline. RNA-seq typically has higher duplication rates due to transcript abundance variation.
4. **RNA-specific QC** — Qualimap RNA-seq module calculates transcript genomic origin (what fraction of reads map to exonic, intronic, or intergenic regions), 5'-to-3' coverage bias (which can reveal RNA degradation), and overall alignment statistics.
5. **Report aggregation** — MultiQC consolidates everything into a single HTML report.

**RSeQC QC battery** (`rseqc_suite.sh`):

A standalone post-alignment diagnostic suite that runs eight specialized RNA-seq quality modules:

| Module | What it measures |
|--------|------------------|
| `read_distribution.py` | Percentage of reads mapping to CDS, 5' UTR, 3' UTR, introns, and intergenic regions |
| `geneBody_coverage.py` | 5'-to-3' coverage uniformity across all transcripts (detects RNA degradation) |
| `junction_saturation.py` | Whether sequencing depth is sufficient to discover all splice junctions |
| `inner_distance.py` | Fragment insert size distribution between paired-end reads |
| `infer_experiment.py` | Automatically detects library strandedness protocol |
| `clipping_profile.py` | Position-dependent adapter and quality clipping rates |
| `bam_stat.py` | General BAM-level alignment statistics |
| `tin.py` | Transcript Integrity Number — a per-gene RNA quality metric |

---

### Project 3 — Coverage Panorama

**Directory**: `03_coverage_panorama/`

A custom Python visualization tool that transforms raw sequencing depth data into publication-ready genome-wide coverage art. It reads `mosdepth` regions BED files and renders each chromosome as a colour-coded horizontal strip, similar to a cytogenetic karyotype painted with coverage intensity.

**How it works:**

- The colour scale maps coverage depth onto a continuous gradient: white (zero coverage), blue (normal/baseline depth), and red (excess coverage or amplification).
- Each chromosome is drawn proportionally to its genomic length, allowing visual identification of large copy number variations (CNVs), target enrichment dropouts, or centromeric gaps.
- A multi-sample comparison mode places samples side by side, making it straightforward to spot differences between tumour and normal, or between sequencing batches.

The companion script `karyotype_painter.py` generates realistic synthetic BED datasets (uniform 30X WGS, CNV-bearing, tumour at 65% purity, and uneven RNA-seq profiles) so the tool can be demonstrated and tested without real sequencing data.

---

### Project 4 — QC Sentinel Dashboard

**Directory**: `04_qc_sentinel/`

A unified quality control monitoring application that scans a results directory, parses alignment and duplication and coverage log files, and renders a colour-coded pass/warn/fail checklist against established consortium guidelines.

**Key features:**

- **Auto-discovery**: Recursively scans for `*.flagstat.txt`, `*.markdup_metrics.txt`, and `*.mosdepth.summary.txt` files — no manual file listing required.
- **Threshold enforcement**: Applies GATK Best Practices and ENCODE consortium thresholds to every metric. For example, mapping rate below 75% triggers a `FAIL`, duplication rate above 35% triggers a `FAIL`, and mean coverage below 15X triggers a `FAIL`.
- **Cohort mode**: When multiple samples are present, displays a side-by-side comparison table so outlier samples are immediately visible.
- **Export formats**: Generates self-contained HTML reports (shareable via email or embedded in lab notebooks) and machine-readable JSON summaries for programmatic downstream analysis.
- **Terminal display**: Uses the Python `Rich` library to render styled, colour-coded tables directly in the terminal.

**Thresholds applied:**

| Metric | PASS | WARN | FAIL |
|--------|------|------|------|
| Mapping rate | >=90% | >=75% | <75% |
| Properly paired | >=85% | >=70% | <70% |
| Duplication rate | <=20% | <=35% | >35% |
| Mean coverage (WGS) | >=25X | >=15X | <15X |

---

### Project 5 — Aligner Benchmark Harness

**Directory**: `05_aligner_benchmark/`

A benchmarking tool that runs multiple alignment algorithms on the same dataset under identical conditions and produces an objective comparison of accuracy, speed, and resource consumption.

**How it works:**

1. **Subsampling**: Extracts a user-defined number of read pairs from the input FASTQs using `seqtk sample`, ensuring all aligners process identical input.
2. **Alignment**: Runs each available aligner (BWA-MEM2, minimap2 in short-read mode, and optionally HISAT2 if an index is provided) with the same thread count and reference genome.
3. **Profiling**: Wraps each aligner invocation with `/usr/bin/time -v` to capture wall-clock runtime (seconds) and peak resident set size (MaxRSS in MB).
4. **Scoring**: Computes a weighted composite score for each aligner:

```
Score = 0.40 x (alignment_rate / best_alignment_rate)
      + 0.25 x (speed / best_speed)
      + 0.20 x (properly_paired_pct / 100)
      + 0.15 x (memory_efficiency / best_memory_efficiency)
```

5. **Reporting**: Generates a four-panel bar chart comparing alignment rate, runtime, memory, and composite score, alongside a self-contained HTML report that highlights the top-scoring aligner.

**Interpreting results**: BWA-MEM2 is the gold standard for DNA-seq variant calling (GATK-compatible). Minimap2 is extremely fast and excels at long-read alignment. HISAT2 is splice-aware and only meaningful for RNA-seq data. The composite score is a general-purpose heuristic — the "best" aligner depends on the downstream application.

---

## Installation and Setup

### 1. Create Conda Environment

Install all core aligners, samtools, and Python libraries by importing the environment specification:

```bash
conda env create -f environment.yml
conda activate ngs-mapping
```

### 2. Prepare Reference and Test Dataset

Run the automated preparation script. It downloads the human Chromosome 22 reference genome from UCSC, retrieves gene annotations (GTF) and transcript database files (BED12), and simulates 50,000 paired-end sequencing reads (150 bp, Q40 quality, excluding centromeric assembly gaps):

```bash
python download_and_prep.py
```

This creates:
- `reference/chr22.fa` — the Chromosome 22 reference FASTA
- `reference/chr22.gtf` — RefGene gene annotations filtered for chr22
- `reference/chr22.bed` — transcript models in BED12 format (required by RSeQC)
- `data/chr22_sample_R1.fastq.gz` and `data/chr22_sample_R2.fastq.gz` — simulated paired-end reads

---

## Running the Orchestration Scripts

Two master scripts are provided to run all tools sequentially on the Chromosome 22 dataset:

### 1. Execute the Alignment Pipelines

Builds reference indexes (BWA-MEM2 and HISAT2) and runs the WGS and RNA-seq alignment pipelines, followed by the full RSeQC QC battery:

```bash
bash run_pipelines.sh
```

### 2. Execute Downstream Reports and Benchmarks

Once the alignment run completes, generates the coverage panorama visualization, the QC Sentinel cohort dashboards (for both WGS and RNA-seq), and the multi-aligner benchmarking comparison:

```bash
bash run_downstream.sh
```

### Output Locations

| Output | Directory | Description |
|--------|-----------|-------------|
| WGS BAM and QC | `results/` | Aligned BAM, flagstat, coverage, MultiQC report |
| RNA-seq BAM and QC | `rnaseq_results/` | Aligned BAM, RSeQC metrics, Qualimap, MultiQC report |
| Benchmark comparison | `benchmark_results/` | Per-aligner BAMs, timing logs, comparison plots, HTML report |

---

## Technology Stack

| Category | Tools |
|----------|-------|
| DNA alignment | BWA-MEM2 |
| RNA alignment | HISAT2 |
| Duplicate marking | Picard MarkDuplicates |
| Coverage analysis | mosdepth |
| Pre-alignment QC | fastp |
| RNA-seq QC | Qualimap, RSeQC (8 modules) |
| Report aggregation | MultiQC |
| Benchmarking | seqtk, GNU time, matplotlib |
| QC dashboard | Python Rich, custom parsers |
| Workflow engine | Nextflow DSL2 (optional) |
| Environment | Conda |

## License

This project is released under the MIT License. See `LICENSE` for details.
