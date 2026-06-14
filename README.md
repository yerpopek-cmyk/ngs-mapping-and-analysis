# NGS Mapping and Analysis Suite

The NGS Mapping and Analysis Suite is a production-grade bioinformatics workspace that implements standard genomic alignment workflows, quality control batteries, and benchmarking harnesses. Designed for whole-genome sequencing (WGS), whole-exome sequencing (WES), and RNA-seq datasets, it provides a structured repository for processing sequencing reads, validating alignment quality, and comparing tool performance.

---

## What This Suite Is About

Aligning high-throughput sequencing data is a critical first step in genomic analysis. This suite wraps best-practice bioinformatics tools into automated pipelines and adds robust downstream utilities to visualize coverage, aggregate quality control metrics, and benchmark alignment algorithms. 

To facilitate testing and development without massive storage requirements, the workspace includes an automated data downloader and read simulator that retrieves human **Chromosome 22** and simulates paired-end reads to run the entire suite locally.

---

## Repository Components

### 1. WGS/WES Alignment Pipeline (`01_wgs_pipeline/`)
This pipeline follows GATK Best Practices to align DNA-seq reads to a reference genome. It features:
- **BWA-MEM2**: The state-of-the-art DNA aligner, utilizing AVX2/AVX-512 acceleration for rapid mapping.
- **Picard MarkDuplicates**: Identifies and flags PCR and optical duplicates to prevent false variant calling.
- **mosdepth**: Calculates genome-wide or target-region coverage depth with high speed.
- **MultiQC**: Aggregates flagstat, markdup, and coverage logs into a single HTML report.

### 2. RNA-seq Alignment Suite (`02_rnaseq_pipeline/`)
A splice-aware alignment workflow designed to map transcriptome sequencing reads across exon junctions. It includes:
- **HISAT2**: A fast, splice-aware aligner that maps reads using hierarchical indexing.
- **Qualimap**: Computes specific RNA-seq metrics, such as transcript genomic origin, 5'-3' bias, and coverage.
- **RSeQC QC Battery**: A comprehensive quality control suite that calculates:
  - **Read Distribution**: Identifies if reads map to CDS, 5' UTR, 3' UTR, introns, or intergenic regions.
  - **Strandedness**: Infers library preparation strandedness (e.g., dUTP / strand-specific reverse).
  - **Inner Distance**: Estimates the insert size between paired-end reads.
  - **Junction Saturation**: Checks if sequencing depth was sufficient to discover splice junctions.

### 3. Coverage Panorama (`03_coverage_panorama/`)
A custom Python visualization tool that reads mosdepth regions files and generates karyotype-style coverage art. It maps coverage depth onto color-coded horizontal bars representing chromosomes, allowing visual identification of copy number variations (CNVs) or target enrichment issues.

### 4. QC Sentinel Dashboard (`04_qc_sentinel/`)
A unified quality control dashboard that recursively scans a results directory, parses alignment/duplication/coverage files, and renders a color-coded pass/warn/fail checklist in the terminal (using the Python `Rich` library) based on GATK and ENCODE consortium guidelines. It also exports self-contained HTML reports to share with collaborators.

### 5. Aligner Benchmark Harness (`05_aligner_benchmark/`)
A benchmarking tool that subsamples input reads and runs them through **BWA-MEM2**, **minimap2** (short-read mode), and **HISAT2** under identical resources. It measures alignment accuracy, runtime, and peak memory usage (MaxRSS) using the GNU `time` utility, and outputs a weighted composite score and comparison plots.

---

## Installation & Setup

### 1. Create Conda Environment
Install all core aligners, samtools, and python libraries by importing the environment specification:

```bash
conda env create -f environment.yml
conda activate ngs-mapping
```

### 2. Prepare Reference & Test Dataset
Run the custom Python preparation script. This script downloads the human Chromosome 22 reference genome from UCSC, retrieves gene annotations (GTF) and database transcripts (BED), and simulates 50,000 paired-end sequencing reads (excluding centromeres/assembly gaps):

```bash
python download_and_prep.py
```

---

## Running the Orchestration Scripts

For convenience, two main orchestrator scripts are provided to run all tools sequentially:

### 1. Execute the Alignment Pipelines
To build reference indexes and run the DNA-seq WGS and RNA-seq alignment pipelines on the Chromosome 22 dataset:

```bash
bash run_pipelines.sh
```

### 2. Execute Downstream Reports and Benchmarks
Once the alignment run completes, generate the coverage panorama bar chart, the QC Sentinel cohort dashboard, and the multi-aligner benchmarking comparison:

```bash
bash run_downstream.sh
```

All final reports, plots, and BAM files will be written to the `results/`, `rnaseq_results/`, and `benchmark_results/` directories.
