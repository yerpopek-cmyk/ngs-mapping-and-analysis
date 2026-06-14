# NGS Mapping and Analysis Suite

A production-grade pipeline suite and bioinformatics utility package for Whole-Genome Sequencing (WGS) and RNA-seq alignment, quality control, and comparison.

## Repository Structure

- `01_wgs_pipeline/`: Whole-genome/exome alignment pipeline using BWA-MEM2, Picard MarkDuplicates, and mosdepth.
- `02_rnaseq_pipeline/`: Splice-aware RNA-seq alignment suite featuring HISAT2 and the RSeQC QC battery.
- `03_coverage_panorama/`: Visualization utility that generates genome-wide, karyotype-style coverage plots.
- `04_qc_sentinel/`: Recursive log aggregator that compiles Picard, samtools, and mosdepth metrics into a unified terminal dashboard and HTML report.
- `05_aligner_benchmark/`: Standardized harness to benchmark alignment rate, insert size, memory usage, and runtime across BWA-MEM2, minimap2, and HISAT2.

## Installation and Setup

### 1. Conda Environment
To install the required dependencies (aligners, samtools, python libraries), create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate ngs-mapping
```

### 2. Prepare Reference Genome & Simulated Dataset
Run the preparation script to download the human Chromosome 22 reference genome from UCSC and simulate 50,000 paired-end sequencing reads (excluding centromeric/gap regions):

```bash
python download_and_prep.py
```

## Running the Pipelines

### Alignment Pipelines (WGS and RNA-seq)
To execute the alignment pipelines on the Chromosome 22 dataset:

```bash
bash run_pipelines.sh
```

This script automates:
- Building BWA-MEM2 and HISAT2 indexes.
- Running WGS alignment (Picard duplicate marking and mosdepth coverage).
- Running RNA-seq splice-aware alignment (Picard duplicate marking, Qualimap RNA-seq metrics, and the RSeQC battery).

### Downstream Analysis & Benchmark
To run the downstream reports, visualizations, and benchmarking comparisons:

```bash
bash run_downstream.sh
```

This script generates:
- Chromosome 22 coverage karyotype art (`results/chr22_panorama.png`).
- Unified QC Sentinel HTML and terminal dashboard reports.
- Multi-aligner benchmarks comparing HISAT2, BWA-MEM2, and minimap2.

## Pipeline Configurations

The individual scripts inside the pipeline folders support custom inputs and multi-threading options. Refer to the documentation in each project subdirectory for more details:
- [WGS Alignment pipeline](./01_wgs_pipeline/README.md)
- [RNA-seq Alignment suite](./02_rnaseq_pipeline/README.md)
- [Coverage Panorama Visualization](./03_coverage_panorama/README.md)
- [QC Sentinel Terminal Dashboard](./04_qc_sentinel/README.md)
- [Aligner Benchmark Comparison](./05_aligner_benchmark/README.md)
