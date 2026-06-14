# NGS Mapping and Analysis Suite

A collection of NGS alignment pipelines and quality control tools for whole-genome sequencing (WGS) and RNA-seq.

## Projects

- **01 WGS Alignment Pipeline**: Single-sample alignment using BWA-MEM2, duplicate marking with Picard, and coverage with mosdepth.
- **02 RNA-seq Alignment Suite**: Splice-aware alignment using HISAT2 and RSeQC transcript quality control.
- **03 Coverage Panorama**: Visualization tool to render genome-wide coverage as a chromosome paint map.
- **04 QC Sentinel**: Unified terminal dashboard and HTML report generator for NGS alignment metrics.
- **05 Aligner Benchmark**: Benchmarking harness comparing BWA-MEM2, minimap2, and HISAT2.

## Quick Start

### 1. Setup Environment
Ensure you have Conda or Mamba installed, then create and activate the environment:
```bash
conda env create -f environment.yml
conda activate ngs-mapping
```

### 2. Download Reference and Test Data
Prepare the human Chromosome 22 reference genome and simulate paired-end sequencing reads:
```bash
python download_and_prep.py
```

### 3. Run Pipelines
Execute the WGS and RNA-seq alignment pipelines:
```bash
bash run_pipelines.sh
```

### 4. Run Downstream Analysis
Generate coverage panorama art, QC dashboards, and aligner benchmarks:
```bash
bash run_downstream.sh
```

## Requirements
- Linux or macOS (WSL2 supported on Windows)
- Conda package manager
- Minimum 16GB RAM recommended for Chromosome 22 runs (32GB+ for full human genome)
