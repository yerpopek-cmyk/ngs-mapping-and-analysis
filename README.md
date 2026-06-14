# 🧬 NGS Mapping and Analysis Suite

> A curated collection of production-grade NGS alignment pipelines and innovative bioinformatics tools for the bioinformatics community.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-brightgreen.svg)](https://nextflow.io/)

---

## 📦 Projects

| # | Project | Type | Tools | Description |
|---|---------|------|-------|-------------|
| 01 | [WGS Alignment Pipeline](./01_wgs_pipeline/) | Pipeline | BWA-MEM2, Picard, mosdepth | Production-grade whole-genome alignment |
| 02 | [RNA-seq Alignment Suite](./02_rnaseq_pipeline/) | Pipeline | HISAT2, STAR, RSeQC | Splice-aware alignment + full QC battery |
| 03 | [Coverage Panorama](./03_coverage_panorama/) | 🎨 Tool | Python, mosdepth | Karyotype-style genome coverage art |
| 04 | [QC Sentinel](./04_qc_sentinel/) | 🖥️ Tool | Python Rich | Unified terminal QC dashboard |
| 05 | [Aligner Benchmark](./05_aligner_benchmark/) | 📊 Tool | Bash + Python | Multi-aligner comparison harness |

---

## 🚀 Quick Start

```bash
# 1. Clone the repo
git clone https://github.com/yerpopek-cmyk/ngs-mapping-and-analysis.git
cd ngs-mapping-and-analysis

# 2. Set up the conda environment
conda env create -f environment.yml
conda activate ngs-mapping

# 3. Run any project
bash 01_wgs_pipeline/wgs_align.sh --help
python 03_coverage_panorama/coverage_panorama.py --help
python 04_qc_sentinel/qc_sentinel.py --help
```

---

## 🧬 Skill Basis

All pipelines follow the **Bioinformatics Mapping and Alignment Skill v2.0.0**,
implementing GATK Best Practices and ENCODE consortium guidelines. See
[skills/skillmapping.md](./skills/skillmapping.md) for the full reference.

---

## 📋 Requirements

- Linux / macOS (WSL2 supported)
- Conda or Mamba (recommended)
- 32GB RAM minimum for human genome alignment
- 8+ CPU cores recommended

**Core tools** (installed via `environment.yml`):
`bwa-mem2 ≥ 2.2.1`, `samtools ≥ 1.15`, `picard ≥ 2.27`, `mosdepth ≥ 0.3.3`,
`hisat2 ≥ 2.2.1`, `STAR ≥ 2.7.10a`, `RSeQC ≥ 3.0.1`, `Qualimap ≥ 2.2.2`,
`preseq ≥ 3.1.2`, `fastp ≥ 0.23`, `MultiQC ≥ 1.14`

---

## 📁 Repository Structure

```
ngs-mapping-and-analysis/
├── README.md
├── LICENSE
├── .gitignore
├── environment.yml
│
├── 01_wgs_pipeline/           # Whole-genome alignment (Bash + Nextflow)
│   ├── wgs_align.sh
│   ├── nextflow/
│   │   ├── main.nf
│   │   └── nextflow.config
│   └── README.md
│
├── 02_rnaseq_pipeline/        # RNA-seq splice-aware alignment
│   ├── rnaseq_hisat2.sh
│   ├── rnaseq_star.sh
│   ├── rseqc_suite.sh
│   └── README.md
│
├── 03_coverage_panorama/      # 🎨 Genome coverage visualization art
│   ├── coverage_panorama.py
│   ├── karyotype_painter.py
│   ├── requirements.txt
│   └── README.md
│
├── 04_qc_sentinel/            # 🖥️ Terminal QC dashboard
│   ├── qc_sentinel.py
│   ├── parsers/
│   │   ├── __init__.py
│   │   ├── flagstat_parser.py
│   │   ├── picard_parser.py
│   │   └── mosdepth_parser.py
│   ├── requirements.txt
│   └── README.md
│
├── 05_aligner_benchmark/      # 📊 Multi-aligner benchmarking
│   ├── run_benchmark.sh
│   ├── compare_aligners.py
│   ├── requirements.txt
│   └── README.md
│
└── skills/
    └── skillmapping.md        # Mapping skill reference
```

---

## 🤝 Contributing

Pull requests welcome. For major changes, open an issue first.

## 📄 License

[MIT](./LICENSE) — GPTomics / yerpopek-cmyk
