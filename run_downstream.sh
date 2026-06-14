#!/usr/bin/env bash
# =============================================================================
#  run_downstream.sh — Run Projects 03, 04, 05 Downstream Pipelines
# =============================================================================
set -euo pipefail

# Enable Conda in current script environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate wgs_align_env
export PATH="$PATH:/home/yer_kanat/miniconda3/envs/QC_fastq/bin:/home/yer_kanat/miniconda3/envs/ucsc/bin"

echo "=================================================="
echo "  NGS Mapping & Analysis Suite Downstream Runner"
echo "  Running Projects 03, 04, and 05 on chr22"
echo "=================================================="

echo -e "
>>> PROJECT 03: Generating Coverage Panorama (coverage_panorama.py)..."
python 03_coverage_panorama/coverage_panorama.py \
  -s results/qc/chr22_sample.regions.bed.gz \
  -n "Chr22 WGS Sample (50k reads)" \
  -o results/chr22_panorama.png

echo -e "
>>> PROJECT 04: Running QC Sentinel Dashboard..."
# Single sample report (WGS)
python 04_qc_sentinel/qc_sentinel.py \
  --dir results/ \
  --html results/wgs_qc_report.html \
  --json results/wgs_qc_summary.json

# Single sample report (RNA-seq)
python 04_qc_sentinel/qc_sentinel.py \
  --dir rnaseq_results/ \
  --html rnaseq_results/rnaseq_qc_report.html \
  --json rnaseq_results/rnaseq_qc_summary.json

echo -e "
>>> PROJECT 05: Running Aligner Benchmark (run_benchmark.sh)..."
bash 05_aligner_benchmark/run_benchmark.sh \
  -1 data/chr22_sample_R1.fastq.gz \
  -2 data/chr22_sample_R2.fastq.gz \
  -r reference/chr22.fa \
  -x reference/chr22_hisat2_index \
  -o benchmark_results \
  -n 50000 \
  -t 8

echo -e "
=================================================="
echo "  [COMPLETE] Downstream runs finished successfully!"
echo "  Project 03 Panorama:    results/chr22_panorama.png"
echo "  Project 04 WGS Report:  results/wgs_qc_report.html"
echo "  Project 04 RNA Report:  rnaseq_results/rnaseq_qc_report.html"
echo "  Project 05 Report:      benchmark_results/benchmark_report.html"
echo "=================================================="
