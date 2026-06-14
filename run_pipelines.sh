#!/usr/bin/env bash
# =============================================================================
#  run_pipelines.sh — NGS Mapping & Analysis Suite Orchestrator
#  Runs WGS alignment and RNA-seq alignment on chr22 dataset sequentially.
# =============================================================================
set -euo pipefail

# Enable Conda in current script environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate wgs_align_env
export PATH="$PATH:/home/yer_kanat/miniconda3/envs/QC_fastq/bin:/home/yer_kanat/miniconda3/envs/ucsc/bin"

echo "=================================================="
echo "  NGS Mapping & Analysis Suite Orchestrator"
echo "  Running WGS and RNA-seq pipelines on chr22"
echo "=================================================="

echo -e "
>>> STEP 1: Running WGS Pipeline (wgs_align.sh)..."
bash 01_wgs_pipeline/wgs_align.sh \
  -s chr22_sample \
  -1 data/chr22_sample_R1.fastq.gz \
  -2 data/chr22_sample_R2.fastq.gz \
  -r reference/chr22.fa \
  -o results/ \
  -t 8

echo -e "
>>> STEP 2: Building HISAT2 index for chr22..."
if [ ! -f reference/chr22_hisat2_index.1.ht2 ]; then
  hisat2-build -p 8 reference/chr22.fa reference/chr22_hisat2_index
else
  echo "HISAT2 index already exists, skipping index build."
fi

echo -e "
>>> STEP 3: Running RNA-seq Pipeline (rnaseq_hisat2.sh)..."
bash 02_rnaseq_pipeline/rnaseq_hisat2.sh \
  -s chr22_sample \
  -1 data/chr22_sample_R1.fastq.gz \
  -2 data/chr22_sample_R2.fastq.gz \
  -i reference/chr22_hisat2_index \
  -a reference/chr22.gtf \
  -o rnaseq_results/ \
  -t 8

echo -e "
>>> STEP 4: Running RSeQC QC Battery (rseqc_suite.sh)..."
bash 02_rnaseq_pipeline/rseqc_suite.sh \
  -b rnaseq_results/aligned/chr22_sample.dedup.bam \
  -r reference/chr22.bed \
  -o rnaseq_results/rseqc/ \
  -s chr22_sample

echo -e "
=================================================="
echo "  [COMPLETE] All Chromosome 22 runs finished successfully!"
echo "  WGS MultiQC Report:      results/multiqc/chr22_sample_multiqc_report.html"
echo "  RNA-seq MultiQC Report:  rnaseq_results/multiqc/chr22_sample_multiqc.html"
echo "=================================================="
