#!/usr/bin/env bash
# =============================================================================
#  rseqc_suite.sh — Complete RSeQC Quality Control Battery
#
#  Runs all major RSeQC modules on a de-duplicated BAM file.
#  Requires a BED12 gene annotation file (converted from GTF).
#
#  Usage:
#    bash rseqc_suite.sh -b sample.dedup.bam -r annotation.bed \
#                        -o output_dir [-s SAMPLE] [-t THREADS]
#
#  Generate annotation BED from GTF:
#    gtfToGenePred -genePredExt -ignoreGroupsWithoutExons genes.gtf genes.genePred
#    genePredToBed genes.genePred genes.bed
# =============================================================================
set -euo pipefail
trap 'printf "\n[ERROR] RSeQC suite failed at line %s\n" "$LINENO"' ERR

THREADS=4
OUTDIR="rseqc_results"
SAMPLE=""

while getopts ":b:r:o:s:t:h" opt; do
  case $opt in
    b) BAM="$OPTARG" ;;
    r) BED="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    s) SAMPLE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) grep "^#  " "$0" | sed 's/^#  //'; exit 0 ;;
    :) printf "[ERROR] -%s requires an argument\n" "$OPTARG"; exit 1 ;;
  esac
done

[[ -z "${BAM:-}" || -z "${BED:-}" ]] && {
  printf "[ERROR] -b (BAM) and -r (BED) are required.\n"; exit 1
}

[[ -z "$SAMPLE" ]] && SAMPLE=$(basename "$BAM" .bam)

mkdir -p "$OUTDIR"

TOTAL=8
step() { printf "\n[STEP %s/%s] %s\n" "$1" "$TOTAL" "$2"; }

printf "═%.0s" {1..65}; printf "\n"
printf "  RSeQC Quality Control Suite — %s\n" "$SAMPLE"
printf "  BAM: %s\n" "$BAM"
printf "  BED: %s\n" "$BED"
printf "═%.0s" {1..65}; printf "\n"

# ── Ensure BAM is indexed ──────────────────────────────────────────────────────
[[ -f "${BAM}.bai" ]] || { printf "  Indexing BAM...\n"; samtools index "$BAM"; }

# ── Step 1: Read distribution ──────────────────────────────────────────────────
step 1 "Read distribution (genomic feature assignment)"
read_distribution.py \
  -r "$BED" \
  -i "$BAM" \
  > "${OUTDIR}/${SAMPLE}_read_distribution.txt"
# Parse and display summary
printf "  Feature assignments:\n"
awk '/CDS_Exons|5.UTR|3.UTR|Introns|Intergenic/ {printf "    %-25s %s%%\n", $1, $3}' \
  "${OUTDIR}/${SAMPLE}_read_distribution.txt"

# ── Step 2: Gene body coverage ────────────────────────────────────────────────
step 2 "Gene body coverage (5' to 3' uniformity)"
geneBody_coverage.py \
  -r "$BED" \
  -i "$BAM" \
  -o "${OUTDIR}/${SAMPLE}_genebody" \
  2>/dev/null
printf "  [✓] Curve: %s/%s_genebody.geneBodyCoverage.curves.pdf\n" "$OUTDIR" "$SAMPLE"

# ── Step 3: Junction saturation ───────────────────────────────────────────────
step 3 "Junction saturation (splice site discovery)"
junction_saturation.py \
  -r "$BED" \
  -i "$BAM" \
  -o "${OUTDIR}/${SAMPLE}_junction" \
  2>/dev/null
printf "  [✓] Plot: %s/%s_junction.junctionSaturation_plot.pdf\n" "$OUTDIR" "$SAMPLE"

# ── Step 4: Inner distance (fragment size) ────────────────────────────────────
step 4 "Inner distance distribution (insert size)"
inner_distance.py \
  -r "$BED" \
  -i "$BAM" \
  -o "${OUTDIR}/${SAMPLE}_inner_dist" \
  2>/dev/null
MEAN_DIST=$(awk 'NR>1 && NF>1 {sum+=$1*$2; total+=$2} END {printf "%.0fbp", sum/total}' \
  "${OUTDIR}/${SAMPLE}_inner_dist.inner_distance_freq.txt" 2>/dev/null || echo "N/A")
printf "  [✓] Mean inner distance: %s\n" "$MEAN_DIST"

# ── Step 5: Infer strandedness ────────────────────────────────────────────────
step 5 "Inferring library strandedness"
infer_experiment.py \
  -r "$BED" \
  -i "$BAM" \
  > "${OUTDIR}/${SAMPLE}_strandedness.txt" \
  2>/dev/null
printf "  Strandedness inference:\n"
grep -E "\"1\+\+|\"1\-\-|Fraction" "${OUTDIR}/${SAMPLE}_strandedness.txt" | \
  awk '{printf "    %s\n", $0}' || true

# ── Step 6: Clipping profile ──────────────────────────────────────────────────
step 6 "Read clipping profile (adapter/quality trimming check)"
clipping_profile.py \
  -i "$BAM" \
  -s "PE" \
  -o "${OUTDIR}/${SAMPLE}_clipping" \
  2>/dev/null
printf "  [✓] Clipping profile plot generated\n"

# ── Step 7: BAM statistics ────────────────────────────────────────────────────
step 7 "BAM statistics (bam_stat.py)"
bam_stat.py \
  -i "$BAM" \
  > "${OUTDIR}/${SAMPLE}_bam_stat.txt" \
  2>/dev/null
printf "  Key stats:\n"
grep -E "^(Total|Mapped|Proper|Read1|Read2|Singleton)" \
  "${OUTDIR}/${SAMPLE}_bam_stat.txt" | \
  awk '{printf "    %-30s %s\n", $1, $NF}' || true

# ── Step 8: mRNA contamination / rRNA check ──────────────────────────────────
step 8 "Tin score (RNA degradation index per transcript)"
tin.py \
  -i "$BAM" \
  -r "$BED" \
  > "${OUTDIR}/${SAMPLE}_tin.txt" \
  2>/dev/null || printf "  [WARN] TIN calculation skipped (may need >1000 transcripts)\n"

printf "\n[RSeQC COMPLETE] %s\n" "$SAMPLE"
printf "  All outputs in: %s/\n" "$OUTDIR"
printf "\n  Quick summary of key files:\n"
for f in read_distribution strandedness bam_stat; do
  printf "    %s/%s_%s.txt\n" "$OUTDIR" "$SAMPLE" "$f"
done
printf "\n  Run 'multiqc %s/' to aggregate into HTML report.\n" "$OUTDIR"
