#!/usr/bin/env bash
# =============================================================================
#  wgs_align.sh — Production WGS/WES Alignment Pipeline
#  Follows GATK Best Practices: BWA-MEM2 → Picard MarkDuplicates → mosdepth
#
#  Usage:
#    bash wgs_align.sh -s SAMPLE -1 R1.fastq.gz -2 R2.fastq.gz \
#                      -r reference.fa -o output_dir [OPTIONS]
#
#  Options:
#    -s  Sample name (required)
#    -1  Read 1 FASTQ (required)
#    -2  Read 2 FASTQ (required)
#    -r  Reference genome FASTA (required, must be indexed)
#    -o  Output directory [default: results]
#    -t  Threads [default: 8]
#    -m  Memory per sort thread [default: 4G]
#    -b  Bed file for WES target regions [optional]
#    -h  Show this help
# =============================================================================
set -euo pipefail
trap 'printf "\n[ERROR] Pipeline failed at line %s — see logs/\n" "$LINENO"' ERR

# ── defaults ──────────────────────────────────────────────────────────────────
THREADS=8
MEM_PER_THREAD="4G"
OUTDIR="results"
BEDFILE=""

usage() {
  grep "^#  " "$0" | sed 's/^#  //'
  exit 1
}

# ── argument parsing ──────────────────────────────────────────────────────────
while getopts ":s:1:2:r:o:t:m:b:h" opt; do
  case $opt in
    s) SAMPLE="$OPTARG" ;;
    1) R1="$OPTARG" ;;
    2) R2="$OPTARG" ;;
    r) REF="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    m) MEM_PER_THREAD="$OPTARG" ;;
    b) BEDFILE="$OPTARG" ;;
    h) usage ;;
    :) printf "[ERROR] Option -%s requires an argument.\n" "$OPTARG"; exit 1 ;;
    \?) printf "[ERROR] Unknown option: -%s\n" "$OPTARG"; exit 1 ;;
  esac
done

[[ -z "${SAMPLE:-}" || -z "${R1:-}" || -z "${R2:-}" || -z "${REF:-}" ]] && {
  printf "[ERROR] -s, -1, -2, -r are required.\n"
  usage
}

# ── setup output directories ──────────────────────────────────────────────────
ALIGNED=$(realpath -m "${OUTDIR}/aligned")
QC=$(realpath -m      "${OUTDIR}/qc")
LOGS=$(realpath -m    "${OUTDIR}/logs")
mkdir -p "$ALIGNED" "$QC" "$LOGS"

TOTAL_STEPS=7
step() { printf "\n[STEP %s/%s] %s\n" "$1" "$TOTAL_STEPS" "$2"; }

printf "═%.0s" {1..65}; printf "\n"
printf "  WGS Alignment Pipeline — %s\n" "$SAMPLE"
printf "  Reference: %s\n" "$REF"
printf "  Threads:   %s\n" "$THREADS"
printf "  Output:    %s\n" "$OUTDIR"
printf "═%.0s" {1..65}; printf "\n"

# ── Step 1: Validate reference index ─────────────────────────────────────────
step 1 "Validating reference genome index"
for ext in .bwt.2bit.64 .0123 .fai; do
  [[ -f "${REF}${ext}" ]] || {
    printf "  [WARN] Missing index: %s%s\n" "$REF" "$ext"
    if [[ "$ext" == ".fai" ]]; then
      printf "  Running: samtools faidx %s\n" "$REF"
      samtools faidx "$REF"
    else
      printf "  Running: bwa-mem2 index %s\n" "$REF"
      bwa-mem2 index "$REF" 2>"${LOGS}/bwa_index.log"
    fi
  }
done
printf "  [✓] Reference index OK\n"

# ── Step 2: FASTQ quality pre-check ──────────────────────────────────────────
step 2 "Pre-alignment FASTQ QC (fastp)"
fastp \
  --in1 "$R1" --in2 "$R2" \
  --json "${QC}/${SAMPLE}_fastp.json" \
  --html "${QC}/${SAMPLE}_fastp.html" \
  --thread "$THREADS" \
  --disable_adapter_trimming \
  --qualified_quality_phred 20 \
  --n_base_limit 5 \
  --low_complexity_filter \
  --report_title "${SAMPLE} pre-alignment QC" \
  2>"${LOGS}/${SAMPLE}_fastp.log"
printf "  [✓] fastp report: %s/%s_fastp.html\n" "$QC" "$SAMPLE"

# ── Step 3: Align with BWA-MEM2 ──────────────────────────────────────────────
step 3 "Aligning reads with BWA-MEM2"
RG_ID="${SAMPLE}_L001"
RG_LINE="@RG\tID:${RG_ID}\tSM:${SAMPLE}\tPU:flowcell.lane1\tLB:${SAMPLE}_lib\tPL:ILLUMINA"

bwa-mem2 mem \
  -t "$THREADS" \
  -M \
  -Y \
  -K 100000000 \
  -R "$RG_LINE" \
  "$REF" "$R1" "$R2" \
  2> >(tee "${LOGS}/${SAMPLE}_bwa.log" >&2) | \
  samtools sort \
    -@ "$((THREADS / 2))" \
    -m "$MEM_PER_THREAD" \
    -o "${ALIGNED}/${SAMPLE}.sorted.bam" \
    -T "${ALIGNED}/${SAMPLE}_sort_tmp" -

printf "  [✓] Sorted BAM: %s/%s.sorted.bam\n" "$ALIGNED" "$SAMPLE"

# ── Step 4: Mark duplicates (Picard) ─────────────────────────────────────────
step 4 "Marking PCR and optical duplicates (Picard)"
picard MarkDuplicates \
  I="${ALIGNED}/${SAMPLE}.sorted.bam" \
  O="${ALIGNED}/${SAMPLE}.markdup.bam" \
  M="${QC}/${SAMPLE}.markdup_metrics.txt" \
  REMOVE_DUPLICATES=false \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=LENIENT \
  CREATE_INDEX=true \
  2>"${LOGS}/${SAMPLE}_picard.log"

# Report duplication rate
DUP_RATE=$(awk 'found {printf "%.1f%%", $9*100; exit} /^LIBRARY/ {found=1}' \
  "${QC}/${SAMPLE}.markdup_metrics.txt")
printf "  [✓] Duplication rate: %s\n" "$DUP_RATE"
printf "  [✓] Metrics: %s/%s.markdup_metrics.txt\n" "$QC" "$SAMPLE"

# ── Step 5: Alignment statistics ─────────────────────────────────────────────
step 5 "Collecting alignment statistics (samtools)"
samtools flagstat -@ "$THREADS" \
  "${ALIGNED}/${SAMPLE}.markdup.bam" > "${QC}/${SAMPLE}.flagstat.txt"

samtools stats -@ "$THREADS" \
  "${ALIGNED}/${SAMPLE}.markdup.bam" | \
  grep "^SN" > "${QC}/${SAMPLE}.stats.txt"

# Extract and display key metrics
MAPPED=$(grep "mapped (" "${QC}/${SAMPLE}.flagstat.txt" | head -1 | grep -o '[0-9.]*%' || echo "N/A")
printf "  [✓] Mapped reads: %s\n" "$MAPPED"

# ── Step 6: Coverage with mosdepth ───────────────────────────────────────────
step 6 "Calculating coverage depth (mosdepth)"
MOSDEPTH_ARGS="-t $THREADS -n --fast-mode -F 1024 -b 1000"

if [[ -n "$BEDFILE" ]]; then
  MOSDEPTH_ARGS="$MOSDEPTH_ARGS --by $BEDFILE"
  printf "  [INFO] WES mode — using target BED: %s\n" "$BEDFILE"
fi

# shellcheck disable=SC2086
mosdepth $MOSDEPTH_ARGS \
  "${QC}/${SAMPLE}" \
  "${ALIGNED}/${SAMPLE}.markdup.bam" \
  2>"${LOGS}/${SAMPLE}_mosdepth.log"

MEAN_COV=$( (zcat "${QC}/${SAMPLE}.mosdepth.summary.txt.gz" 2>/dev/null || \
  cat "${QC}/${SAMPLE}.mosdepth.summary.txt" 2>/dev/null) | \
  awk 'NR>1 && $1=="total" {printf "%.1fX", $4}')
printf "  [✓] Mean coverage: %s\n" "$MEAN_COV"

# ── Step 7: MultiQC aggregate report ─────────────────────────────────────────
step 7 "Generating MultiQC aggregate report"
multiqc "$OUTDIR" \
  -n "${SAMPLE}_multiqc_report" \
  -o "${OUTDIR}/multiqc" \
  -f -z \
  2>"${LOGS}/${SAMPLE}_multiqc.log"
printf "  [✓] Report: %s/multiqc/%s_multiqc_report.html\n" "$OUTDIR" "$SAMPLE"

# ── Summary ───────────────────────────────────────────────────────────────────
printf "\n"
printf "═%.0s" {1..65}; printf "\n"
printf "  [PIPELINE COMPLETE] %s\n" "$SAMPLE"
printf "  Steps completed: %s/%s\n" "$TOTAL_STEPS" "$TOTAL_STEPS"
printf "\n"
printf "  OUTPUT FILES:\n"
printf "    BAM:     %s/%s.markdup.bam\n"           "$ALIGNED" "$SAMPLE"
printf "    flagstat:%s/%s.flagstat.txt\n"           "$QC" "$SAMPLE"
printf "    Markdup: %s/%s.markdup_metrics.txt\n"   "$QC" "$SAMPLE"
printf "    Coverage:%s/%s.mosdepth.summary.txt\n"  "$QC" "$SAMPLE"
printf "    Report:  %s/multiqc/%s_multiqc_report.html\n" "$OUTDIR" "$SAMPLE"
printf "═%.0s" {1..65}; printf "\n"
