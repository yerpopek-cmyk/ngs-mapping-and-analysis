#!/usr/bin/env bash
# =============================================================================
#  rnaseq_hisat2.sh — Splice-Aware RNA-seq Alignment Pipeline (HISAT2)
#
#  Usage:
#    bash rnaseq_hisat2.sh -s SAMPLE -1 R1.fastq.gz -2 R2.fastq.gz \
#                          -i HISAT2_INDEX -o output_dir [OPTIONS]
#
#  Options:
#    -s  Sample name (required)
#    -1  Read 1 FASTQ (required)
#    -2  Read 2 FASTQ (required)
#    -i  HISAT2 index prefix (required)
#    -a  GTF annotation file (required for Qualimap)
#    -o  Output directory [default: results]
#    -t  Threads [default: 8]
#    -S  Strandedness: RF|FR|unstranded [default: RF]
#    -h  Show this help
# =============================================================================
set -euo pipefail
trap 'printf "\n[ERROR] Pipeline halted at line %s\n" "$LINENO"' ERR

# ── defaults ──────────────────────────────────────────────────────────────────
THREADS=8
OUTDIR="results"
STRAND="RF"
GTF=""

while getopts ":s:1:2:i:a:o:t:S:h" opt; do
  case $opt in
    s) SAMPLE="$OPTARG" ;;
    1) R1="$OPTARG" ;;
    2) R2="$OPTARG" ;;
    i) INDEX="$OPTARG" ;;
    a) GTF="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    S) STRAND="$OPTARG" ;;
    h) grep "^#  " "$0" | sed 's/^#  //'; exit 0 ;;
    :) printf "[ERROR] -%s requires an argument\n" "$OPTARG"; exit 1 ;;
    \?) printf "[ERROR] Unknown option: -%s\n" "$OPTARG"; exit 1 ;;
  esac
done

[[ -z "${SAMPLE:-}" || -z "${R1:-}" || -z "${R2:-}" || -z "${INDEX:-}" ]] && {
  printf "[ERROR] -s, -1, -2, -i are required.\n"; exit 1
}

# ── directories ───────────────────────────────────────────────────────────────
ALIGNED=$(realpath -m "${OUTDIR}/aligned")
QC=$(realpath -m      "${OUTDIR}/qc")
LOGS=$(realpath -m    "${OUTDIR}/logs")
mkdir -p "$ALIGNED" "$QC" "$LOGS"

TOTAL_STEPS=7
step() { printf "\n[STEP %s/%s] %s\n" "$1" "$TOTAL_STEPS" "$2"; }

printf "═%.0s" {1..65}; printf "\n"
printf "  RNA-seq HISAT2 Pipeline — %s  [%s strand]\n" "$SAMPLE" "$STRAND"
printf "═%.0s" {1..65}; printf "\n"

# ── Strandedness mapping ───────────────────────────────────────────────────────
case "$STRAND" in
  RF) HISAT2_STRAND="--rna-strandness RF" ; QUALIMAP_STRAND="strand-specific-reverse" ;;
  FR) HISAT2_STRAND="--rna-strandness FR" ; QUALIMAP_STRAND="strand-specific-forward" ;;
  *)  HISAT2_STRAND=""                    ; QUALIMAP_STRAND="non-strand-specific" ;;
esac

# ── Step 1: FASTQ preprocessing ───────────────────────────────────────────────
step 1 "Preprocessing FASTQ with fastp"
fastp \
  --in1 "$R1" --in2 "$R2" \
  --out1 "${OUTDIR}/trimmed_R1.fastq.gz" \
  --out2 "${OUTDIR}/trimmed_R2.fastq.gz" \
  --json "${QC}/${SAMPLE}_fastp.json" \
  --html "${QC}/${SAMPLE}_fastp.html" \
  --thread "$THREADS" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 36 \
  --low_complexity_filter \
  2>"${LOGS}/${SAMPLE}_fastp.log"
printf "  [✓] Trimmed reads ready\n"

# ── Step 2: Align with HISAT2 ─────────────────────────────────────────────────
step 2 "Aligning with HISAT2 (splice-aware)"
# shellcheck disable=SC2086
hisat2 -p "$THREADS" \
  -x "$INDEX" \
  -1 "${OUTDIR}/trimmed_R1.fastq.gz" \
  -2 "${OUTDIR}/trimmed_R2.fastq.gz" \
  $HISAT2_STRAND \
  --rg-id "$SAMPLE" \
  --rg "SM:${SAMPLE}" \
  --rg "PL:ILLUMINA" \
  --novel-splicesite-outfile "${QC}/${SAMPLE}_novel_splicesites.txt" \
  --summary-file "${QC}/${SAMPLE}_hisat2_summary.txt" \
  --dta \
  2>"${LOGS}/${SAMPLE}_hisat2.log" | \
  samtools sort -@ "$((THREADS / 2))" -m 4G \
    -o "${ALIGNED}/${SAMPLE}.sorted.bam" -

ALIGN_RATE=$(grep "overall alignment rate" "${QC}/${SAMPLE}_hisat2_summary.txt" | \
  awk '{print $1}')
printf "  [✓] Overall alignment rate: %s\n" "$ALIGN_RATE"

# ── Step 3: Mark duplicates ────────────────────────────────────────────────────
step 3 "Marking duplicates (keep for RNA-seq)"
picard MarkDuplicates \
  I="${ALIGNED}/${SAMPLE}.sorted.bam" \
  O="${ALIGNED}/${SAMPLE}.markdup.bam" \
  M="${QC}/${SAMPLE}.markdup_metrics.txt" \
  REMOVE_DUPLICATES=false \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=LENIENT \
  CREATE_INDEX=true \
  2>"${LOGS}/${SAMPLE}_picard.log"
printf "  [✓] Marked BAM: %s/%s.markdup.bam\n" "$ALIGNED" "$SAMPLE"

# ── Step 4: Create dedup BAM for RSeQC ────────────────────────────────────────
step 4 "Creating duplicate-removed BAM for RSeQC"
samtools view -@ "$THREADS" -h -b -F 1024 \
  "${ALIGNED}/${SAMPLE}.markdup.bam" > "${ALIGNED}/${SAMPLE}.dedup.bam"
samtools index "${ALIGNED}/${SAMPLE}.dedup.bam"
printf "  [✓] Dedup BAM ready for RSeQC\n"

# ── Step 5: samtools QC ────────────────────────────────────────────────────────
step 5 "samtools flagstat and stats"
samtools flagstat -@ "$THREADS" \
  "${ALIGNED}/${SAMPLE}.markdup.bam" > "${QC}/${SAMPLE}.flagstat.txt"
printf "  [✓] flagstat written\n"

# ── Step 6: Qualimap RNA-seq report ───────────────────────────────────────────
step 6 "Qualimap RNA-seq QC report"
if [[ -n "$GTF" ]]; then
  qualimap rnaseq \
    -bam "${ALIGNED}/${SAMPLE}.markdup.bam" \
    -gtf "$GTF" \
    -outdir "${QC}/qualimap_${SAMPLE}" \
    --java-mem-size=8G \
    -p "$QUALIMAP_STRAND" \
    -pe \
    2>"${LOGS}/${SAMPLE}_qualimap.log"
  printf "  [✓] Qualimap report: %s/qualimap_%s/\n" "$QC" "$SAMPLE"
else
  printf "  [SKIP] No GTF provided; skipping Qualimap rnaseq\n"
fi

# ── Step 7: MultiQC ───────────────────────────────────────────────────────────
step 7 "MultiQC aggregate report"
multiqc "$OUTDIR" -n "${SAMPLE}_multiqc" -o "${OUTDIR}/multiqc" -f -z \
  2>"${LOGS}/${SAMPLE}_multiqc.log"
printf "  [✓] Report: %s/multiqc/%s_multiqc.html\n" "$OUTDIR" "$SAMPLE"

printf "\n[PIPELINE COMPLETE] %s\n" "$SAMPLE"
printf "  Dedup BAM:  %s/%s.dedup.bam\n"  "$ALIGNED" "$SAMPLE"
printf "  MultiQC:    %s/multiqc/%s_multiqc.html\n" "$OUTDIR" "$SAMPLE"
