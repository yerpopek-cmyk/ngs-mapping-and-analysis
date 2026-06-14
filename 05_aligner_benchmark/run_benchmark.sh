#!/usr/bin/env bash
# =============================================================================
#  run_benchmark.sh — Multi-Aligner Benchmarking Harness
#
#  Runs BWA-MEM2, minimap2, and HISAT2 on the same FASTQ input and
#  collects runtime, memory, and alignment quality metrics for comparison.
#
#  Usage:
#    bash run_benchmark.sh -1 R1.fastq.gz -2 R2.fastq.gz \
#                          -r reference.fa \
#                          -x hisat2_index \
#                          -o benchmark_results [OPTIONS]
#
#  Options:
#    -1  Read 1 FASTQ (required)
#    -2  Read 2 FASTQ (required)
#    -r  Reference FASTA for BWA-MEM2 + minimap2 (required)
#    -x  HISAT2 index prefix (optional, skip if not provided)
#    -o  Output directory [default: benchmark_results]
#    -t  Threads [default: 8]
#    -n  Number of reads to subsample for speed [default: 2000000]
#    -h  Help
# =============================================================================
set -euo pipefail
trap 'printf "\n[ERROR] Benchmark halted at line %s\n" "$LINENO"' ERR

THREADS=8
OUTDIR="benchmark_results"
N_READS=2000000
HISAT2_INDEX=""

while getopts ":1:2:r:x:o:t:n:h" opt; do
  case $opt in
    1) R1="$OPTARG" ;; 2) R2="$OPTARG" ;;
    r) REF="$OPTARG" ;; x) HISAT2_INDEX="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;; t) THREADS="$OPTARG" ;;
    n) N_READS="$OPTARG" ;; h) grep "^#  " "$0" | sed 's/^#  //'; exit 0 ;;
    :) printf "[ERROR] -%s requires an argument\n" "$OPTARG"; exit 1 ;;
  esac
done

[[ -z "${R1:-}" || -z "${R2:-}" || -z "${REF:-}" ]] && {
  printf "[ERROR] -1, -2, -r are required.\n"; exit 1
}

mkdir -p "${OUTDIR}"/{bam,logs,metrics,subsample}
METRICS="${OUTDIR}/metrics/benchmark_raw.tsv"
printf "aligner\truntime_sec\tmax_rss_mb\talignment_rate\tproperly_paired_pct\tmapped_reads\ttotal_reads\n" \
  > "$METRICS"

printf "═%.0s" {1..65}; printf "\n"
printf "  Aligner Benchmark Harness\n"
printf "  Subsampling to %'d reads\n" "$N_READS"
printf "  Threads: %s\n" "$THREADS"
printf "═%.0s" {1..65}; printf "\n"

# ── Helper: time + memory wrapper ─────────────────────────────────────────────
# Usage: time_cmd OUTPUT_PREFIX cmd args...
time_cmd() {
  local prefix="$1"; shift
  local time_file="${OUTDIR}/logs/${prefix}_time.txt"
  # GNU time format: elapsed_sec rss_kb
  /usr/bin/time -f "%e %M" -o "$time_file" "$@"
  local elapsed; elapsed=$(awk '{print $1}' "$time_file")
  local rss_mb;  rss_mb=$(awk '{printf "%.1f", $2/1024}' "$time_file")
  echo "$elapsed $rss_mb"
}

# ── Helper: parse alignment rate from flagstat ────────────────────────────────
parse_flagstat() {
  local bam="$1"
  samtools flagstat "$bam" > "${OUTDIR}/metrics/$(basename "$bam" .bam).flagstat.txt"
  local flagstat="${OUTDIR}/metrics/$(basename "$bam" .bam).flagstat.txt"
  local total;  total=$(grep "in total"   "$flagstat" | awk '{print $1}')
  local mapped; mapped=$(grep "^[0-9]* + [0-9]* mapped" "$flagstat" | awk '{print $1}')
  local paired; paired=$(grep "properly paired" "$flagstat" | grep -o '[0-9.]*%' | tr -d '%' || echo "0.0")
  local paired_pct; paired_pct="${paired}"
  local rate; rate=$(echo "$mapped $total" | awk '{printf "%.2f", $1/$2*100}')
  echo "$rate $paired_pct $mapped $total"
}

# ── Step 1: Subsample reads ───────────────────────────────────────────────────
printf "\n[STEP 1/5] Subsampling to %'d reads with seqtk\n" "$N_READS"
seqtk sample -s42 "$R1" "$N_READS" | gzip > "${OUTDIR}/subsample/sub_R1.fastq.gz"
seqtk sample -s42 "$R2" "$N_READS" | gzip > "${OUTDIR}/subsample/sub_R2.fastq.gz"
SUB_R1="${OUTDIR}/subsample/sub_R1.fastq.gz"
SUB_R2="${OUTDIR}/subsample/sub_R2.fastq.gz"
printf "  [✓] Subsampled FASTQs ready\n"

# ── Step 2: BWA-MEM2 ──────────────────────────────────────────────────────────
printf "\n[STEP 2/5] Running BWA-MEM2\n"
ALIGNER="bwa-mem2"
read -r elapsed rss < <(
  time_cmd bwa_mem2 \
    sh -c "bwa-mem2 mem -t $THREADS -M -Y -R '@RG\tID:bench\tSM:bench\tPL:ILLUMINA' $REF $SUB_R1 $SUB_R2 2>${OUTDIR}/logs/bwa_mem2.log | samtools sort -@ 4 -m 2G -o ${OUTDIR}/bam/bwa_mem2.bam -"
)
samtools index "${OUTDIR}/bam/bwa_mem2.bam"
read -r rate paired mapped total < <(parse_flagstat "${OUTDIR}/bam/bwa_mem2.bam")
printf "  Runtime: %ss  MaxRSS: %sMB  Mapped: %s%%\n" "$elapsed" "$rss" "$rate"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "bwa-mem2" "$elapsed" "$rss" "$rate" "$paired" "$mapped" "$total" >> "$METRICS"

# ── Step 3: minimap2 (short-read mode) ───────────────────────────────────────
printf "\n[STEP 3/5] Running minimap2 (sr mode)\n"
read -r elapsed rss < <(
  time_cmd minimap2 \
    sh -c "minimap2 -ax sr -t $THREADS -R '@RG\tID:bench\tSM:bench\tPL:ILLUMINA' $REF $SUB_R1 $SUB_R2 2>${OUTDIR}/logs/minimap2.log | samtools sort -@ 4 -m 2G -o ${OUTDIR}/bam/minimap2.bam -"
)
samtools index "${OUTDIR}/bam/minimap2.bam"
read -r rate paired mapped total < <(parse_flagstat "${OUTDIR}/bam/minimap2.bam")
printf "  Runtime: %ss  MaxRSS: %sMB  Mapped: %s%%\n" "$elapsed" "$rss" "$rate"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "minimap2" "$elapsed" "$rss" "$rate" "$paired" "$mapped" "$total" >> "$METRICS"

# ── Step 4: HISAT2 (optional) ─────────────────────────────────────────────────
if [[ -n "$HISAT2_INDEX" ]]; then
  printf "\n[STEP 4/5] Running HISAT2 (splice-aware)\n"
  read -r elapsed rss < <(
    time_cmd hisat2 \
      sh -c "hisat2 -p $THREADS -x $HISAT2_INDEX -1 $SUB_R1 -2 $SUB_R2 --rg-id bench --rg 'SM:bench' --summary-file ${OUTDIR}/logs/hisat2_summary.txt 2>${OUTDIR}/logs/hisat2.log | samtools sort -@ 4 -m 2G -o ${OUTDIR}/bam/hisat2.bam -"
  )
  samtools index "${OUTDIR}/bam/hisat2.bam"
  read -r rate paired mapped total < <(parse_flagstat "${OUTDIR}/bam/hisat2.bam")
  printf "  Runtime: %ss  MaxRSS: %sMB  Mapped: %s%%\n" "$elapsed" "$rss" "$rate"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "hisat2" "$elapsed" "$rss" "$rate" "$paired" "$mapped" "$total" >> "$METRICS"
else
  printf "\n[STEP 4/5] Skipping HISAT2 (no -x index provided)\n"
fi

# ── Step 5: Generate comparison report ───────────────────────────────────────
printf "\n[STEP 5/5] Generating comparison report\n"
python "$(dirname "$0")/compare_aligners.py" \
  --metrics "$METRICS" \
  --output "${OUTDIR}/benchmark_report.html" \
  --plot   "${OUTDIR}/benchmark_plots.png" \
  --reads  "$N_READS"

printf "\n[BENCHMARK COMPLETE]\n"
printf "  Raw metrics: %s\n" "$METRICS"
printf "  Report:      %s/benchmark_report.html\n" "$OUTDIR"
printf "  Plots:       %s/benchmark_plots.png\n" "$OUTDIR"
