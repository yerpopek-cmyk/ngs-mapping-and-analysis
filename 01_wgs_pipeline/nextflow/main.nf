#!/usr/bin/env nextflow
// =============================================================================
//  main.nf — NGS Mapping and Analysis: WGS Alignment Pipeline (DSL2)
//  Tools: BWA-MEM2 → samtools sort → Picard MarkDuplicates → mosdepth → MultiQC
// =============================================================================
nextflow.enable.dsl=2

// ── Parameters ────────────────────────────────────────────────────────────────
params.reads        = "data/*_{R1,R2}*.fastq.gz"
params.reference    = "reference/genome.fa"
params.outdir       = "results"
params.threads      = 8
params.remove_dups  = false
params.bed          = ""   // Optional WES target BED

// (Moved help message check inside workflow block)

// ── Process: BWA-MEM2 Index ───────────────────────────────────────────────────
process BWA_INDEX {
  tag "index"
  container "quay.io/biocontainers/bwa-mem2:2.2.1--he4a0461_0"
  publishDir "${params.outdir}/reference", mode: 'copy'

  input:
  path reference

  output:
  tuple path(reference), path("${reference}.*"), emit: index

  script:
  """
  bwa-mem2 index ${reference}
  """
}

// ── Process: BWA-MEM2 Align ───────────────────────────────────────────────────
process BWA_MEM2 {
  tag "${sample_id}"
  container "quay.io/biocontainers/mulled-v2-ad9b139c2084ba65780287f7a7596b6e4e082c97:6ab1da9c968f237bf320c9044d030cfd80e14a1e-0"
  cpus params.threads
  memory { 8.GB * task.attempt }
  errorStrategy { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
  maxRetries 2

  input:
  tuple val(sample_id), path(reads)
  tuple path(reference), path(index_files)

  output:
  tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam

  script:
  def rg = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}_lib"
  """
  bwa-mem2 mem \\
    -t ${task.cpus} \\
    -M -Y -K 100000000 \\
    -R "${rg}" \\
    ${reference} \\
    ${reads[0]} ${reads[1]} | \\
  samtools sort \\
    -@ ${task.cpus} \\
    -m 4G \\
    -o ${sample_id}.sorted.bam -
  """
}

// ── Process: Picard MarkDuplicates ────────────────────────────────────────────
process PICARD_MARKDUP {
  tag "${sample_id}"
  container "broadinstitute/picard:2.27.5"
  cpus 4
  memory "2 GB"

  input:
  tuple val(sample_id), path(bam)

  output:
  tuple val(sample_id), path("${sample_id}.markdup.bam"), emit: bam
  path "${sample_id}.markdup.bai",                        emit: bai
  path "${sample_id}.markdup_metrics.txt",                emit: metrics

  script:
  def remove = params.remove_dups ? "true" : "false"
  """
  picard MarkDuplicates \\
    I=${bam} \\
    O=${sample_id}.markdup.bam \\
    M=${sample_id}.markdup_metrics.txt \\
    REMOVE_DUPLICATES=${remove} \\
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \\
    ASSUME_SORTED=true \\
    VALIDATION_STRINGENCY=LENIENT \\
    CREATE_INDEX=true
  """
}

// ── Process: samtools flagstat ─────────────────────────────────────────────────
process FLAGSTAT {
  tag "${sample_id}"
  container "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
  publishDir "${params.outdir}/qc", mode: 'copy'

  input:
  tuple val(sample_id), path(bam)

  output:
  path "${sample_id}.flagstat.txt"

  script:
  """
  samtools flagstat -@ ${task.cpus} ${bam} > ${sample_id}.flagstat.txt
  """
}

// ── Process: mosdepth coverage ────────────────────────────────────────────────
process MOSDEPTH {
  tag "${sample_id}"
  container "quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2"
  cpus 4
  publishDir "${params.outdir}/qc", mode: 'copy'

  input:
  tuple val(sample_id), path(bam)
  path bai

  output:
  path "${sample_id}.*"

  script:
  def bed_arg = params.bed ? "--by ${params.bed}" : "--by 1000"
  """
  mosdepth \\
    -t ${task.cpus} \\
    -n --fast-mode -F 1024 \\
    ${bed_arg} \\
    ${sample_id} \\
    ${bam}
  """
}

// ── Process: MultiQC ──────────────────────────────────────────────────────────
process MULTIQC {
  container "quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"
  publishDir "${params.outdir}/multiqc", mode: 'copy'

  input:
  path '*'

  output:
  path "multiqc_report.html"
  path "multiqc_report_data.zip"

  script:
  """
  multiqc . -n multiqc_report.html -f -z
  """
}

// ── Workflow ───────────────────────────────────────────────────────────────────
workflow {
  // Help message
  if (params.help) {
    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║          WGS Alignment Pipeline (Nextflow DSL2)      ║
    ╠══════════════════════════════════════════════════════╣
    ║  --reads       Input FASTQ glob    [${params.reads}]
    ║  --reference   Reference FASTA     [${params.reference}]
    ║  --outdir      Output directory    [${params.outdir}]
    ║  --threads     CPUs per task       [${params.threads}]
    ║  --remove_dups Remove duplicates   [${params.remove_dups}]
    ║  --bed         WES target BED      [optional]
    ╚══════════════════════════════════════════════════════╝
    """.stripIndent()
    exit 0
  }
  // Input channels
  reads_ch = Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .ifEmpty { error "No reads found matching: ${params.reads}" }

  reference_ch = Channel.fromPath(params.reference, checkIfExists: true)

  // Pipeline
  BWA_INDEX(reference_ch)
  BWA_MEM2(reads_ch, BWA_INDEX.out.index)
  PICARD_MARKDUP(BWA_MEM2.out.bam)

  FLAGSTAT(PICARD_MARKDUP.out.bam)
  MOSDEPTH(PICARD_MARKDUP.out.bam, PICARD_MARKDUP.out.bai)

  // Collect all QC for MultiQC
  qc_files = FLAGSTAT.out
    .mix(MOSDEPTH.out)
    .mix(PICARD_MARKDUP.out.metrics)
    .collect()

  MULTIQC(qc_files)

  // Emit final BAMs
  PICARD_MARKDUP.out.bam
    .map { sample_id, bam -> "${sample_id}: ${bam}" }
    .view { "  [OUTPUT] $it" }
}
