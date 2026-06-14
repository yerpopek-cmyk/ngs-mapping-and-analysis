```yaml
---
name: Bioinformatics Mapping and Alignment
description: Comprehensive skill for high-throughput sequencing read mapping, alignment quality control, duplicate marking, coverage analysis, and library complexity assessment. Supports DNA-seq (WGS/WES), RNA-seq, and specialized alignment workflows using BWA-MEM2, HISAT2, STAR, samtools, Picard, Mosdepth, Qualimap, RSeQC, and Preseq.
author: GPTomics
repository: https://github.com/GPTomics/bioSkills.git
version: 2.0.0
license: MIT
tags:
  - bioinformatics
  - genomics
  - transcriptomics
  - alignment
  - mapping
  - NGS
  - QC
  - BAM
  - SAM
  - FASTQ
  - reference-genome
  - variant-calling
  - expression-analysis
---

# Bioinformatics Read Mapping and Alignment Suite

## 🧬 Overview

This skill provides expert guidance for mapping high-throughput sequencing reads (DNA-seq, RNA-seq, ChIP-seq, ATAC-seq) to reference genomes, performing quality control on alignments, manipulating SAM/BAM files, marking and removing duplicate reads, assessing library complexity, calculating coverage metrics, and generating comprehensive QC reports. The skill integrates best practices from industry standards (GATK best practices, ENCODE consortium guidelines) and is optimized for both small-scale (single sample) and large-scale (population genomics) projects.

## 🎯 When to Use This Skill

**Trigger conditions:**
- You have raw FASTQ files and need to align them to a reference genome
- You need to convert between SAM/BAM/CRAM formats
- You require duplicate marking (PCR/optical) for WGS, WES, or RNA-seq
- You need coverage calculations (depth, breadth, uniformity)
- You require library complexity estimation (saturation curves, unique molecules)
- You need RNA-seq splice-aware alignment and transcriptome quantification
- You want comprehensive alignment QC reports (GC bias, insert size, mapping statistics)
- You work with ChIP-seq or ATAC-seq and need alignment with fragment length analysis
- You need to filter, subset, or manipulate BAM files for downstream analysis (variant calling, differential expression)

**Primary use cases:**
- Whole Genome Sequencing (WGS) alignment at 30X depth
- Whole Exome Sequencing (WES) alignment from capture data
- Bulk RNA-seq alignment and splice junction discovery
- Single-cell RNA-seq alignment (with cell barcode handling)
- ChIP-seq alignment for transcription factor binding analysis
- ATAC-seq alignment for chromatin accessibility
- Metagenomic alignment against multi-species reference databases
- Long-read alignment (PacBio, Oxford Nanopore) using minimap2

## 🔧 Primary Tools and Versions

| Tool | Minimum Version | Purpose | Best For |
|------|----------------|---------|----------|
| **bwa-mem2** | >= 2.2.1 | BWT-based short read aligner | DNA-seq (WGS/WES) with high speed |
| **bwa** (original) | 0.7.17 | Classic BWA algorithm | Legacy workflows, compatibility |
| **minimap2** | >= 2.24 | Long-read and splice-aware aligner | PacBio, ONT, Iso-seq |
| **hisat2** | >= 2.2.1 | Hierarchical indexing for splice-aware | RNA-seq with large genomes |
| **STAR** | >= 2.7.10a | Ultrafast splice-aware aligner | High-sensitivity RNA-seq |
| **samtools** | >= 1.15 | SAM/BAM/CRAM manipulation, sorting, indexing | All workflows |
| **picard** | >= 2.27.0 | MarkDuplicates, CollectMetrics, ReorderSam | Duplicate marking, QC metrics |
| **mosdepth** | >= 0.3.3 | Fast coverage calculation | Quick depth/coverage analysis |
| **qualimap** | >= 2.2.2 | Comprehensive QC reporting | HTML reports for DNA/RNA-seq |
| **RSeQC** | >= 3.0.1 | RNA-seq specific QC metrics | Gene body coverage, read distribution |
| **preseq** | >= 3.1.2 | Library complexity estimation | Saturation curves, unique molecules |
| **bedtools** | >= 2.30.0 | Interval operations, coverage | BED manipulation, intersection |
| **sambamba** | >= 0.8.2 | Parallel duplicate marking | Large datasets (faster than Picard) |
| **fastp** | >= 0.23.0 | FASTQ preprocessing | Adapter trimming, quality filtering |

## 📋 Core Workflows

### Workflow 1: DNA-Seq Alignment (WGS/WES) - Production Pipeline

This is the standard GATK-recommended pipeline for germline variant calling.

```bash
#!/bin/bash
# Input: sample_R1.fastq.gz, sample_R2.fastq.gz
# Output: sample.dedup.bam, sample.dedup.bai, sample.markduplicates.metrics.txt

REFERENCE="reference/hg38.fa"
SAMPLE_NAME="sample"
RG_ID="${SAMPLE_NAME}_flowcell_lane"
RG_PU="flowcell.lane"
RG_LB="${SAMPLE_NAME}_library"
RG_PL="ILLUMINA"

# Step 1: Align with BWA-MEM2 (faster than original BWA)
bwa-mem2 index -p ${REFERENCE%.fa} ${REFERENCE}  # Only once per reference

bwa-mem2 mem -t 32 -M -Y -K 100000000 \
  -R "@RG\tID:${RG_ID}\tSM:${SAMPLE_NAME}\tPU:${RG_PU}\tLB:${RG_LB}\tPL:${RG_PL}" \
  ${REFERENCE%.fa} \
  ${SAMPLE_NAME}_R1.fastq.gz \
  ${SAMPLE_NAME}_R2.fastq.gz | \
  samtools view -bS - | \
  samtools sort -@ 8 -m 4G -o ${SAMPLE_NAME}.sorted.bam -

# Step 2: Mark duplicates (PCR and optical)
picard MarkDuplicates \
  I=${SAMPLE_NAME}.sorted.bam \
  O=${SAMPLE_NAME}.dedup.bam \
  M=${SAMPLE_NAME}.markduplicates.metrics.txt \
  REMOVE_DUPLICATES=false \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=LENIENT \
  CREATE_INDEX=true

# Step 3: Index the deduplicated BAM
samtools index ${SAMPLE_NAME}.dedup.bam

# Step 4: Basic QC metrics (optional but recommended)
samtools flagstat ${SAMPLE_NAME}.dedup.bam > ${SAMPLE_NAME}.flagstat.txt
samtools stats ${SAMPLE_NAME}.dedup.bam | grep -E "^(SN|COV)" > ${SAMPLE_NAME}.stats.txt

# Step 5: Calculate coverage with mosdepth (much faster than bedtools)
mosdepth -t 8 -n --fast-mode -F 1024 \
  -b 1000 \
  ${SAMPLE_NAME}.coverage \
  ${SAMPLE_NAME}.dedup.bam

# Step 6: Extract coverage statistics (mean depth, % > 20X)
awk -F'\t' 'NR>1 {print $4}' ${SAMPLE_NAME}.coverage.regions.bed.gz | \
  awk '{sum+=$1; count++} END {print "Mean coverage:", sum/count}'
```

**Rationale for parameters:**
- `-M`: Marks split hits as secondary (required for Picard MarkDuplicates)
- `-Y`: Uses soft-clipping for supplementary alignments (preserves sequence for SV detection)
- `-K 100000000`: Batch size (100M) for improved performance on large genomes
- `-t 32`: Number of threads for alignment (scale to available cores)
- `OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500`: Default for Illumina (adjust to sequencer type)

### Workflow 2: RNA-Seq Alignment with Splice-Aware Tools

Two complementary approaches: **HISAT2** (faster, lower memory) vs **STAR** (more sensitive, higher memory).

**Option A: HISAT2 (Recommended for most RNA-seq projects)**

```bash
#!/bin/bash
# Build HISAT2 index (once per genome + transcriptome)
hisat2-build -p 8 \
  --ss reference/splicesites.txt \
  --exon reference/exons.txt \
  reference/hg38.fa \
  reference/hg38_transcriptome

# Align paired-end RNA-seq reads
hisat2 -p 16 \
  -x reference/hg38_transcriptome \
  -1 ${SAMPLE}_R1.fastq.gz \
  -2 ${SAMPLE}_R2.fastq.gz \
  --known-splicesite-infile reference/splicesites.txt \
  --novel-splicesite-outfile ${SAMPLE}_novel_splicesites.txt \
  --rna-strandness RF \  # RF = stranded library (dUTP method)
  --summary-file ${SAMPLE}.hisat2_summary.txt \
  2> ${SAMPLE}.hisat2.log | \
  samtools view -bS - | \
  samtools sort -@ 8 -m 4G -o ${SAMPLE}.hisat2.sorted.bam -

# Mark duplicates (for RNA-seq, marking but not removing recommended)
picard MarkDuplicates \
  I=${SAMPLE}.hisat2.sorted.bam \
  O=${SAMPLE}.hisat2.marked.bam \
  M=${SAMPLE}.hisat2.duplicate_metrics.txt \
  REMOVE_DUPLICATES=false

# Remove duplicates physically for RSeQC compatibility (required!)
samtools view -h -b -F 1024 ${SAMPLE}.hisat2.marked.bam > ${SAMPLE}.hisat2.dedup.bam
samtools index ${SAMPLE}.hisat2.dedup.bam
```

**Option B: STAR (Higher sensitivity, more accurate for novel junctions)**

```bash
#!/bin/bash
# Build STAR index (requires substantial memory: ~30GB for human)
STAR --runMode genomeGenerate \
  --genomeDir reference/STAR_hg38 \
  --genomeFastaFiles reference/hg38.fa \
  --sjdbGTFfile reference/genes.gtf \
  --sjdbOverhang 100 \
  --runThreadN 16

# Align reads with STAR
STAR --genomeDir reference/STAR_hg38 \
  --readFilesIn ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz \
  --readFilesCommand zcat \
  --runThreadN 32 \
  --outFileNamePrefix ${SAMPLE}_STAR_ \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMcompression 6 \
  --outSAMattrRGline "ID:${SAMPLE} SM:${SAMPLE}" \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverLmax 0.04 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --chimSegmentMin 15 \
  --chimOutType WithinBAM \
  --quantMode GeneCounts

# STAR outputs: ${SAMPLE}_STAR_Aligned.sortedByCoord.out.bam
# Convert to deduped BAM for RSeQC
samtools view -h -b -F 1024 ${SAMPLE}_STAR_Aligned.sortedByCoord.out.bam > ${SAMPLE}_STAR.dedup.bam
samtools index ${SAMPLE}_STAR.dedup.bam
```

**Strandedness parameter guide:**
- `RF`: Read 1 (first in pair) matches transcript's reverse strand (dUTP, Illumina stranded)
- `FR`: Read 1 matches transcript's forward strand (standard stranded)
- `--fr/--rf` for HISAT2; `--outSAMstrandField intronMotif` for STAR

### Workflow 3: RNA-Seq Quality Control with RSeQC

RSeQC provides detailed RNA-seq metrics that are missing from standard BAM QC tools.

```bash
#!/bin/bash
# Prerequisites: Convert GTF to BED format (RSeQC requires BED)
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons annotation.gtf annotation.genePred
genePredToBed annotation.genePred annotation.bed

# 1. Read distribution (what regions do reads map to?)
read_distribution.py \
  -r annotation.bed \
  -i ${SAMPLE}.dedup.bam \
  > ${SAMPLE}.read_distribution.txt

# Output includes: CDS_Exons, 5'UTR_Exons, 3'UTR_Exons, Introns, Intergenic

# 2. Gene body coverage (5' to 3' bias)
geneBody_coverage.py \
  -r annotation.bed \
  -i ${SAMPLE}.dedup.bam \
  -o ${SAMPLE}_genebody

# 3. Insert size distribution (fragment length)
collect_insert_size_metrics.py \
  -i ${SAMPLE}.dedup.bam \
  -o ${SAMPLE}_insert_size.txt

# 4. Junction saturation (how many known vs novel junctions detected)
junction_saturation.py \
  -r annotation.bed \
  -i ${SAMPLE}.dedup.bam \
  -o ${SAMPLE}_junction

# 5. Inner distance (fragment length between pairs)
inner_distance.py \
  -r annotation.bed \
  -i ${SAMPLE}.dedup.bam \
  -o ${SAMPLE}_inner_distance
```

### Workflow 4: Comprehensive Qualimap Reports

Qualimap generates interactive HTML reports with coverage heatmaps, GC bias plots, and insert size distributions.

```bash
#!/bin/bash
# DNA-seq (WGS/WES) QC report
qualimap bamqc \
  -bam ${SAMPLE}.dedup.bam \
  -gff annotation.gff \
  -outdir qualimap_${SAMPLE}_bamqc \
  --java-mem-size=8G \
  -c \
  -nt 8

# RNA-seq QC report (requires GTF, not GFF)
qualimap rnaseq \
  -bam ${SAMPLE}.hisat2.marked.bam \
  -gtf annotation.gtf \
  -outdir qualimap_${SAMPLE}_rnaseq \
  --java-mem-size=8G \
  -p strand-specific-reverse \  # Adjust strandedness
  -a quant \
  -pe

# Optional: ChIP-seq QC report
qualimap chipseq \
  -bam ${SAMPLE}.dedup.bam \
  -gff peaks.bed \
  -outdir qualimap_${SAMPLE}_chipseq \
  --java-mem-size=8G
```

### Workflow 5: Library Complexity and Saturation Analysis

Use preseq to extrapolate library complexity and determine if deeper sequencing would yield more unique molecules.

```bash
#!/bin/bash
# For paired-end libraries
preseq lc_extrap \
  -pe \
  -seg_len 100000 \
  -B ${SAMPLE}.sorted.bam \
  -o ${SAMPLE}_complexity_curve.txt \
  -v \
  -extrap 5e9

# Generate complexity metrics summary
preseq c_curve \
  -B ${SAMPLE}.sorted.bam \
  -o ${SAMPLE}_complexity_metrics.txt

# Plot using R (optional)
Rscript -e '
data <- read.table("sample_complexity_curve.txt", header=TRUE)
pdf("complexity_plot.pdf")
plot(data$TOTAL_READS, data$EXPECTED_DISTINCT, type="l", 
     xlab="Total Reads Sequenced", ylab="Expected Distinct Reads",
     main="Library Complexity Curve")
dev.off()
'

# Expected interpretation:
# - Steep curve (linear increase) → Complex library, more sequencing beneficial
# - Plateauing curve → Library complexity exhausted, deeper sequencing unnecessary
```

### Workflow 6: ChIP-seq Alignment and Fragment Length Analysis

ChIP-seq requires additional steps for fragment length estimation and cross-correlation analysis.

```bash
#!/bin/bash
# Align with BWA-MEM (same as DNA-seq)
bwa-mem2 mem -t 16 -M -Y \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
  reference/hg38 \
  ${SAMPLE}_R1.fastq.gz \
  ${SAMPLE}_R2.fastq.gz | \
  samtools sort -@ 8 -o ${SAMPLE}.sorted.bam -

# Mark duplicates (recommended for ChIP-seq)
picard MarkDuplicates \
  I=${SAMPLE}.sorted.bam \
  O=${SAMPLE}.dedup.bam \
  M=${SAMPLE}.duplicates.txt \
  REMOVE_DUPLICATES=true \  # Remove duplicates for ChIP-seq
  ASSUME_SORTED=true

# Filter to keep only high-quality properly paired reads (MAPQ >= 30)
samtools view -h -q 30 -f 0x2 -b ${SAMPLE}.dedup.bam > ${SAMPLE}.filtered.bam
samtools index ${SAMPLE}.filtered.bam

# Cross-correlation analysis for fragment length (ENCODE quality metric)
# Requires hotspot or SPP package
Rscript compute_spp.R -i=${SAMPLE}.filtered.bam -savp=${SAMPLE}.spp.pdf \
  -frag=${SAMPLE}.fragment_length.txt -peak=peaks.bed

# Calculate NSC (Normalized Strand Cross-correlation) and RSC (Relative Strand Cross-correlation)
# NSC > 1.05, RSC > 0.8 indicates good ChIP signal
```

### Workflow 7: ATAC-seq Alignment and Fragment Size Filtering

ATAC-seq requires special handling for mitochondrial reads and nucleosome-free fragments.

```bash
#!/bin/bash
# ATAC-seq alignment (default BWA-MEM2)
bwa-mem2 mem -t 16 -M -Y \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
  reference/hg38 \
  ${SAMPLE}_R1.fastq.gz \
  ${SAMPLE}_R2.fastq.gz | \
  samtools view -bS -f 0x2 -q 30 | \
  samtools sort -@ 8 -o ${SAMPLE}.bam -

# Remove duplicates
picard MarkDuplicates \
  I=${SAMPLE}.bam \
  O=${SAMPLE}.dedup.bam \
  M=${SAMPLE}.duplicates.txt \
  REMOVE_DUPLICATES=true

# Remove mitochondrial reads (optional, unless studying mito)
samtools view -h ${SAMPLE}.dedup.bam | \
  awk '$3 != "chrM" && $3 != "MT"' | \
  samtools view -bS - > ${SAMPLE}.nuclear.bam

# Filter by fragment length (nucleosome-free: <120bp, mononucleosome: 180-247bp)
samtools view -h ${SAMPLE}.nuclear.bam | \
  awk 'function abs(x){return ((x < 0.0) ? -x : x)} 
       $9 != 0 { fraglen = abs($9); 
       if (fraglen <= 120 || (fraglen >= 180 && fraglen <= 247)) print $0 }' | \
  samtools view -bS - > ${SAMPLE}.fragment_filtered.bam

# Generate fragment length histogram for QC
samtools view ${SAMPLE}.nuclear.bam | \
  awk '$9 != 0 {print abs($9)}' | \
  sort -n | \
  uniq -c > ${SAMPLE}.fragment_lengths.txt
```

### Workflow 8: Long-Read Alignment (PacBio / Oxford Nanopore)

For Iso-seq, cDNA-seq, or genomic long-reads using minimap2.

```bash
#!/bin/bash
# Align PacBio HiFi reads (accuracy >99%)
minimap2 -ax map-hifi \
  -t 32 \
  reference/hg38.fa \
  ${SAMPLE}.fastq.gz | \
  samtools view -bS - | \
  samtools sort -@ 8 -o ${SAMPLE}.hifi.bam -

# Align ONT reads (lower accuracy, more permissive)
minimap2 -ax map-ont \
  -t 32 \
  -L \  # Use long-read splice mode for RNA
  reference/hg38.fa \
  ${SAMPLE}.fastq.gz | \
  samtools view -bS - | \
  samtools sort -@ 8 -o ${SAMPLE}.ont.bam -

# Index and flagstat
samtools index ${SAMPLE}.hifi.bam
samtools flagstat ${SAMPLE}.hifi.bam > ${SAMPLE}.hifi.flagstat.txt

# Convert to CRAM for space savings (10x compression vs BAM)
samtools view -C -T reference/hg38.fa -o ${SAMPLE}.hifi.cram ${SAMPLE}.hifi.bam
```

## 🔬 Advanced Techniques and Best Practices

### 1. Handling High-Copy Number Regions and Repetitive DNA

For genomes with high repeat content (plants, some fungi), use alternative aligners or filtering:

```bash
# BWA-MEM with lower mismatch penalty (-B) and higher gap penalty (-E)
bwa-mem2 mem -B 4 -E 2 -O 6,6 -L 5,5 -t 32 ... 

# Filter multimapping reads (NH:i tag > 1)
samtools view -h -b -e 'NH:i==1' ${SAMPLE}.bam > ${SAMPLE}.unique.bam
```

### 2. Parallel Processing for Large Cohorts (100+ samples)

```bash
#!/bin/bash
# Using GNU parallel for bulk processing
export -f align_sample
parallel -j 10 align_sample {} ::: sample_names.txt

# Or using Snakemake workflow management
snakefile -s alignment.smk --cores 32 --use-conda
```

### 3. Memory-Efficient Processing of Large BAMs

```bash
# Use sambamba (Rust-based, 2x faster than samtools)
sambamba view -h -f bam -t 8 -F "mapping_quality >= 30" input.bam > filtered.bam

# Use CRAM format (reference-based compression)
samtools view -C -T reference.fa -o output.cram input.bam
```

### 4. Quality Score Calibration for RNA-seq

Some aligners need quality score adjustments for RNA-seq data with high mismatch rates:

```bash
# STAR option for higher mismatch tolerance (e.g., for poor quality libraries)
STAR --outFilterMismatchNoverLmax 0.1 --outFilterMultimapNmax 100 ...
```

## 📊 Quality Control Metrics and Interpretation

### Key Metrics to Monitor:

| Metric | Good Range | Action if Out of Range |
|--------|-----------|------------------------|
| **Alignment rate** | >90% for DNA-seq, >70% for RNA-seq | Check adapter contamination, reference genome version |
| **Properly paired rate** | >85% | Check fragment size distribution, library prep |
| **Duplication rate** | <20% for WGS, <40% for WES, <30% for RNA-seq | Over-amplification during PCR, consider deeper sequencing if complex |
| **Mean coverage depth** | 30-40X for WGS, 100-150X for WES, 30-50M reads for RNA-seq | Insufficient sequencing depth, need more lanes |
| **Coverage uniformity** (fold 80) | <2.0 for WGS, <3.0 for WES | GC bias, capture inefficiency |
| **GC bias** | Pearson correlation <0.1 | Library prep issues, PCR bias |
| **Insert size** (mean) | 300-500bp for short-read | Fragmentation or size selection issues |
| **rRNA contamination** (RNA-seq) | <5% | Incomplete rRNA depletion |
| **Chimerics** (RNA-seq) | <5% | Incomplete cDNA synthesis or alignment artifacts |
| **Fraction of reads in targets** (WES) | >60% | Off-target capture |

### Alignment Summary with MultiQC

Aggregate metrics across multiple samples using MultiQC:

```bash
# Collect all QC outputs into one report
multiqc . -n multiqc_report.html -f -z
```

## 🚨 Common Issues and Troubleshooting

### Issue 1: Low alignment rate (<70%)
**Solution:**
```bash
# Check FASTQ quality and contamination
fastqc ${SAMPLE}_R1.fastq.gz -o fastqc_reports/
multiqc fastqc_reports/

# Trim adapters and low-quality bases
fastp -i ${SAMPLE}_R1.fastq.gz -I ${SAMPLE}_R2.fastq.gz \
  -o trimmed_R1.fastq.gz -O trimmed_R2.fastq.gz \
  -q 20 -u 30 -n 5 -l 36

# Use more permissive alignment parameters
bwa-mem2 mem -B 4 -O 6,6 -L 5,5 -t 32 ...
```

### Issue 2: High duplication rate (>50%)
**Solution:**
```bash
# Check if due to low complexity library
preseq lc_extrap -pe -B input.bam -o complexity.txt

# If complexity is low, re-sequence deeper
# If complexity is adequate, use stricter duplicate removal
picard MarkDuplicates REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=100

# For RNA-seq, consider UMIs (unique molecular identifiers)
umi_tools dedup -I input.bam -S dedup.bam --method unique
```

### Issue 3: Coverage is uneven (fold80 > 2.5)
**Solution:**
```bash
# Compute GC bias
qualimap bamqc -bam dedup.bam --java-mem-size=8G

# Apply GC correction (for CNV analysis)
# Use CNVkit or other GC-aware tools
```

### Issue 4: RNA-seq has high intronic mapping (>30%)
**Solution:**
```bash
# DNA contamination? Run FastQ-Screen
fastq_screen --conf fastq_screen.conf --aligner bowtie2 --top 100000 --threads 8 sample.fastq

# Or use HISAT2 with stricter parameters
hisat2 --max-intronlen 100000 --min-intronlen 20 --no-temp-splicesite ...
```

## 🎓 Learning Resources and Reference

- **GATK Best Practices**: https://gatk.broadinstitute.org/hc/en-us/articles/360035531912
- **ENCODE Consortium Guidelines**: https://www.encodeproject.org/about/experiment-guidelines/
- **SAM Format Specification**: https://samtools.github.io/hts-specs/SAMv1.pdf
- **RSeQC Documentation**: http://rseqc.sourceforge.net/
- **Qualimap Tutorial**: http://qualimap.conesalab.org/tutorial.shtml

## 📝 Example Prompts for Claude

Use these prompts to trigger this skill effectively:

> "Align my paired-end WGS reads (FASTQ) to the human genome (hg38) using BWA-MEM2, mark duplicates with Picard, and calculate coverage with mosdepth. Generate a Qualimap QC report."

> "I have RNA-seq data from a stranded dUTP library. Use HISAT2 to align to the mouse genome (mm10), mark duplicates (but don't remove), then run RSeQC's geneBody_coverage and read_distribution. Also, convert the BED for annotation."

> "Perform ChIP-seq alignment for transcription factor H3K27ac using BWA-MEM, remove duplicates (REMOVE_DUPLICATES=true), filter for MAPQ >= 30, and compute fragment length distribution and cross-correlation metrics. Assess library complexity with preseq."

> "I have Oxford Nanopore long reads for a bacterial genome (E. coli). Use minimap2 with map-ont preset, convert output to sorted BAM, flagstat, and optionally convert to CRAM for space efficiency."

> "Run a complete RNA-seq QC suite on my STAR-aligned BAM files: gene body coverage, read distribution, insert size metrics, junction saturation, and generate an HTML report with Qualimap rnaseq."

## 🧪 Validation and Testing

To validate your alignment pipeline:

```bash
# Test with a subsampled dataset (100k reads)
seqtk sample -s42 ${SAMPLE}_R1.fastq.gz 100000 > test_R1.fastq
seqtk sample -s42 ${SAMPLE}_R2.fastq.gz 100000 > test_R2.fastq

# Run full alignment pipeline
bash align_pipeline.sh test

# Compare with expected metrics
# For human: expected alignment rate ~98% (DNA), ~85% (RNA)
```

## 📌 Version Control and Reproducibility

Document all tool versions for reproducibility:

```bash
# Capture environment
conda list --export > conda_environment.txt
pip freeze > pip_requirements.txt

# Record git commit of pipeline
git rev-parse HEAD > pipeline_version.txt

# Docker/Singularity for complete reproducibility
docker pull quay.io/biocontainers/bwa-mem2:2.2.1--he4a0461_0
```

---

**Skill Version:** 2.0.0  
**Last Updated:** 2024  
**Maintainer:** GPTomics  
**License:** MIT  
**Contributions Welcome:** Submit PRs to https://github.com/GPTomics/bioSkills.git
```