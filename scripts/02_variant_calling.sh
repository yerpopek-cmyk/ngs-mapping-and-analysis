#!/bin/bash

# ==============================================================================
# Скрипт 02: Картирование ридов, Variant Calling и Интерпретация
# ------------------------------------------------------------------------------
# От FASTQ файлов к аннотированному VCF. Включает контроль качества (MAPQ), 
# нормализацию инделов и фильтрацию ложноположительных артефактов.
# ==============================================================================

set -e

REF="hg38.fa" # Замените на путь к вашему референсу
READS1="sample_R1.fastq.gz"
READS2="sample_R2.fastq.gz"
SAMPLE="S1"
THREADS=10

echo "🗺️ Шаг 1: Индексация референса и картирование (BWA-MEM)"
samtools faidx $REF
# BWA-MEM генерирует выравнивание. Передаем в samtools для сортировки.
# MAPQ = -10 * log10(P(ошибка картирования))
bwa mem -t $THREADS $REF $READS1 $READS2 | samtools sort -@ $THREADS -o ${SAMPLE}.bam
samtools index -@ $THREADS ${SAMPLE}.bam

# Опционально: Маркировка ПЦР-дубликатов с помощью Picard
# picard MarkDuplicates -I ${SAMPLE}.bam -O ${SAMPLE}_markdup.bam -M metrics.txt
# samtools index ${SAMPLE}_markdup.bam

echo "🔍 Шаг 2: Вызов вариантов (Variant Calling) с помощью FreeBayes"
# FreeBayes использует байесовские гаплотипные модели. 
# QUAL вычисляется как вероятность того, что локус полиморфен.
freebayes -f $REF -b ${SAMPLE}.bam \
  --min-alternate-fraction 0.2 \
  --min-base-quality 20 \
  > ${SAMPLE}.raw.vcf

echo "🧹 Шаг 3: Нормализация VCF (bcftools norm)"
# Декомпозиция мультиаллельных вариантов (-m -any) и 
# Left-alignment (сдвиг индела максимально влево для унификации)
bcftools norm -f $REF -m -any ${SAMPLE}.raw.vcf -O z -o ${SAMPLE}.norm.vcf.gz
tabix -p vcf ${SAMPLE}.norm.vcf.gz

echo "🛡️ Шаг 4: Жесткая фильтрация (Hard Filtering)"
# Фильтруем по глубине (DP) и качеству (QUAL).
# Требуем наличие чтений для реф (RO) и альт (AO) аллелей.
bcftools filter -i 'QUAL>30 && INFO/DP>10 && FMT/AO>0 && FMT/RO>0' \
  ${SAMPLE}.norm.vcf.gz -O z -o ${SAMPLE}.filt.vcf.gz
tabix -p vcf ${SAMPLE}.filt.vcf.gz

echo "📈 Шаг 5: Быстрая QC статистика"
# Соотношение Ti/Tv для человека ожидается ~2.0 - 2.1
bcftools stats ${SAMPLE}.filt.vcf.gz > ${SAMPLE}.stats.txt
grep "number of SNPs:" ${SAMPLE}.stats.txt
grep "TSTV" ${SAMPLE}.stats.txt

echo "🏷️ Шаг 6: Функциональная аннотация"
# Оценка эффекта мутаций (VEP-подобная аналитика). 
# Узнаем является ли мутация missense, frameshift, stop_gained и т.д.
# Для этого метода требуется файл аннотации генов: genes.gff3.gz
if [ -f "genes.gff3.gz" ]; then
    bcftools csq -f $REF -g genes.gff3.gz ${SAMPLE}.filt.vcf.gz -O v -o ${SAMPLE}.annotated.vcf
else
    echo "Файл genes.gff3.gz не найден, пропускаем bcftools csq. Воспользуйтесь snpEff:"
    echo "java -Xmx4g -jar snpEff.jar GRCh38.99 ${SAMPLE}.filt.vcf.gz > ${SAMPLE}.annotated.vcf"
fi

echo "✅ Пайплайн Variant Calling завершен!"
