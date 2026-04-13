``
🧬 NGS Mapping & Bioinformatics Analysis Pipeline

Добро пожаловать в репозиторий NGS Mapping & Analysis. Этот проект объединяет ключевые этапы биоинформатического анализа данных секвенирования следующего поколения (NGS): от сырых прочтений до аннотации геномов, поиска мутаций (Variant Calling) и филогенетических реконструкций.

📂 Структура репозитория

В репозитории файлы распределены по их назначению для обеспечения порядка и воспроизводимости:

docs/: Содержит теоретическую базу. В файле theory_and_formulas.md собраны все ключевые математические формулы (E-value, QUAL, алгоритм Neighbor-Joining, скоринг Prodigal) и принципы работы инструментов.

scripts/: Директория с исполняемыми bash-скриптами. Каждый скрипт снабжен подробными комментариями.

01_annotation.sh — пайплайн структурной и функциональной аннотации генома.

02_variant_calling.sh — пайплайн картирования (mapping), вызова вариантов и их фильтрации.

03_phylogenetics.sh — пайплайн множественного выравнивания для филогенетики.

LICENSE: Правила использования кода (MIT License).

🛠️ Используемые инструменты (Tools Stack)

Этап анализа

Инструменты

Назначение

Mapping & QC

BWA-MEM, Samtools, Picard

Картирование ридов на референс, сортировка, маркировка дубликатов

Variant Calling

FreeBayes, bcftools, vt

Поиск SNV и Indel, байесовский вывод, нормализация VCF

Annotation (Genomic)

Prokka, Prodigal

Поиск ORF, РНК, структурная аннотация бактерий/вирусов

Annotation (Variants)

snpEff, bcftools csq

Оценка функционального эффекта мутаций (VEP-аналоги)

Phylogenetics

MAFFT, MEGA

Множественное выравнивание (FFT-алгоритм) и построение NJ деревьев

Workflow

Nextflow

Развертывание конвейеров (опционально для функциональной аннотации)

🚀 Как использовать (Quick Start)

Клонируйте репозиторий:

git clone [https://github.com/your-username/ngs-mapping-and-analysis.git](https://github.com/your-username/ngs-mapping-and-analysis.git)
cd ngs-mapping-and-analysis


Убедитесь, что у вас установлены необходимые инструменты (рекомендуется использовать conda):

conda create -n ngs_env -c bioconda bwa samtools freebayes bcftools prokka mafft snpeff
conda activate ngs_env


Перейдите в папку со скриптами и запустите нужный пайплайн:

cd scripts
bash 02_variant_calling.sh


🧠 Логика пайплайна (Evidence Ladder)

Наш анализ строится на концептуальной "лестнице доказательств":
Reads → QC → Assembly → Annotation → Downstream (Variants / Pathways / Comparative)

(Подробную математическую базу и теорию смотрите в docs/theory_and_formulas.md)
``
