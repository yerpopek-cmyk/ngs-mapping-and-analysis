# 🎨 Project 3: Coverage Panorama

> Transform raw sequencing depth data into **publication-ready genome-wide coverage art**.

## What it does

`coverage_panorama.py` reads [mosdepth](https://github.com/brentp/mosdepth) regions BED files
and renders each chromosome as a colour-coded horizontal strip — like a genome karyotype painted
with coverage depth. The colour scale runs from white (no coverage) through blue (normal depth)
to red (excess / amplification).

`karyotype_painter.py` generates realistic synthetic BED data for demos and testing — no real
sequencing data needed to explore the tool.

## Output examples

| Mode | Description |
|------|-------------|
| Single sample | Karyotype panorama with per-chromosome mean depth |
| Compare mode | Side-by-side comparison of multiple samples |
| CNV simulation | Gain/loss events visible as colour shifts |

## Quick start

```bash
pip install -r requirements.txt

# 1. Generate synthetic demo data (no BAM needed)
python karyotype_painter.py --mode all

# 2. Render single sample
python coverage_panorama.py \
  -s wgs_normal.bed.gz \
  -n "Normal WGS 30X" \
  -o normal_panorama.png

# 3. Multi-sample comparison (compare CNV vs normal)
python coverage_panorama.py \
  -s wgs_normal.bed.gz wgs_cnv.bed.gz wgs_tumour.bed.gz \
  -n "Normal" "CNV" "Tumour 65%" \
  --mode compare \
  -o comparison.png
```

## Using with real data

```bash
# 1. Run mosdepth on your BAM
mosdepth -t 8 -n --fast-mode -b 1000 sample sample.markdup.bam

# 2. Render the panorama
python coverage_panorama.py \
  -s sample.regions.bed.gz \
  -n "Sample ID — WGS 35X" \
  -o sample_panorama.png --dpi 300
```

## Parameters

| Flag | Description | Default |
|------|-------------|---------|
| `--sample / -s` | mosdepth BED file(s) | required |
| `--name / -n`   | Sample label(s)       | filename |
| `--output / -o` | Output file (.png/.svg/.pdf) | coverage_panorama.png |
| `--mode`        | `single` or `compare` | single |
| `--dpi`         | Output resolution     | 180 |
| `--no-summary`  | Skip terminal table   | off |

## Simulation modes (`karyotype_painter.py`)

| Mode | Description |
|------|-------------|
| `wgs_normal`   | Uniform 30X coverage with mild GC bias |
| `cnv`          | chr8 gain (1.6×), chr17 loss (0.5×) |
| `tumour`       | Multi-region CNV at 65% purity |
| `rna_uneven`   | Gene-density-biased RNA-seq coverage |
| `all`          | Generate all four datasets |
