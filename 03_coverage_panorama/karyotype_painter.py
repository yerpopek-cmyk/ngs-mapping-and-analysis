#!/usr/bin/env python3
"""
karyotype_painter.py — Synthetic Coverage Data Generator
=========================================================
Generates realistic synthetic mosdepth-style BED files for testing
and demonstrating coverage_panorama.py without real sequencing data.

Simulates:
  • Normal WGS coverage (Gaussian around 30X)
  • WGS with focal amplification (chromosome arm gain)
  • WGS with chromosome loss (deletion / LOH)
  • RNA-seq non-uniform coverage (gene-density biased)
  • Tumour heterogeneity (mixed copy-number states)

Usage:
  python karyotype_painter.py --mode wgs_normal   --out normal.bed.gz
  python karyotype_painter.py --mode cnv           --out cnv.bed.gz
  python karyotype_painter.py --mode rna_uneven    --out rna.bed.gz
  python karyotype_painter.py --mode tumour        --out tumour.bed.gz

  # Then visualise
  python coverage_panorama.py -s normal.bed.gz cnv.bed.gz \\
    -n "Normal WGS" "CNV sample" --mode compare -o comparison.png
"""

import argparse
import gzip
import random
import math
from pathlib import Path

import numpy as np

# hg38 chromosome sizes in base pairs
HG38_SIZES_BP = {
    "chr1":  249_250_621, "chr2":  243_199_373, "chr3":  198_022_430,
    "chr4":  191_154_276, "chr5":  180_915_260, "chr6":  171_115_067,
    "chr7":  159_138_663, "chr8":  146_364_022, "chr9":  141_213_431,
    "chr10": 135_534_747, "chr11": 135_006_516, "chr12": 133_851_895,
    "chr13": 115_169_878, "chr14": 107_349_540, "chr15": 102_531_392,
    "chr16":  90_354_753, "chr17":  81_195_210, "chr18":  78_077_248,
    "chr19":  59_128_983, "chr20":  63_025_520, "chr21":  48_129_895,
    "chr22":  51_304_566, "chrX":  155_270_560, "chrY":   59_373_566,
}

BIN_SIZE = 1_000  # 1 kb bins (matches mosdepth default)


import pandas as pd

import subprocess

def write_bed_gz(df: pd.DataFrame, out_path: str):
    """Write dataframe to a gzipped BED file."""
    df.to_csv(out_path, sep="\t", header=False, index=False, float_format="%.2f", compression="gzip")
    print(f"  [OK] Written: {out_path}  ({len(df):,} bins)")


def simulate_wgs_normal(
    mean_depth: float = 30.0,
    dispersion: float = 0.08,
    pct_low_complex: float = 0.03,
) -> pd.DataFrame:
    """Simulate uniform WGS coverage with mild GC bias."""
    rng = np.random.default_rng(42)
    dfs = []
    for chrom, size in HG38_SIZES_BP.items():
        num_bins = math.ceil(size / BIN_SIZE)
        starts = np.arange(0, size, BIN_SIZE)
        ends = np.minimum(starts + BIN_SIZE, size)
        
        gc_pos = starts / size
        gc_factors = 1.0 + 0.05 * np.sin(gc_pos * np.pi * 8)
        
        depths = rng.negative_binomial(
            n=mean_depth / dispersion,
            p=1 / (1 + dispersion),
            size=num_bins
        ) * gc_factors
        
        dropouts = rng.random(size=num_bins) < pct_low_complex
        dropout_mult = rng.uniform(0.0, 0.2, size=num_bins)
        depths = np.where(dropouts, depths * dropout_mult, depths)
        depths = np.maximum(depths, 0)
        
        dfs.append(pd.DataFrame({
            "chrom": chrom,
            "start": starts,
            "end": ends,
            "depth": depths
        }))
    return pd.concat(dfs, ignore_index=True)


def simulate_cnv(
    mean_depth: float = 30.0,
    gain_chrom: str = "chr8",
    loss_chrom: str = "chr17",
    gain_factor: float = 1.6,
    loss_factor: float = 0.5,
) -> pd.DataFrame:
    """Simulate WGS with one gained and one lost chromosome."""
    rng = np.random.default_rng(99)
    dfs = []
    for chrom, size in HG38_SIZES_BP.items():
        factor = 1.0
        if chrom == gain_chrom:
            factor = gain_factor
        elif chrom == loss_chrom:
            factor = loss_factor

        num_bins = math.ceil(size / BIN_SIZE)
        starts = np.arange(0, size, BIN_SIZE)
        ends = np.minimum(starts + BIN_SIZE, size)
        
        depths = rng.poisson(mean_depth * factor, size=num_bins) * rng.uniform(0.92, 1.08, size=num_bins)
        depths = np.maximum(depths, 0)
        
        dfs.append(pd.DataFrame({
            "chrom": chrom,
            "start": starts,
            "end": ends,
            "depth": depths
        }))
    return pd.concat(dfs, ignore_index=True)


def simulate_tumour_heterogeneity(
    mean_depth: float = 60.0,
    purity: float = 0.65,
) -> pd.DataFrame:
    """Simulate a tumour sample with clonal copy-number aberrations."""
    rng = np.random.default_rng(7)

    cnv_regions = [
        ("chr3",  0.0,  0.5,  3),   # chr3p gain
        ("chr8",  0.0,  1.0,  4),   # chr8 amplification
        ("chr13", 0.0,  1.0,  1),   # chr13 deletion
        ("chr17", 0.4,  0.8,  5),   # focal amplification
        ("chr22", 0.0,  0.3,  1),   # chr22q loss
    ]

    dfs = []
    for chrom, size in HG38_SIZES_BP.items():
        num_bins = math.ceil(size / BIN_SIZE)
        starts = np.arange(0, size, BIN_SIZE)
        ends = np.minimum(starts + BIN_SIZE, size)
        
        fracs = starts / size
        cn = np.full(num_bins, 2)
        for c, sf, ef, cn_val in cnv_regions:
            if chrom == c:
                mask = (fracs >= sf) & (fracs <= ef)
                cn[mask] = cn_val
                
        effective_cn = purity * cn + (1 - purity) * 2
        n_param = np.maximum(mean_depth * effective_cn / 2, 1)
        depths = rng.negative_binomial(n=n_param, p=0.5)
        depths = np.maximum(depths, 0)
        
        dfs.append(pd.DataFrame({
            "chrom": chrom,
            "start": starts,
            "end": ends,
            "depth": depths
        }))
    return pd.concat(dfs, ignore_index=True)


def simulate_rna_coverage(
    mean_depth: float = 25.0,
    gene_density_bias: float = 0.4,
) -> pd.DataFrame:
    """Simulate RNA-seq coverage — higher in gene-dense regions."""
    rng = np.random.default_rng(13)

    dense_regions = {
        "chr19": [(0.1, 0.9)],
        "chr17": [(0.3, 0.7)],
        "chr11": [(0.5, 0.7)],
        "chr6":  [(0.2, 0.5)],
    }

    dfs = []
    for chrom, size in HG38_SIZES_BP.items():
        num_bins = math.ceil(size / BIN_SIZE)
        starts = np.arange(0, size, BIN_SIZE)
        ends = np.minimum(starts + BIN_SIZE, size)
        
        fracs = starts / size
        factors = np.ones(num_bins)
        for start_f, end_f in dense_regions.get(chrom, []):
            mask = (fracs >= start_f) & (fracs <= end_f)
            factors[mask] += gene_density_bias
            
        n_param = np.maximum(mean_depth * factors, 1).astype(int)
        depths = rng.negative_binomial(n=n_param, p=0.5) * rng.uniform(0.8, 1.2, size=num_bins)
        
        zero_mask = rng.random(size=num_bins) < 0.15
        depths[zero_mask] = 0
        depths = np.maximum(depths, 0)
        
        dfs.append(pd.DataFrame({
            "chrom": chrom,
            "start": starts,
            "end": ends,
            "depth": depths
        }))
    return pd.concat(dfs, ignore_index=True)


def main():
    parser = argparse.ArgumentParser(
        description="Synthetic coverage data generator for Coverage Panorama demos",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--mode",
        choices=["wgs_normal", "cnv", "tumour", "rna_uneven", "all"],
        default="all",
        help="Simulation mode [default: all]",
    )
    parser.add_argument(
        "--out", "-o", default=None,
        help="Output BED.gz path (auto-named if not set)",
    )
    parser.add_argument(
        "--depth", "-d", type=float, default=30.0,
        help="Target mean coverage depth [default: 30]",
    )
    args = parser.parse_args()

    modes = {
        "wgs_normal": ("wgs_normal.bed.gz",    lambda: simulate_wgs_normal(args.depth)),
        "cnv":        ("wgs_cnv.bed.gz",        lambda: simulate_cnv(args.depth)),
        "tumour":     ("wgs_tumour.bed.gz",     lambda: simulate_tumour_heterogeneity(args.depth * 2)),
        "rna_uneven": ("rnaseq_uneven.bed.gz",  lambda: simulate_rna_coverage(args.depth)),
    }

    to_run = list(modes.keys()) if args.mode == "all" else [args.mode]

    print(f"\n  Karyotype Painter — generating synthetic coverage data")
    print(f"  Bin size: {BIN_SIZE:,} bp | Target depth: {args.depth}X\n")

    for mode in to_run:
        filename, generator = modes[mode]
        out = args.out if (args.out and len(to_run) == 1) else filename
        print(f"  Simulating: {mode}")
        df = generator()
        write_bed_gz(df, out)

    if args.mode == "all":
        print("\n  All datasets ready! Now run:")
        print("    python coverage_panorama.py \\")
        print("      -s wgs_normal.bed.gz wgs_cnv.bed.gz wgs_tumour.bed.gz rnaseq_uneven.bed.gz \\")
        print("      -n 'Normal WGS' 'CNV sample' 'Tumour 65%' 'RNA-seq' \\")
        print("      --mode compare -o panorama_comparison.png")


if __name__ == "__main__":
    import sys
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass
    main()
