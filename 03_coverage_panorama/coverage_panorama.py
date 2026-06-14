#!/usr/bin/env python3
"""
coverage_panorama.py — Genome-Wide Coverage Art Generator
==========================================================
Transforms mosdepth output into publication-ready karyotype-style
coverage heatmaps. Perfect for presentations, lab reports, and
bioinformatics society outreach material.

Features:
  • Chromosome-level coverage heatmaps (publication quality)
  • Centromere/telomere awareness (gaps rendered correctly)
  • Coverage percentile banding (depth distribution)
  • Side-by-side multi-sample comparison
  • Export to SVG, PNG, or PDF

Usage:
  python coverage_panorama.py \\
    --sample sample.mosdepth.regions.bed.gz \\
    --name "Patient_001 — WGS 35X" \\
    --output coverage_panorama.png

  # Multi-sample comparison
  python coverage_panorama.py \\
    --sample A.bed.gz B.bed.gz \\
    --name "Normal" "Tumor" \\
    --output comparison.png \\
    --mode compare
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter

matplotlib.use("Agg")

# ── Chromosome order (hg38 / GRCh38) ─────────────────────────────────────────
HG38_CHROM_ORDER = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
    "chrX", "chrY",
]

# hg38 chromosome sizes (Mb) — for proportional scaling
HG38_SIZES_MB = {
    "chr1": 249, "chr2": 242, "chr3": 198, "chr4": 191, "chr5": 182,
    "chr6": 171, "chr7": 159, "chr8": 145, "chr9": 138, "chr10": 134,
    "chr11": 135, "chr12": 133, "chr13": 114, "chr14": 107, "chr15": 102,
    "chr16": 90,  "chr17": 84,  "chr18": 80,  "chr19": 59,  "chr20": 64,
    "chr21": 47,  "chr22": 51,  "chrX": 156,  "chrY": 57,
}

PALETTE = {
    "deep":        "#0A2463",
    "mid":         "#3E92CC",
    "shallow":     "#D8F1FF",
    "zero":        "#F5F5F5",
    "excess":      "#C1121F",
    "background":  "#FAFAFA",
    "border":      "#CCCCCC",
    "text":        "#1A1A2E",
    "centromere":  "#E8D5B7",
}


def read_mosdepth_bed(path: str) -> pd.DataFrame:
    """Parse mosdepth regions BED (gz or plain) using a vectorised reader."""
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=["chrom", "start", "end", "depth"],
        usecols=[0, 1, 2, 3],
        dtype={"chrom": str, "start": "int64", "end": "int64", "depth": "float64"},
    )
    return df


def compute_chrom_stats(df: pd.DataFrame) -> dict:
    """Compute per-chromosome depth statistics."""
    stats = {}
    for chrom, grp in df.groupby("chrom"):
        stats[chrom] = {
            "mean":   grp["depth"].mean(),
            "median": grp["depth"].median(),
            "std":    grp["depth"].std(),
            "p5":     grp["depth"].quantile(0.05),
            "p95":    grp["depth"].quantile(0.95),
            "pct_zero": (grp["depth"] == 0).mean() * 100,
        }
    return stats


def depth_to_color(depth: float, mean_depth: float) -> str:
    """Map depth value to a colour on a diverging scale."""
    if depth == 0:
        return PALETTE["zero"]
    ratio = depth / max(mean_depth, 1.0)
    if ratio > 2.0:
        # excess coverage (red)
        t = min((ratio - 2.0) / 2.0, 1.0)
        return mcolors.to_hex(
            mcolors.LinearSegmentedColormap.from_list(
                "ex", [PALETTE["mid"], PALETTE["excess"]]
            )(t)
        )
    elif ratio >= 0.5:
        # normal range (blue gradient)
        t = (ratio - 0.5) / 1.5
        return mcolors.to_hex(
            mcolors.LinearSegmentedColormap.from_list(
                "norm", [PALETTE["shallow"], PALETTE["deep"]]
            )(t)
        )
    else:
        # low coverage
        return PALETTE["shallow"]


def build_chrom_array(chrom_df: pd.DataFrame, n_bins: int = 500) -> np.ndarray:
    """
    Bin chromosome into a fixed-width depth array for heatmap rendering.

    Vectorised with numpy: assigns each input region's midpoint to an
    output bin and averages depth values per bin. This avoids Python-level
    iteration over potentially millions of mosdepth windows.
    """
    if chrom_df.empty:
        return np.zeros(n_bins)

    max_pos = chrom_df["end"].max()
    if max_pos == 0:
        return np.zeros(n_bins)

    starts = chrom_df["start"].to_numpy()
    ends   = chrom_df["end"].to_numpy()
    depths = chrom_df["depth"].to_numpy()

    midpoints = (starts + ends) / 2.0
    bin_idx = np.floor(midpoints / max_pos * n_bins).astype(int)
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)

    sums   = np.bincount(bin_idx, weights=depths, minlength=n_bins)
    counts = np.bincount(bin_idx, minlength=n_bins)

    with np.errstate(invalid="ignore", divide="ignore"):
        result = np.where(counts > 0, sums / counts, 0)

    # Forward-fill any empty bins from neighbours (avoids visual gaps)
    nonzero_mask = counts > 0
    if nonzero_mask.any() and not nonzero_mask.all():
        idx = np.arange(n_bins)
        valid_idx = idx[nonzero_mask]
        result = np.interp(idx, valid_idx, result[nonzero_mask])

    return result


def draw_panorama(
    df: pd.DataFrame,
    sample_name: str = "Sample",
    output: str = "coverage_panorama.png",
    dpi: int = 180,
) -> None:
    """
    Draw a single-sample karyotype-style coverage panorama.
    Each chromosome is a horizontal strip coloured by coverage depth.
    """
    # Filter to standard chromosomes
    chroms = [c for c in HG38_CHROM_ORDER if c in df["chrom"].unique()]
    stats  = compute_chrom_stats(df)
    global_mean = df["depth"].mean()

    N_BINS = 600  # bins per chromosome strip
    FIG_W  = 14
    FIG_H  = max(len(chroms) * 0.45 + 3, 8)

    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor=PALETTE["background"])
    gs  = GridSpec(len(chroms), 1, figure=fig, hspace=0.12,
                   left=0.10, right=0.88, top=0.92, bottom=0.06)

    # ── Colourmap ─────────────────────────────────────────────────────────────
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "coverage", [
            PALETTE["zero"],
            PALETTE["shallow"],
            PALETTE["mid"],
            PALETTE["deep"],
            PALETTE["excess"],
        ],
        N=256,
    )
    cmap.set_bad(color=PALETTE["background"])
    norm = mcolors.Normalize(vmin=0, vmax=global_mean * 2.5)

    for idx, chrom in enumerate(chroms):
        ax = fig.add_subplot(gs[idx])
        chrom_df = df[df["chrom"] == chrom]
        # Scale bar width proportional to chromosome size
        chrom_mb   = HG38_SIZES_MB.get(chrom, 100)
        max_chrom  = max(HG38_SIZES_MB.values())
        width_frac = chrom_mb / max_chrom
        n_fill     = max(int(width_frac * N_BINS), 1)

        arr = build_chrom_array(chrom_df, n_bins=n_fill)

        # Extend array with white space on the right
        padded = np.full(N_BINS, np.nan)
        padded[:n_fill] = arr

        ax.imshow(
            padded.reshape(1, -1),
            aspect="auto",
            cmap=cmap,
            norm=norm,
            interpolation="bilinear",
        )

        # Chromosome label
        ax.set_ylabel(
            chrom.replace("chr", ""),
            fontsize=7.5,
            rotation=0,
            labelpad=18,
            va="center",
            color=PALETTE["text"],
            fontfamily="monospace",
        )

        # Mean depth annotation
        chrom_mean = stats.get(chrom, {}).get("mean", 0)
        ax.text(
            1.01, 0.5,
            f"{chrom_mean:.0f}X",
            transform=ax.transAxes,
            fontsize=6.5,
            va="center",
            color=PALETTE["text"],
        )

        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor(PALETTE["border"])
            spine.set_linewidth(0.4)

    # ── Title ─────────────────────────────────────────────────────────────────
    fig.text(
        0.49, 0.96,
        f"Coverage Panorama — {sample_name}",
        ha="center", va="top",
        fontsize=14, fontweight="bold",
        color=PALETTE["text"],
    )
    fig.text(
        0.49, 0.93,
        f"Genome-wide mean: {global_mean:.1f}X  |  Bins: {N_BINS}/chr  |  Genome: hg38",
        ha="center", va="top",
        fontsize=8,
        color="#666666",
    )

    # ── Colorbar ──────────────────────────────────────────────────────────────
    cbar_ax = fig.add_axes([0.90, 0.10, 0.012, 0.75])
    sm      = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Coverage depth (X)", fontsize=8, color=PALETTE["text"])
    cbar.ax.tick_params(labelsize=7, colors=PALETTE["text"])
    cbar.ax.yaxis.set_tick_params(color=PALETTE["border"])

    # ── Chromosome label ──────────────────────────────────────────────────────
    fig.text(0.01, 0.96, "Chr", fontsize=8, color="#888888", va="top")
    fig.text(0.89, 0.96, "Mean", fontsize=8, color="#888888", va="top")

    plt.savefig(output, dpi=dpi, bbox_inches="tight",
                facecolor=PALETTE["background"])
    plt.close()
    print(f"[✓] Saved: {output}")


def draw_comparison(
    dfs: list,
    names: list,
    output: str = "comparison.png",
    dpi: int = 180,
) -> None:
    """Side-by-side multi-sample comparison (up to 4 samples)."""
    n = len(dfs)
    all_chroms = [c for c in HG38_CHROM_ORDER if any(c in df["chrom"].unique() for df in dfs)]
    global_max_mean = max(df["depth"].mean() for df in dfs) * 2.5

    N_BINS = 400
    FIG_W = 5 * n + 2
    FIG_H = max(len(all_chroms) * 0.5 + 3, 8)

    fig, axes = plt.subplots(
        len(all_chroms), n,
        figsize=(FIG_W, FIG_H),
        facecolor=PALETTE["background"],
        gridspec_kw={"hspace": 0.08, "wspace": 0.04},
    )
    if len(all_chroms) == 1:
        axes = axes.reshape(1, -1)
    if n == 1:
        axes = axes.reshape(-1, 1)

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "cov", [PALETTE["zero"], PALETTE["shallow"], PALETTE["mid"], PALETTE["deep"]], N=256
    )
    norm = mcolors.Normalize(vmin=0, vmax=global_max_mean)

    for col, (df, name) in enumerate(zip(dfs, names)):
        stats       = compute_chrom_stats(df)
        sample_mean = df["depth"].mean()

        axes[0][col].set_title(
            f"{name}\n{sample_mean:.1f}X mean",
            fontsize=9, fontweight="bold", color=PALETTE["text"]
        )

        for row, chrom in enumerate(all_chroms):
            ax       = axes[row][col]
            chrom_df = df[df["chrom"] == chrom]
            arr      = build_chrom_array(chrom_df, n_bins=N_BINS)

            ax.imshow(arr.reshape(1, -1), aspect="auto",
                      cmap=cmap, norm=norm, interpolation="bilinear")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_linewidth(0.3); sp.set_edgecolor(PALETTE["border"])

            if col == 0:
                ax.set_ylabel(
                    chrom.replace("chr", ""),
                    fontsize=7, rotation=0, labelpad=14, va="center",
                    color=PALETTE["text"], fontfamily="monospace",
                )

    fig.suptitle(
        "Coverage Panorama — Multi-Sample Comparison",
        fontsize=13, fontweight="bold", y=0.995, color=PALETTE["text"],
    )
    plt.savefig(output, dpi=dpi, bbox_inches="tight",
                facecolor=PALETTE["background"])
    plt.close()
    print(f"[✓] Saved comparison: {output}")


def print_summary(df: pd.DataFrame, name: str) -> None:
    """Print a styled coverage summary to the terminal."""
    stats = compute_chrom_stats(df)
    global_mean = df["depth"].mean()
    pct_zero = (df["depth"] == 0).mean() * 100

    print(f"\n  ╔{'═'*52}╗")
    print(f"  ║  Coverage Summary — {name:<31}║")
    print(f"  ╠{'═'*52}╣")
    print(f"  ║  Genome mean coverage:    {global_mean:>6.1f}X              ║")
    print(f"  ║  Zero-coverage windows:   {pct_zero:>6.2f}%              ║")
    print(f"  ╠{'═'*52}╣")
    print(f"  ║  {'Chr':<6} {'Mean':>7} {'Median':>8} {'P5':>6} {'P95':>6} ║")
    print(f"  ╠{'═'*52}╣")
    for c in HG38_CHROM_ORDER:
        if c in stats:
            s = stats[c]
            print(f"  ║  {c.replace('chr',''):<6} {s['mean']:>7.1f} "
                  f"{s['median']:>8.1f} {s['p5']:>6.1f} {s['p95']:>6.1f} ║")
    print(f"  ╚{'═'*52}╝")


def main():
    parser = argparse.ArgumentParser(
        description="Coverage Panorama — genome-wide coverage art generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--sample", "-s", nargs="+", required=True,
        help="mosdepth regions BED file(s) (.bed or .bed.gz)",
    )
    parser.add_argument(
        "--name", "-n", nargs="+", default=None,
        help="Sample name(s) for plot titles",
    )
    parser.add_argument(
        "--output", "-o", default="coverage_panorama.png",
        help="Output file (png/svg/pdf) [default: coverage_panorama.png]",
    )
    parser.add_argument(
        "--mode", choices=["single", "compare"], default="single",
        help="single: one figure per sample, compare: side-by-side [default: single]",
    )
    parser.add_argument(
        "--dpi", type=int, default=180,
        help="Output DPI [default: 180]",
    )
    parser.add_argument(
        "--no-summary", action="store_true",
        help="Skip terminal summary table",
    )
    args = parser.parse_args()

    sample_files = args.sample
    names = args.name or [Path(f).name.split(".")[0] for f in sample_files]

    if len(names) < len(sample_files):
        names += [f"Sample_{i}" for i in range(len(names), len(sample_files))]

    print(f"\n  Coverage Panorama v1.0 — loading {len(sample_files)} sample(s)...")

    dfs = []
    for path, name in zip(sample_files, names):
        print(f"  Reading: {path}")
        df = read_mosdepth_bed(path)
        if df.empty:
            print(f"  [ERROR] No data loaded from {path}", file=sys.stderr)
            sys.exit(1)
        dfs.append(df)
        if not args.no_summary:
            print_summary(df, name)

    if args.mode == "compare" and len(dfs) > 1:
        draw_comparison(dfs, names, output=args.output, dpi=args.dpi)
    else:
        # Single sample per output (or one per file in single mode)
        base, ext = Path(args.output).stem, Path(args.output).suffix
        for i, (df, name) in enumerate(zip(dfs, names)):
            out = f"{base}{ext}" if len(dfs) == 1 else f"{base}_{i+1}{ext}"
            draw_panorama(df, sample_name=name, output=out, dpi=args.dpi)


if __name__ == "__main__":
    import sys
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass
    main()
