#!/usr/bin/env python3
"""
compare_aligners.py — Multi-Aligner Benchmark Comparison
=========================================================
Reads the TSV produced by run_benchmark.sh, generates comparison
bar plots, and writes an HTML report.

Usage:
  python compare_aligners.py \\
    --metrics benchmark_results/metrics/benchmark_raw.tsv \\
    --output  benchmark_results/benchmark_report.html \\
    --plot    benchmark_results/benchmark_plots.png \\
    --reads   2000000
"""

import argparse
import base64
import json
from io import BytesIO
from pathlib import Path

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

matplotlib.use("Agg")

ALIGNER_COLOURS = {
    "bwa-mem2":  "#3E92CC",
    "minimap2":  "#E8A838",
    "hisat2":    "#2D9E4E",
    "bwa":       "#9B59B6",
}
DEFAULT_COLOUR = "#999999"


def load_metrics(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()
    # Ensure numeric columns
    for col in ["runtime_sec", "max_rss_mb", "alignment_rate",
                "properly_paired_pct", "mapped_reads", "total_reads"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def make_plots(df: pd.DataFrame, n_reads: int = 0) -> str:
    """Generate comparison plots and return base64-encoded PNG."""
    aligners = df["aligner"].tolist()
    colours  = [ALIGNER_COLOURS.get(a, DEFAULT_COLOUR) for a in aligners]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), facecolor="#FAFAFA")
    fig.suptitle(
        f"Multi-Aligner Benchmark" + (f"  ({n_reads:,} reads)" if n_reads else ""),
        fontsize=14, fontweight="bold", color="#1A1A2E", y=1.01,
    )

    def bar(ax, values, title, ylabel, colour_list=None, threshold=None):
        bars = ax.bar(aligners, values,
                      color=colour_list or colours, width=0.5,
                      edgecolor="white", linewidth=1.5)
        ax.set_title(title, fontweight="bold", fontsize=10, color="#1A1A2E")
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_ylim(0, max(values) * 1.2)
        ax.set_facecolor("#F5F5F5")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Value labels on bars
        for b, v in zip(bars, values):
            ax.text(
                b.get_x() + b.get_width() / 2.0,
                b.get_height() + max(values) * 0.01,
                f"{v:.1f}",
                ha="center", va="bottom", fontsize=9, color="#333333",
            )

        if threshold:
            ax.axhline(threshold, color="#D63031", linestyle="--",
                       linewidth=1.2, alpha=0.7, label=f"Threshold: {threshold}")
            ax.legend(fontsize=8)

    # Plot 1: Alignment rate
    bar(axes[0][0], df["alignment_rate"].tolist(),
        "Alignment Rate", "% reads mapped",
        threshold=90.0)

    # Plot 2: Runtime
    bar(axes[0][1], df["runtime_sec"].tolist(),
        "Runtime", "seconds  (lower = faster)",
        colour_list=[ALIGNER_COLOURS.get(a, DEFAULT_COLOUR) for a in aligners])

    # Plot 3: Memory usage
    bar(axes[1][0], df["max_rss_mb"].tolist(),
        "Peak Memory (MaxRSS)", "MB  (lower = better)")

    # Plot 4: Properly paired
    pp_vals = df["properly_paired_pct"].tolist() if "properly_paired_pct" in df.columns else [0] * len(df)
    bar(axes[1][1], pp_vals,
        "Properly Paired Rate", "% reads properly paired",
        threshold=85.0)

    plt.tight_layout()

    buf = BytesIO()
    plt.savefig(buf, format="png", dpi=150, bbox_inches="tight",
                facecolor="#FAFAFA")
    buf.seek(0)
    img_b64 = base64.b64encode(buf.read()).decode("ascii")
    plt.close()
    return img_b64


def speed_score(df: pd.DataFrame) -> pd.DataFrame:
    """Compute a composite benchmark score (higher = better overall)."""
    df = df.copy()
    # Normalise metrics (0–1), higher is always better
    df["score_align"]   = df["alignment_rate"]       / df["alignment_rate"].max()
    df["score_speed"]   = (1 / df["runtime_sec"])    / (1 / df["runtime_sec"]).max()
    df["score_memory"]  = (1 / df["max_rss_mb"])     / (1 / df["max_rss_mb"]).max()
    df["score_paired"]  = df.get("properly_paired_pct", pd.Series([1]*len(df))) / 100

    # Weighted composite
    df["composite"] = (
        df["score_align"]  * 0.40 +
        df["score_speed"]  * 0.25 +
        df["score_memory"] * 0.15 +
        df["score_paired"] * 0.20
    ) * 100
    return df.sort_values("composite", ascending=False)


def make_html(df: pd.DataFrame, img_b64: str, n_reads: int) -> str:
    scored = speed_score(df)

    table_rows = ""
    for _, row in scored.iterrows():
        aligner = row["aligner"]
        colour  = ALIGNER_COLOURS.get(aligner, DEFAULT_COLOUR)
        top     = "⭐" if _ == scored.index[0] else ""
        table_rows += f"""
        <tr>
          <td><b style="color:{colour}">{aligner}</b> {top}</td>
          <td>{row['alignment_rate']:.2f}%</td>
          <td>{row.get('properly_paired_pct', 'N/A')}</td>
          <td>{row['runtime_sec']:.1f}s</td>
          <td>{row['max_rss_mb']:.0f} MB</td>
          <td><b>{row['composite']:.1f}</b></td>
        </tr>"""

    recommendation = scored.iloc[0]["aligner"]

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Aligner Benchmark Report</title>
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ font-family: 'Segoe UI', Roboto, sans-serif; background: #f8f9fa;
            color: #212529; padding: 2rem; line-height: 1.6; }}
    h1 {{ color: #0d6efd; margin-bottom: 0.25rem; }}
    .subtitle {{ color: #6c757d; margin-bottom: 2rem; font-size: 0.9rem; }}
    .card {{ background: white; border-radius: 10px; padding: 1.5rem;
             margin-bottom: 1.5rem; box-shadow: 0 2px 8px rgba(0,0,0,0.07); }}
    h2 {{ color: #343a40; margin-bottom: 1rem; font-size: 1.1rem; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th {{ background: #0d6efd; color: white; padding: 10px 14px; text-align: left;
          font-weight: 600; }}
    td {{ padding: 10px 14px; border-bottom: 1px solid #dee2e6; }}
    tr:hover {{ background: #f1f3f5; }}
    .badge {{ display: inline-block; padding: 4px 10px; border-radius: 20px;
              font-size: 0.8rem; font-weight: 600; background: #d1ecf1;
              color: #0c5460; }}
    .rec {{ background: #d4edda; color: #155724; border-left: 4px solid #28a745;
            padding: 1rem 1.5rem; border-radius: 4px; margin-top: 1rem; }}
    img {{ max-width: 100%; border-radius: 8px; }}
    .score-note {{ font-size: 0.8rem; color: #6c757d; margin-top: 0.5rem; }}
  </style>
</head>
<body>
  <h1>🧬 Multi-Aligner Benchmark Report</h1>
  <p class="subtitle">
    {n_reads:,} reads benchmarked | {len(df)} aligners tested
  </p>

  <div class="card">
    <h2>📊 Performance Comparison</h2>
    <img src="data:image/png;base64,{img_b64}" alt="Benchmark plots">
  </div>

  <div class="card">
    <h2>📋 Results Table</h2>
    <table>
      <thead>
        <tr>
          <th>Aligner</th>
          <th>Alignment Rate</th>
          <th>Properly Paired</th>
          <th>Runtime</th>
          <th>Peak Memory</th>
          <th>Composite Score</th>
        </tr>
      </thead>
      <tbody>{table_rows}</tbody>
    </table>
    <p class="score-note">
      Composite score: 40% alignment rate + 25% speed + 20% properly paired + 15% memory.
    </p>

    <div class="rec">
      <b>Recommendation:</b> <code>{recommendation}</code> has the best overall
      composite score for this dataset. Consider your primary use case:
      BWA-MEM2 for DNA-seq variant calling; HISAT2/STAR for RNA-seq;
      minimap2 for long reads or rapid alignment checks.
    </div>
  </div>
</body>
</html>"""


def main():
    parser = argparse.ArgumentParser(description="Multi-aligner comparison and report")
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--output",  default="benchmark_report.html")
    parser.add_argument("--plot",    default="benchmark_plots.png")
    parser.add_argument("--reads",   type=int, default=0)
    args = parser.parse_args()

    df  = load_metrics(args.metrics)
    print(f"  Loaded {len(df)} aligner results")

    img_b64 = make_plots(df, n_reads=args.reads)

    # Also save standalone plot
    with open(args.plot, "wb") as fh:
        fh.write(base64.b64decode(img_b64))
    print(f"  [OK] Plot saved: {args.plot}")

    html = make_html(df, img_b64, n_reads=args.reads)
    Path(args.output).write_text(html, encoding="utf-8")
    print(f"  [OK] Report saved: {args.output}")

    # Print summary to terminal
    scored = speed_score(df)
    print("\n  ┌─ Composite Scores ─────────────────────┐")
    for _, row in scored.iterrows():
        bar = "█" * int(row["composite"] / 5)
        print(f"  │  {row['aligner']:<12} {bar:<20} {row['composite']:.1f}")
    print("  └────────────────────────────────────────┘")


if __name__ == "__main__":
    import sys
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass
    main()
