#!/usr/bin/env python3
"""
qc_sentinel.py — NGS QC Sentinel: Unified Quality Control Dashboard
====================================================================
Aggregates samtools flagstat, Picard MarkDuplicates metrics, mosdepth
summary, and HISAT2/STAR alignment logs into a single Rich terminal
dashboard with colour-coded pass/warn/fail thresholds.

Features:
  • Auto-discovers QC files in a results directory
  • Colour-coded metrics (green/yellow/red) per GATK/ENCODE thresholds
  • Multi-sample table for cohort-level overview
  • JSON export for downstream processing
  • HTML report generation

Usage:
  # Single sample
  python qc_sentinel.py --dir results/

  # Multi-sample cohort
  python qc_sentinel.py --dir results/ --multi

  # Export JSON
  python qc_sentinel.py --dir results/ --json qc_summary.json

  # Export HTML report
  python qc_sentinel.py --dir results/ --html qc_report.html
"""

import argparse
import json
import re
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

from rich.console import Console
from rich.table import Table
from rich import box
from rich.panel import Panel
from rich.columns import Columns
from rich.text import Text
from rich.progress import track
from rich.rule import Rule

# ── Thresholds (GATK Best Practices + ENCODE) ─────────────────────────────────
THRESHOLDS = {
    "mapping_rate": {
        "wgs":  {"pass": 90.0, "warn": 75.0},
        "rna":  {"pass": 70.0, "warn": 55.0},
    },
    "properly_paired": {"pass": 85.0, "warn": 70.0},
    "duplication_rate": {
        "wgs":  {"pass": 20.0, "warn": 35.0},   # lower is better
        "wes":  {"pass": 40.0, "warn": 60.0},
        "rna":  {"pass": 30.0, "warn": 50.0},
    },
    "mean_coverage_wgs": {"pass": 25.0, "warn": 15.0},
    "pct_callable_wgs":  {"pass": 85.0, "warn": 70.0},
}

console = Console(highlight=False)


# ─────────────────────────────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class FlagstatMetrics:
    total_reads:       int   = 0
    mapped_reads:      int   = 0
    mapped_pct:        float = 0.0
    properly_paired:   int   = 0
    properly_paired_pct: float = 0.0
    duplicates:        int   = 0
    secondary:         int   = 0
    supplementary:     int   = 0


@dataclass
class PicardDupMetrics:
    estimated_library_size: int   = 0
    duplication_rate:       float = 0.0
    optical_dup_rate:       float = 0.0
    pcr_dup_rate:           float = 0.0
    read_pairs_examined:    int   = 0


@dataclass
class MosdepthMetrics:
    mean_coverage:    float = 0.0
    median_coverage:  float = 0.0
    pct_1x:           float = 0.0
    pct_10x:          float = 0.0
    pct_20x:          float = 0.0
    pct_30x:          float = 0.0


@dataclass
class SampleQC:
    sample_name: str = ""
    flagstat:    FlagstatMetrics  = field(default_factory=FlagstatMetrics)
    picard:      PicardDupMetrics = field(default_factory=PicardDupMetrics)
    mosdepth:    MosdepthMetrics  = field(default_factory=MosdepthMetrics)
    source_files: dict = field(default_factory=dict)


# ─────────────────────────────────────────────────────────────────────────────
# Parsers
# ─────────────────────────────────────────────────────────────────────────────

def parse_flagstat(path: Path) -> FlagstatMetrics:
    """Parse samtools flagstat output."""
    m   = FlagstatMetrics()
    txt = path.read_text(encoding="utf-8")

    def extract(pattern: str, default=0):
        match = re.search(pattern, txt, re.MULTILINE)
        return int(match.group(1)) if match else default

    def extract_pct(pattern: str) -> float:
        match = re.search(pattern, txt, re.MULTILINE)
        return float(match.group(1)) if match else 0.0

    m.total_reads     = extract(r"^(\d+) \+ \d+ in total")
    m.mapped_reads    = extract(r"^(\d+) \+ \d+ mapped")
    m.properly_paired = extract(r"^(\d+) \+ \d+ properly paired")
    m.duplicates      = extract(r"^(\d+) \+ \d+ duplicates")
    m.secondary       = extract(r"^(\d+) \+ \d+ secondary")
    m.supplementary   = extract(r"^(\d+) \+ \d+ supplementary")

    if m.total_reads > 0:
        m.mapped_pct        = m.mapped_reads / m.total_reads * 100
        m.properly_paired_pct = m.properly_paired / m.total_reads * 100

    return m


def parse_picard_metrics(path: Path) -> PicardDupMetrics:
    """Parse Picard MarkDuplicates metrics file."""
    m   = PicardDupMetrics()
    txt = path.read_text(encoding="utf-8")

    # Picard metrics table has a header row followed by data row
    lines = [ln for ln in txt.splitlines() if ln and not ln.startswith("#")]
    header_idx = next((i for i, l in enumerate(lines) if l.startswith("LIBRARY")), None)
    if header_idx is None:
        return m

    try:
        headers = lines[header_idx].split("\t")
        values  = lines[header_idx + 1].split("\t")
        row     = dict(zip(headers, values))

        m.estimated_library_size = int(float(row.get("ESTIMATED_LIBRARY_SIZE", 0) or 0))
        m.duplication_rate       = float(row.get("PERCENT_DUPLICATION", 0) or 0) * 100
        m.read_pairs_examined    = int(float(row.get("READ_PAIRS_EXAMINED", 0) or 0))

        total_dup = float(row.get("READ_PAIR_DUPLICATES", 0) or 0)
        optical   = float(row.get("READ_PAIR_OPTICAL_DUPLICATES", 0) or 0)
        examined  = max(m.read_pairs_examined, 1)

        m.optical_dup_rate = optical / examined * 100
        m.pcr_dup_rate     = (total_dup - optical) / examined * 100
    except (IndexError, ValueError, KeyError):
        pass

    return m


def parse_mosdepth_summary(path: Path) -> MosdepthMetrics:
    """Parse mosdepth summary file (mosdepth.summary.txt)."""
    m   = MosdepthMetrics()
    txt = path.read_text(encoding="utf-8")

    for line in txt.splitlines():
        if line.startswith("total") or line.startswith("genome"):
            parts = line.split("\t")
            if len(parts) >= 4:
                try:
                    m.mean_coverage   = float(parts[3])
                    m.median_coverage = float(parts[4]) if len(parts) > 4 else 0.0
                except ValueError:
                    pass

    thresh_path = Path(str(path).replace("summary.txt", "thresholds.bed.gz"))
    if thresh_path.exists():
        import gzip
        with gzip.open(thresh_path, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("chrom"):
                    continue
                parts = line.strip().split("\t")
                # thresholds file has columns: chrom start end 1X 10X 20X 30X
                if len(parts) >= 7 and parts[0] in ("total", "genome"):
                    try:
                        total = float(parts[3]) if float(parts[3]) > 0 else 1
                        m.pct_1x  = float(parts[4]) / total * 100
                        m.pct_10x = float(parts[5]) / total * 100
                        m.pct_20x = float(parts[6]) / total * 100
                        m.pct_30x = float(parts[7]) / total * 100 if len(parts) > 7 else 0
                    except (ValueError, IndexError):
                        pass
    return m


def discover_qc_files(base_dir: Path, sample_name: str = None) -> SampleQC:
    """
    Auto-discover QC files in a directory tree.
    Looks for *.flagstat.txt, *.markdup_metrics.txt, *.mosdepth.summary.txt
    """
    qc = SampleQC()

    # Try to find files matching common naming patterns
    flagstat_files = list(base_dir.rglob("*.flagstat.txt"))
    picard_files   = list(base_dir.rglob("*.markdup_metrics.txt"))
    mosdepth_files = list(base_dir.rglob("*.mosdepth.summary.txt"))

    if not flagstat_files and not picard_files:
        return qc

    # Use first match or filter by sample name
    def pick(files, name=None):
        if name:
            hits = [f for f in files if name in f.name]
            return hits[0] if hits else (files[0] if files else None)
        return files[0] if files else None

    fs = pick(flagstat_files, sample_name)
    ps = pick(picard_files, sample_name)
    ms = pick(mosdepth_files, sample_name)

    if fs:
        qc.flagstat    = parse_flagstat(fs)
        qc.source_files["flagstat"] = str(fs)
    if ps:
        qc.picard      = parse_picard_metrics(ps)
        qc.source_files["picard"] = str(ps)
    if ms:
        qc.mosdepth    = parse_mosdepth_summary(ms)
        qc.source_files["mosdepth"] = str(ms)

    if fs:
        qc.sample_name = sample_name or fs.name.replace(".flagstat.txt", "")

    return qc


# ─────────────────────────────────────────────────────────────────────────────
# Threshold evaluation
# ─────────────────────────────────────────────────────────────────────────────

def status(value: float, thresholds: dict, lower_is_better: bool = False) -> tuple[str, str]:
    """
    Returns (status_label, rich_style) for a metric value.
    status: PASS / WARN / FAIL
    """
    p, w = thresholds["pass"], thresholds["warn"]
    if lower_is_better:
        if value <= p:      return "PASS", "bold green"
        elif value <= w:    return "WARN", "bold yellow"
        else:               return "FAIL", "bold red"
    else:
        if value >= p:      return "PASS", "bold green"
        elif value >= w:    return "WARN", "bold yellow"
        else:               return "FAIL", "bold red"


def fmt_status(value: float, thresholds: dict, lower_is_better: bool = False,
               fmt: str = ".1f", unit: str = "") -> Text:
    label, style = status(value, thresholds, lower_is_better)
    t = Text()
    t.append(f"{value:{fmt}}{unit}  ", style="white")
    t.append(f"[{label}]", style=style)
    return t


# ─────────────────────────────────────────────────────────────────────────────
# Display functions
# ─────────────────────────────────────────────────────────────────────────────

def display_single(qc: SampleQC) -> None:
    console.print()
    console.print(Rule(f"[bold cyan]QC Sentinel — {qc.sample_name}[/bold cyan]"))

    # ── Alignment metrics ────────────────────────────────────────────────────
    align_table = Table(
        title="[bold]Alignment Metrics[/bold]  (samtools flagstat)",
        box=box.ROUNDED,
        show_header=True,
        header_style="bold blue",
        min_width=60,
    )
    align_table.add_column("Metric",   style="cyan", min_width=28)
    align_table.add_column("Value",    justify="right")
    align_table.add_column("Status",   justify="center")

    f = qc.flagstat
    align_table.add_row(
        "Total reads",
        f"{f.total_reads:,}",
        Text("—", style="dim"),
    )
    align_table.add_row(
        "Mapped reads",
        f"{f.mapped_reads:,}",
        fmt_status(f.mapped_pct, THRESHOLDS["mapping_rate"]["wgs"], fmt=".1f", unit="%"),
    )
    align_table.add_row(
        "Properly paired",
        f"{f.properly_paired:,}",
        fmt_status(f.properly_paired_pct, THRESHOLDS["properly_paired"], fmt=".1f", unit="%"),
    )
    align_table.add_row(
        "Duplicate reads",
        f"{f.duplicates:,}",
        Text("see below", style="dim"),
    )
    align_table.add_row(
        "Secondary alignments",
        f"{f.secondary:,}",
        Text("—", style="dim"),
    )

    # ── Duplication metrics ──────────────────────────────────────────────────
    dup_table = Table(
        title="[bold]Duplication Metrics[/bold]  (Picard MarkDuplicates)",
        box=box.ROUNDED,
        show_header=True,
        header_style="bold blue",
        min_width=60,
    )
    dup_table.add_column("Metric",   style="cyan", min_width=28)
    dup_table.add_column("Value",    justify="right")
    dup_table.add_column("Status",   justify="center")

    p = qc.picard
    dup_table.add_row(
        "Duplication rate",
        f"{p.duplication_rate:.2f}%",
        fmt_status(p.duplication_rate, THRESHOLDS["duplication_rate"]["wgs"],
                   lower_is_better=True, fmt=".2f", unit="%"),
    )
    dup_table.add_row(
        "PCR duplicates",
        f"{p.pcr_dup_rate:.2f}%",
        Text("—", style="dim"),
    )
    dup_table.add_row(
        "Optical duplicates",
        f"{p.optical_dup_rate:.2f}%",
        Text("—", style="dim"),
    )
    dup_table.add_row(
        "Estimated library size",
        f"{p.estimated_library_size:,}" if p.estimated_library_size > 0 else "N/A",
        Text("—", style="dim"),
    )

    # ── Coverage metrics ─────────────────────────────────────────────────────
    cov_table = Table(
        title="[bold]Coverage Metrics[/bold]  (mosdepth)",
        box=box.ROUNDED,
        show_header=True,
        header_style="bold blue",
        min_width=60,
    )
    cov_table.add_column("Metric",   style="cyan", min_width=28)
    cov_table.add_column("Value",    justify="right")
    cov_table.add_column("Status",   justify="center")

    md = qc.mosdepth
    cov_table.add_row(
        "Mean coverage",
        f"{md.mean_coverage:.1f}X",
        fmt_status(md.mean_coverage, THRESHOLDS["mean_coverage_wgs"], fmt=".1f", unit="X"),
    )
    cov_table.add_row(
        "Median coverage",
        f"{md.median_coverage:.1f}X",
        Text("—", style="dim"),
    )
    if md.pct_20x > 0:
        cov_table.add_row(
            "≥20X callable bases",
            f"{md.pct_20x:.1f}%",
            fmt_status(md.pct_20x, THRESHOLDS["pct_callable_wgs"], fmt=".1f", unit="%"),
        )
    if md.pct_30x > 0:
        cov_table.add_row(
            "≥30X callable bases",
            f"{md.pct_30x:.1f}%",
            Text("—", style="dim"),
        )

    console.print(Columns([align_table, dup_table]))
    console.print()
    console.print(cov_table)

    # ── Source files ─────────────────────────────────────────────────────────
    if qc.source_files:
        src = Table(box=box.SIMPLE, show_header=False, min_width=60)
        src.add_column("Type", style="dim")
        src.add_column("Path", style="dim italic")
        for k, v in qc.source_files.items():
            src.add_row(k, v)
        console.print(Panel(src, title="[dim]Source files[/dim]", border_style="dim"))


def display_multi(samples: list[SampleQC]) -> None:
    """Multi-sample overview table."""
    console.print()
    console.print(Rule("[bold cyan]QC Sentinel — Cohort Overview[/bold cyan]"))

    t = Table(box=box.ROUNDED, header_style="bold blue", min_width=100)
    t.add_column("Sample",          style="cyan bold",  min_width=20)
    t.add_column("Total reads",     justify="right",    min_width=14)
    t.add_column("Mapped %",        justify="center",   min_width=12)
    t.add_column("Paired %",        justify="center",   min_width=12)
    t.add_column("Dup rate",        justify="center",   min_width=12)
    t.add_column("Mean cov",        justify="center",   min_width=10)
    t.add_column("Library size",    justify="right",    min_width=14)

    for qc in samples:
        f  = qc.flagstat
        p  = qc.picard
        md = qc.mosdepth

        _, map_style  = status(f.mapped_pct, THRESHOLDS["mapping_rate"]["wgs"])
        _, pair_style = status(f.properly_paired_pct, THRESHOLDS["properly_paired"])
        _, dup_style  = status(p.duplication_rate, THRESHOLDS["duplication_rate"]["wgs"], lower_is_better=True)
        _, cov_style  = status(md.mean_coverage, THRESHOLDS["mean_coverage_wgs"])

        t.add_row(
            qc.sample_name,
            f"{f.total_reads:,}",
            Text(f"{f.mapped_pct:.1f}%",          style=map_style),
            Text(f"{f.properly_paired_pct:.1f}%",  style=pair_style),
            Text(f"{p.duplication_rate:.1f}%",     style=dup_style),
            Text(f"{md.mean_coverage:.1f}X",        style=cov_style),
            f"{p.estimated_library_size:,}" if p.estimated_library_size > 0 else "—",
        )

    console.print(t)
    console.print()
    console.print(
        "[dim]Legend: [bold green]PASS[/bold green] / "
        "[bold yellow]WARN[/bold yellow] / [bold red]FAIL[/bold red][/dim]"
    )


def export_json(samples: list, path: str) -> None:
    data = [asdict(s) for s in samples]
    Path(path).write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")
    console.print(f"[green][✓] JSON exported:[/green] {path}")


def export_html(samples: list, path: str) -> None:
    """Generate a minimal self-contained HTML summary."""
    rows = ""
    for qc in samples:
        f, p, md = qc.flagstat, qc.picard, qc.mosdepth
        st_map,  _ = status(f.mapped_pct,         THRESHOLDS["mapping_rate"]["wgs"])
        st_pair, _ = status(f.properly_paired_pct,THRESHOLDS["properly_paired"])
        st_dup,  _ = status(p.duplication_rate,   THRESHOLDS["duplication_rate"]["wgs"], True)
        st_cov,  _ = status(md.mean_coverage,      THRESHOLDS["mean_coverage_wgs"])
        colour  = {"PASS": "#2d9e4e", "WARN": "#e8a838", "FAIL": "#d63031"}
        rows += f"""
        <tr>
          <td><b>{qc.sample_name}</b></td>
          <td>{f.total_reads:,}</td>
          <td style="color:{colour[st_map]}">{f.mapped_pct:.1f}%</td>
          <td style="color:{colour[st_pair]}">{f.properly_paired_pct:.1f}%</td>
          <td style="color:{colour[st_dup]}">{p.duplication_rate:.1f}%</td>
          <td style="color:{colour[st_cov]}">{md.mean_coverage:.1f}X</td>
          <td>{p.estimated_library_size:,}</td>
        </tr>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>QC Sentinel Report</title>
  <style>
    body {{ font-family: 'Segoe UI', sans-serif; background: #f8f9fa; color: #212529; padding: 2rem; }}
    h1 {{ color: #0d6efd; }} table {{ border-collapse: collapse; width: 100%; }}
    th {{ background: #0d6efd; color: white; padding: 10px 14px; text-align: left; }}
    td {{ padding: 8px 14px; border-bottom: 1px solid #dee2e6; }}
    tr:hover {{ background: #f1f3f5; }}
    .badge {{ display: inline-block; padding: 2px 8px; border-radius: 4px; font-size: 0.75rem; }}
  </style>
</head>
<body>
  <h1>🧬 QC Sentinel Report</h1>
  <p>Generated for {len(samples)} sample(s). Thresholds: GATK Best Practices / ENCODE.</p>
  <table>
    <thead>
      <tr>
        <th>Sample</th><th>Total reads</th><th>Mapped %</th>
        <th>Properly paired %</th><th>Dup rate</th>
        <th>Mean coverage</th><th>Library size</th>
      </tr>
    </thead>
    <tbody>{rows}</tbody>
  </table>
</body>
</html>"""
    Path(path).write_text(html, encoding="utf-8")
    console.print(f"[green][✓] HTML report:[/green] {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="QC Sentinel — unified NGS QC terminal dashboard",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--dir",    "-d", required=True, help="Results directory to scan")
    parser.add_argument("--sample", "-s", nargs="+",     help="Sample name(s) to include")
    parser.add_argument("--multi",  "-m", action="store_true", help="Multi-sample cohort view")
    parser.add_argument("--json",         default=None,  help="Export JSON to this path")
    parser.add_argument("--html",         default=None,  help="Export HTML report to this path")
    args = parser.parse_args()

    base = Path(args.dir)
    if not base.exists():
        console.print(f"[red][ERROR] Directory not found: {base}[/red]")
        sys.exit(1)

    # Discover samples
    if args.sample:
        names = args.sample
    else:
        # Auto-discover from flagstat files
        names = sorted(
            set(f.name.replace(".flagstat.txt", "")
                for f in base.rglob("*.flagstat.txt"))
        ) or [""]

    samples = []
    for name in names:
        qc = discover_qc_files(base, name if name else None)
        if qc.sample_name or qc.flagstat.total_reads > 0:
            samples.append(qc)

    if not samples:
        console.print("[yellow][WARN] No QC files found. Check your --dir path.[/yellow]")
        sys.exit(0)

    if args.multi or len(samples) > 1:
        display_multi(samples)
    else:
        display_single(samples[0])

    if args.json:
        export_json(samples, args.json)
    if args.html:
        export_html(samples, args.html)


if __name__ == "__main__":
    import sys
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass
    main()
