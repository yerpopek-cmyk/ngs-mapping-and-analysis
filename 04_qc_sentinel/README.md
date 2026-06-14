# 🖥️ Project 4: QC Sentinel

> One command. Every metric. Colour-coded pass/warn/fail against GATK + ENCODE thresholds.

## What it does

`qc_sentinel.py` scans a results directory for samtools flagstat, Picard MarkDuplicates metrics,
and mosdepth summary files, then renders a unified Rich terminal dashboard — no more opening
five different files to understand whether your sample passed QC.

## Features

- **Auto-discovery**: scans `results/` for `*.flagstat.txt`, `*.markdup_metrics.txt`, `*.mosdepth.summary.txt`
- **Threshold-aware**: green/yellow/red per GATK Best Practices and ENCODE guidelines
- **Cohort mode**: multi-sample side-by-side comparison table
- **Export**: JSON (machine-readable) and HTML (shareable report)

## Quick start

```bash
pip install -r requirements.txt

# Single sample
python qc_sentinel.py --dir results/

# Cohort overview
python qc_sentinel.py --dir results/ --multi

# Export HTML report to share with collaborators
python qc_sentinel.py --dir results/ --html qc_report.html --json qc_summary.json
```

## Expected file structure in `--dir`

QC Sentinel searches recursively, so any of these layouts work:

```
results/
├── qc/
│   ├── sample1.flagstat.txt
│   ├── sample1.markdup_metrics.txt
│   └── sample1.mosdepth.summary.txt
└── ...
```

## Thresholds

| Metric | PASS | WARN | FAIL |
|--------|------|------|------|
| Mapping rate (WGS) | ≥90% | ≥75% | <75% |
| Properly paired | ≥85% | ≥70% | <70% |
| Duplication rate (WGS) | ≤20% | ≤35% | >35% |
| Mean coverage (WGS) | ≥25X | ≥15X | <15X |
| Callable at 20X | ≥85% | ≥70% | <70% |

Sources: [GATK Best Practices](https://gatk.broadinstitute.org) and [ENCODE QC metrics](https://www.encodeproject.org/about/experiment-guidelines/).

## Output examples

### Terminal dashboard
```
══════════════════════════════════════════════════
  QC Sentinel — sample1
══════════════════════════════════════════════════
 Alignment Metrics (samtools flagstat)
┌──────────────────────────┬──────────┬──────────┐
│ Metric                   │ Value    │ Status   │
├──────────────────────────┼──────────┼──────────┤
│ Total reads              │ 98,450,200│   —     │
│ Mapped reads             │ 97,103,847│ 98.6%   PASS │
│ Properly paired          │ 95,200,100│ 96.7%   PASS │
└──────────────────────────┴──────────┴──────────┘
```

### HTML report
A self-contained HTML file with colour-coded table, shareable via email or web.
