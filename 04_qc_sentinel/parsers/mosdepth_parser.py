"""
mosdepth_parser.py — Parse mosdepth summary and thresholds output

mosdepth produces:
  - {prefix}.mosdepth.summary.txt        per-chromosome + total depth stats
  - {prefix}.mosdepth.global.dist.txt    cumulative distribution of depth
  - {prefix}.thresholds.bed.gz           (optional) bases at depth thresholds
  - {prefix}.regions.bed.gz              per-window depth (if -b/--by used)

summary.txt format:
    chrom	length	bases	mean	min	max
    chr1	248956422	7468692660	30.00	0	245
    ...
    total	3088269832	92648094960	30.00	0	312
"""

import gzip
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class MosdepthSummary:
    mean_coverage:     float = 0.0
    min_coverage:      float = 0.0
    max_coverage:      float = 0.0
    total_bases:       int = 0
    genome_length:     int = 0
    per_chrom:         dict = field(default_factory=dict)


@dataclass
class MosdepthThresholds:
    pct_1x:  float = 0.0
    pct_5x:  float = 0.0
    pct_10x: float = 0.0
    pct_20x: float = 0.0
    pct_30x: float = 0.0
    pct_50x: float = 0.0


def parse_mosdepth_summary(path: str | Path) -> MosdepthSummary:
    """
    Parse a {prefix}.mosdepth.summary.txt file.

    Args:
        path: Path to the mosdepth summary.txt file

    Returns:
        MosdepthSummary with overall and per-chromosome stats.
        The 'total' row (genome-wide, includes all reads) is preferred;
        falls back to a region-restricted summary if absent.
    """
    path = Path(path)
    result = MosdepthSummary()

    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return result

    header = lines[0].split("\t")
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) != len(header):
            continue
        row = dict(zip(header, parts))
        chrom = row.get("chrom", "")

        try:
            mean_val = float(row.get("mean", 0))
            length   = int(row.get("length", 0))
            bases    = int(row.get("bases", 0))
            min_val  = float(row.get("min", 0))
            max_val  = float(row.get("max", 0))
        except ValueError:
            continue

        if chrom in ("total", "total_region"):
            result.mean_coverage = mean_val
            result.genome_length = length
            result.total_bases   = bases
            result.min_coverage  = min_val
            result.max_coverage  = max_val
        else:
            result.per_chrom[chrom] = {
                "mean": mean_val, "length": length,
                "bases": bases, "min": min_val, "max": max_val,
            }

    return result


def parse_mosdepth_thresholds(path: str | Path) -> MosdepthThresholds:
    """
    Parse a {prefix}.thresholds.bed.gz file to compute % of genome
    covered at standard depth thresholds (requires mosdepth run with
    --thresholds 1,5,10,20,30,50).

    Returns:
        MosdepthThresholds with percentages of genome at each depth cutoff.
    """
    path = Path(path)
    result = MosdepthThresholds()

    if not path.exists():
        return result

    opener = gzip.open if str(path).endswith(".gz") else open

    totals = {}     # threshold_col -> sum of bases meeting threshold
    total_len = 0
    header = None

    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#chrom") or line.startswith("chrom"):
                header = line.strip().lstrip("#").split("\t")
                for col in header[4:]:
                    totals[col] = 0
                continue
            if header is None:
                continue

            parts = line.strip().split("\t")
            start, end = int(parts[1]), int(parts[2])
            region_len = end - start
            total_len += region_len

            for i, col in enumerate(header[4:], start=4):
                if i < len(parts):
                    try:
                        totals[col] += int(parts[i])
                    except ValueError:
                        pass

    if total_len == 0:
        return result

    def pct(col_name: str) -> float:
        for key in totals:
            if col_name in key:
                return totals[key] / total_len * 100
        return 0.0

    result.pct_1x  = pct("1X")
    result.pct_5x  = pct("5X")
    result.pct_10x = pct("10X")
    result.pct_20x = pct("20X")
    result.pct_30x = pct("30X")
    result.pct_50x = pct("50X")

    return result


if __name__ == "__main__":
    import sys, json
    from dataclasses import asdict

    if len(sys.argv) != 2:
        print("Usage: python mosdepth_parser.py <prefix.mosdepth.summary.txt>")
        sys.exit(1)

    result = parse_mosdepth_summary(sys.argv[1])
    print(json.dumps(asdict(result), indent=2))
