"""
flagstat_parser.py — Parse samtools flagstat output

samtools flagstat output format (modern, samtools >= 1.13):

    98450200 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    0 + 0 supplementary
    0 + 0 duplicates
    97103847 + 0 mapped (98.63% : N/A)
    98450200 + 0 paired in sequencing
    49225100 + 0 read1
    49225100 + 0 read2
    95200100 + 0 properly paired (96.66% : N/A)
    96800000 + 0 with itself and mate mapped
    303847 + 0 singletons (0.31% : N/A)
    0 + 0 with mate mapped to a different chr
    0 + 0 with mate mapped to a different chr (mapQ>=5)
"""

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass
class FlagstatMetrics:
    total_reads:         int = 0
    secondary:           int = 0
    supplementary:       int = 0
    duplicates:          int = 0
    mapped_reads:        int = 0
    mapped_pct:          float = 0.0
    paired_in_sequencing: int = 0
    read1:               int = 0
    read2:               int = 0
    properly_paired:     int = 0
    properly_paired_pct: float = 0.0
    singletons:          int = 0
    singletons_pct:      float = 0.0


def parse_flagstat(path: str | Path) -> FlagstatMetrics:
    """
    Parse a samtools flagstat output file.

    Args:
        path: Path to the flagstat .txt file

    Returns:
        FlagstatMetrics dataclass with parsed values
    """
    path = Path(path)
    txt = path.read_text(encoding="utf-8")
    m = FlagstatMetrics()

    def _int(pattern: str) -> int:
        match = re.search(pattern, txt, re.MULTILINE)
        return int(match.group(1)) if match else 0

    def _pct(pattern: str) -> float:
        match = re.search(pattern, txt, re.MULTILINE)
        return float(match.group(1)) if match else 0.0

    m.total_reads          = _int(r"^(\d+) \+ \d+ in total")
    m.secondary            = _int(r"^(\d+) \+ \d+ secondary")
    m.supplementary        = _int(r"^(\d+) \+ \d+ supplementary")
    m.duplicates           = _int(r"^(\d+) \+ \d+ duplicates")
    m.mapped_reads         = _int(r"^(\d+) \+ \d+ mapped \(")
    m.mapped_pct           = _pct(r"mapped \(([\d.]+)%")
    m.paired_in_sequencing = _int(r"^(\d+) \+ \d+ paired in sequencing")
    m.read1                = _int(r"^(\d+) \+ \d+ read1")
    m.read2                = _int(r"^(\d+) \+ \d+ read2")
    m.properly_paired      = _int(r"^(\d+) \+ \d+ properly paired")
    m.properly_paired_pct  = _pct(r"properly paired \(([\d.]+)%")
    m.singletons           = _int(r"^(\d+) \+ \d+ singletons")
    m.singletons_pct       = _pct(r"singletons \(([\d.]+)%")

    return m


if __name__ == "__main__":
    import sys, json
    from dataclasses import asdict

    if len(sys.argv) != 2:
        print("Usage: python flagstat_parser.py <flagstat.txt>")
        sys.exit(1)

    result = parse_flagstat(sys.argv[1])
    print(json.dumps(asdict(result), indent=2))
