"""
picard_parser.py — Parse Picard MarkDuplicates metrics files

Picard metrics files have the structure:

    ## htsjdk.samtools.metrics.StringHeader
    # picard.sam.markduplicates.MarkDuplicates ...
    ## METRICS CLASS	picard.sam.DuplicationMetrics
    LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	...	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
    sample_lib	1234	49000000	...	0.123456	38500000

    ## HISTOGRAM ...
"""

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass
class PicardDupMetrics:
    library:                       str = ""
    unpaired_reads_examined:       int = 0
    read_pairs_examined:           int = 0
    secondary_or_supplementary:    int = 0
    unmapped_reads:                int = 0
    unpaired_read_duplicates:      int = 0
    read_pair_duplicates:          int = 0
    read_pair_optical_duplicates:  int = 0
    percent_duplication:           float = 0.0
    estimated_library_size:        int = 0

    @property
    def duplication_rate_pct(self) -> float:
        return self.percent_duplication * 100

    @property
    def pcr_duplicate_rate_pct(self) -> float:
        """PCR duplicates (excluding optical) as % of read pairs examined."""
        if self.read_pairs_examined == 0:
            return 0.0
        pcr_dups = self.read_pair_duplicates - self.read_pair_optical_duplicates
        return max(pcr_dups, 0) / self.read_pairs_examined * 100

    @property
    def optical_duplicate_rate_pct(self) -> float:
        if self.read_pairs_examined == 0:
            return 0.0
        return self.read_pair_optical_duplicates / self.read_pairs_examined * 100


def parse_picard_metrics(path: str | Path) -> PicardDupMetrics:
    """
    Parse a Picard MarkDuplicates metrics text file.

    Args:
        path: Path to the *.markdup_metrics.txt / *.duplicate_metrics.txt file

    Returns:
        PicardDupMetrics dataclass with parsed values.
        Returns a zero-valued instance if the METRICS table is not found.
    """
    path = Path(path)
    lines = [ln.rstrip("\n") for ln in path.read_text(encoding="utf-8").splitlines()]

    # Find the metrics header/data pair: a line starting with "LIBRARY"
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("LIBRARY"):
            header_idx = i
            break

    if header_idx is None or header_idx + 1 >= len(lines):
        return PicardDupMetrics()

    headers = lines[header_idx].split("\t")
    values  = lines[header_idx + 1].split("\t")
    row = dict(zip(headers, values))

    def _int(key: str) -> int:
        try:
            return int(float(row.get(key, 0) or 0))
        except ValueError:
            return 0

    def _float(key: str) -> float:
        try:
            return float(row.get(key, 0) or 0)
        except ValueError:
            return 0.0

    return PicardDupMetrics(
        library=row.get("LIBRARY", ""),
        unpaired_reads_examined=_int("UNPAIRED_READS_EXAMINED"),
        read_pairs_examined=_int("READ_PAIRS_EXAMINED"),
        secondary_or_supplementary=_int("SECONDARY_OR_SUPPLEMENTARY_RDS"),
        unmapped_reads=_int("UNMAPPED_READS"),
        unpaired_read_duplicates=_int("UNPAIRED_READ_DUPLICATES"),
        read_pair_duplicates=_int("READ_PAIR_DUPLICATES"),
        read_pair_optical_duplicates=_int("READ_PAIR_OPTICAL_DUPLICATES"),
        percent_duplication=_float("PERCENT_DUPLICATION"),
        estimated_library_size=_int("ESTIMATED_LIBRARY_SIZE"),
    )


if __name__ == "__main__":
    import sys, json
    from dataclasses import asdict

    if len(sys.argv) != 2:
        print("Usage: python picard_parser.py <markdup_metrics.txt>")
        sys.exit(1)

    result = parse_picard_metrics(sys.argv[1])
    out = asdict(result)
    out["duplication_rate_pct"] = result.duplication_rate_pct
    out["pcr_duplicate_rate_pct"] = result.pcr_duplicate_rate_pct
    out["optical_duplicate_rate_pct"] = result.optical_duplicate_rate_pct
    print(json.dumps(out, indent=2))
