"""
parsers — Standalone QC file parsers used by qc_sentinel.py

These are exposed separately so other scripts can import individual
parsers without pulling in the Rich dashboard dependencies.
"""

from .flagstat_parser import parse_flagstat
from .picard_parser import parse_picard_metrics
from .mosdepth_parser import parse_mosdepth_summary

__all__ = [
    "parse_flagstat",
    "parse_picard_metrics",
    "parse_mosdepth_summary",
]
