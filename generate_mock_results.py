#!/usr/bin/env python3
import gzip
from pathlib import Path

def write_flagstat(path: Path, total_reads: int, mapped_pct: float, properly_paired_pct: float):
    mapped_reads = int(total_reads * mapped_pct / 100)
    properly_paired = int(total_reads * properly_paired_pct / 100)
    singletons = int(total_reads * 0.0031)
    
    content = f"""{total_reads} + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
{mapped_reads} + 0 mapped ({mapped_pct:.2f}% : N/A)
{total_reads} + 0 paired in sequencing
{total_reads // 2} + 0 read1
{total_reads // 2} + 0 read2
{properly_paired} + 0 properly paired ({properly_paired_pct:.2f}% : N/A)
{int(total_reads * 0.98)} + 0 with itself and mate mapped
{singletons} + 0 singletons (0.31% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
"""
    path.write_text(content)

def write_markdup(path: Path, lib_name: str, read_pairs: int, dup_rate: float):
    pair_duplicates = int(read_pairs * dup_rate)
    optical_duplicates = int(pair_duplicates * 0.1)
    est_lib_size = int(read_pairs * (1 - dup_rate))
    
    content = f"""## METRICS CLASS\tpicard.sam.DuplicationMetrics
LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE
{lib_name}\t0\t{read_pairs}\t0\t1000\t0\t{pair_duplicates}\t{optical_duplicates}\t{dup_rate:.6f}\t{est_lib_size}
"""
    path.write_text(content)

def write_mosdepth(prefix: Path, mean_cov: float, pct_1x: float, pct_10x: float, pct_20x: float, pct_30x: float):
    # Summary file
    summary_content = f"""chrom\tlength\tbases\tmean\tmin\tmax
chr1\t248956422\t{int(248956422 * mean_cov)}\t{mean_cov:.2f}\t0\t245
total\t3088269832\t{int(3088269832 * mean_cov)}\t{mean_cov:.2f}\t0\t312
"""
    prefix.with_suffix(".mosdepth.summary.txt").write_text(summary_content)
    
    # Thresholds file
    # Format: chrom start end 1X 10X 20X 30X
    # Let's generate a gzipped thresholds bed file
    thresholds_path = prefix.with_suffix(".mosdepth.thresholds.bed.gz")
    total_len = 3088269832
    b_1x = int(total_len * pct_1x / 100)
    b_10x = int(total_len * pct_10x / 100)
    b_20x = int(total_len * pct_20x / 100)
    b_30x = int(total_len * pct_30x / 100)
    
    with gzip.open(thresholds_path, "wt") as fh:
        fh.write("#chrom\tstart\tend\t1X\t10X\t20X\t30X\n")
        fh.write(f"total\t0\t{total_len}\t{b_1x}\t{b_10x}\t{b_20x}\t{b_30x}\n")

def write_benchmark_tsv(path: Path):
    content = """aligner\truntime_sec\tmax_rss_mb\talignment_rate\tproperly_paired_pct\tmapped_reads\ttotal_reads
bwa-mem2\t120.5\t2400.0\t98.63\t96.66\t97103847\t98450200
minimap2\t45.2\t4800.0\t97.12\t94.20\t95614820\t98450200
hisat2\t85.7\t1800.0\t92.40\t88.10\t90967980\t98450200
"""
    path.write_text(content)

def main():
    out_dir = Path("mock_results")
    out_dir.mkdir(exist_ok=True)
    qc_dir = out_dir / "qc"
    qc_dir.mkdir(exist_ok=True)
    metrics_dir = out_dir / "metrics"
    metrics_dir.mkdir(exist_ok=True)
    
    # Write Sample 1 (Normal) - Pass QC
    write_flagstat(qc_dir / "sample_normal.flagstat.txt", 98450200, 98.63, 96.66)
    write_markdup(qc_dir / "sample_normal.markdup_metrics.txt", "normal_lib", 49225100, 0.100000)
    write_mosdepth(qc_dir / "sample_normal", 30.5, 99.8, 98.5, 95.2, 88.6)
    
    # Write Sample 2 (Tumour) - Warning QC (slightly higher duplication, lower coverage)
    write_flagstat(qc_dir / "sample_tumour.flagstat.txt", 92450100, 95.12, 91.20)
    write_markdup(qc_dir / "sample_tumour.markdup_metrics.txt", "tumour_lib", 46225050, 0.280000)
    write_mosdepth(qc_dir / "sample_tumour", 21.2, 98.5, 91.2, 78.4, 52.1)

    # Write Aligner Benchmark TSV
    write_benchmark_tsv(metrics_dir / "benchmark_raw.tsv")
    
    print("[OK] Mock results and benchmark data generated under mock_results/")

if __name__ == "__main__":
    import sys
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass
    main()
