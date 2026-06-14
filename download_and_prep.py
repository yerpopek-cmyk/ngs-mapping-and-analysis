#!/usr/bin/env python3
import os
import gzip
import random
import urllib.request
from pathlib import Path

# Paths
ROOT_DIR = Path(__file__).parent
REF_DIR = ROOT_DIR / "reference"
DATA_DIR = ROOT_DIR / "data"

# URLs
CHR22_FA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
REFGENE_GTF_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz"
REFGENE_TXT_URL = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"

def download_file(url: str, dest_path: Path):
    print(f"Downloading {url} to {dest_path}...")
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Use urllib to stream download
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req) as response, open(dest_path, 'wb') as out_file:
        meta = response.info()
        file_size = int(meta.get("Content-Length", 0))
        print(f"File size: {file_size / (1024*1024):.2f} MB")
        
        downloaded = 0
        block_size = 8192
        while True:
            buffer = response.read(block_size)
            if not buffer:
                break
            downloaded += len(buffer)
            out_file.write(buffer)
            # Simple progress update every 5MB
            if downloaded % (5 * 1024 * 1024) < block_size:
                print(f"  Downloaded: {downloaded / (1024*1024):.2f} MB / {file_size / (1024*1024):.2f} MB")
    print(f"[OK] Download complete: {dest_path.name}")

def extract_chr22_fa(gz_path: Path, dest_path: Path):
    print(f"Decompressing {gz_path.name} to {dest_path.name}...")
    with gzip.open(gz_path, 'rt', encoding='utf-8') as f_in, open(dest_path, 'w', encoding='utf-8') as f_out:
        for line in f_in:
            f_out.write(line)
    print(f"[OK] Decompressed reference genome: {dest_path.name}")

def filter_chr22_gtf(gz_path: Path, dest_path: Path):
    print(f"Filtering {gz_path.name} for chr22 and saving to {dest_path.name}...")
    count = 0
    with gzip.open(gz_path, 'rt', encoding='utf-8') as f_in, open(dest_path, 'w', encoding='utf-8') as f_out:
        for line in f_in:
            # GTF lines start with the chromosome name
            if line.startswith("chr22\t"):
                f_out.write(line)
                count += 1
            elif line.startswith("#"):
                f_out.write(line)
    print(f"[OK] Extracted {count} features to {dest_path.name}")

def convert_refgene_to_bed12(gz_path: Path, dest_path: Path):
    print(f"Converting {gz_path.name} to BED12 and saving to {dest_path.name}...")
    count = 0
    with gzip.open(gz_path, 'rt', encoding='utf-8') as f_in, open(dest_path, 'w', encoding='utf-8') as f_out:
        for line in f_in:
            parts = line.strip().split('\t')
            if len(parts) < 15:
                continue
            
            chrom = parts[2]
            if chrom != "chr22":
                continue
                
            transcript_id = parts[1]
            strand = parts[3]
            txStart = int(parts[4])
            txEnd = int(parts[5])
            cdsStart = int(parts[6])
            cdsEnd = int(parts[7])
            exonCount = int(parts[8])
            
            starts_str = parts[9].strip(",").split(",")
            ends_str = parts[10].strip(",").split(",")
            gene_name = parts[12]
            
            if len(starts_str) != exonCount or len(ends_str) != exonCount:
                continue
                
            exon_starts = [int(x) for x in starts_str]
            exon_ends = [int(x) for x in ends_str]
            
            block_sizes = []
            block_starts = []
            for s, e in zip(exon_starts, exon_ends):
                block_sizes.append(str(e - s))
                block_starts.append(str(s - txStart))
                
            block_sizes_str = ",".join(block_sizes) + ","
            block_starts_str = ",".join(block_starts) + ","
            
            bed_name = f"{transcript_id}_{gene_name}"
            
            # BED12 Columns: chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
            bed_line = f"{chrom}\t{txStart}\t{txEnd}\t{bed_name}\t0\t{strand}\t{cdsStart}\t{cdsEnd}\t0\t{exonCount}\t{block_sizes_str}\t{block_starts_str}\n"
            f_out.write(bed_line)
            count += 1
            
    print(f"[OK] Converted {count} transcripts to BED12 at {dest_path.name}")

def simulate_reads_chr22(fasta_path: Path, out_r1: Path, out_r2: Path, num_pairs: int = 50000, read_len: int = 150):
    print(f"Loading reference {fasta_path.name} to memory for simulation...")
    with open(fasta_path, 'r', encoding='utf-8') as f:
        lines = f.read().splitlines()
    seq = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
    seq_len = len(seq)
    print(f"Loaded sequence of length {seq_len:,} bp")
    
    print(f"Simulating {num_pairs:,} read pairs...")
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    
    def reverse_complement(s: str) -> str:
        mapping = str.maketrans("ATCG", "TAGC")
        return s.translate(mapping)[::-1]
        
    random.seed(42) # Replicability
    
    with gzip.open(out_r1, 'wt', compresslevel=3) as f1, gzip.open(out_r2, 'wt', compresslevel=3) as f2:
        simulated = 0
        attempts = 0
        max_attempts = num_pairs * 10
        
        while simulated < num_pairs and attempts < max_attempts:
            attempts += 1
            frag_len = int(random.normalvariate(250, 50))
            frag_len = max(read_len + 10, min(frag_len, 450))
            
            start = random.randint(0, seq_len - frag_len - 1)
            frag = seq[start : start + frag_len]
            
            # Exclude sequences with high N count (e.g. > 2%)
            if frag.count("N") > frag_len * 0.02:
                continue
                
            r1_seq = frag[:read_len]
            r2_seq_raw = frag[-read_len:]
            r2_seq = reverse_complement(r2_seq_raw)
            
            # Sanger quality scores (Q40 = "I")
            qual = "I" * read_len
            
            # Write FASTQ
            f1.write(f"@SIM_CHR22_{simulated}/1\n{r1_seq}\n+\n{qual}\n")
            f2.write(f"@SIM_CHR22_{simulated}/2\n{r2_seq}\n+\n{qual}\n")
            
            simulated += 1
            if simulated % 10000 == 0:
                print(f"  Simulated {simulated:,} pairs...")
                
    print(f"[OK] Simulated {simulated:,} paired-end reads in fastq.gz format.")

def main():
    print("=== CHROMOSOME 22 DATASET PREPARATION ===")
    
    # 1. Download
    chr22_gz = REF_DIR / "chr22.fa.gz"
    refgene_gtf_gz = REF_DIR / "hg38.refGene.gtf.gz"
    refgene_txt_gz = REF_DIR / "refGene.txt.gz"
    
    if not chr22_gz.exists():
        download_file(CHR22_FA_URL, chr22_gz)
    if not refgene_gtf_gz.exists():
        download_file(REFGENE_GTF_URL, refgene_gtf_gz)
    if not refgene_txt_gz.exists():
        download_file(REFGENE_TXT_URL, refgene_txt_gz)
        
    # 2. Extract and format
    chr22_fa = REF_DIR / "chr22.fa"
    chr22_gtf = REF_DIR / "chr22.gtf"
    chr22_bed = REF_DIR / "chr22.bed"
    
    if not chr22_fa.exists():
        extract_chr22_fa(chr22_gz, chr22_fa)
    if not chr22_gtf.exists():
        filter_chr22_gtf(refgene_gtf_gz, chr22_gtf)
    if not chr22_bed.exists():
        convert_refgene_to_bed12(refgene_txt_gz, chr22_bed)
        
    # 3. Simulate Reads
    r1_fastq = DATA_DIR / "chr22_sample_R1.fastq.gz"
    r2_fastq = DATA_DIR / "chr22_sample_R2.fastq.gz"
    
    if not r1_fastq.exists() or not r2_fastq.exists():
        simulate_reads_chr22(chr22_fa, r1_fastq, r2_fastq)
        
    print("\n[COMPLETE] All Chromosome 22 files prepared successfully!")
    print(f"Reference:  {chr22_fa}")
    print(f"GTF:        {chr22_gtf}")
    print(f"BED:        {chr22_bed}")
    print(f"Read 1:     {r1_fastq}")
    print(f"Read 2:     {r2_fastq}")

if __name__ == "__main__":
    main()
