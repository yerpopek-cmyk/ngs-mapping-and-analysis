#!/usr/bin/env python3
import gzip
import random
import sys
from pathlib import Path

def reverse_complement(seq: str) -> str:
    mapping = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(mapping)[::-1]

def simulate_reads(fasta_path: Path, out_r1: Path, out_r2: Path, num_pairs: int = 10000, read_len: int = 150):
    # Read FASTA
    lines = fasta_path.read_text(encoding="utf-8").splitlines()
    seq = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
    seq_len = len(seq)
    
    print(f"Loaded reference sequence of length {seq_len} bp")
    
    # Open gzip output files
    with gzip.open(out_r1, "wt") as f1, gzip.open(out_r2, "wt") as f2:
        for idx in range(1, num_pairs + 1):
            # Normal distribution for fragment length (mean=250, std=50)
            frag_len = int(random.normalvariate(250, 50))
            frag_len = max(read_len + 10, min(frag_len, 450)) # bounds
            
            # Pick a random starting position on the circular/linear genome
            if seq_len <= frag_len:
                # Genome is shorter than fragment length (common for PhiX under some fragment size distributions)
                # Just adjust fragment length to genome length
                frag_len = seq_len - 1
            
            start = random.randint(0, seq_len - frag_len - 1)
            
            # Extract reads
            r1_seq = seq[start : start + read_len]
            r2_seq_raw = seq[start + frag_len - read_len : start + frag_len]
            r2_seq = reverse_complement(r2_seq_raw)
            
            # Simple sequencing quality scores (Sanger Q40: "I")
            qual = "I" * read_len
            
            # Write FASTQ format
            # Read 1
            f1.write(f"@SIM_{idx}/1\n{r1_seq}\n+\n{qual}\n")
            # Read 2
            f2.write(f"@SIM_{idx}/2\n{r2_seq}\n+\n{qual}\n")
            
    print(f"Simulated {num_pairs} read pairs saved to:\n  - {out_r1}\n  - {out_r2}")

def main():
    fasta_path = Path("phix.fa")
    out_r1 = Path("phix_R1.fastq.gz")
    out_r2 = Path("phix_R2.fastq.gz")
    simulate_reads(fasta_path, out_r1, out_r2)

if __name__ == "__main__":
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass
    main()
