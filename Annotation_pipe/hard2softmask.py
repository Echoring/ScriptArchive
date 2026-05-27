#!/usr/bin/env python3
# Usage: python3 hard2softmask.py genome.fasta hardmask.fasta
# Revert hardmasked genome to softmasked genome.
import sys
import numpy as np
from Bio import SeqIO

def process_with_numpy(genome_file, hardmask_file):
    genome_records = {rec.id: rec for rec in SeqIO.parse(genome_file, "fasta")}
    hardmask_records = {rec.id: rec for rec in SeqIO.parse(hardmask_file, "fasta")}
    
    missing_in_hardmask = set(genome_records.keys()) - set(hardmask_records.keys())
    if missing_in_hardmask:
        print(f"Warning: The following chromosomes are missing in hardmask file: {missing_in_hardmask}", file=sys.stderr)
    
    for chrom_id, genome_rec in genome_records.items():
        if chrom_id not in hardmask_records:
            print(f"Warning: Skipping {chrom_id} as it's not in hardmask file", file=sys.stderr)
            continue
            
        hardmask_rec = hardmask_records[chrom_id]
        
        genome_str = str(genome_rec.seq).upper()
        mask_str = str(hardmask_rec.seq).upper()
        
        if len(genome_str) != len(mask_str):
            print(f"Error: Length mismatch for {chrom_id} ({len(genome_str)} vs {len(mask_str)})", file=sys.stderr)
            continue
            
        genome_arr = np.frombuffer(bytearray(genome_str.encode('ascii')), dtype='u1')
        mask_arr = np.frombuffer(mask_str.encode('ascii'), dtype='u1')
        
        genome_arr[mask_arr == ord('N')] += 32
        
        print(f">{genome_rec.id}")
        print(genome_arr.tobytes().decode('ascii'))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 hard2softmask.py genome.fasta hardmask.fasta", file=sys.stderr)
        sys.exit(1)
    process_with_numpy(sys.argv[1], sys.argv[2])