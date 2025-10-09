#!/usr/bin/env python3
"""
extract_all_training_promoters.py
Extract promoter sequences for all regulators based on training gene lists.
"""

from Bio import SeqIO
import glob
import os

# Input multi-strain promoter FASTA
fasta_file = "all_promoters.fasta"

# Get all *_training_genes.txt files
training_files = glob.glob("*_training_genes.txt")

# Output directory
output_dir = "training_promoters"
os.makedirs(output_dir, exist_ok=True)

for tf_file in training_files:
    regulator = tf_file.split("_training_genes.txt")[0]
    output_file = os.path.join(output_dir, f"{regulator}_promoters.fasta")
    
    # Read gene names
    with open(tf_file) as f:
        genes = set(line.strip() for line in f if line.strip())
    
    # Extract matching sequences
    count = 0
    with open(output_file, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            header_parts = record.id.split("|")
            gene_name = header_parts[1]  # second part is gene name
            if gene_name in genes:
                SeqIO.write(record, out, "fasta")
                count += 1

    print(f"[{regulator}] Extracted {count} promoters to {output_file}")

print("✅ All promoter sequences extracted!")
