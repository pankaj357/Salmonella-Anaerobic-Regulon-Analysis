#!/usr/bin/env python3
"""
Extract promoter sequences from genome FASTA using GFF annotation.

Outputs:
1. FASTA of promoters
2. CSV mapping promoters to gene coordinates

Author: Pankaj
"""

import sys
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def extract_promoters(genome_fasta, gff_file, output_prefix, upstream=300, downstream=50):
    # Load genome
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    
    promoters = []
    mapping = []

    with open(gff_file) as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attributes = parts
            if feature.lower() != "gene":
                continue
            start = int(start)
            end = int(end)
            # Parse gene name or locus_tag
            attr_dict = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, value = attr.split("=",1)
                    attr_dict[key.strip()] = value.strip()
                elif " " in attr:
                    key, value = attr.split(" ",1)
                    attr_dict[key.strip()] = value.strip()
            gene_name = attr_dict.get("gene", attr_dict.get("Name", attr_dict.get("locus_tag","NA")))
            locus_tag = attr_dict.get("locus_tag", "NA")

            # Extract promoter sequence
            seq = genome[chrom].seq
            if strand == "+":
                p_start = max(0, start - upstream -1)
                p_end = min(len(seq), start + downstream -1)
                promoter_seq = seq[p_start:p_end]
            else:
                p_start = max(0, end - downstream -1)
                p_end = min(len(seq), end + upstream -1)
                promoter_seq = seq[p_start:p_end].reverse_complement()

            header = f">{locus_tag}|{gene_name}|{chrom}:{start}-{end}|{strand}"
            promoters.append((header, str(promoter_seq)))
            mapping.append([locus_tag, gene_name, chrom, start, end, strand, p_start+1, p_end])

    # Write FASTA
    fasta_file = f"{output_prefix}_promoters.fasta"
    with open(fasta_file, "w") as f:
        for header, seq in promoters:
            f.write(f"{header}\n{seq}\n")

    # Write CSV
    df = pd.DataFrame(mapping, columns=["locus_tag","gene_name","chrom","start","end","strand","promoter_start","promoter_end"])
    csv_file = f"{output_prefix}_promoters.csv"
    df.to_csv(csv_file, index=False)

    print(f"Promoters written to: {fasta_file}")
    print(f"Mapping CSV written to: {csv_file}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python extract_promoters.py genome.fasta genome.gff output_prefix")
        sys.exit(1)
    genome_fasta = sys.argv[1]
    gff_file = sys.argv[2]
    output_prefix = sys.argv[3]
    extract_promoters(genome_fasta, gff_file, output_prefix)
