#!/usr/bin/env python3
import sys
from Bio import Entrez

Entrez.email = "test@example.com"

def download_genome(accession, output_file):
    print(f"Downloading {accession} to {output_file}...")
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        with open(output_file, "w") as f:
            f.write(handle.read())
        handle.close()
        print("Download complete.")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        # Default for backward compatibility/testing
        accession = "NC_000913.3"
        output_file = "ecoli.fasta"
        print(f"Usage: python download_genome.py <accession> <output_file>")
        print(f"Using defaults: {accession} -> {output_file}")
    else:
        accession = sys.argv[1]
        output_file = sys.argv[2]

    download_genome(accession, output_file)
