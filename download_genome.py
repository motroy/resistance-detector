from Bio import Entrez
import sys

Entrez.email = "test@example.com"
accession = "NC_000913.3"
output_file = "test/ecoli.fasta"

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
