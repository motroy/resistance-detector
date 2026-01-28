#!/usr/bin/env python3
import sys
from Bio import Entrez
import urllib.request
import gzip
import shutil
import os

Entrez.email = "tool_test@example.com"

def download_assembly(accession):
    print(f"Searching for {accession}...")
    handle = Entrez.esearch(db="assembly", term=accession)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        print("Assembly not found.")
        sys.exit(1)

    uid = record["IdList"][0]
    print(f"Found UID: {uid}")

    handle = Entrez.esummary(db="assembly", id=uid, report="full")
    record = Entrez.read(handle)
    handle.close()

    try:
        summary = record['DocumentSummarySet']['DocumentSummary'][0]
        url = summary.get("FtpPath_RefSeq", "")
        if not url:
            url = summary.get("FtpPath_GenBank", "")
    except Exception as e:
        print(f"Error parsing summary: {e}")
        sys.exit(1)

    if not url:
        print("No FTP path found.")
        sys.exit(1)

    print(f"FTP URL: {url}")

    label = url.split("/")[-1]
    fasta_url = f"{url}/{label}_genomic.fna.gz"

    print(f"Downloading {fasta_url}...")
    filename = f"{accession}.fasta"

    try:
        with urllib.request.urlopen(fasta_url) as response:
            with gzip.GzipFile(fileobj=response) as uncompressed:
                with open(filename, 'wb') as out_file:
                    shutil.copyfileobj(uncompressed, out_file)
    except Exception as e:
        print(f"Download failed: {e}")
        sys.exit(1)

    print(f"Saved to {filename}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        accession = sys.argv[1]
    else:
        print("Usage: python download_assembly.py <accession>")
        sys.exit(1)
    download_assembly(accession)
