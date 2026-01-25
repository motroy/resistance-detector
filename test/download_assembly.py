#!/usr/bin/env python3
import sys
import os
import urllib.request
import gzip
import shutil
from Bio import Entrez

Entrez.email = "test@example.com"

def get_assembly_ftp(accession):
    """Get FTP link for assembly"""
    print(f"Searching for {accession}...")
    try:
        handle = Entrez.esearch(db="assembly", term=accession)
        record = Entrez.read(handle)
        handle.close()

        if not record['IdList']:
            print("Assembly not found")
            return None

        id = record['IdList'][0]
        handle = Entrez.esummary(db="assembly", id=id)
        summary = Entrez.read(handle)
        handle.close()

        ftp_path = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        return ftp_path

    except Exception as e:
        print(f"Error finding assembly: {e}")
        return None

def download_assembly(accession, output_file):
    ftp_path = get_assembly_ftp(accession)
    if not ftp_path:
        sys.exit(1)

    print(f"Found FTP path: {ftp_path}")

    file_name = ftp_path.split('/')[-1] + "_genomic.fna.gz"
    url = f"{ftp_path}/{file_name}"

    print(f"Downloading {url}...")

    try:
        urllib.request.urlretrieve(url, "temp.gz")

        print(f"Decompressing to {output_file}...")
        with gzip.open("temp.gz", 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove("temp.gz")
        print("Done.")

    except Exception as e:
        print(f"Error downloading/extracting: {e}")
        if os.path.exists("temp.gz"):
            os.remove("temp.gz")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python download_assembly.py <accession> <output_file>")
        sys.exit(1)

    download_assembly(sys.argv[1], sys.argv[2])
