#!/bin/bash
set -e

# BLAST+
if ! command -v blastn &> /dev/null
then
    echo "Installing BLAST+..."
    wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
    tar -xzf ncbi-blast-2.16.0+-x64-linux.tar.gz
    export PATH=$PWD/ncbi-blast-2.16.0+/bin:$PATH
fi

# GAMMA (Gene Allele Mutation Microbial Assessment)
if ! command -v GAMMA.py &> /dev/null
then
    echo "Installing GAMMA..."
    pip install gamma-amr
fi

# GAMMA_DB_Maker (database preparation tool for GAMMA)
if [ ! -f GAMMA_DB_Maker.py ]
then
    echo "Downloading GAMMA_DB_Maker..."
    wget -q https://raw.githubusercontent.com/rastanton/GAMMA_DB_Maker/main/GAMMA_DB_Maker.py
fi

# SeqKit
if ! command -v seqkit &> /dev/null
then
    echo "Installing SeqKit..."
    wget -q https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz
    tar -xzf seqkit_linux_amd64.tar.gz
    export PATH=$PWD:$PATH
fi

echo "Dependencies check complete."
