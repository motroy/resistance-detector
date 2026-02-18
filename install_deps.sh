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

# Miniprot
if ! command -v miniprot &> /dev/null
then
    echo "Installing Miniprot..."
    git clone https://github.com/lh3/miniprot
    cd miniprot
    make
    export PATH=$PWD:$PATH
    cd ..
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
