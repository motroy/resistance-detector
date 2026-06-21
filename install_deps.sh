#!/bin/bash
set -e

BIN_DIR="${1:-$PWD}"
mkdir -p "$BIN_DIR"

# BLAST+
if ! command -v blastn &> /dev/null
then
    echo "Installing BLAST+..."
    wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
    tar -xzf ncbi-blast-2.16.0+-x64-linux.tar.gz
    export PATH=$PWD/ncbi-blast-2.16.0+/bin:$PATH
fi

# BLAT (required by GAMMA for protein-level alignment)
# Note: the file at .../linux.x86_64/blat is itself a directory containing the
# binary (also named "blat") - .../linux.x86_64/blat/blat is the real download.
if ! command -v blat &> /dev/null
then
    echo "Installing BLAT..."
    wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat -O "$BIN_DIR/blat"
    chmod +x "$BIN_DIR/blat"
    export PATH="$BIN_DIR:$PATH"
    if ! command -v blat &> /dev/null
    then
        echo "WARNING: BLAT install failed (UCSC mirror may be unreachable)." >&2
        echo "  Install manually, e.g. via conda: conda install -c bioconda ucsc-blat" >&2
    fi
fi

# GAMMA (Gene Allele Mutation Microbial Assessment)
# Not a PyPI package - install by fetching the script directly from GitHub.
if ! command -v GAMMA.py &> /dev/null
then
    echo "Installing GAMMA..."
    pip install -q unidecode biopython
    wget -q https://raw.githubusercontent.com/rastanton/GAMMA/master/GAMMA.py -O "$BIN_DIR/GAMMA.py"
    chmod +x "$BIN_DIR/GAMMA.py"
    export PATH="$BIN_DIR:$PATH"
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
echo "Make sure $BIN_DIR is on your PATH in future shells, e.g.:"
echo "  export PATH=\"$BIN_DIR:\$PATH\""
