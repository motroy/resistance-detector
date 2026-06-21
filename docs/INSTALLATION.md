# Installation

## Prerequisites

- **NCBI BLAST+**
- **GAMMA** (Gene Allele Mutation Microbial Assessment)
- **BLAT** (required by GAMMA for its protein-level alignment search)
- **GAMMA_DB_Maker** (for preparing nucleotide databases for GAMMA)
- **SeqKit**
- **Python 3** with **Biopython**

## Option 1: pip (Recommended)

```bash
pip install fos-cazavi
```

This installs the `fos-cazavi` command and bundles the reference data files (`example_database_deduplicated.fasta`, `primers.tsv`). System tools (BLAST+, GAMMA, seqkit) must still be installed separately.

## Option 2: Conda

```bash
conda env create -f environment.yaml
conda activate resistance_detector
pip install -e .
```

## Option 3: Pixi

```bash
pixi install
```

## Option 4: Manual Installation

```bash
# Install system dependencies (BLAST+, BLAT, GAMMA, SeqKit, GAMMA_DB_Maker)
# into ./bin (or pass a different directory as the first argument)
bash install_deps.sh ./bin
export PATH="$PWD/bin:$PATH"

# Install Python dependencies
pip install biopython

# Make script executable (for running directly from the repo)
chmod +x fos-cazavi
```

GAMMA is not distributed on PyPI; `install_deps.sh` fetches `GAMMA.py` directly from its
GitHub repository (rastanton/GAMMA) and `unidecode`/`biopython`, its Python dependencies. BLAT
is fetched from the UCSC binary mirror; if your network blocks that mirror, install it via
conda instead (`conda install -c bioconda ucsc-blat`).
