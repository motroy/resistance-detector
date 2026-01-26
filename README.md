# FOS-CAZAVI Resistance Detector

A CLI tool for detecting fosfomycin (FOS) and ceftazidime-avibactam (CAZAVI) resistance genes and mutations from bacterial genome assemblies.

## Features

- **Modular CLI**: Separate commands for database creation, acquired gene detection, and mutation analysis.
- **Gene Detection**: Identifies resistance genes (fosA variants, blaKPC, blaOXA-48, etc.) using BLAST+.
- **Mutation Detection**: Detects known resistance mutations (D179Y, V240G, T243M, etc.) and analyzes protein mutations using miniprot.
- **Amplicon Detection**: Uses seqkit amplicon to find PCR products from primer pairs and checks them for resistance mutations.
- **Sequence Extraction**: Outputs detected gene sequences to FASTA.
- **Multiple Output Formats**: TSV results, human-readable summary, and raw BLAST/Miniprot output.
- **Logging**: Comprehensive logging of commands, parameters, and tool versions.

## Repository Structure

```
resistance-detector/
├── fos_cazavi/              # Python package
│   ├── __init__.py
│   ├── cli.py               # CLI entry point
│   ├── db.py                # Database creation module
│   ├── acquired.py          # BLAST-based detection module
│   ├── mutations.py         # Miniprot/Amplicon detection module
│   └── utils.py             # Shared utilities
├── fos-cazavi               # Executable script
├── create_test_genomes.py   # Test genome generator
├── batch_analysis.sh        # Batch processing script
├── Snakefile                # Snakemake workflow
├── config.yaml              # Snakemake configuration
├── environment.yaml         # Conda environment
├── data/
│   ├── example_database.fasta   # Reference nucleotide database
│   ├── cazavi_proteins.fasta    # CAZAVI resistance proteins for miniprot
│   └── primers.tsv              # Primer sequences for amplicon detection
├── example_results/         # Example outputs
└── README.md
```

## Installation

### Prerequisites

- **NCBI BLAST+**
- **Miniprot**
- **SeqKit**
- **Python 3** with **Biopython**

### Option 1: Conda (Recommended)

```bash
conda env create -f environment.yaml
conda activate resistance_detector
```

### Option 2: Manual Installation

```bash
# Install system dependencies
sudo apt-get install ncbi-blast+

# Install Miniprot and SeqKit (e.g., via brew or download binaries)
# brew install miniprot seqkit

# Install Python dependencies
pip install biopython

# Make script executable
chmod +x fos-cazavi
```

## Quick Start

### 1. Create Reference Database

Download sequences from NCBI (AMRfinderPlus) and build the database:

```bash
./fos-cazavi create-db -e your.email@example.com -o resistance_db
```

### 2. Run Full Analysis

Detect acquired genes, mutations, and amplicons:

```bash
./fos-cazavi fos-cazavi-all \
    -a your_assembly.fasta \
    -d resistance_db.fasta \
    -o results \
    --proteins data/cazavi_proteins.fasta \
    --primers data/primers.tsv
```

## Usage

### Main Command: `fos-cazavi`

The tool is divided into subcommands:

#### `create-db`
Creates the reference database.

```bash
./fos-cazavi create-db -e <email> -o <output_prefix>
```

#### `fos-cazavi-acquired`
Detects acquired resistance genes using BLAST.

```bash
./fos-cazavi fos-cazavi-acquired \
    -a <assembly> \
    -d <database> \
    -o <output_prefix> \
    [--min_id 90] [--min_cov 80]
```

#### `fos-cazavi-mutations`
Detects mutations using Miniprot (protein alignment) and SeqKit (amplicon detection).

```bash
./fos-cazavi fos-cazavi-mutations \
    -a <assembly> \
    -o <output_prefix> \
    --proteins <proteins.fasta> \
    --primers <primers.tsv>
```

#### `fos-cazavi-all`
Runs the complete pipeline (acquired + mutations).

```bash
./fos-cazavi fos-cazavi-all \
    -a <assembly> \
    -d <database> \
    -o <output_prefix> \
    --proteins <proteins.fasta> \
    --primers <primers.tsv>
```

## Output Files

| File | Description |
|------|-------------|
| `*_results.tsv` | Tab-delimited gene detection results (BLAST) |
| `*_genes.fasta` | FASTA sequences of detected genes |
| `*_summary.txt` | Human-readable summary of all findings |
| `*_analysis.log` | Log of command, parameters, and tool versions |
| `*_blast.txt` | Raw BLAST output |
| `*_amplicons.tsv` | Amplicon detection results (with --primers) |
| `*_protein_mutations.tsv` | Protein mutation results (with --proteins) |
| `*_miniprot.paf` | Raw miniprot PAF output (with --proteins) |

## Example Results

See `example_results/` folder for full file content.

### Summary Output (`results_summary.txt`)

```
======================================================================
FOS-CAZAVI Resistance Detection Summary
======================================================================

Assembly: test_genomes/ecoli_multi_resistance.fasta
Total genes detected: 3
Method: BLAST+

FOSFOMYCIN RESISTANCE GENES:
--------------------------------------------------
  fosA3: 100.00% identity, 100.00% coverage

CEFTAZIDIME-AVIBACTAM RESISTANCE (KPC):
--------------------------------------------------
  blaKPC-3: 99.52% identity, 100.10% coverage
    Mutations: D179Y/N,V240Q,T243M

CEFTAZIDIME-AVIBACTAM RESISTANCE (OXA):
--------------------------------------------------
  blaOXA-48: 100.00% identity, 100.00% coverage
    Mutations: P68D,Y211A
```

### BLAST Results (`results_results.tsv`)

```tsv
Contig	Gene	Identity%	Coverage%	Mutations	Method
contig_plasmid1_fosA3	fosA3	100.00	100.00	-	BLAST
contig_plasmid2_blaKPC3	blaKPC-3	99.52	100.10	D179Y/N,V240Q,T243M	BLAST
contig_plasmid3_blaOXA48	blaOXA-48	100.00	100.00	P68D,Y211A	BLAST
```

## Snakemake Workflow

The repository includes a Snakemake workflow for batch processing.

```bash
# Edit config.yaml with your samples
snakemake --cores 4
```

## License

MIT License
