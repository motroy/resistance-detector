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
├── fos_cazavi/              # Python package (installed via pip)
│   ├── __init__.py
│   ├── cli.py               # CLI entry point
│   ├── db.py                # Database creation module
│   ├── acquired.py          # BLAST-based detection module
│   ├── mutations.py         # Miniprot/Amplicon detection module
│   ├── utils.py             # Shared utilities
│   └── data/                # Bundled reference data
│       ├── example_database.fasta   # Example nucleotide database
│       ├── cazavi_proteins.fasta    # CAZAVI resistance proteins for miniprot
│       └── primers.tsv              # Primer sequences for amplicon detection
├── fos-cazavi               # Executable script (for repo use without install)
├── pyproject.toml           # PyPI package configuration
├── create_test_genomes.py   # Test genome generator
├── batch_analysis.sh        # Batch processing script
├── Snakefile                # Snakemake workflow
├── config.yaml              # Snakemake configuration
├── environment.yaml         # Conda environment
├── scripts/
│   ├── download_assembly.py      # Download assemblies from NCBI
│   └── simulate_esbl_genomes.py  # Download real ESBL genomes and simulate resistance
├── tests/                   # Unit test suite (194 tests)
│   ├── conftest.py
│   ├── test_blast_parsing.py
│   ├── test_genome_creation.py
│   ├── test_miniprot_parsing.py
│   └── test_mutation_detection.py
├── example_results/         # Example outputs
├── LICENSE
└── README.md
```

## Installation

### Prerequisites

- **NCBI BLAST+**
- **Miniprot**
- **SeqKit**
- **Python 3** with **Biopython**

### Option 1: pip (Recommended)

```bash
pip install fos-cazavi
```

This installs the `fos-cazavi` command and bundles the reference data files (`cazavi_proteins.fasta`, `primers.tsv`). System tools (BLAST+, miniprot, seqkit) must still be installed separately.

### Option 2: Conda

```bash
conda env create -f environment.yaml
conda activate resistance_detector
pip install -e .
```

### Option 3: Pixi

```bash
pixi install
```

### Option 4: Manual Installation

```bash
# Install system dependencies
bash install_deps.sh

# Install Python dependencies
pip install biopython

# Make script executable (for running directly from the repo)
chmod +x fos-cazavi
```

## Quick Start

### 1. Create Reference Database

Download sequences from NCBI (AMRfinderPlus) and build the database:

```bash
fos-cazavi create-db -e your.email@example.com -o resistance_db
```

### 2. Run Full Analysis

Detect acquired genes, mutations, and amplicons (bundled proteins and primers are used by default):

```bash
fos-cazavi fos-cazavi-all \
    -a your_assembly.fasta \
    -d resistance_db.fasta \
    -o results
```

## Usage

### Main Command: `fos-cazavi`

The tool is divided into subcommands:

#### `create-db`
Creates the reference database.

```bash
fos-cazavi create-db -e <email> -o <output_prefix>
```

#### `fos-cazavi-acquired`
Detects acquired resistance genes using BLAST.

```bash
fos-cazavi fos-cazavi-acquired \
    -a <assembly> \
    -d <database> \
    -o <output_prefix> \
    [--min_id 90] [--min_cov 80]
```

#### `fos-cazavi-mutations`
Detects mutations using Miniprot (protein alignment) and SeqKit (amplicon detection).
`--proteins` and `--primers` default to the bundled reference files.

```bash
fos-cazavi fos-cazavi-mutations \
    -a <assembly> \
    -o <output_prefix> \
    [--proteins <proteins.fasta>] \
    [--primers <primers.tsv>]
```

#### `fos-cazavi-all`
Runs the complete pipeline (acquired + mutations).
`--proteins` and `--primers` default to the bundled reference files.

```bash
fos-cazavi fos-cazavi-all \
    -a <assembly> \
    -d <database> \
    -o <output_prefix> \
    [--proteins <proteins.fasta>] \
    [--primers <primers.tsv>]
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

### Summary Output (`*_summary.txt`)

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

### BLAST Results (`*_results.tsv`)

```tsv
Contig	Gene	Identity%	Coverage%	Mutations	Method
contig_plasmid1_fosA3	fosA3	100.00	100.00	-	BLAST
contig_plasmid2_blaKPC3	blaKPC-3	99.52	100.10	D179Y/N,V240Q,T243M	BLAST
contig_plasmid3_blaOXA48	blaOXA-48	100.00	100.00	P68D,Y211A	BLAST
```

## Testing

### Running Tests

Install dependencies and run the full test suite:

```bash
pip install "fos-cazavi[dev]"
python -m pytest tests/ -v
```

### Test Results

```
============================= test session starts ==============================
platform linux -- Python 3.11.14, pytest-9.0.2
collected 194 items

194 passed in 8.85s
```

**All 194 tests pass.** No external tools (BLAST, miniprot, seqkit) are required to run the tests — all tool-dependent logic is tested via synthetic inputs.

### Test Coverage Summary

| Test Module | Tests | What Is Covered |
|---|---|---|
| `test_blast_parsing.py` | 22 | BLAST output parsing: identity/coverage filtering, gene name extraction from various ID formats, multi-hit handling, edge cases |
| `test_genome_creation.py` | 74 | Reference database loading (30 genes verified), `introduce_mutation()` utility, synthetic contig construction for all FOS and CAZAVI genome scenarios |
| `test_miniprot_parsing.py` | 20 | Miniprot CS-tag parsing for KPC (D179Y, V240G, T243M), OXA-48 (P68A, Y211S), CMY-178 (N70T), porins (OmpK35/36), AcrB, and CS-tag edge cases |
| `test_mutation_detection.py` | 78 | `detect_mutations()` for all 27 gene families: wildtype no-call verification plus every documented resistance mutation across FOS and CAZAVI pathways |

### Genes and Mutations Covered by Tests

**Fosfomycin (FOS) — Plasmidic:**

| Gene | Mutations Tested |
|---|---|
| fosA3, fosA4, fosA5, fosA11 | K90E, H119Q (+ wildtype) |
| fosAKP | I91V (+ wildtype) |

**Fosfomycin (FOS) — Chromosomal:**

| Gene | Mutations Tested |
|---|---|
| murA | D369N, L370I |
| uhpB | G469R, H350Y, H350Q |
| uhpC | F384L |
| uhpA | D54N, R139C, R139H |
| uhpT | G55D, W198\*, E258\*, W350\* |
| glpT | E44\*, W88\*, G90D, W234\*, R362C, R362H |
| cyaA | G463D, G463\* |
| ptsI | H191Y, H191Q |
| galU | R282V |
| lon | Q558\* |

**Ceftazidime-Avibactam (CAZAVI) — Plasmidic:**

| Gene | Mutations Tested |
|---|---|
| blaKPC-2, -3, -31, -190 | D179Y, V240G, T243M (incl. double mutants) |
| blaOXA-48 | P68A, Y211S (incl. double mutant) |
| blaCMY-178 | N70T |
| blaSHV-12 | G238S, E240K |

**Ceftazidime-Avibactam (CAZAVI) — Chromosomal:**

| Gene | Mutations Tested |
|---|---|
| ompK36 | G134D, G134\*, D135\*, G213D |
| ompK35 | G134D, D135\*, D181G, D181\* |
| acrB | G617D, G617N, F626L, A628T, A628V |
| mexR | W69G, W69\*, A75V |
| nalD | Q153\*, L174R |
| ftsI | A333V, A333T, Y350C, Y350S, S357N |
| envZ | G244S, G244D, T324I, T324A |

## Testing with Real Genomes

### Recommended Public Assemblies

The following NCBI assemblies are well-characterized and can be used to validate the full pipeline end-to-end (requires BLAST+, miniprot, seqkit):

| Accession | Organism | Resistance Genes | Use |
|---|---|---|---|
| GCA_000281535.1 | *K. pneumoniae* KPNIH1 | blaKPC-3 | CAZAVI positive control |
| GCF_000016305.1 | *K. pneumoniae* MGH 78578 | none | Negative control |
| GCF_000005845.2 | *E. coli* K-12 MG1655 | murA, uhpT, glpT (chromosomal) | Chromosomal FOS targets |

GCA_000281535.1 (*K. pneumoniae* KPNIH1) is from the 2011 NIH clinical outbreak and is confirmed to carry blaKPC-3. It is widely used as a reference for KPC carbapenemase studies.

### Downloading an Assembly

Use the provided script to download any NCBI assembly:

```bash
python scripts/download_assembly.py GCA_000281535.1
# Saves to GCA_000281535.1.fasta
```

### Running the Full Pipeline on a Real Genome

```bash
# 1. Download a KPC-carrying K. pneumoniae assembly
python scripts/download_assembly.py GCA_000281535.1

# 2. Create BLAST database (fetches real sequences from NCBI)
fos-cazavi create-db -e your.email@example.com -o resistance_db
makeblastdb -in resistance_db.fasta -dbtype nucl -out resistance_db

# 3. Run full analysis (bundled proteins and primers used by default)
fos-cazavi fos-cazavi-all \
    -a GCA_000281535.1.fasta \
    -d resistance_db \
    -o kpnih1_results

# 4. View summary
cat kpnih1_results_summary.txt
```

Expected output for KPNIH1 will include blaKPC-3 at ≥99% identity with CAZAVI resistance mutations reported.

### Generating Synthetic Test Genomes

If you do not have a real assembly, generate synthetic test genomes that embed all resistance scenarios:

```bash
python create_test_genomes.py  # produces test_genomes/ directory
```

This produces FASTA files in `test_genomes/` covering negative controls, single-gene plasmidic/chromosomal scenarios, and a multi-resistance genome (fosA3 + blaKPC-3 + blaOXA-48).

## Snakemake Workflow

The repository includes a Snakemake workflow for batch processing.

```bash
# Edit config.yaml with your samples
snakemake --cores 4
```

## License

MIT License
