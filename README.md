# FOS-CAZAVI Resistance Detector

A CLI tool for detecting fosfomycin (FOS) and ceftazidime-avibactam (CAZAVI) resistance genes and mutations from bacterial genome assemblies.

## Features

- **Modular CLI**: Separate commands for database creation, acquired gene detection, and mutation analysis.
- **Gene Detection**: Identifies resistance genes (fosA variants, blaKPC, blaOXA-48, etc.) using BLAST+.
- **Mutation Detection**: Detects known resistance mutations (D179Y, V240G, T243M, etc.) using GAMMA (protein-level gene alignment via translated nucleotide CDS) and SeqKit amplicon analysis.
- **Indel-aware mutation calling**: Both the BLAST and SeqKit mutation-detection paths walk a real gapped alignment (BLAST's own `qseq`/`sseq` for the BLAST path, a fresh `Bio.Align.PairwiseAligner` alignment against the reference gene for the SeqKit amplicon path) rather than doing a fixed-position lookup into the translated query. This means insertions/deletions relative to the reference gene are reported as e.g. `L166del` instead of producing fabricated, frame-shifted point-substitution calls at every downstream tracked position.
- **Amplicon Detection**: Uses seqkit amplicon to find PCR products from primer pairs and checks them for resistance mutations.
- **Sequence Extraction**: Outputs detected gene sequences to FASTA.
- **Multiple Output Formats**: TSV results, human-readable summary, and raw BLAST/GAMMA output.
- **Logging**: Comprehensive logging of commands, parameters, and tool versions.

## Repository Structure

```
resistance-detector/
├── fos_cazavi/              # Python package (installed via pip)
│   ├── __init__.py
│   ├── cli.py               # CLI entry point
│   ├── db.py                # Database creation module
│   ├── acquired.py          # BLAST-based detection module
│   ├── mutations.py         # GAMMA/Amplicon detection module
│   ├── utils.py             # Shared utilities
│   └── data/                # Bundled reference data
│       ├── example_database.fasta              # Raw nucleotide CDS sequences
│       ├── example_database_deduplicated.fasta # GAMMA-ready database (GAMMA_DB_Maker output)
│       ├── example_database_mutations.tsv      # Mutation definitions TSV
│       └── primers.tsv                         # Primer sequences for amplicon detection
├── fos-cazavi               # Executable script (for repo use without install)
├── pyproject.toml           # PyPI package configuration
├── GAMMA_DB_Maker.py        # GAMMA database preparation tool
├── create_test_genomes.py   # Test genome generator
├── batch_analysis.sh        # Batch processing script
├── Snakefile                # Snakemake workflow
├── config.yaml              # Snakemake configuration
├── environment.yaml         # Conda environment
├── scripts/
│   ├── download_assembly.py      # Download assemblies from NCBI
│   └── simulate_esbl_genomes.py  # Download real ESBL genomes and simulate resistance
├── tests/                   # Unit test suite (285 tests)
│   ├── conftest.py
│   ├── test_blast_parsing.py
│   ├── test_dual_detection.py
│   ├── test_gamma_parsing.py
│   ├── test_genome_creation.py
│   └── test_mutation_detection.py
├── example_results/         # Example outputs
├── bioproject_tests/        # Real-genome validation against published bioprojects
│   ├── PRJNA595047_test/        # 4 K. pneumoniae assemblies (NCBI),
│   │                            #   incl. comparison against a published in vitro selection study
│   ├── PRJNA741867_test_results/ # 6 K. pneumoniae ST307 assemblies (NCBI),
│   │                            #   incl. comparison against a published clinical CAZ/AVI-resistance study
│   └── PRJNA1086695_test/       # 2 myloasm-assembled isolates
├── LICENSE
└── README.md
```

## Installation

### Prerequisites

- **NCBI BLAST+**
- **GAMMA** (Gene Allele Mutation Microbial Assessment)
- **BLAT** (required by GAMMA for its protein-level alignment search)
- **GAMMA_DB_Maker** (for preparing nucleotide databases for GAMMA)
- **SeqKit**
- **Python 3** with **Biopython**

### Option 1: pip (Recommended)

```bash
pip install fos-cazavi
```

This installs the `fos-cazavi` command and bundles the reference data files (`example_database_deduplicated.fasta`, `primers.tsv`). System tools (BLAST+, GAMMA, seqkit) must still be installed separately.

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

## Quick Start

### 1. Create Reference Database

Download sequences from NCBI (AMRfinderPlus) and build the database. The `create-db` command fetches nucleotide CDS sequences and then runs GAMMA_DB_Maker to prepare a properly formatted GAMMA database:

```bash
fos-cazavi create-db -e your.email@example.com -o . -p resistance_db
```

This produces `resistance_db.fasta` (raw nucleotide CDS) and `resistance_db_deduplicated.fasta` (GAMMA-ready database) in the specified directory. Use `resistance_db_deduplicated.fasta` with `--genes`.

### 2. Run Full Analysis

Detect acquired genes, mutations, and amplicons (bundled genes and primers are used by default):

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
Creates the reference database and prepares it for GAMMA with GAMMA_DB_Maker.

```bash
fos-cazavi create-db -e <email> -o <output_dir> [-p <prefix>]
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
Detects mutations using GAMMA (protein-level alignment of nucleotide CDS) and SeqKit (amplicon detection).
`--genes` and `--primers` default to the bundled reference files.

```bash
fos-cazavi fos-cazavi-mutations \
    -a <assembly> \
    -o <output_prefix> \
    [--genes <genes.fasta>] \
    [--primers <primers.tsv>]
```

#### `fos-cazavi-all`
Runs the complete pipeline (acquired + mutations).
`--genes` and `--primers` default to the bundled reference files.

```bash
fos-cazavi fos-cazavi-all \
    -a <assembly> \
    -d <database> \
    -o <output_prefix> \
    [--genes <genes.fasta>] \
    [--primers <primers.tsv>]
```

## Output Files

| File | Description |
|------|-------------|
| `*_results.tsv` | Tab-delimited gene detection results (BLAST), including a `Copy_Number` column |
| `*_genes.fasta` | FASTA sequences of detected genes |
| `*_summary.txt` | Human-readable summary of all findings, including per-gene copy number |
| `*_summary.json` | Machine-parsable summary: per-gene copy number and loci, plus mutations/gene-alignments/amplicons |
| `*_summary.tsv` | Machine-parsable, one-row-per-gene summary with `Copy_Number`, loci, and max identity/coverage |
| `*_all_results.tsv` | Combined TSV of all detections (acquired genes, mutations, gene alignments, amplicons), including a `Copy_Number` column |
| `*_analysis.log` | Log of command, parameters, and tool versions |
| `*_blast.txt` | Raw BLAST output |
| `*_amplicons.tsv` | Amplicon detection results (with --primers) |
| `*_protein_mutations.tsv` | Protein mutation results (GAMMA) |
| `*_gamma_prefix.gamma` | Raw GAMMA output |
| `*_unified_mutations.tsv` | Unified dual-method mutation report |

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
  fosA3 (copy number: 1): 100.00% identity, 100.00% coverage

CEFTAZIDIME-AVIBACTAM RESISTANCE (KPC):
--------------------------------------------------
  blaKPC-3 (copy number: 1): 99.52% identity, 100.10% coverage
    Mutations: D179Y/N,V240Q,T243M

CEFTAZIDIME-AVIBACTAM RESISTANCE (OXA):
--------------------------------------------------
  blaOXA-48 (copy number: 1): 100.00% identity, 100.00% coverage
    Mutations: P68D,Y211A
```

### BLAST Results (`*_results.tsv`)

```tsv
Contig	Gene	Identity%	Coverage%	Mutations	Method	Copy_Number
contig_plasmid1_fosA3	fosA3	100.00	100.00	-	BLAST	1
contig_plasmid2_blaKPC3	blaKPC-3	99.52	100.10	D179Y/N,V240Q,T243M	BLAST	1
contig_plasmid3_blaOXA48	blaOXA-48	100.00	100.00	P68D,Y211A	BLAST	1
```

`Copy_Number` is the number of distinct genomic loci where that gene was detected (after redundancy filtering). A gene detected on two different contigs/loci will show `Copy_Number: 2` on both rows.

### Machine-Readable Summary (`*_summary.tsv` / `*_summary.json`)

`*_summary.tsv` aggregates results to one row per gene, making copy number easy to script against:

```tsv
Sample	Gene	Copy_Number	Loci	Max_Identity%	Max_Coverage%	Mutations
sample.fasta	fosA3	2	contig1:1-576;contig2:50-626	100.00	100.00	-
sample.fasta	blaKPC-3	1	contig3:2000-3000	99.52	100.10	D179Y/N,V240Q,T243M
```

`*_summary.json` contains the same per-gene aggregation (with full per-locus detail) plus mutation, gene-alignment, and amplicon results, for easy parsing by downstream scripts/pipelines.

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
collected 285 items

285 passed in 8.37s
```

**All 285 tests pass.** No external tools (BLAST, GAMMA, seqkit) are required to run the tests — all tool-dependent logic is tested via synthetic inputs.

### Test Coverage Summary

| Test Module | Tests | What Is Covered |
|---|---|---|
| `test_blast_parsing.py` | 23 | BLAST output parsing: identity/coverage filtering, gene name extraction from various ID formats, multi-hit handling, edge cases |
| `test_genome_creation.py` | 53 | Reference database loading (30 genes verified), `introduce_mutation()` utility, synthetic contig construction for all FOS and CAZAVI genome scenarios |
| `test_gamma_parsing.py` | 62 | GAMMA Codon_Changes field parsing for KPC (D179Y, V240G, T243M), OXA-48 (P68A, Y211S), CMY-178 (N70T), porins (OmpK35/36), AcrB, gene name normalization, and fosA variant names (fosA3/4/5/7/11) |
| `test_mutation_detection.py` | 105 | `detect_mutations()` for all gene families: wildtype no-call verification plus every documented resistance mutation across FOS and CAZAVI pathways |
| `test_dual_detection.py` | 42 | Dual-method (GAMMA + SeqKit amplicon) confidence scoring: 100% when both methods agree, 50% when only one detects a mutation |

### Validation Against Real, Published Genomes

Beyond the synthetic unit-test suite, the pipeline has been run end-to-end (BLAST+GAMMA+SeqKit, all external tools) against real NCBI assemblies tied to published studies, with results cross-checked against each paper's reported findings:

| Folder | Genomes | Paper compared against |
|---|---|---|
| `bioproject_tests/PRJNA595047_test/` | 4 *K. pneumoniae* assemblies | Pariona et al. 2024 (doi:10.1128/spectrum.01173-24) — *in vitro* meropenem-selected blaKPC reversion/Ω-loop deletion mutants |
| `bioproject_tests/PRJNA741867_test_results/` | 6 *K. pneumoniae* ST307 assemblies | Hernández-García et al. 2022 (JCM 60:e02245-21) — clinical ceftazidime-avibactam-selected blaKPC-46/-66/-92 X-loop variants |
| `bioproject_tests/PRJNA1086695_test/` | 2 myloasm-assembled isolates | — (assembly + detection validation only) |

These real-genome runs caught and validated the fix for an indel-misrepresentation bug: the BLAST and SeqKit mutation callers used to do a fixed-position lookup into the translated query with no sequence alignment, so any real insertion/deletion relative to the reference gene shifted every downstream "known mutation position" and produced fabricated point-substitution calls. See `bioproject_tests/PRJNA595047_test/RESULTS_SUMMARY.md` and `bioproject_tests/PRJNA741867_test_results/COMPARISON_TO_PAPER.md` for the full writeups, including before/after comparisons.

### Genes and Mutations Covered by Tests

**Fosfomycin (FOS) — Plasmidic:**

| Gene | Mutations Tested |
|---|---|
| fosA3, fosA4, fosA5, fosA7, fosA11 | K90E, H119Q (+ wildtype) |
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

The following NCBI assemblies are well-characterized and can be used to validate the full pipeline end-to-end (requires BLAST+, GAMMA, seqkit):

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

# 2. Create and prepare the GAMMA database (fetches real sequences from NCBI,
#    then runs GAMMA_DB_Maker to validate reading frames and deduplicate)
fos-cazavi create-db -e your.email@example.com -o resistance_db
makeblastdb -in resistance_db.fasta -dbtype nucl -out resistance_db

# 3. Run full analysis (bundled genes and primers used by default)
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
