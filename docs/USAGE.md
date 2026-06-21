# Usage

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

## Main Command: `fos-cazavi`

The tool is divided into subcommands:

### `create-db`
Creates the reference database and prepares it for GAMMA with GAMMA_DB_Maker.

```bash
fos-cazavi create-db -e <email> -o <output_dir> [-p <prefix>]
```

### `fos-cazavi-acquired`
Detects acquired resistance genes using BLAST.

```bash
fos-cazavi fos-cazavi-acquired \
    -a <assembly> \
    -d <database> \
    -o <output_prefix> \
    [--min_id 90] [--min_cov 80]
```

### `fos-cazavi-mutations`
Detects mutations using GAMMA (protein-level alignment of nucleotide CDS) and SeqKit (amplicon detection).
`--genes` and `--primers` default to the bundled reference files.

```bash
fos-cazavi fos-cazavi-mutations \
    -a <assembly> \
    -o <output_prefix> \
    [--genes <genes.fasta>] \
    [--primers <primers.tsv>]
```

### `fos-cazavi-all`
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

See [OUTPUT_FILES.md](OUTPUT_FILES.md) for a description of the files each command produces.

## Snakemake Workflow

The repository includes a Snakemake workflow for batch processing.

```bash
# Edit config.yaml with your samples
snakemake --cores 4
```
