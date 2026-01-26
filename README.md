# FOS-CAZAVI Resistance Detector

A BLAST-based tool for detecting fosfomycin (FOS) and ceftazidime-avibactam (CAZAVI) resistance genes and mutations from bacterial genome assemblies.

## Features

- **Gene Detection**: Identifies resistance genes (fosA variants, blaKPC, blaOXA-48, etc.)
- **Mutation Detection**: Detects known resistance mutations (D179Y, V240G, T243M, etc.)
- **Sequence Extraction**: Outputs detected gene sequences to FASTA
- **Multiple Output Formats**: TSV results, human-readable summary, and raw BLAST output

## Repository Structure

```
resistance-detector/
├── resistance_detector.py      # Main detection tool
├── create_reference_database.py # Database builder
├── batch_analysis.sh           # Batch processing script
├── Snakefile                   # Snakemake workflow
├── config.yaml                 # Snakemake configuration
├── environment.yaml            # Conda environment
├── data/
│   ├── example_database.fasta  # Example reference database
│   └── primers.tsv             # Primer sequences for validation
└── README.md
```

## Installation

### Option 1: Conda (Recommended)

```bash
conda env create -f environment.yaml
conda activate resistance_detector
```

### Option 2: Manual Installation

```bash
# Install BLAST+
sudo apt-get install ncbi-blast+  # Ubuntu/Debian
# OR
brew install blast               # macOS

# Install Python dependencies
pip install biopython
```

## Quick Start

### 1. Run with Example Database

```bash
python resistance_detector.py \
    -a your_assembly.fasta \
    -d data/example_database.fasta \
    -o sample1
```

### 2. Build Complete Database (Recommended)

```bash
# Create database from NCBI
python create_reference_database.py \
    -e your.email@example.com \
    -o resistance_db.fasta

# Run detection
python resistance_detector.py \
    -a your_assembly.fasta \
    -d resistance_db.fasta \
    -o sample1
```

### 3. Check Results

```bash
cat sample1_summary.txt      # Human-readable summary
cat sample1_results.tsv      # Tab-delimited results
head sample1_genes.fasta     # Detected gene sequences
```

## Usage

### Main Tool: resistance_detector.py

```
usage: resistance_detector.py [-h] -a ASSEMBLY -d DATABASE -o OUTPUT
                             [--min_id MIN_ID] [--min_cov MIN_COV]

Required arguments:
  -a, --assembly    Input assembly file (FASTA)
  -d, --database    Resistance gene database (FASTA)
  -o, --output      Output prefix

Optional arguments:
  --min_id          Minimum percent identity (default: 90)
  --min_cov         Minimum percent coverage (default: 80)
```

### Database Builder: create_reference_database.py

```
usage: create_reference_database.py [-h] -e EMAIL [-o OUTPUT]

Required arguments:
  -e, --email       Email for NCBI Entrez

Optional arguments:
  -o, --output      Output file (default: resistance_genes.fasta)
```

### Batch Processing

```bash
./batch_analysis.sh -d resistance_db.fasta -i assemblies/ -o results/
```

## Output Files

| File | Description |
|------|-------------|
| `*_results.tsv` | Tab-delimited results (Contig, Gene, Identity%, Coverage%, Mutations) |
| `*_genes.fasta` | FASTA sequences of detected genes |
| `*_summary.txt` | Human-readable summary |
| `*_blast.txt` | Raw BLAST output |

## Detected Resistance Mechanisms

### Fosfomycin (FOS)

**Plasmid-mediated:**
- fosA3, fosA4, fosA5, fosA7, fosA11 (glutathione S-transferases)

**Chromosomal mutations:**
- uhpB, uhpC, uhpT, uhpA, glpT (transporters)
- galU, lon, cyaA, ptsI (regulatory genes)

### Ceftazidime-Avibactam (CAZAVI)

**KPC variants:**
- blaKPC-2, blaKPC-3, blaKPC-31, blaKPC-190
- Key mutations: D179Y, V240G, T243M

**OXA-48 variants:**
- blaOXA-48, blaOXA-181, blaOXA-232
- Key mutations: P68A, Y211S

**Chromosomal:**
- ompK35, ompK36 (porins)
- acrB, envZ, ftsI

## Snakemake Workflow

```bash
# Edit config.yaml with your samples
snakemake --cores 4
```

## Citation

If you use this tool, please cite:

- FOS resistance: PMID 33193186, PMID 28928846
- KPC mutations: PMID 28031275
- OXA-48 mutations: PMID 30700621

## License

MIT License
