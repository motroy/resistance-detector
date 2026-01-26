# FOS-CAZAVI Resistance Detector

A BLAST-based tool for detecting fosfomycin (FOS) and ceftazidime-avibactam (CAZAVI) resistance genes and mutations from bacterial genome assemblies.

## Features

- **Gene Detection**: Identifies resistance genes (fosA variants, blaKPC, blaOXA-48, etc.)
- **Mutation Detection**: Detects known resistance mutations (D179Y, V240G, T243M, etc.)
- **Amplicon Detection**: Uses seqkit amplicon to find PCR products from primer pairs
- **Sequence Extraction**: Outputs detected gene sequences to FASTA
- **Multiple Output Formats**: TSV results, human-readable summary, and raw BLAST output

## Repository Structure

```
resistance-detector/
├── resistance_detector.py       # Main detection tool
├── create_reference_database.py # Database builder
├── create_test_genomes.py       # Test genome generator
├── batch_analysis.sh            # Batch processing script
├── Snakefile                    # Snakemake workflow
├── config.yaml                  # Snakemake configuration
├── environment.yaml             # Conda environment
├── data/
│   ├── example_database.fasta   # Reference database
│   └── primers.tsv              # Primer sequences for amplicon detection
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
# Install BLAST+ and seqkit
sudo apt-get install ncbi-blast+  # Ubuntu/Debian
# OR
brew install blast seqkit         # macOS

# Install seqkit (if not via apt/brew)
conda install -c bioconda seqkit

# Install Python dependencies
pip install biopython
```

## Quick Start

```bash
# Basic detection (gene + mutation detection)
python resistance_detector.py \
    -a your_assembly.fasta \
    -d data/example_database.fasta \
    -o results

# With amplicon detection using primers
python resistance_detector.py \
    -a your_assembly.fasta \
    -d data/example_database.fasta \
    --primers data/primers.tsv \
    -o results

# View results
cat results_summary.txt
```

## Amplicon Detection with Primers

The tool uses **seqkit amplicon** to detect PCR products from primer pairs. This is useful for:
- Validating the presence of mutation sites
- Confirming gene regions targeted by diagnostic primers
- Checking if resistance-associated regions are amplifiable

### Primers File Format

The primers file (`data/primers.tsv`) contains primer pairs for detecting FOS resistance mutations:

| Primer | Nucleotide_sequence | Purpose | Mutation | Gene | Pair_ID |
|--------|---------------------|---------|----------|------|---------|
| uhpB_F | ACTGGGCGTCAGTAACGACG | Verification | - | uhpB | uhpB_ver |
| uhpB_R | ATGGCGCATCGGCAGGCGCT | Verification | - | uhpB | uhpB_ver |
| Mut3/uhpB_G469R/xbaI-F1p | CCGTCTAGAGGTGGATTTATTGCTCTCGCTG | Mutagenesis | uhpB G469R | uhpB | Mut3_uhpB_1 |

**Included primer pairs target:**
- uhpB G469R mutation (fosfomycin resistance)
- uhpC F384L mutation (fosfomycin resistance)
- galU R282V mutation (fosfomycin resistance)
- lon Q558* mutation (fosfomycin resistance)

### Running with Amplicon Detection

```bash
$ python resistance_detector.py \
    -a sample.fasta \
    -d data/example_database.fasta \
    --primers data/primers.tsv \
    -o sample_results

======================================================================
FOS-CAZAVI Resistance Detector
======================================================================
Running BLAST search (min_id=90.0%, min_cov=80.0%)...
Found 2 gene hits passing thresholds
Analyzing hits and detecting mutations...
Loading primers from data/primers.tsv...
Loaded 43 primers
Detecting amplicons with seqkit...
Found 5 amplicons
Checking for mutations within amplicons...
======================================================================
Analysis complete!
Results written to sample_results_*
======================================================================
```

### Amplicon Output Files

When using `--primers`, additional output files are generated:

| File | Description |
|------|-------------|
| `*_amplicons.tsv` | Amplicon detection results with coordinates |
| `*_seqkit_primers.tsv` | Primer pairs formatted for seqkit |

### Example Amplicon Results

```
AMPLICON DETECTION:
--------------------------------------------------
  Amplicon: uhpB_ver (1250 bp)
    Location: contig_5:125000-126250
    Primers: uhpB_F -> uhpB_R
    Mutations/Genes in amplicon:
      - uhpB (G469R detected)
```

## Example Usage and Results

### Running the Detector

```bash
$ python resistance_detector.py \
    -a sample_klebsiella.fasta \
    -d data/example_database.fasta \
    -o kp_sample

======================================================================
FOS-CAZAVI Resistance Detector
======================================================================
Running BLAST search (min_id=90.0%, min_cov=80.0%)...
Found 3 gene hits passing thresholds
Analyzing hits and detecting mutations...
Writing results to kp_sample_results.tsv...
Detected 3 resistance genes
Writing 3 gene sequences to kp_sample_genes.fasta...
Writing summary to kp_sample_summary.txt...
======================================================================
Analysis complete!
Results written to kp_sample_*
======================================================================
```

### Example Summary Output

```
======================================================================
FOS-CAZAVI Resistance Detection Summary
======================================================================

Assembly: sample_klebsiella.fasta
Total genes detected: 3

FOSFOMYCIN RESISTANCE GENES:
--------------------------------------------------
  fosA3: 100.00% identity, 100.00% coverage

CEFTAZIDIME-AVIBACTAM RESISTANCE (KPC):
--------------------------------------------------
  blaKPC-3: 99.52% identity, 100.10% coverage
    Mutations: D179Y/N,T243M

CEFTAZIDIME-AVIBACTAM RESISTANCE (OXA):
--------------------------------------------------
  blaOXA-48: 100.00% identity, 100.00% coverage
```

### Example TSV Output

| Contig | Gene | Identity | Coverage | Mutations |
|--------|------|----------|----------|-----------|
| contig_plasmid1 | fosA3 | 100.00 | 100.00 | - |
| contig_plasmid2 | blaKPC-3 | 99.52 | 100.10 | D179Y/N,T243M |
| contig_plasmid3 | blaOXA-48 | 100.00 | 100.00 | - |

### Testing with Synthetic Genomes

Generate test genomes with known mutations to validate the tool:

```bash
# Create test genomes
python create_test_genomes.py

# Run detection with amplicon analysis
python resistance_detector.py \
    -a test_genomes/ecoli_chromosomal_fos.fasta \
    -d data/example_database.fasta \
    --primers data/primers.tsv \
    -o test_chromosomal
```

**Validation Results:**

| Test Genome | Expected | Detected | Status |
|-------------|----------|----------|--------|
| ecoli_negative_control | None | 0 genes | PASS |
| ecoli_fosA3_positive | fosA3 | fosA3 (100%) | PASS |
| ecoli_blaKPC3_D179Y | blaKPC-3 + D179Y | blaKPC-3 (99.7%) + D179Y | PASS |
| ecoli_multi_resistance | 3 genes + mutations | fosA3, blaKPC-3 (D179Y,T243M), blaOXA-48 | PASS |

## Usage Reference

### Main Tool: resistance_detector.py

```
usage: resistance_detector.py [-h] -a ASSEMBLY -d DATABASE -o OUTPUT
                             [--min_id MIN_ID] [--min_cov MIN_COV]
                             [--primers PRIMERS]

Required arguments:
  -a, --assembly    Input assembly file (FASTA)
  -d, --database    Resistance gene database (FASTA)
  -o, --output      Output prefix

Optional arguments:
  --min_id          Minimum percent identity (default: 90)
  --min_cov         Minimum percent coverage (default: 80)
  --primers         Primer file for amplicon detection (TSV format)
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
| `*_amplicons.tsv` | Amplicon detection results (when using --primers) |

## Detected Resistance Mechanisms

### Fosfomycin (FOS)

**Plasmid-mediated:**
- fosA3, fosA4, fosA5, fosA7, fosA11 (glutathione S-transferases)

**Chromosomal mutations:**
- uhpB (G469R), uhpC (F384L), uhpT, uhpA, glpT (transporters)
- galU (R282V), lon (Q558*), cyaA, ptsI (regulatory genes)

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
- seqkit: doi:10.1371/journal.pone.0163962

## License

MIT License
