# FOS-CAZAVI Resistance Detector

A BLAST-based tool for detecting fosfomycin (FOS) and ceftazidime-avibactam (CAZAVI) resistance genes and mutations from bacterial genome assemblies.

## Features

- **Gene Detection**: Identifies resistance genes (fosA variants, blaKPC, blaOXA-48, etc.)
- **Mutation Detection**: Detects known resistance mutations (D179Y, V240G, T243M, etc.)
- **Protein Mutation Analysis**: Uses miniprot for porin/efflux pump mutation detection
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
│   ├── example_database.fasta   # Reference nucleotide database
│   ├── cazavi_proteins.fasta    # CAZAVI resistance proteins for miniprot
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
# Install BLAST+, seqkit, and miniprot
sudo apt-get install ncbi-blast+  # Ubuntu/Debian
# OR
brew install blast seqkit miniprot  # macOS

# Install via conda if not in apt/brew
conda install -c bioconda seqkit miniprot

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

# Full analysis with protein mutations and amplicon detection
python resistance_detector.py \
    -a your_assembly.fasta \
    -d data/example_database.fasta \
    --proteins data/cazavi_proteins.fasta \
    --primers data/primers.tsv \
    -o results

# View results
cat results_summary.txt
```

## Protein Mutation Analysis with Miniprot

The tool uses **miniprot** for protein-to-genome alignment to detect mutations in CAZAVI resistance-associated proteins. This method is based on [Xiong et al. 2025](https://doi.org/10.3389/fcimb.2025.1645042).

### Reference Proteins

The `data/cazavi_proteins.fasta` contains reference sequences for:

| Protein | Accession | Function | Resistance Association |
|---------|-----------|----------|------------------------|
| AcrA | WP_002892072.1 | Efflux pump membrane fusion | CAZAVI resistance |
| AcrB | WP_002892069.1 | Efflux RND transporter | CAZAVI resistance |
| OmpK35 | CAA09665.1 | Outer membrane porin | CAZAVI resistance |
| OmpK36 | ADG56549.1 | Outer membrane porin | CAZAVI resistance |

### Running with Protein Analysis

```bash
$ python resistance_detector.py \
    -a klebsiella_sample.fasta \
    -d data/example_database.fasta \
    --proteins data/cazavi_proteins.fasta \
    -o kp_results

======================================================================
FOS-CAZAVI Resistance Detector
======================================================================
Running BLAST search (min_id=90.0%, min_cov=80.0%)...
Found 2 gene hits passing thresholds
Analyzing hits and detecting mutations...
Running miniprot for protein mutation analysis...
Found 4 protein alignments
  OmpK35_CAA09665.1: 3 mutations detected
  OmpK36_ADG56549.1: 2 mutations detected
======================================================================
Analysis complete!
Results written to kp_results_*
======================================================================
```

### Miniprot Output

The `*_protein_mutations.tsv` file contains:

| Column | Description |
|--------|-------------|
| Protein | Reference protein name |
| Contig | Assembly contig with match |
| Start/End | Genomic coordinates |
| Identity | Percent amino acid identity |
| Coverage | Percent of protein aligned |
| Mutations | Detected amino acid changes (e.g., G134D, D135N) |
| CS_Tag | Raw miniprot cs tag for verification |

### Example Protein Mutation Results

```
PROTEIN MUTATION ANALYSIS (miniprot):
--------------------------------------------------
Reference: doi.org/10.3389/fcimb.2025.1645042

Outer Membrane Porins:
  OmpK35_CAA09665.1: 98.5% identity, 100.0% coverage
    Location: contig_3:125000-126200
    Mutations: G134D, D135N, ins136_GD

  OmpK36_ADG56549.1: 95.2% identity, 98.5% coverage
    Location: contig_3:128000-129100
    Mutations: D135*, del213_2aa

Efflux Pump Components:
  acrB_WP_002892069.1: 99.8% identity, 100.0% coverage
    Location: contig_5:450000-453200
    No mutations detected
```

### CS Tag Interpretation

Miniprot's cs tag encodes mutations:
- `:N` - N identical amino acids
- `*XY` - Substitution (ref X → query Y)
- `+AA` - Insertion of amino acids
- `-nt` - Deletion of nucleotides
- Frameshift and intron operators for complex events

## Amplicon Detection with Primers

The tool uses **seqkit amplicon** to detect PCR products from primer pairs.

### Running with Amplicon Detection

```bash
python resistance_detector.py \
    -a sample.fasta \
    -d data/example_database.fasta \
    --primers data/primers.tsv \
    -o sample_results
```

### Primers File Format

| Primer | Nucleotide_sequence | Purpose | Mutation | Gene | Pair_ID |
|--------|---------------------|---------|----------|------|---------|
| uhpB_F | ACTGGGCGTCAGTAACGACG | Verification | - | uhpB | uhpB_ver |
| uhpB_R | ATGGCGCATCGGCAGGCGCT | Verification | - | uhpB | uhpB_ver |

**Included primer pairs target:**
- uhpB G469R, uhpC F384L, galU R282V, lon Q558* (fosfomycin resistance)

## Example Usage and Results

### Full Analysis

```bash
$ python resistance_detector.py \
    -a klebsiella_sample.fasta \
    -d data/example_database.fasta \
    --proteins data/cazavi_proteins.fasta \
    --primers data/primers.tsv \
    -o kp_full

======================================================================
FOS-CAZAVI Resistance Detector
======================================================================
Running BLAST search (min_id=90.0%, min_cov=80.0%)...
Found 3 gene hits passing thresholds
Analyzing hits and detecting mutations...
Loading primers from data/primers.tsv...
Loaded 42 primers
Detecting amplicons with seqkit...
Found 2 amplicons
Running miniprot for protein mutation analysis...
Found 4 protein alignments
  OmpK36_ADG56549.1: 2 mutations detected
======================================================================
Analysis complete!
Results written to kp_full_*
======================================================================
```

### Example Summary Output

```
======================================================================
FOS-CAZAVI Resistance Detection Summary
======================================================================

Assembly: klebsiella_sample.fasta
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

PROTEIN MUTATION ANALYSIS (miniprot):
--------------------------------------------------
Reference: doi.org/10.3389/fcimb.2025.1645042

Outer Membrane Porins:
  OmpK36_ADG56549.1: 95.2% identity, 98.5% coverage
    Location: contig_3:128000-129100
    Mutations: G134D, D135N
```

## Usage Reference

### Main Tool: resistance_detector.py

```
usage: resistance_detector.py [-h] -a ASSEMBLY -d DATABASE -o OUTPUT
                             [--min_id MIN_ID] [--min_cov MIN_COV]
                             [--primers PRIMERS] [--proteins PROTEINS]

Required arguments:
  -a, --assembly    Input assembly file (FASTA)
  -d, --database    Resistance gene database (FASTA)
  -o, --output      Output prefix

Optional arguments:
  --min_id          Minimum percent identity (default: 90)
  --min_cov         Minimum percent coverage (default: 80)
  --primers         Primer file for amplicon detection (TSV)
  --proteins        Protein sequences for miniprot analysis (FASTA)
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
| `*_results.tsv` | Tab-delimited gene detection results |
| `*_genes.fasta` | FASTA sequences of detected genes |
| `*_summary.txt` | Human-readable summary |
| `*_blast.txt` | Raw BLAST output |
| `*_amplicons.tsv` | Amplicon detection results (with --primers) |
| `*_protein_mutations.tsv` | Protein mutation results (with --proteins) |
| `*_miniprot.paf` | Raw miniprot PAF output (with --proteins) |

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

**Porin mutations (detected via miniprot):**
- OmpK35: Truncations, insertions, loop mutations
- OmpK36: G134D, D135N, GD/TD loop insertions

**Efflux pump mutations (detected via miniprot):**
- AcrA, AcrB: Overexpression-associated mutations

## Snakemake Workflow

```bash
# Edit config.yaml with your samples
snakemake --cores 4
```

## Citation

If you use this tool, please cite:

- Miniprot method: [Xiong et al. 2025](https://doi.org/10.3389/fcimb.2025.1645042)
- miniprot: [Li H. 2023](https://doi.org/10.1093/bioinformatics/btad014)
- seqkit: [Shen et al. 2016](https://doi.org/10.1371/journal.pone.0163962)
- FOS resistance: PMID 33193186, PMID 28928846
- KPC mutations: PMID 28031275
- OXA-48 mutations: PMID 30700621

## License

MIT License
