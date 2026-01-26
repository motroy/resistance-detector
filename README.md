# FOS-CAZAVI Resistance Detector

A tool for detecting fosfomycin (FOS) and ceftazidime-avibactam (CAZAVI) resistance genes and mutations from whole genome sequencing data.

## Features

- **Gene Detection**: BLAST-based detection of resistance genes (fosA variants, blaKPC, blaOXA-48, etc.)
- **Mutation Detection**: Identifies known resistance mutations (D179Y, V240G, etc.)
- **Sequence Extraction**: Outputs detected genes to FASTA file
- **Comprehensive Reporting**: TSV results, human-readable summary, and raw BLAST output

## Installation

### Requirements

- Python 3.6+
- BLAST+ toolkit
- Biopython

### Install dependencies

```bash
# Ubuntu/Debian
sudo apt-get install ncbi-blast+
pip install biopython

# Conda
conda install -c bioconda blast biopython
```

## Quick Start

### 1. Create Reference Database

```bash
# Download and create database from NCBI
python create_reference_database.py -e your.email@example.com -o resistance_db.fasta

# OR use a pre-built database (see Database section below)
```

### 2. Run Detection

```bash
python resistance_detector.py \
    -a assembly.fasta \
    -d resistance_db.fasta \
    -o sample1
```

### 3. Check Results

```bash
# View summary
cat sample1_summary.txt

# View detailed results
cat sample1_results.tsv

# View detected gene sequences
head sample1_genes.fasta
```

## Usage

### Main Tool: resistance_detector.py

```bash
usage: resistance_detector.py [-h] -a ASSEMBLY -d DATABASE -o OUTPUT
                             [--min_id MIN_ID] [--min_cov MIN_COV] [-v]

Detect fosfomycin and ceftazidime-avibactam resistance genes and mutations

required arguments:
  -a, --assembly    Input assembly file (FASTA format)
  -d, --database    Resistance gene database (FASTA format)
  -o, --output      Output prefix for result files

optional arguments:
  --min_id         Minimum percent identity (default: 90)
  --min_cov        Minimum percent coverage (default: 80)
  -v, --version    Show version

examples:
  resistance_detector.py -a assembly.fasta -d resistance_db.fasta -o sample1
  resistance_detector.py -a assembly.fasta -d resistance_db.fasta -o sample1 --min_id 95
```

### Database Builder: create_reference_database.py

```bash
usage: create_reference_database.py [-h] -e EMAIL [-o OUTPUT]

Create reference database for FOS-CAZAVI resistance detection

required arguments:
  -e, --email      Email for NCBI Entrez (required by NCBI)

optional arguments:
  -o, --output     Output database file (default: resistance_genes.fasta)

example:
  create_reference_database.py -e your.email@example.com -o resistance_db.fasta
```

## Output Files

The tool generates four output files:

1. **`<prefix>_results.tsv`** - Tab-delimited results table
   ```
   Contig    Gene      Identity%    Coverage%    Mutations
   contig1   blaKPC-3  99.85        100.00       D179Y,T243M
   contig2   fosA3     100.00       99.50        -
   ```

2. **`<prefix>_genes.fasta`** - FASTA file with detected gene sequences
   ```
   >contig1_blaKPC-3 identity=99.85% coverage=100.00% mutations=D179Y,T243M
   ATGTCACTGTATCGCCTTCTCCTTATTGCTATTG...
   >contig2_fosA3 identity=100.00% coverage=99.50% mutations=-
   ATGAACATTGTGAAAATTATTGGGCACCAGTCTG...
   ```

3. **`<prefix>_summary.txt`** - Human-readable summary
4. **`<prefix>_blast.txt`** - Raw BLAST output (for debugging)

## Detected Resistance Mechanisms

### Fosfomycin Resistance

**Plasmid-mediated:**
- fosA3, fosA4, fosA5, fosA11 (glutathione S-transferases)

**Chromosomal mutations:**
- uhpB (G469R), uhpC (F384L) - transporters
- galU (R282V), lon (Q558*) - secondary mutations
- uhpT, glpT, uhpA, cyaA, ptsI - transporter/regulatory genes

### Ceftazidime-Avibactam Resistance

**KPC variants:**
- blaKPC-2, blaKPC-3, blaKPC-31, blaKPC-190
- Key mutations: D179Y, V240G, T243M (Ω-loop and 240-loop)

**OXA-48 variants:**
- blaOXA-48
- Key mutations: P68A, Y211S

**Chromosomal changes:**
- ompK36, acrB, envZ, ftsI (porins/efflux)

## Database Options

### Option 1: Auto-build from NCBI (Recommended)

```bash
python create_reference_database.py -e your.email@example.com -o resistance_db.fasta
```

### Option 2: Use Pre-built Databases

**CARD Database:**
```bash
wget https://card.mcmaster.ca/latest/data
tar -xvf data
# Extract relevant genes
```

**ResFinder Database:**
```bash
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
cat resfinder_db/*.fsa > resistance_db.fasta
```

### Option 3: Manual Curation

Create a FASTA file with your reference sequences:

```fasta
>fosA3_reference
ATGAACATTGTGAAAATTATTGGGCACCAGTCTGGCGCTGGCAAAACCACGCTGCTGAAC...

>blaKPC-3_reference
ATGTCACTGTATCGCCTTCTCCTTATTGCTATTG...

>blaOXA-48_reference  
ATGTCAAAGAATTATTTTCAGCGATTTCGTATATATTTTTTATG...
```

## Advanced Usage

### Batch Processing

```bash
# Process multiple assemblies
for assembly in assemblies/*.fasta; do
    sample=$(basename $assembly .fasta)
    python resistance_detector.py \
        -a $assembly \
        -d resistance_db.fasta \
        -o results/$sample
done

# Combine results
echo -e "Sample\tGene\tIdentity\tCoverage\tMutations" > all_results.tsv
for result in results/*_results.tsv; do
    sample=$(basename $result _results.tsv)
    tail -n +2 $result | awk -v s=$sample '{print s"\t"$0}' >> all_results.tsv
done
```

### Adjust Thresholds

```bash
# Strict detection (high identity/coverage)
python resistance_detector.py -a assembly.fasta -d resistance_db.fasta -o sample1 \
    --min_id 95 --min_cov 95

# Permissive detection (lower thresholds)
python resistance_detector.py -a assembly.fasta -d resistance_db.fasta -o sample1 \
    --min_id 80 --min_cov 70
```

### Integration with Other Tools

```bash
# Run after assembly with SPAdes
spades.py -1 R1.fastq -2 R2.fastq -o assembly_out
python resistance_detector.py -a assembly_out/contigs.fasta -d resistance_db.fasta -o sample1

# Run after assembly with MEGAHIT
megahit -1 R1.fastq -2 R2.fastq -o assembly_out
python resistance_detector.py -a assembly_out/final.contigs.fa -d resistance_db.fasta -o sample1
```

## Troubleshooting

### BLAST not found
```bash
# Check if BLAST is installed
blastn -version

# If not, install:
conda install -c bioconda blast
# OR
sudo apt-get install ncbi-blast+
```

### Biopython not found
```bash
pip install biopython
# OR
conda install -c conda-forge biopython
```

### No genes detected
- Check database file exists and is properly formatted
- Try lowering thresholds (--min_id 80 --min_cov 70)
- Verify assembly quality
- Check raw BLAST output: `<prefix>_blast.txt`

### Database creation fails
- Verify email is valid
- Check internet connection
- Some sequences may need manual addition

## Citation

If you use this tool, please cite the relevant papers for the resistance mechanisms:

**Fosfomycin:**
- PMID: 33193186 (uhpB, uhpC mutations)
- PMID: 28928846 (fosA variants)

**Ceftazidime-Avibactam:**
- PMID: 28031275 (KPC mutations)
- PMID: 30700621 (OXA-48 mutations)

## Author

Motroy

## License

MIT License - Free to use and modify
