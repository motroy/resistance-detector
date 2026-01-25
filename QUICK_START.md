# Quick Start Guide - FOS-CAZAVI Resistance Detector

## Installation (5 minutes)

```bash
# 1. Install dependencies
conda create -n resistance_detector python=3.9 blast biopython -c bioconda -c conda-forge
conda activate resistance_detector

# OR with pip
pip install biopython
# Install BLAST separately (system package manager)

# 2. Download the tools (or use the files you already have)
# All files are in this directory
```

## Quick Test (2 minutes)

```bash
# Run the example test
./test_example.sh

# This will:
# - Create a test assembly with resistance genes
# - Run the detector
# - Show results
```

## Basic Usage

### Option 1: Using the Example Database (Quick Test)

```bash
# Run on your assembly with the example database
python resistance_detector.py \
    -a your_assembly.fasta \
    -d example_database.fasta \
    -o your_sample
```

### Option 2: Build Complete Database (Recommended for Real Analysis)

```bash
# Step 1: Build database from NCBI
python create_reference_database.py \
    -e your.email@example.com \
    -o resistance_genes_complete.fasta

# Step 2: Run detection
python resistance_detector.py \
    -a your_assembly.fasta \
    -d resistance_genes_complete.fasta \
    -o your_sample
```

### Option 3: Batch Analysis (Multiple Samples)

```bash
# Create directory with assemblies
mkdir assemblies
cp assembly1.fasta assembly2.fasta assembly3.fasta assemblies/

# Run batch analysis
./batch_analysis.sh \
    -d resistance_genes.fasta \
    -i assemblies/ \
    -o results/
```

### Option 4: Snakemake Workflow (For Pipeline Integration)

```bash
# Edit config.yaml to list your assemblies
nano config.yaml

# Run workflow
snakemake --cores 4 --use-conda
```

## Output Files Explained

After running, you get 4 files:

1. **`sample_results.tsv`** - Main results table
   - Columns: Contig, Gene, Identity%, Coverage%, Mutations
   - Import into Excel or R for further analysis

2. **`sample_genes.fasta`** - Extracted gene sequences
   - Use for phylogenetic analysis
   - Confirm with PCR/Sanger sequencing
   - Submit to databases

3. **`sample_summary.txt`** - Human-readable summary
   - Quick overview
   - Share with collaborators

4. **`sample_blast.txt`** - Raw BLAST output
   - For debugging
   - Detailed alignment information

## Common Use Cases

### Use Case 1: Screen Clinical Isolates
```bash
# Quick screen of 10 clinical isolate assemblies
./batch_analysis.sh -d resistance_db.fasta -i clinical_isolates/ -o clinical_results/
```

### Use Case 2: Track Outbreak Strains
```bash
# Detect mutations in outbreak-associated KPC variants
python resistance_detector.py -a outbreak_strain.fasta -d resistance_db.fasta -o outbreak_001
grep "D179Y\|V240G\|T243M" outbreak_001_results.tsv
```

### Use Case 3: Publication-Ready Data
```bash
# High-stringency detection
python resistance_detector.py \
    -a assembly.fasta \
    -d resistance_db.fasta \
    -o sample \
    --min_id 95 \
    --min_cov 95

# Extract genes for phylogenetic tree
cat sample_genes.fasta
```

## Integrating with Your Existing Pipelines

### Integration with SPAdes Assembly Pipeline
```bash
# Add to your assembly script
spades.py -1 R1.fq -2 R2.fq -o spades_out
python resistance_detector.py \
    -a spades_out/contigs.fasta \
    -d resistance_db.fasta \
    -o ${sample_id}_resistance
```

### Integration with Prokka Annotation
```bash
# Run after Prokka
prokka --outdir prokka_out assembly.fasta
python resistance_detector.py \
    -a assembly.fasta \
    -d resistance_db.fasta \
    -o ${sample_id}_resistance

# Combine results
paste prokka_out/results.tsv ${sample_id}_resistance_results.tsv > combined.tsv
```

### Integration with Your Snakemake Pipelines
```python
# Add to your Snakefile
rule resistance_detection:
    input:
        assembly = "assembly/{sample}.fasta",
        db = "databases/resistance_genes.fasta"
    output:
        results = "resistance/{sample}_results.tsv",
        genes = "resistance/{sample}_genes.fasta"
    params:
        prefix = "resistance/{sample}"
    shell:
        "python resistance_detector.py -a {input.assembly} -d {input.db} -o {params.prefix}"
```

## Troubleshooting

**Problem: No genes detected**
```bash
# Solution 1: Lower thresholds
python resistance_detector.py -a assembly.fasta -d db.fasta -o sample --min_id 85 --min_cov 75

# Solution 2: Check assembly quality
# Assembly too fragmented? Try long-read assembly or hybrid approach
```

**Problem: BLAST not found**
```bash
# Install BLAST
conda install -c bioconda blast
# OR
sudo apt-get install ncbi-blast+
```

**Problem: Mutations not detected**
```bash
# Check the raw gene sequence
grep ">" sample_genes.fasta
# Manually inspect in Artemis/Geneious
# The gene might be truncated or have novel mutations
```

## Next Steps

1. ✅ Run test example
2. ✅ Build your database
3. ✅ Analyze your samples
4. ⭐ Validate interesting hits with PCR/Sanger
5. ⭐ Correlate with phenotypic resistance data
6. ⭐ Include in your publications

## Support

- Check README.md for detailed documentation
- Review the code - it's well-commented
- File issues or email: motroy@post.bgu.ac.il

## Citation

If you use this tool, please cite the resistance mechanism papers:
- PMID: 33193186 (FOS resistance)
- PMID: 28928846 (fosA variants)
- PMID: 28031275 (KPC mutations)
- PMID: 30700621 (OXA-48 mutations)
