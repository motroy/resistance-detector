"""
Snakemake workflow for FOS-CAZAVI resistance detection
Integrates with bacterial genomics pipelines

Author: Motroy
"""

import os
from pathlib import Path

# Configuration
configfile: "config.yaml"

# Input assemblies
ASSEMBLIES = config.get("assemblies", [])
SAMPLES = [Path(a).stem for a in ASSEMBLIES]

# Parameters
RESISTANCE_DB = config.get("resistance_db", "resistance_genes.fasta")
MIN_IDENTITY = config.get("min_identity", 90)
MIN_COVERAGE = config.get("min_coverage", 80)
OUTPUT_DIR = config.get("output_dir", "resistance_results")

# Rules
rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "{sample}_results.tsv"), sample=SAMPLES),
        expand(os.path.join(OUTPUT_DIR, "{sample}_genes.fasta"), sample=SAMPLES),
        os.path.join(OUTPUT_DIR, "combined_results.tsv"),
        os.path.join(OUTPUT_DIR, "summary_statistics.txt")

rule detect_resistance:
    input:
        assembly = lambda wildcards: [a for a in ASSEMBLIES if Path(a).stem == wildcards.sample][0],
        database = RESISTANCE_DB
    output:
        results = os.path.join(OUTPUT_DIR, "{sample}_results.tsv"),
        genes = os.path.join(OUTPUT_DIR, "{sample}_genes.fasta"),
        summary = os.path.join(OUTPUT_DIR, "{sample}_summary.txt"),
        blast = os.path.join(OUTPUT_DIR, "{sample}_blast.txt")
    params:
        prefix = os.path.join(OUTPUT_DIR, "{sample}"),
        min_id = MIN_IDENTITY,
        min_cov = MIN_COVERAGE
    log:
        os.path.join(OUTPUT_DIR, "logs", "{sample}.log")
    conda:
        "envs/resistance_detector.yaml"
    shell:
        """
        ./fos-cazavi fos-cazavi-acquired \
            -a {input.assembly} \
            -d {input.database} \
            -o {params.prefix} \
            --min_id {params.min_id} \
            --min_cov {params.min_cov} \
            2>&1 | tee {log}
        """

rule combine_results:
    input:
        expand(os.path.join(OUTPUT_DIR, "{sample}_results.tsv"), sample=SAMPLES)
    output:
        os.path.join(OUTPUT_DIR, "combined_results.tsv")
    run:
        with open(output[0], 'w') as outf:
            outf.write("Sample\tContig\tGene\tIdentity%\tCoverage%\tMutations\n")
            
            for result_file in input:
                sample = Path(result_file).stem.replace("_results", "")
                with open(result_file) as inf:
                    next(inf)  # Skip header
                    for line in inf:
                        outf.write(f"{sample}\t{line}")

rule summarize_results:
    input:
        os.path.join(OUTPUT_DIR, "combined_results.tsv")
    output:
        os.path.join(OUTPUT_DIR, "summary_statistics.txt")
    run:
        import pandas as pd
        from collections import Counter
        
        df = pd.read_csv(input[0], sep='\t')
        
        with open(output[0], 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("FOS-CAZAVI Resistance Detection Summary\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Total samples analyzed: {len(SAMPLES)}\n")
            f.write(f"Samples with resistance genes: {df['Sample'].nunique()}\n\n")
            
            # Gene frequency
            f.write("Gene frequency:\n")
            gene_counts = df['Gene'].value_counts()
            for gene, count in gene_counts.items():
                f.write(f"  {gene:<20} {count}\n")
            
            # Mutation frequency
            f.write("\nMutation frequency:\n")
            mutations = df[df['Mutations'] != '-']['Mutations'].str.split(',').explode()
            mut_counts = mutations.value_counts()
            for mut, count in mut_counts.items():
                f.write(f"  {mut:<20} {count}\n")
