# Validation Against Real, Published Genomes

Beyond the synthetic unit-test suite, the pipeline has been run end-to-end (BLAST+GAMMA+SeqKit, all external tools) against real NCBI assemblies tied to published studies, with results — including the genotype-derived FOS/CAZ-AVI phenotype predictions (see [OUTPUT_FILES.md](OUTPUT_FILES.md#predicted-phenotypes)) — cross-checked against each paper's reported findings:

| Folder | Genomes | Paper compared against |
|---|---|---|
| `bioproject_tests/PRJNA595047_test/` | 4 *K. pneumoniae* assemblies | Pariona et al. 2024 (doi:10.1128/spectrum.01173-24) — *in vitro* meropenem-selected blaKPC reversion/Ω-loop deletion mutants |
| `bioproject_tests/PRJNA741867_test_results/` | 6 *K. pneumoniae* ST307 assemblies | Hernández-García et al. 2022 (JCM 60:e02245-21) — clinical ceftazidime-avibactam-selected blaKPC-46/-66/-92 X-loop variants |
| `bioproject_tests/PRJNA1086695_test/` | 2 myloasm-assembled isolates | — (assembly + detection validation only) |

For PRJNA741867, the genotype-predicted CAZ/AVI phenotype matches the paper's reported susceptibility for all 6 isolates: KPC-3 (native, no tracked mutation) isolates A-1/B-1/C-1 predict Susceptible, while the KPC-46/-66/-92 X-loop variant isolates A-2/B-2/C-2 predict Resistant. For PRJNA595047, the `novelKPC-MUT1/2` Ω-loop deletion isolates predict Resistant while the wildtype `KPC2-MUT1/2` isolates predict Susceptible.

These real-genome runs caught and validated the fix for an indel-misrepresentation bug: the BLAST and SeqKit mutation callers used to do a fixed-position lookup into the translated query with no sequence alignment, so any real insertion/deletion relative to the reference gene shifted every downstream "known mutation position" and produced fabricated point-substitution calls. See `bioproject_tests/PRJNA595047_test/RESULTS_SUMMARY.md` and `bioproject_tests/PRJNA741867_test_results/COMPARISON_TO_PAPER.md` for the full writeups, including before/after comparisons.

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
