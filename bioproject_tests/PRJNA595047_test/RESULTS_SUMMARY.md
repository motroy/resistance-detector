# PRJNA595047 — FOS-CAZAVI Resistance Detection Results

Four *Klebsiella pneumoniae* assemblies (see
[paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11448024/)), from an in vitro
selection experiment in which ceftazidime-avibactam-resistant mutants were
derived from a KPC-2-producing parent.

## Command used

```bash
fos-cazavi fos-cazavi-all \
    -a <accession>_genomic.fna \
    -o <accession> \
    -d fos_cazavi/data/example_database.fasta \
    --organism Klebsiella_pneumoniae
```

Reference data: AMRFinderPlus 2026-08-07.1.

## Results

| Accession | Strain (per `data_summary.tsv`) | blaKPC call | Changes (Ambler) | Predicted CAZ/AVI |
|---|---|---|---|---|
| GCA_038433235.1 | novelKPC-MUT1 | novel blaKPC variant | R164_I173del | **Resistant** |
| GCA_038433245.1 | novelKPC-MUT2 | novel blaKPC variant | W165_E168del | **Resistant** |
| GCA_038433285.1 | KPC2-MUT2 | blaKPC-2 | none | Susceptible |
| GCA_038433305.1 | KPC2-MUT1 | blaKPC-2 | none | Susceptible |

The two strains the study names `novelKPC-MUT*` each carry an in-frame deletion
inside the Omega loop (Ambler 164-179) of KPC: ten residues in MUT1, four in
MUT2. Neither deletion matches a named NCBI allele, so both are reported as
*novel blaKPC variants* with their changes listed, and called Resistant on the
mechanism — an in-frame Omega-loop deletion is the best-described route to
avibactam escape. The two `KPC2-MUT*` strains carry unmodified KPC-2 and are
called Susceptible.

All four carry CTX-M-15 and SHV-12 (neither of which survives avibactam) and
the intrinsic chromosomal `fosAKP`; none carries an acquired fosA-family
enzyme, so all four are predicted fosfomycin-susceptible.

## Why this result differs from the earlier committed run

The earlier pipeline had no way to name a KPC variant that is not in its
reference set, and its phenotype rule keyed on a short list of literal mutation
labels, so novel Omega-loop deletions did not produce a Resistant call. Variants
are now called against the canonical KPC-2 protein with indels reported as
indels, and the phenotype rule works from the structural mechanism as well as
from named alleles. See [../../docs/METHODS.md](../../docs/METHODS.md).
