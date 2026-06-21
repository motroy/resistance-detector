# Output Files

| File | Description |
|------|-------------|
| `*_results.tsv` | Tab-delimited gene detection results (BLAST), including a `Copy_Number` column |
| `*_genes.fasta` | FASTA sequences of detected genes |
| `*_summary.txt` | Human-readable summary of all findings: predicted FOS/CAZ-AVI phenotypes, per-gene copy number, mutations |
| `*_summary.json` | Machine-parsable summary: predicted phenotypes, per-gene copy number and loci, plus mutations/gene-alignments/amplicons |
| `*_summary.tsv` | Machine-parsable, one-row-per-gene summary with `Copy_Number`, loci, and max identity/coverage (predicted phenotypes are included as leading `#`-prefixed comment lines) |
| `*_all_results.tsv` | Combined TSV of all detections (acquired genes, mutations, gene alignments, amplicons), including a `Copy_Number` column |
| `*_analysis.log` | Log of command, parameters, and tool versions |
| `*_blast.txt` | Raw BLAST output |
| `*_amplicons.tsv` | Amplicon detection results (with --primers) |
| `*_protein_mutations.tsv` | Protein mutation results (GAMMA) |
| `*_gamma_prefix.gamma` | Raw GAMMA output |
| `*_unified_mutations.tsv` | Unified dual-method mutation report |

## Predicted Phenotypes

Every summary output includes a genotype-derived, conclusive **Susceptible**/**Resistant** call for
each drug, plus the supporting evidence (which gene/mutation triggered the call):

- **Fosfomycin (FOS)**: Resistant if an acquired `fosA`-family enzyme (`fosA`, `fosA3/4/5/7/11`) is
  detected, or if a loss-of-function mutation (premature stop codon or in-frame deletion) is found in
  a fosfomycin uptake/regulatory gene (`murA, uhpT, uhpA, uhpB, uhpC, glpT, cyaA, ptsI, galU, lon`).
  The near-universal chromosomal `fosAKP` I91V variant is intrinsic and is *not* treated as a
  resistance signal.
- **Ceftazidime-Avibactam (CAZ/AVI)**: Resistant if a `blaKPC` hit carries a tracked Omega-loop/X-loop
  substitution (e.g. L166W, E167D, N169D, D179Y/N, V240G, T243M) or an in-frame indel in that region;
  wildtype `blaKPC` (no tracked mutation) is Susceptible.

This is implemented in `fos_cazavi/phenotype.py` and is a **genotype-based prediction only** — not a
substitute for phenotypic antimicrobial susceptibility testing (AST). See
[VALIDATION.md](VALIDATION.md) for cross-checks of these calls against published, phenotypically
characterized genomes.

## Example Results

### Summary Output (`*_summary.txt`)

```
======================================================================
FOS-CAZAVI Resistance Detection Summary
======================================================================

Assembly: test_genomes/ecoli_multi_resistance.fasta

PREDICTED PHENOTYPES (genotype-based):
--------------------------------------------------
  Fosfomycin (FOS): Resistant
    - Acquired fosfomycin-inactivating enzyme fosA3 detected (100.00% identity, 100.00% coverage)
  Ceftazidime-Avibactam (CAZ/AVI): Resistant
    - blaKPC variant blaKPC-3 carries Omega-loop/X-loop resistance marker(s): D179Y/N,V240Q,T243M
  Note: Genotype-based prediction only; not a substitute for phenotypic antimicrobial susceptibility testing (AST).

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

`*_summary.tsv` leads with `#`-prefixed phenotype comment lines, then aggregates results to one row per gene, making copy number easy to script against:

```tsv
# Predicted_Phenotype_Fosfomycin	Resistant	Acquired fosfomycin-inactivating enzyme fosA3 detected (100.00% identity, 100.00% coverage)
# Predicted_Phenotype_Ceftazidime_Avibactam	Resistant	blaKPC variant blaKPC-3 carries Omega-loop/X-loop resistance marker(s): D179Y/N,V240Q,T243M
# Disclaimer	Genotype-based prediction only; not a substitute for phenotypic antimicrobial susceptibility testing (AST).
Sample	Gene	Copy_Number	Loci	Max_Identity%	Max_Coverage%	Mutations
sample.fasta	fosA3	2	contig1:1-576;contig2:50-626	100.00	100.00	-
sample.fasta	blaKPC-3	1	contig3:2000-3000	99.52	100.10	D179Y/N,V240Q,T243M
```

`*_summary.json` contains the same `predicted_phenotypes` block plus the per-gene aggregation (with full per-locus detail), mutation, gene-alignment, and amplicon results, for easy parsing by downstream scripts/pipelines.
