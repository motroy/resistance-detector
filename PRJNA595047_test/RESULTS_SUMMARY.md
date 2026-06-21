# PRJNA595047 — FOS-CAZAVI Resistance Detection Results

Four *Klebsiella pneumoniae* assemblies from BioProject PRJNA1100451 (the assemblies bundled in
`PRJNA595047.ncbi_dataset.zip`; see [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11448024/)) were
run through `fos-cazavi fos-cazavi-all` using the bundled reference database
(`fos_cazavi/data/example_database.fasta` / `example_database_deduplicated.fasta` /
`example_database_mutations.tsv`) and default primers/genes.

## Genomes analyzed

| Accession | Strain (per `data_summary.tsv`) | Size | Contigs N50 |
|---|---|---|---|
| GCA_038433235.1 | novelKPC-MUT1 | 5,590,552 bp | 156,428 bp |
| GCA_038433245.1 | novelKPC-MUT2 | 5,567,630 bp | 141,399 bp |
| GCA_038433285.1 | KPC2-MUT2 | 5,590,586 bp | 129,255 bp |
| GCA_038433305.1 | KPC2-MUT1 | 5,589,487 bp | 141,077 bp |

## Command used

```
unzip PRJNA595047.ncbi_dataset.zip   # produces ncbi_dataset/data/<accession>/<accession>_*_genomic.fna

fos-cazavi fos-cazavi-all \
    -a ncbi_dataset/data/<accession>/<accession>_*_genomic.fna \
    -o <accession> \
    -d fos_cazavi/data/example_database.fasta \
    --genes fos_cazavi/data/example_database_deduplicated.fasta \
    --mutations fos_cazavi/data/example_database_mutations.tsv
```

(Required external tools — BLAST+, GAMMA, SeqKit, BLAT — were installed locally to run the pipeline;
they are not part of this repo and are not part of the committed results.)

## Key findings

All four isolates carry the same chromosomal background (**fosAKP** I91V, **acrB**, **envZ**, **ompK35**,
**ompK36**, **ftsI**, **blaSHV-12**), but differ sharply in **blaKPC-2**, exactly matching the strain
naming scheme (`novelKPC-MUT*` vs `KPC2-MUT*`):

| Gene | GCA_038433235.1 (novelKPC-MUT1) | GCA_038433245.1 (novelKPC-MUT2) | GCA_038433285.1 (KPC2-MUT2) | GCA_038433305.1 (KPC2-MUT1) |
|---|---|---|---|---|
| **blaKPC-2** | 96.60% id / 100% cov — **L166A, E167R, L168D, N169T, D178S** | 98.64% id / 100% cov — **L166S, E167A, L168I, N169P, D178P, V239A, T242Y** | 100.00% id / 100% cov — wild-type | 100.00% id / 100% cov — wild-type |
| fosAKP | 98.81% id / 100% cov — I91V | 98.81% id / 100% cov — I91V | 98.81% id / 100% cov — I91V | 98.81% id / 100% cov — I91V |
| ompK36 | 99.90% id / 91.58% cov — G134L, D135\*/N, G213T | 99.90% id / 91.30% cov — G134N, D135F | 99.90% id / 91.58% cov — G134L, D135\*/N, G213T | 99.90% id / 91.58% cov — G134L, D135\*/N, G213T |
| ompK35 | 99.44% id / 100% cov — D135G, D181R | same | same | same |
| acrB | 99.87% id / 100% cov — G617A, F626A, A628T/V | same | same | same |
| envZ | 99.85% id / 100% cov — G244S/D, T324Q | same | same | same |
| ftsI | 99.43% id / 100% cov — A333P, Y350L, S357Q | same | same | same |
| blaSHV-12 | 99.42% id / 100% cov — G238A, E240G | same | same | same |

**Interpretation:** `fos-cazavi` correctly separates the two strain backgrounds described by the assembly
metadata — the `KPC2-MUT*` isolates carry plain wild-type **blaKPC-2** (no mutations relative to the
reference, consistent with KPC-2 being the "MUT-free" comparator), while the `novelKPC-MUT*` isolates
carry distinct clusters of substitutions in the Ω-loop region of blaKPC (residues 166–178, plus V239A/T242Y
in MUT2) — the region associated with expanded-spectrum cephalosporin/ceftazidime-avibactam resistance in
KPC variants. The mutation calls are corroborated by both detection methods: GAMMA (protein-level alignment)
flags `blaKPC-2` mutations at 50% confidence on its own, and the independent SeqKit amplicon method
(PCR-primer-based extraction + sequence comparison) detects the *same* substitutions from the
`blaKPC_F`/`blaKPC_R` amplicon (852–882 bp), confirming the calls are not assembly/alignment artifacts
(see `*_unified_mutations.tsv` and `*_amplicons.tsv`).

No other acquired resistance genes (besides the chromosomal set above) were detected by BLAST in any of
the four genomes (`*_results.tsv` / `*_all_results.tsv`).

## Files in this folder

- `PRJNA595047.ncbi_dataset.zip` — original NCBI dataset download (unzip to get
  `ncbi_dataset/data/<accession>/<accession>_*_genomic.fna` assemblies plus
  `data_summary.tsv`/`assembly_data_report.jsonl` metadata; not kept unzipped here to avoid duplicating
  the ~22 MB of genome data already in the zip)
- `<accession>_*` — fos-cazavi outputs per genome: `_summary.txt` (human-readable report),
  `_results.tsv`/`_all_results.tsv` (BLAST gene detection), `_unified_mutations.tsv` (cross-method
  mutation confidence), `_genes.fasta` (extracted gene sequences), `_blast.txt` (raw BLAST output),
  `_gamma.gamma`/`_gamma.psl` (raw GAMMA output), `_amplicons.tsv`/`_seqkit_primers.tsv` (SeqKit amplicon
  detection)
