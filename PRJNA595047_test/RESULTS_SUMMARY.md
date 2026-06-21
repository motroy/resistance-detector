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
| **blaKPC-2** | 96.60% id / 100% cov — **R163_I172 deletion** (10 aa) | 98.64% id / 100% cov — **W164_E167 deletion** (4 aa) | 100.00% id / 100% cov — wild-type | 100.00% id / 100% cov — wild-type |
| fosAKP | 98.81% id / 100% cov — I91V | 98.81% id / 100% cov — I91V | 98.81% id / 100% cov — I91V | 98.81% id / 100% cov — I91V |
| ompK36 | 99.90% id / 91.58% cov — G213I (contig-edge truncation, no internal mutation) | 99.90% id / 91.30% cov — G213I (contig-edge truncation) | 99.90% id / 91.58% cov — G213I (contig-edge truncation) | 99.90% id / 91.58% cov — G213I (contig-edge truncation) |
| ompK35 | 99.44% id / 100% cov — D135G, D181R | same | same | same |
| acrB | 99.87% id / 100% cov — G617A, F626A, A628T/V | same | same | same |
| envZ | 99.85% id / 100% cov — G244S/D, T324Q | same | same | same |
| ftsI | 99.43% id / 100% cov — A333P, Y350L, S357Q | same | same | same |
| blaSHV-12 | 99.42% id / 100% cov — G238A, E240G | same | same | same |

**Interpretation:** `fos-cazavi` correctly separates the two strain backgrounds described by the assembly
metadata — the `KPC2-MUT*` isolates carry plain wild-type **blaKPC-2** (no mutations relative to the
reference, consistent with KPC-2 being the "MUT-free" comparator), while the `novelKPC-MUT*` isolates
carry small in-frame deletions in the Ω-loop region of blaKPC (residues 163–172 in MUT1, 164–167 in
MUT2) — the region associated with expanded-spectrum cephalosporin/ceftazidime-avibactam resistance in
KPC variants. This is corroborated by GAMMA's independent protein-level alignment, which classifies both
hits as `Match_Type=Indel` (vs. `Native` for the two wild-type genomes) — see `*_protein_mutations.tsv`.

No other acquired resistance genes (besides the chromosomal set above) were detected by BLAST in any of
the four genomes (`*_results.tsv` / `*_all_results.tsv`).

## Bug found and fixed: indels misreported as point substitutions

The BLAST-based mutation caller (`fos_cazavi/utils.py::detect_mutations()`, used from
`fos_cazavi/acquired.py`) did a **fixed-position lookup** into the translated query protein
(`protein[pos-1]`) with no sequence alignment. Any real insertion/deletion in the query relative to the
reference gene shifts every downstream "known mutation position" by the indel length, so the function
read amino acids from the wrong (frame-shifted) codon and fabricated plausible-looking but false
substitution calls.

This was caught on these exact genomes: GCA_038433235.1 and GCA_038433245.1 (`novelKPC-MUT1/2`) actually
carry small in-frame **deletions** in blaKPC-2 (confirmed independently by GAMMA's `Match_Type=Indel`
classification), but the old code reported them as fake point substitutions —
`L166A,E167R,L168D,N169T,D178S` and `L166S,E167A,L168I,N169P,D178P,V239A,T242Y` respectively. The same
bug also affected `ompK36` in all four genomes (a contig-edge truncation, not an indel within the gene,
but still a length mismatch that the fixed-position lookup couldn't handle), producing spurious
`G134L,D135*/N,G213T`-style calls that differed between genomes only because of slightly different
truncation lengths, even though the genomes are otherwise identical at this locus.

Note that the SeqKit-amplicon cross-check (`fos_cazavi/mutations.py`) calls the *same* underlying
`detect_mutations()` on a translated PCR amplicon, so it shares this exact limitation — its agreement
with the old BLAST calls in `*_unified_mutations.tsv` was **not** independent confirmation, just the same
bug reproduced twice. GAMMA's protein-level alignment (which does proper indel-aware alignment) was the
only one of the three methods that got this right all along, correctly flagging `Match_Type=Indel`.

**Fix:** added `detect_mutations_aligned()` in `fos_cazavi/utils.py`, which walks BLAST's gapped
nucleotide alignment (`qseq`/`sseq`, now requested via `-outfmt`) column-by-column to map each known
mutation's reference codon to its true corresponding query codon (or to a gap, reported as e.g.
`L166del`), instead of assuming the query and reference are the same length. `fos_cazavi/acquired.py`
now uses this function. All 285 existing unit tests pass; the genomes in this folder were re-run with
the fix and the result files in this folder (`*_results.tsv`, `*_all_results.tsv`, `*_summary.txt`,
`*_blast.txt`, `*_amplicons.tsv`) now reflect the corrected output.

## Comparison to Pariona et al. 2024 (DOI 10.1128/spectrum.01173-24)

The paper describes *in vitro* meropenem selection of parental strain 331 (KPC-114-positive, ST11)
yielding "KPC-2 mutant #1/#2" via a `543_548delCTCATC` deletion in *bla*<sub>KPC-114</sub> that reverts it
to full-length, wild-type **blaKPC-2** (Table S1 of the supplement). Two of the genomes in this dataset
are an exact match:

- **GCA_038433305.1** (`KPC2-MUT1`, contig prefix `JBCEBH`) and **GCA_038433285.1** (`KPC2-MUT2`, contig
  prefix `JBCEBI`) are the GenBank assemblies corresponding to the paper's "KPC-2 mutant #1/#2." Both show
  100.00% identity / 100.00% coverage, **no mutations**, exactly matching the paper's reversion-to-wild-type
  finding.

The other two genomes do **not** correspond to anything in this paper's main text or supplementary
tables:

- **GCA_038433235.1** (`novelKPC-MUT1`) and **GCA_038433245.1** (`novelKPC-MUT2`) carry small in-frame
  deletions at residues 163–172 and 164–167 respectively — upstream of and distinct from the
  `543_548delCTCATC`/S181_P182 deletion described in Table S1. These appear to be a different, unrelated
  blaKPC variant/experiment bundled into the same BioProject, not additional reversion mutants from the
  same selection.

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
