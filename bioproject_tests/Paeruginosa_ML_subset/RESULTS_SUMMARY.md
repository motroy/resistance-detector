# *Pseudomonas aeruginosa* subset — scope and robustness test

Twelve complete *P. aeruginosa* genomes, selected from the 1,437-genome
supplementary table of Noman *et al.*, *Machine Learning Techniques for
Antimicrobial Resistance Prediction of Pseudomonas aeruginosa from Whole Genome
Sequence Data*.

## What this set is, and is not

It is **not** a phenotype validation. The source table's per-drug columns are
almost invariant (fosfomycin = 1 in 1,425 of 1,437 rows) with no stated S/I/R
semantics, and 1,054 of the 1,437 rows are labelled *Computational Prediction*
rather than measured AST — using those as ground truth would be circular. The
table also reports **ceftazidime**, not ceftazidime-avibactam; they are
different agents and one cannot validate the other.

What it **is**: a test of the tool on a species it barely supported, and a
gene-level comparison against the paper's own gene calls (genotype vs genotype,
which is a fair comparison). It was chosen to exercise metallo-beta-lactamase
detection, KPC in a non-*Enterobacterales* host, and the intrinsic-gene logic.

## Genomes

| Strain | Accession(s) | Carbapenemase per paper |
|---|---|---|
| ST773 | CP041945 | blaNDM-1 |
| CDN118 | CP054591 | blaVIM-2 |
| PA99 | CP042967 | blaIMP-1 |
| AG1 | CP045739 | blaIMP-18, blaVIM-2 |
| R31 | CP061850, CP061851 | blaKPC-2 |
| P23 | CP065417, CP065418 | blaKPC-2 |
| SE5331 | CP046402 | blaGES-5/-6/-7/-13 |
| SE5352 | CP054843 | blaGES-1 |
| AR_0360 | CP027165 | none |
| AR_0440 | CP029148 | none |
| AR442 | CP029090 | none |
| AR_0095 | CP027538 | none |

## Results

| Strain | Beta-lactamases detected | Predicted FOS | Predicted CAZ/AVI |
|---|---|---|---|
| ST773 | blaPDC-like, **blaNDM-1** | Resistant (intrinsic) | **Resistant** |
| CDN118 | blaPDC-like, **blaVIM-2** | Resistant (intrinsic) | **Resistant** |
| PA99 | blaPDC-like, **blaIMP-1 ×2, blaIMP-like** | Resistant (intrinsic) | **Resistant** |
| AG1 | blaPDC-like, **blaVIM-2, blaIMP-like** | Resistant (intrinsic) | **Resistant** |
| R31 | blaPDC-like, blaKPC-2 | Resistant (intrinsic) | Indeterminate |
| P23 | blaPDC-like, blaKPC-2 | Resistant (intrinsic) | Indeterminate |
| SE5331 | blaPDC-like, blaGES-like | Resistant (intrinsic) | Indeterminate |
| SE5352 | blaPDC-like, blaGES-like | Resistant (intrinsic) | Indeterminate |
| AR_0360 | blaPDC-like | Resistant (intrinsic) | Indeterminate |
| AR_0440 | blaPDC-like | Resistant (intrinsic) | Indeterminate |
| AR442 | blaPDC-like | Resistant (intrinsic) | Indeterminate |
| AR_0095 | blaPDC-like | Resistant (intrinsic) | Indeterminate |

Gene-level concordance with the paper, for genes in this tool's scope: the
chromosomal AmpC (`blaPAO` / `blaPDC`) and chromosomal `fosA` were found in
12/12, and every carbapenemase the paper reports was found — including
`blaIMP-18` in AG1, which an earlier version of the database missed.

`PA99` carries three `blaIMP` copies at genuinely distinct chromosomal loci
(4.66 Mb, 5.28 Mb, 5.74 Mb), which the copy-number reporting shows correctly.

Not detected, and out of this tool's scope: the acquired class D oxacillinases
`blaOXA-2`, `blaOXA-4` and `blaOXA-21`, and the non-beta-lactam genes
(aminoglycoside, sulphonamide, tetracycline, `crpP`, `catB7`, `qacE`).

## What this set changed in the tool

Running it exposed four real problems, all now fixed:

1. **Fosfomycin was reported Susceptible for every *P. aeruginosa*.** The
   species is intrinsically fosfomycin-resistant and has no fosfomycin
   breakpoints; the rule that intrinsic `fosA` is not scored as *acquired*
   resistance had been wrongly letting the species-level call come out
   susceptible. Intrinsic species resistance now takes precedence.
2. **Ceftazidime-avibactam was reported Susceptible when the dominant mechanism
   had not been looked at.** In *P. aeruginosa*, PDC/AmpC derepression and PDC
   variants drive most CAZ/AVI resistance, and neither is assessed here. A
   negative result is now `Indeterminate` with that stated, not `Susceptible`.
3. **Most blaIMP alleles were undetectable.** 62 of the 108 known blaIMP
   alleles are below 90% identity to blaIMP-1, so two reference alleles missed
   most of the family at the default threshold. All alleles of the diverse
   carbapenemase families (NDM, VIM, IMP, SPM, GIM, SIM, GES — 375 sequences)
   are now included.
4. **Allele names were asserted more precisely than the data supports.**
   blaIMP-18 and blaIMP-99 are 99.7% identical, so "closest reference" is not
   evidence of which allele is present. Outside blaKPC — the one family with
   real allele typing — an inexact match is now reported as the family
   (`blaIMP-like`), with the closest reference still shown in the `Gene`
   column.

## Reproducing

```bash
# accessions are in subset_accessions.tsv
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=<acc>&rettype=fasta&retmode=text" > <strain>.fna

fos-cazavi fos-cazavi-all -a <strain>.fna -o <strain> \
    -d fos_cazavi/data/example_database.fasta \
    --organism Pseudomonas_aeruginosa
```

## Limits

Every CAZ/AVI call for *P. aeruginosa* that is not driven by an MBL is
`Indeterminate`, which is honest but not very useful. Making the tool actually
informative for this species would need PDC variant typing and some proxy for
expression — see [../../docs/METHODS.md](../../docs/METHODS.md#9-known-limits).
