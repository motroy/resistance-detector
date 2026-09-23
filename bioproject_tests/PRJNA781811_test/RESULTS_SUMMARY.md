# PRJNA781811 — FOS-CAZAVI Resistance Detection Results

18 *Klebsiella pneumoniae* / *K. variicola* assemblies from BioProject PRJNA781811
(listed in `PRJNA781811.ncbi_datasets.tsv`; see
[Arena et al. 2022, Front. Microbiol. 13:983294](https://doi.org/10.3389/fmicb.2022.983294)),
run through the current pipeline.

## Command used

```bash
datasets download genome accession <accession> --include genome
unzip ncbi_dataset.zip

fos-cazavi fos-cazavi-all \
    -a ncbi_dataset/data/<accession>/<accession>_*_genomic.fna \
    -o <accession> \
    -d fos_cazavi/data/example_database.fasta \
    --organism Klebsiella_pneumoniae
```

Reference data: AMRFinderPlus 2026-08-07.1 (see `fos_cazavi/data/DATA_VERSION.txt`).

## Results

| Accession | Beta-lactamases | fosA-family | Predicted FOS | Predicted CAZ/AVI |
|---|---|---|---|---|
| GCA_027151785.1 | SHV-12 | fosAKP (intrinsic) | **Resistant** | Susceptible |
| GCA_027151795.1 | SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027151835.1 | KPC-3, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027151845.1 | SHV-12 | **fosA3** + intrinsic fosAKP | **Resistant** | Susceptible |
| GCA_027151875.1 | KPC-3, CMY-2, CTX-M-15, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027151985.1 | SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027151995.1 | CTX-M-15, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152005.1 | CTX-M-15, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152065.1 | KPC-3, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152105.1 | CTX-M-15, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152185.1 | SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152205.1 | SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152215.1 | SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152225.1 | KPC-2, CTX-M-15, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152245.1 | CTX-M-15, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152405.1 | SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152445.1 | KPC-3, SHV-12 | fosAKP (intrinsic) | **Resistant** | Susceptible |
| GCA_027152495.1 (*K. variicola*) | SHV-12 | fosA9, no intrinsic copy | **Indeterminate** | Susceptible |

One isolate (GCA_027151845.1) carries a fosA3 **alongside** an intact intrinsic
chromosomal `fosAKP`, which makes the fosA3 unambiguously acquired: predicted
fosfomycin-resistant.

One isolate — **GCA_027152495.1, the sole *K. variicola* genome in this set** —
carries a single fosA-family hit with no intrinsic-named copy recognised, and
that hit is not close enough to the (*K. pneumoniae*-sourced) `fosAKP`
reference to be resolved as the same gene (see "Changes" below). Every
*Klebsiella* has a chromosomal fosA, so this cannot be resolved by sequence
identity alone — the hit may be a divergent copy of the species' own gene
rather than an acquired one, and being a different species from the reference
strain makes that the more likely explanation here. It is reported as
**Indeterminate** with that ambiguity stated, rather than asserted as resistant.

Two isolates — **GCA_027151785.1** and **GCA_027152445.1** — carry an intact
intrinsic `fosAKP` but also a premature stop in a fosfomycin transport gene
(`uhpB` at residue 77/491, and `glpT` at residue 393/448, respectively), each
independently confirmed by GAMMA's own alignment. Both are predicted
fosfomycin-resistant on that basis. Neither gene was checked for this species
before the reference-data fix described below; see
[`ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md`](../ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md#2-a-real-reference-database-gap--now-fixed).
No phenotypic AST data exists for this BioProject, so these two calls cannot
be checked against a measured MIC — they are reported as a genuinely new
genotype finding, not a validated one.

The remaining 14 carry only the intrinsic `fosAKP`, with no loss of function
in any fosfomycin gene, and are predicted susceptible.

All 18 are predicted ceftazidime-avibactam-susceptible. Four carry blaKPC
(KPC-2 or KPC-3) with no Omega-loop, 237–243 or insertion-loop change, so
avibactam is expected to inhibit the enzyme; none carries a
metallo-beta-lactamase.

## Changes from earlier runs of this BioProject

These results have gone through three corrections since first committed:

* **GCA_027152215.1** was originally reported as carrying `fosA5` at 96.19%
  identity and called fosfomycin-resistant. Allele assignment was fixed to
  choose by bitscore across all fosA references rather than from a truncated
  hit list, and the same locus matched intrinsic `fosAKP` better — susceptible.
* Once fixed, **GCA_027151985.1**, **GCA_027152185.1** and **GCA_027152495.1**
  still had their single fosA locus reported as an ambiguous acquired allele
  (`fosA10-like`, `fosA5-like`, `fosA9-like`) rather than the intrinsic gene,
  because the intrinsic reference lost the bitscore race by a hair at the
  identical span — a real hit-selection bug, not a labelling choice. All three
  were reported `Indeterminate`.
* That bug is now fixed (`BlastDetector._prefer_intrinsic_naming()`; see
  [`ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md`](../ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md)
  for the full account, found via real MIC-tested genomes). **GCA_027151985.1**
  and **GCA_027152185.1** now correctly resolve to intrinsic `fosAKP` —
  susceptible. **GCA_027152495.1** does not: it is *K. variicola*, a different
  species from the `fosAKP` reference strain, and its native copy sits outside
  the fix's identity tolerance — it remains genuinely ambiguous and stays
  `Indeterminate`. This is the fix discriminating correctly: real ambiguity
  (different species) stays flagged; spurious ambiguity (ordinary same-species
  strain divergence) gets resolved.
* blaKPC alleles are now named from the observed protein changes against
  KPC-2, rather than from whichever KPC reference happened to win the BLAST
  hit.
* **GCA_027151785.1** and **GCA_027152445.1** now report fosfomycin
  `Resistant`, previously `Susceptible`. The fosfomycin transport genes
  (`uhpT`, `uhpA`, `uhpB`, `uhpC`, `glpT`, `cyaA`, `ptsI`, `galU`, `murA`)
  were *E. coli*-only references and invisible in this species at the default
  identity threshold; a *K. pneumoniae*-specific reference was added for each
  (see [`ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md`](../ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md#2-a-real-reference-database-gap--now-fixed)),
  and both isolates' premature stops (in `uhpB` and `glpT` respectively) are
  now detected and independently confirmed by GAMMA.

See [../../docs/METHODS.md](../../docs/METHODS.md) for how allele assignment
and the phenotype rules work.

## Limits

No phenotypic AST results are available for these isolates in the referenced
study's public data, so this is a genotype-only comparison: it shows the
pipeline's gene and allele assignments on real assemblies, not its accuracy
against measured MICs.
