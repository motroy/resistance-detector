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
| GCA_027151785.1 | SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027151795.1 | SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027151835.1 | KPC-3, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027151845.1 | SHV-12 | **fosA3** + fosAKP | **Resistant** | Susceptible |
| GCA_027151875.1 | KPC-3, CMY-2, CTX-M-15, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027151985.1 | SHV-12 | **fosA10** | **Resistant** | Susceptible |
| GCA_027151995.1 | CTX-M-15, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152005.1 | CTX-M-15, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152065.1 | KPC-3, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152105.1 | CTX-M-15, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152185.1 | SHV-12 | **fosA5** | **Resistant** | Susceptible |
| GCA_027152205.1 | SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152215.1 | SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152225.1 | KPC-2, CTX-M-15, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152245.1 | CTX-M-15, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152405.1 | SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152445.1 | KPC-3, SHV-12 | fosAKP (intrinsic) | Susceptible | Susceptible |
| GCA_027152495.1 | SHV-12 | **fosA5** | **Resistant** | Susceptible |

4 of 18 carry an acquired fosA-family enzyme (fosA3, fosA5 ×2, fosA10) and are
predicted fosfomycin-resistant. The other 14 carry only the intrinsic
chromosomal `fosAKP`, which is present in fosfomycin-susceptible *K.
pneumoniae* and is not scored as acquired resistance.

All 18 are predicted ceftazidime-avibactam-susceptible. Four carry blaKPC
(KPC-2 or KPC-3) with no Omega-loop, 237–243 or insertion-loop change, so
avibactam is expected to inhibit the enzyme; none carries a
metallo-beta-lactamase.

## Changes from the earlier run of this BioProject

These results differ from the ones previously committed here, because the
detection method changed:

* **GCA_027152215.1** was previously reported as carrying `fosA5` at 96.19%
  identity and called fosfomycin-resistant. The same locus now matches the
  intrinsic `fosAKP` reference better (96.19% to fosAKP, chosen by bitscore
  across all fosA references rather than from a truncated hit list), and the
  isolate is called susceptible. The earlier call was a false positive
  produced by assigning an intrinsic gene to an acquired allele.
* **GCA_027151985.1** is now reported as `fosA10` (99.29%), an allele that was
  not in the earlier reference set.
* blaKPC alleles are now named from the observed protein changes against
  KPC-2, rather than from whichever KPC reference happened to win the BLAST
  hit.

See [../../docs/METHODS.md](../../docs/METHODS.md) for how allele assignment
and the phenotype rules work.

## Limits

No phenotypic AST results are available for these isolates in the referenced
study's public data, so this is a genotype-only comparison: it shows the
pipeline's gene and allele assignments on real assemblies, not its accuracy
against measured MICs.
