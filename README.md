# FOS-CAZAVI Resistance Detector

A CLI tool for detecting fosfomycin (FOS) and ceftazidime-avibactam (CAZ/AVI)
resistance genes and mutations in bacterial genome assemblies.

Both drugs are treated as first-class targets: acquired enzymes and the
chromosomal mechanisms for each, with every call carrying the evidence behind
it and the drug it belongs to.

## What it does

- **Types blaKPC alleles.** The observed protein changes are compared with every
  blaKPC allele defined by NCBI, so a hit is reported as `blaKPC-31`,
  `blaKPC-66`, or as a novel variant with its changes listed — not just as
  "blaKPC present".
- **Calls variants properly.** The gene span is recovered from the contig,
  translated, and aligned to the reference protein. Substitutions, in-frame
  indels, insertions, premature stops and frameshifts are each reported as what
  they are, in standardised **Ambler numbering** for class A beta-lactamases.
- **Predicts phenotypes with stated evidence** — `Resistant`,
  `Indeterminate` or `Susceptible`, each with the specific finding behind it.
  Hotspot changes that are not documented give `Indeterminate` rather than a
  guess in either direction.
- **Detects metallo-beta-lactamases** (NDM, VIM, IMP, SPM, GIM, SIM), which
  avibactam does not inhibit, and which therefore make CAZ/AVI inactive
  regardless of any KPC present.
- **Distinguishes acquired from intrinsic.** The chromosomal fosA of *K.
  pneumoniae* and *P. aeruginosa* is never scored as acquired fosfomycin
  resistance, and a lone fosA hit that cannot be told apart from the species'
  own chromosomal copy is reported as Indeterminate rather than asserted.
- **Keeps to its two drugs.** Most curated mutations in these genes were
  described for other antibiotics; each is tagged with the drug it belongs to,
  so a carbapenem or tigecycline mutation is never presented as a FOS or
  CAZ/AVI finding.
- **Cross-checks with a second caller.** GAMMA independently reports codon
  changes; agreement between the two is recorded per change.
- **Says when it cannot tell.** Genes running off a contig boundary are flagged
  rather than silently called.

## Quick start

```bash
pip install fos-cazavi

fos-cazavi fos-cazavi-all \
    -a your_assembly.fasta \
    -o results \
    --organism Klebsiella_pneumoniae
```

The reference data is bundled; `-d` is optional. `--organism` is required for
curated chromosomal point mutations, which are only meaningful against a
species-matched reference — see [docs/USAGE.md](docs/USAGE.md).

Example output:

```
PREDICTED PHENOTYPES (genotype-based):
  Fosfomycin (FOS): Resistant
    - Acquired fosfomycin-modifying enzyme fosA3 (100.00% identity, 100.00% coverage)
  Ceftazidime-Avibactam (CAZ/AVI): Resistant
    - blaKPC-33 on contig_blaKPC: blaKPC-33 is curated by NCBI as "inhibitor-resistant
      extended-spectrum class A beta-lactamase KPC-33" (subclass CEPHALOSPORIN)
    - blaKPC-33 on contig_blaKPC: D179Y: documented ceftazidime-avibactam resistance substitution
  Note: Genotype-based prediction only; not a substitute for phenotypic AST.
```

## Validation

Against **measured MICs** in ten carbapenem-resistant *E. coli*: fosfomycin
10/10 correct, ceftazidime-avibactam 9/10 resistant with the tenth reported
`Indeterminate` rather than wrong. On six *K. pneumoniae* ST307 assemblies from
a study of CAZ/AVI resistance emerging on therapy, it reproduced all six
reported phenotypes and assigned the exact blaKPC allele (KPC-46, KPC-66,
KPC-92) in each resistant isolate. Five further sets are included. See
[docs/VALIDATION.md](docs/VALIDATION.md), and
[docs/METHODS.md](docs/METHODS.md#9-known-limits) for what the method cannot do.

## Documentation

- [Methods — how every call is made, and its limits](docs/METHODS.md)
- [Installation](docs/INSTALLATION.md)
- [Usage / CLI reference](docs/USAGE.md)
- [Output files](docs/OUTPUT_FILES.md)
- [Testing](docs/TESTING.md)
- [Validation against published genomes](docs/VALIDATION.md)

## Repository layout

```
resistance-detector/
├── fos_cazavi/
│   ├── cli.py            # CLI entry point and report writing
│   ├── acquired.py       # BLAST detection and per-locus variant calling
│   ├── variants.py       # Protein-level variant caller, Ambler numbering
│   ├── betalactamase.py  # blaKPC allele typing, CAZ/AVI marker rules
│   ├── phenotype.py      # Genotype-to-phenotype prediction
│   ├── references.py     # Loading and re-validating the reference data
│   ├── mutations.py      # GAMMA cross-check, seqkit amplicon mapping
│   ├── build_data.py     # Rebuilds all reference data from AMRFinderPlus
│   ├── db.py             # `create-db`, a thin wrapper over build_data
│   ├── utils.py          # Logging, dependency checks, primer loading
│   └── data/             # Bundled reference data (see METHODS.md)
├── create_test_genomes.py  # Synthetic genomes + their expected results
├── tests/                  # Unit and end-to-end tests
├── docs/
├── bioproject_tests/       # Real-genome validation runs
└── example_results/
```

## Requirements

Python ≥3.8, Biopython, and BLAST+ (`blastn`, `makeblastdb`). `GAMMA.py` (with
`blat`) and `seqkit` are optional; without them the cross-check and amplicon
mapping are skipped.

## Caveat

These are genotypic predictions. They are not a substitute for phenotypic
antimicrobial susceptibility testing.

## License

MIT License
