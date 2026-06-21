# FOS-CAZAVI Resistance Detector

A CLI tool for detecting fosfomycin (FOS) and ceftazidime-avibactam (CAZAVI) resistance genes and mutations from bacterial genome assemblies.

## Features

- **Modular CLI**: Separate commands for database creation, acquired gene detection, and mutation analysis.
- **Gene Detection**: Identifies resistance genes (fosA variants, blaKPC, blaOXA-48, etc.) using BLAST+.
- **Mutation Detection**: Detects known resistance mutations (D179Y, V240G, T243M, etc.) using GAMMA (protein-level gene alignment via translated nucleotide CDS) and SeqKit amplicon analysis.
- **Indel-aware mutation calling**: Both the BLAST and SeqKit mutation-detection paths walk a real gapped alignment rather than doing a fixed-position lookup, so insertions/deletions relative to the reference gene are reported as e.g. `L166del` instead of producing fabricated, frame-shifted point-substitution calls.
- **Amplicon Detection**: Uses seqkit amplicon to find PCR products from primer pairs and checks them for resistance mutations.
- **Copy-Number Reporting**: Annotates every detected gene with the number of distinct genomic loci it was found at.
- **Predicted Phenotypes**: Derives a conclusive, genotype-based Susceptible/Resistant call for both fosfomycin and ceftazidime-avibactam, with supporting evidence. See [docs/OUTPUT_FILES.md](docs/OUTPUT_FILES.md#predicted-phenotypes).
- **Multiple Output Formats**: Human-readable summary plus machine-parsable TSV/JSON.
- **Logging**: Comprehensive logging of commands, parameters, and tool versions.

## Documentation

- [Installation](docs/INSTALLATION.md)
- [Usage / CLI Reference](docs/USAGE.md)
- [Output Files Reference](docs/OUTPUT_FILES.md)
- [Testing](docs/TESTING.md)
- [Validation Against Real, Published Genomes](docs/VALIDATION.md)

## Repository Structure

```
resistance-detector/
├── fos_cazavi/              # Python package (installed via pip)
│   ├── __init__.py
│   ├── cli.py               # CLI entry point
│   ├── db.py                # Database creation module
│   ├── acquired.py          # BLAST-based detection module
│   ├── mutations.py         # GAMMA/Amplicon detection module
│   ├── phenotype.py         # Genotype-to-phenotype prediction (FOS, CAZ/AVI)
│   ├── utils.py             # Shared utilities
│   └── data/                # Bundled reference data
│       ├── example_database.fasta              # Raw nucleotide CDS sequences
│       ├── example_database_deduplicated.fasta # GAMMA-ready database (GAMMA_DB_Maker output)
│       ├── example_database_mutations.tsv      # Mutation definitions TSV
│       └── primers.tsv                         # Primer sequences for amplicon detection
├── fos-cazavi               # Executable script (for repo use without install)
├── pyproject.toml           # PyPI package configuration
├── GAMMA_DB_Maker.py        # GAMMA database preparation tool
├── create_test_genomes.py   # Test genome generator
├── batch_analysis.sh        # Batch processing script
├── Snakefile                # Snakemake workflow
├── config.yaml              # Snakemake configuration
├── environment.yaml         # Conda environment
├── scripts/
│   ├── download_assembly.py      # Download assemblies from NCBI
│   └── simulate_esbl_genomes.py  # Download real ESBL genomes and simulate resistance
├── tests/                   # Unit test suite (305 tests)
├── docs/                    # Documentation (installation, usage, outputs, testing, validation)
├── example_results/         # Example outputs
├── bioproject_tests/        # Real-genome validation against published bioprojects
│   ├── PRJNA595047_test/        # 4 K. pneumoniae assemblies (NCBI),
│   │                            #   incl. comparison against a published in vitro selection study
│   ├── PRJNA741867_test_results/ # 6 K. pneumoniae ST307 assemblies (NCBI),
│   │                            #   incl. comparison against a published clinical CAZ/AVI-resistance study
│   └── PRJNA1086695_test/       # 2 myloasm-assembled isolates
├── LICENSE
└── README.md
```

## Quick Start

```bash
pip install fos-cazavi

fos-cazavi fos-cazavi-all \
    -a your_assembly.fasta \
    -d resistance_db.fasta \
    -o results
```

See [docs/INSTALLATION.md](docs/INSTALLATION.md) and [docs/USAGE.md](docs/USAGE.md) for full setup and CLI details.

## License

MIT License
