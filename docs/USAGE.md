# Usage

## Quick start

The reference data ships with the package, so a run needs only an assembly:

```bash
fos-cazavi fos-cazavi-all \
    -a your_assembly.fasta \
    -o results \
    --organism Klebsiella_pneumoniae
```

## `--organism` matters

Chromosomal point mutations are only meaningful against a species-matched
reference: without it, ordinary between-species sequence differences are
indistinguishable from resistance mutations. So `--organism` is required for
those calls, and the tool says so on stderr when it is missing.

Supported values (AMRFinderPlus taxgroup names):

* `Escherichia`
* `Klebsiella_pneumoniae`
* `Pseudomonas_aeruginosa`

Without `--organism` the tool still reports, in full:

* acquired resistance genes and their alleles,
* blaKPC typing and the ceftazidime-avibactam call,
* loss of function in chromosomal genes (premature stops, frameshifts,
  truncations), which is species independent.

## Subcommands

### `fos-cazavi-all`

The full pipeline: gene detection and variant calling, the GAMMA cross-check,
amplicon mapping, and all summary outputs.

```bash
fos-cazavi fos-cazavi-all \
    -a <assembly> \
    -o <output_prefix> \
    [--organism <organism>] \
    [-d <database>] [--genes <fasta>] [--primers <tsv>] \
    [--min_id 90] [--min_cov 80]
```

### `fos-cazavi-acquired`

Gene detection and variant calling only (BLAST). No GAMMA, no amplicons.

```bash
fos-cazavi fos-cazavi-acquired -a <assembly> -o <prefix> [--organism <organism>]
```

### `fos-cazavi-mutations`

The GAMMA and `seqkit` analyses on their own, without the BLAST caller.

```bash
fos-cazavi fos-cazavi-mutations -a <assembly> -o <prefix> [--genes <fasta>] [--primers <tsv>]
```

### `batch`

Analyse a collection and get one table for the whole thing.

```bash
fos-cazavi batch \
    -i assemblies/ \
    -o results/ \
    --organism Klebsiella_pneumoniae \
    -j 8
```

`-i` takes files, directories, or a mix. Each sample gets the usual per-sample
outputs in `-o`, plus two combined tables:

| File | Contents |
|---|---|
| `batch_combined_summary.tsv` | One row per sample: both phenotype calls, the genes found, loss-of-function findings, QC columns and the full evidence |
| `batch_combined_genes.tsv` | One row per detected gene copy across all samples |

A sample that fails is reported with its traceback and does not stop the run;
the exit status is non-zero if any sample failed.

### `combine`

Build the same two tables from results you already have, including runs done
separately:

```bash
fos-cazavi combine -i results_run1/ results_run2/ -o all
```

It reads `*_summary.json`, so it works on any output directory.

### `create-db`

Rebuilds the reference data from the current AMRFinderPlus release:

```bash
fos-cazavi create-db -e your.email@example.com -o <output_dir>
```

This runs `python3 -m fos_cazavi.build_data`, which downloads the AMRFinderPlus
files, rebuilds the sequence database, the reference proteins, the validated
point-mutation table and the blaKPC allele table, and reports anything it had
to drop. (`-e/--email` is kept for backwards compatibility and is not used.)

## Options

| Option | Default | Meaning |
|---|---|---|
| `-a, --assembly` | required | Input assembly (FASTA) |
| `-o, --output` | required | Output prefix |
| `--organism` | none | Sample organism; enables curated chromosomal point mutations |
| `-d, --database` | bundled | Nucleotide reference database |
| `--genes` | bundled | Nucleotide CDS database for GAMMA |
| `--primers` | bundled | Primer definitions for amplicon mapping |
| `--mutations` | bundled | Point-mutation definitions |
| `--min_id` | 90 | Minimum percent identity |
| `--min_cov` | 80 | Minimum percent coverage of the reference gene |
| `-t, --threads` | 1 | Threads handed to `blastn` and `seqkit` within one sample |
| `-j, --jobs` | cores − 1 | (`batch` only) Assemblies analysed in parallel |

## Parallelism

Two knobs, and they trade off:

* **`--jobs`** runs that many assemblies at once, each in its own process.
* **`--threads`** is passed to the external tools within a single assembly
  (`blastn -num_threads`, `seqkit -j`).

For a collection, `--jobs` is almost always the better lever: one assembly takes
a couple of seconds, and BLAST scales poorly across threads on a query that
small. On ten *E. coli* genomes, `-j 4` cut a 15.9 s serial run to 6.1 s.

Total CPU use is roughly `jobs × threads`, so keep the product at or below your
core count. Results are identical either way — parallelism is a performance
knob, and the test suite asserts that serial and parallel runs produce
byte-identical tables.

## Required external tools

`blastn` and `makeblastdb` are required. `GAMMA.py` (with `blat`) and `seqkit`
are optional: without them the cross-check and amplicon mapping are skipped and
the run continues, with a warning.

See [METHODS.md](METHODS.md) for how the calls are made and what they mean.
