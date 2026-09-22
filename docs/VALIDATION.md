# Validation against real, published genomes

The pipeline was run end-to-end (BLAST + GAMMA + seqkit) on real NCBI
assemblies from published studies, and the calls compared with what those
studies reported. Every result folder under `bioproject_tests/` was regenerated
with the current code and the current reference data (AMRFinderPlus
2026-08-07.1).

| Folder | Genomes | Study | Outcome |
|---|---|---|---|
| `PRJNA741867_test_results/` | 6 *K. pneumoniae* ST307 | Clinical ceftazidime-avibactam-selected KPC variants | **6/6 concordant**, exact allele assignment for all three resistant isolates |
| `PRJNA595047_test/` | 4 *K. pneumoniae* | In vitro selection of KPC Omega-loop deletion mutants | **4/4 concordant** with the study's own strain naming |
| `PRJNA1086695_test/` | 2 long-read assemblies | Assembly + detection | blaKPC-179 identified in one isolate |
| `PRJNA781811_test/` | 18 *K. pneumoniae* / *K. variicola* | Bacteraemia isolate collection | Genotype-only comparison; 1 unambiguous acquired fosA, 3 ambiguous (Indeterminate) |

## PRJNA741867 — the clearest test

Three patients, each with a susceptible baseline isolate and a
ceftazidime-avibactam-resistant isolate that emerged on therapy.

| Sample | Paper | Called allele | Changes (Ambler) | Predicted |
|---|---|---|---|---|
| A-1 | Susceptible | blaKPC-3 | H274Y | Susceptible ✅ |
| A-2 | Resistant | blaKPC-46 | L169P; H274Y | Resistant ✅ |
| B-1 | Susceptible | blaKPC-3 | H274Y | Susceptible ✅ |
| B-2 | Resistant | blaKPC-66 | E166_L167del; H274Y | Resistant ✅ |
| C-1 | Susceptible | blaKPC-3 | H274Y | Susceptible ✅ |
| C-2 | Resistant | blaKPC-92 | E168D; L169_N170del; H274Y | Resistant ✅ |

The discriminating detail: all six carry `H274Y`, which is simply what makes a
KPC a KPC-3. It is reported among the protein changes but is not treated as a
resistance marker, which is why the three baselines come out susceptible.

## What the current version changed

Re-running the previously committed BioProject results with the corrected
pipeline changed several calls. Each change is a correction:

* **A false-positive fosfomycin call removed.** In PRJNA781811,
  GCA_027152215.1 was previously reported as carrying acquired `fosA5` at
  96.19% identity and called resistant. The locus in fact matches the intrinsic
  chromosomal `fosAKP` better; allele assignment now ranks every reference by
  bitscore instead of choosing from a truncated hit list, and the isolate is
  called susceptible.
* **A false-negative ceftazidime-avibactam call fixed.** In PRJNA1086695,
  SRR28296939 carries an insertion immediately after Omega-loop residue D179,
  which matches NCBI allele KPC-179 — curated as an inhibitor-resistant
  extended-spectrum enzyme. The earlier version could not represent insertions
  and reported a wild-type blaKPC with a susceptible call.
* **Novel variants are now called on mechanism.** In PRJNA595047, the two
  strains the study itself names `novelKPC-MUT1/2` carry Omega-loop deletions
  that match no named allele. They are now reported as novel blaKPC variants
  with their changes listed, and called resistant because an in-frame
  Omega-loop deletion is a documented avibactam-escape mechanism.
* **Spurious chromosomal "mutations" gone.** Calls such as `G213I` in ompK36 or
  `A333P` in ftsI were *K. pneumoniae*-versus-*E. coli* sequence differences
  scored against a mismatched reference. Chromosomal point mutations are now
  numbered against an organism-matched reference protein and only reported for
  curated positions.

## Reproducing this

```bash
datasets download genome accession <accession> --include genome
unzip ncbi_dataset.zip

fos-cazavi fos-cazavi-all \
    -a ncbi_dataset/data/<accession>/<accession>_*_genomic.fna \
    -o <accession> \
    -d fos_cazavi/data/example_database.fasta \
    --organism Klebsiella_pneumoniae
```

Each folder's `RESULTS_SUMMARY.md` or `COMPARISON_TO_PAPER.md` has the full
per-genome detail.

## The fosfomycin side

The fosfomycin half has no paired MIC data in any of these BioProjects, so it is
validated by mechanism rather than against measured phenotypes:

* Acquired enzyme detection was exercised on the 18-genome PRJNA781811 set,
  where it separates one unambiguous acquired fosA3 from 14 intrinsic-only
  isolates and three that cannot be resolved by sequence alone.
* Loss-of-function detection (nonsense, frameshift, truncation in the uptake and
  regulatory genes) and the curated fosfomycin mutations are covered by the
  synthetic scenarios, which declare their expected result up front.
* Two curated mutations sitting in fosfomycin genes but belonging to *other*
  drugs — `cyaA_S352T` (fosmidomycin) and `galU_R101C` (cephalosporin) — have
  explicit regression tests asserting they do **not** produce a fosfomycin call.

## How far this goes

These are genotype-to-published-phenotype comparisons on a few dozen genomes,
most of them *K. pneumoniae*. The ceftazidime-avibactam side is validated
against reported phenotypes; the fosfomycin side is not, because no isolate here
has a measured fosfomycin MIC. There is also no *E. coli* or *P. aeruginosa*
validation set. The synthetic scenarios in `create_test_genomes.py` cover the
logic paths the real genomes do not, but a synthetic genome only tests that the
code does what it was designed to do — not that the design matches biology.
Treat the tool accordingly, and see [METHODS.md](METHODS.md#9-known-limits).
