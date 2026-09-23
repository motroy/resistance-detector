# Validation against real, published genomes

The pipeline was run end-to-end (BLAST + GAMMA + seqkit) on real NCBI
assemblies from published studies, and the calls compared with what those
studies reported. Every result folder under `bioproject_tests/` was regenerated
with the current code and the current reference data (AMRFinderPlus
2026-08-07.1).

| Folder | Genomes | Study | Outcome |
|---|---|---|---|
| `CREC_fosA3_China/` | 10 *E. coli* | fosA3 in carbapenem-resistant *E. coli* (Zhang *et al.* 2025) | **Measured MICs for both drugs**: fosfomycin 10/10, CAZ/AVI 9/10 resistant + 1 indeterminate, 0 wrong |
| `ESKAPE_fos_GOLD_Kpneumoniae/` | 21 *K. pneumoniae* | Gold-standard ESKAPE fosfomycin AST set | **Measured MICs**: specificity 10/10, sensitivity 2/11 — found and fixed a hit-selection bug and a reference-database gap (9 transport genes were undetectable in this species; now fixed) |
| `ESKAPE_fos_GOLD_Paeruginosa/` | 15 *P. aeruginosa* | Gold-standard ESKAPE fosfomycin AST set | **Measured MICs**; falsified the tool's prior "always Resistant" rule on 10/15 isolates and drove the fix |
| `ESKAPE_fos_Kpneumoniae_round2/` | 24 *K. pneumoniae* | Same source, disjoint accessions | **Measured MICs**: specificity 6/6, sensitivity 0/18 even with the fix applied — replicates that the remaining gap is not a detection problem (combined with round 1: 2/29) |
| `PRJNA741867_test_results/` | 6 *K. pneumoniae* ST307 | Clinical ceftazidime-avibactam-selected KPC variants | **6/6 concordant**, exact allele assignment for all three resistant isolates |
| `PRJNA595047_test/` | 4 *K. pneumoniae* | In vitro selection of KPC Omega-loop deletion mutants | **4/4 concordant** with the study's own strain naming |
| `PRJNA1086695_test/` | 2 long-read assemblies | Assembly + detection | blaKPC-179 identified in one isolate |
| `PRJNA781811_test/` | 18 *K. pneumoniae* / *K. variicola* | Bacteraemia isolate collection | **Measured MICs for both drugs** (paper's Table 3): fosfomycin specificity 8/8 (+1 correctly Indeterminate), sensitivity 3/9; CAZ/AVI 17/18, with the one miss showing three genomes sharing an identical genotype split 2 susceptible/1 resistant |
| `Paeruginosa_ML_subset/` | 12 *P. aeruginosa* | ML AMR-prediction dataset (Noman *et al.*) | Scope/robustness test on a new species; gene-level concordance, not phenotype |
| `Kleborate_cross_check/` | 63 *K. pneumoniae* / *K. variicola* (reuses the sets above) | Independent third-party tool comparison, not a published study | 62/63 concordant on acquired fosA-gene presence and CAZ/AVI-relevant beta-lactamase families; found and fixed two real gaps (blaCTX-M group-9 variants, blaVEB) |

## CREC fosA3 — the first set with measured MICs

Ten carbapenem-resistant *E. coli* with broth MICs for **both** drugs (Zhang
*et al.*, J Glob Antimicrob Resist 42 (2025) 80–87, Table 1; accessions from
Table S4). At the time this was the only measured-MIC set here; the two
ESKAPE-GOLD sets below (*K. pneumoniae*, *P. aeruginosa*) added fosfomycin
MICs, including real susceptible isolates. Every other set compares genotype
with a study's reported genotype or narrative phenotype.

* **Fosfomycin 10/10 correct.** All carry `fosA3` (FOS MIC 256–>256). *E. coli*
  has no intrinsic chromosomal fosA, so there is none of the ambiguity that
  makes lone fosA hits uncertain in *Klebsiella*. This is the first test of the
  fosfomycin side against real MICs.
* **Ceftazidime-avibactam 9/10 resistant, 1 indeterminate, 0 wrong.** Nine carry
  an NDM metallo-beta-lactamase.
* **E2257** is CAZ/AVI resistant (>128) with *no carbapenemase* anywhere in the
  assembly. It carries `blaCMY-2`, `blaCTX-M-15` and a premature stop at residue
  257 of OmpF — AmpC plus lost permeability. It is reported `Indeterminate` with
  that mechanism named, rather than susceptible.

Finding E2257 changed the tool: *E. coli* `ompC`/`ompF` were not in the database
at all, and the porin-amplification rule was gated on blaKPC alone. See
[../bioproject_tests/CREC_fosA3_China/RESULTS_SUMMARY.md](../bioproject_tests/CREC_fosA3_China/RESULTS_SUMMARY.md).

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

Every set, end to end, with one command. Assemblies are fetched from NCBI on
first run (accessions are in each set's `accessions.tsv`) and cached:

```bash
scripts/run_validation.sh -j 6
```

That analyses all 52 assemblies and writes
`bioproject_tests/all_validation_combined_summary.tsv` — one row per sample
across every set — plus the per-gene table beside it. It takes about 35 seconds
on six cores once the assemblies are cached.

### Provenance

The committed results were produced by a run logged with
[dochist](https://github.com/motroy/dochist-docs), which records each command
and checksums every artifact it produces:

* [`docs/PROVENANCE.md`](PROVENANCE.md) — FAIR compliance report for the run:
  7 commands, 112 artifacts with SHA-256 checksums, environment snapshot.
* [`scripts/rerun_validation.sh`](../scripts/rerun_validation.sh) — the
  reproduction script `dochist extract` derived from that session.
* [`demo/validation-run.cast`](../demo/validation-run.cast) — the run itself,
  replayable with `asciinema play`, and rendered as a GIF in
  [`demo/`](../demo/README.md).

Each folder's `RESULTS_SUMMARY.md` or `COMPARISON_TO_PAPER.md` has the full
per-genome detail.

## The *P. aeruginosa* subset

Twelve complete genomes from a 1,437-genome machine-learning AMR dataset. This
is explicitly **not** a phenotype validation — that table's per-drug labels are
near-invariant and mostly computational predictions, and it reports ceftazidime
rather than ceftazidime-avibactam. It is a scope and robustness test, and a
gene-level comparison against the paper's own gene calls.

It found four real defects, all since fixed: fosfomycin reported Susceptible for
an intrinsically resistant species; CAZ/AVI reported Susceptible without having
looked at the dominant mechanism (PDC/AmpC); most blaIMP alleles undetectable
because the family is far more diverse than two references cover; and allele
names asserted more precisely than 99.7%-identical references can support. See
[../bioproject_tests/Paeruginosa_ML_subset/RESULTS_SUMMARY.md](../bioproject_tests/Paeruginosa_ML_subset/RESULTS_SUMMARY.md).

## The fosfomycin side

The fosfomycin half is now validated against measured MICs across two species,
including — for the first time — real **susceptible** isolates, closing what
was previously this tool's biggest untested gap:

* **CREC (*E. coli*), 10 resistant isolates:** 10/10 correct (all carry
  `fosA3`).
* **ESKAPE-GOLD (*K. pneumoniae*), 21 isolates, 10 susceptible:**
  **specificity 10/10** — no false Resistant call anywhere. Sensitivity on the
  resistant/intermediate isolates started at 0/11, which was a real and
  significant finding, not noise: it traced to (a) a genuine detection bug
  that this run found and fixed (four genomes had their single, native `fosA`
  copy mis-named as an ambiguous acquired allele) and (b) a genuine
  reference-database gap — the fosfomycin transport genes this tool checks
  were *E. coli*-only references, so *K. pneumoniae*'s own orthologs sat below
  the detection threshold and those genes were not actually being checked for
  this species at all. Both are now fixed: a K. pneumoniae-specific reference
  for each of the 9 transport/regulatory genes was added
  (`fos_cazavi/references.py::GENE_ALIASES`), and re-running this set found
  two real, GAMMA-confirmed premature stops in `uhpB` (isolates KP_R_02,
  KP_R_05) — sensitivity is now **2/11**. See
  [`ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md`](../bioproject_tests/ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md).
* **ESKAPE-GOLD (*K. pneumoniae*) round 2, 24 more isolates, disjoint
  accessions, 6 susceptible:** a replication run, weighted toward
  resistant/intermediate isolates specifically to test whether 0/11 was
  sampling noise. It wasn't at the time: **specificity 6/6**, and with the
  transport-gene fix now applied, sensitivity is still **0/18** — this
  round's 18 non-susceptible isolates carry no coding-sequence change in any
  of these 9 genes, consistent with the literature that much fosfomycin
  resistance in Enterobacterales is promoter/IS-element-driven rather than a
  coding change a CDS-level tool can see. Combined sensitivity across both
  rounds: **2/29**. See
  [`ESKAPE_fos_Kpneumoniae_round2/RESULTS_SUMMARY.md`](../bioproject_tests/ESKAPE_fos_Kpneumoniae_round2/RESULTS_SUMMARY.md).
* **ESKAPE-GOLD (*P. aeruginosa*), 15 isolates, 10 susceptible, 1 resistant:**
  this run is what caught and fixed an outright wrong assumption — the tool
  previously asserted fosfomycin `Resistant` for every *P. aeruginosa* isolate
  unconditionally, and this real data contradicted that on all 10 susceptible
  isolates on the first run. See
  [`ESKAPE_fos_GOLD_Paeruginosa/RESULTS_SUMMARY.md`](../bioproject_tests/ESKAPE_fos_GOLD_Paeruginosa/RESULTS_SUMMARY.md).
* **PRJNA781811, 18 isolates, measured MICs for both drugs** (Arena *et al.*
  2022, Table 3 — previously a genotype-only comparison, upgraded once the
  paper's own susceptibility table was supplied): fosfomycin specificity
  8/8 (+1 correctly `Indeterminate` for the sole *K. variicola* genome),
  sensitivity 3/9 (`fosA3`, `uhpB`, `glpT` — the same three mechanisms
  already validated in the ESKAPE-GOLD sets). Ceftazidime-avibactam 17/18,
  with one real miss (GMR140) that turned out to be genuinely informative:
  two other isolates in the same set share its exact detectable genotype
  (`blaKPC-3` H274Y, `blaSHV-like` S235G/K236E, the same `ompK36` insertion,
  the same `ompK35` truncation) and are measured susceptible — concrete,
  real-world confirmation that `ompK35` truncation is correctly *not*
  scored as evidence (it appears in 4 of these 18 genomes; only one is
  actually resistant) and a real instance of the documented "expression is
  not measured" limitation. See
  [`PRJNA781811_test/RESULTS_SUMMARY.md`](../bioproject_tests/PRJNA781811_test/RESULTS_SUMMARY.md).
* Two curated mutations sitting in fosfomycin genes but belonging to *other*
  drugs — `cyaA_S352T` (fosmidomycin) and `galU_R101C` (cephalosporin) — have
  explicit regression tests asserting they do **not** produce a fosfomycin call.

## A second opinion from an unrelated tool

Every validation above compares this tool against a study's own reported
genotype or a measured MIC. As a different kind of check, the 63
*K. pneumoniae*/*K. variicola* genomes across the ESKAPE-GOLD fosfomycin sets
and PRJNA781811 (all three now have measured phenotypes) were also run through
[Kleborate](https://github.com/klebgenomics/Kleborate), the
community-standard *K. pneumoniae* genotyping tool, independently built and
maintained, using a different reference database (CARD) entirely.

Kleborate's acquired-gene screen does cover both drugs relevantly: fosfomycin
(`Fcyn_acquired` — acquired fosA-family enzymes, though notably *no*
chromosomal transport-gene tracking at all) and ceftazidime-avibactam
(`Bla_Carb_acquired`/`Bla_ESBL_acquired` — KPC, OXA-48-like, NDM, VIM, IMP,
CTX-M, SHV-ESBL, VEB and more), though it reports gene presence rather than
predicting a phenotype the way this tool does.

Result: **62/63 concordant** on acquired-fosA-gene presence (the one
disagreement is this tool's own already-flagged ambiguous *K. variicola*
isolate — Kleborate independently agreeing there's no confident acquired
call corroborates that `Indeterminate`, it doesn't contradict it) and,
after this cross-check found and this tool fixed two real gaps
(`blaCTX-M` group-9 variants, `blaVEB` — both previously entirely absent),
**62/63 concordant** on CAZ/AVI-relevant beta-lactamase families too. See
[`Kleborate_cross_check/RESULTS_SUMMARY.md`](../bioproject_tests/Kleborate_cross_check/RESULTS_SUMMARY.md)
for the full methodology, including how Kleborate was installed (it is not
a dependency of this project) and a CLI bug that had to be worked around to
get complete results.

## How far this goes

About 112 genomes across four species (*E. coli*, *K. pneumoniae*,
*K. variicola*, *P. aeruginosa*), of which 88 have measured MICs for at least
one of the two drugs, spanning both resistant and susceptible isolates.

**Specificity is well tested for fosfomycin**: 34 real, lab-confirmed
susceptible isolates (24 *K. pneumoniae*, 10 *P. aeruginosa*) produce zero
false Resistant calls, plus one further *K. variicola* isolate that
correctly reads `Indeterminate` rather than falsely `Resistant` — 24/24
correctly `Susceptible` for *K. pneumoniae*, 10/10 honestly `Indeterminate`
(never `Resistant`) for *P. aeruginosa*, where no clinical breakpoint exists
to be susceptible *against*.
**Sensitivity for fosfomycin in *K. pneumoniae* is now a three-times-replicated
measurement with the reference-database gap closed**: 5/38
resistant/intermediate isolates called correctly across three independent
real-MIC samples (2/11 and 0/18 in the two ESKAPE-GOLD rounds, 3/9 in
PRJNA781811 — the same three mechanisms, `fosA3`/`uhpB`/`glpT`, account for
every correct call across all three). The replication says the remaining gap
is not sampling noise: most fosfomycin resistance in this data is not
explained by a coding-sequence change in the 9 transport/regulatory genes
this tool checks, consistent with the literature on promoter/IS-element-driven
`uhpT` regulation, which is outside what a CDS-level tool can see. See
[METHODS.md](METHODS.md#9-known-limits).

**Ceftazidime-avibactam now has two independent real-MIC tests**: CREC
(9/10 resistant + 1 indeterminate, 0 wrong) and PRJNA781811 (17/18, with the
one miss shown above to be a genuine expression-level effect no gene-content
caller could resolve, not a rule this tool got wrong). Susceptible-isolate
CAZ/AVI testing beyond these two sets is still limited — most other sets are
fosfomycin-only or, for the ESKAPE-GOLD sets, contain measured phenotypes
for fosfomycin only.

A single master lookup table,
[`bioproject_tests/GOLD_FOSFOMYCIN_GENOME_PHENOTYPES.tsv`](../bioproject_tests/GOLD_FOSFOMYCIN_GENOME_PHENOTYPES.tsv),
records every genome drawn from the two source CSVs across all three
fosfomycin-GOLD runs above (60 rows: accession, measured phenotype and
testing method, predicted FOS/CAZ-AVI phenotype, fosA locus, concordance),
so future work can look up a genome's result without re-running the
pipeline or re-deriving it from each set's `RESULTS_SUMMARY.md`.
