# ESKAPE fosfomycin GOLD — *K. pneumoniae* subset

21 real *K. pneumoniae* clinical isolates with **measured fosfomycin
susceptibility phenotypes** (broth/agar dilution or disk diffusion against
EUCAST or CLSI breakpoints), drawn from a curated gold-standard ESKAPE
fosfomycin AST dataset (`source: BV-BRC`) and a broader BV-BRC AMR extract, as
two CSVs of laboratory-confirmed phenotypes with NCBI assembly accessions.
Selected as a stratified sample (fixed seed) of 10 Susceptible, 8 Resistant, 3
Intermediate, to test both sensitivity and specificity — something the
CREC set already did for *E. coli*, and this now does for *K. pneumoniae*.

## Method

```bash
fos-cazavi batch -i <genomes>/ -o bioproject_tests/ESKAPE_fos_GOLD_Kpneumoniae \
    --organism Klebsiella_pneumoniae
```

`ground_truth.tsv` carries the full source record for every isolate (measured
phenotype, testing method, standard, PMID, strain, source file). Reference
data: AMRFinderPlus 2026-08-07.1.

## Results (current)

| Sample | Lab phenotype (FOS) | Predicted FOS | fosA locus | Notes |
|---|---|---|---|---|
| KP_S_01–10 | Susceptible | **Susceptible** ×10 | fosAKP (intrinsic) | ✅ |
| KP_I_01–03 | Intermediate | Susceptible ×3 | fosAKP (intrinsic) | see below |
| KP_R_02, KP_R_05 | **Resistant** | **Resistant** ×2 | fosAKP (intrinsic) | uhpB premature stop, GAMMA-confirmed ✅ |
| KP_R_01, KP_R_03, KP_R_04, KP_R_06, KP_R_07, KP_R_08 | **Resistant** | Susceptible ×6 | fosAKP (intrinsic) | see below |

**Specificity: 10/10.** Every genuinely Susceptible isolate is correctly
called Susceptible, with no false Resistant/Indeterminate call anywhere in
the set — unchanged by everything below.

**Sensitivity: 2/11 (Resistant + Intermediate), up from 0/11.** This set went
through two real fixes since first committed; both are described in full
below, but the short version: a hit-selection bug that hid the genome's own
`fosAKP` behind an ambiguous acquired-allele name (fixed, see "1" below), and
a reference-database gap that made 9 fosfomycin transport/regulatory genes
undetectable in this species at all (fixed, see "2" below). Fixing the second
one is what moved sensitivity from 0/11 to 2/11: `KP_R_02` and `KP_R_05` — both
lab-confirmed Resistant by broth dilution — carry a premature stop at the same
residue (353 of 491) in `uhpB`, the sensor histidine kinase that induces
`uhpT` transcription. Both calls are independently confirmed by GAMMA's own
alignment ("truncation at codon 353"), agreeing exactly with the BLAST-based
call. The remaining 9 non-susceptible isolates carry no acquired fosA-family
enzyme and no coding-sequence loss of function in any of the 9 transport
genes — see "Is 2/11 the ceiling?" below for why that is not surprising.

### 1. A real detection bug, found and fixed here

Four of these 21 genomes (`KP_R_01`, `KP_R_02`, `KP_S_01`, `KP_S_10`) initially
had their single fosA copy reported as `fosA5-like` or `fosA10-like` — an
**ambiguous, ostensibly-acquired** call — rather than as the intrinsic
`fosAKP`. Investigation showed the intrinsic reference *did* hit the identical
genomic span in every case, at 94.9–99.1% identity, but lost the bitscore race
to a family member (`fosA5`, `fosA10`) by less than 2 percentage points,
because `fosAKP` is defined from a single reference strain and other lineages'
native copy is sometimes, by chance, a point or two closer to a different
catalogued allele than to that one reference. Since each genome carries
exactly one fosA locus, it can only be the native chromosomal gene — there is
nothing else it could be.

This was a real bug in hit selection, not a labelling nuance: it would have
reported four genuinely wild-type isolates as ambiguous/possibly-resistant.
Fixed in `BlastDetector._prefer_intrinsic_naming()` — when an acquired-named
hit is within 3 percentage points of a competing intrinsic-gene hit at the
same span, the intrinsic name and reference win. See
[`fos_cazavi/acquired.py`](../../fos_cazavi/acquired.py) and
[`tests/test_acquired.py`](../../tests/test_acquired.py) for the fix and its
regression tests. All four now correctly report `fosAKP`.

### 2. A real reference-database gap — now fixed

The fosfomycin transport/regulatory genes this tool checks for loss of
function and curated point mutations (`uhpT`, `uhpA`, `uhpB`, `uhpC`, `glpT`,
`cyaA`, `ptsI`, `galU`, `murA`) were **all sourced from *E. coli* K-12** in the
bundled nucleotide database — none had a *K. pneumoniae*-specific reference.
Direct BLAST search confirmed *K. pneumoniae*'s own orthologs are present,
full-length (97–100% coverage), but only **84–89% nucleotide identity** to the
*E. coli* reference — below this tool's 90% default detection threshold. The
genes were not absent; they were invisible to the detector for this species.

**The fix**: a second, *K. pneumoniae*-specific reference for each of the 9
genes was added, sourced from *K. pneumoniae* subsp. *pneumoniae* HS11286
(RefSeq `NC_016845.1`, a complete, closed reference genome), located by
identifying the corresponding annotated locus in that genome's own feature
table (e.g. `uhpT`/`uhpC`/`uhpB`/`uhpA` as the adjacent four-gene operon
matching the product descriptions "hexose phosphate transport protein",
"regulatory protein UhpC", "two-component regulatory system sensor histidine
kinase" and "DNA-binding transcriptional activator UhpA"). Each new reference
matches the *E. coli* CDS at 77–89% nucleotide identity — the same range the
diagnostic BLAST search had already found — confirming these are the correct
orthologs, and matches real clinical *K. pneumoniae* assemblies (tested across
Resistant, Intermediate and Susceptible isolates from this set) at
**98.6–100% identity**, comfortably above the 90% threshold.

Because a BLAST nucleotide database needs a unique sequence ID per entry, the
new reference is stored under its own name (`uhpT_Kpn`, `uhpB_Kpn`, ...)
rather than replacing the *E. coli* entry, and mapped back to the canonical
gene name (`uhpT`, `uhpB`, ...) everywhere else in the tool via
[`fos_cazavi/references.py::GENE_ALIASES`](../../fos_cazavi/references.py) —
the Gene column, `FOS_TRANSPORT_GENES` membership, loss-of-function reporting
and GAMMA cross-checking all see one name regardless of which reference
actually matched. See
[`tests/test_references.py`](../../tests/test_references.py) and
[`tests/test_pipeline.py::TestKlebsiellaPneumoniaeFosfomycinTransportGenes`](../../tests/test_pipeline.py)
for the regression tests, including one asserting a locus is never
double-counted now that two differently-named references can match it.

### Is 2/11 the ceiling?

Likely close to it for this tool's method, and that is an expected, literature
consistent result rather than a disappointing one. Fosfomycin resistance in
Enterobacterales that is *not* explained by an acquired `fosA`-family enzyme
is a well-documented, partly-unexplained area in the literature — much of it
driven by promoter or regulatory changes (including IS-element insertions
disrupting `uhpT` expression) rather than coding-sequence loss of function,
and some fraction remains mechanistically unexplained even with full genome
sequencing. A tool limited to gene-content and coding-sequence variant calling
was never going to catch all of that; `uhpB`'s premature stop is exactly the
kind of mechanism it *can* catch, which is why fixing the detection gap
recovered it. The remaining 9 isolates' resistance likely lies upstream of
what a CDS-level caller can see.

**Update:** a second, independent 24-genome sample from the same source data
(disjoint accessions, weighted toward resistant/intermediate isolates) found
one more case of the same kind (a `glpT` premature stop, in a genotype-only
isolate from a different validation set) but no new cases within its own 18
non-susceptible isolates — 0/18, for a combined **2/29** across both rounds.
See
[`ESKAPE_fos_Kpneumoniae_round2/RESULTS_SUMMARY.md`](../ESKAPE_fos_Kpneumoniae_round2/RESULTS_SUMMARY.md).

## What this set changed

* Fixed a real intrinsic/acquired fosA mis-naming bug (above), verified on 4
  real genomes and locked in with 9 unit tests.
* Identified, and then fixed, a significant reference-database gap:
  fosfomycin transport-gene detection was *E. coli*-only and silently
  ineffective for *K. pneumoniae* at the default identity threshold. Fixing it
  recovered 2 real, GAMMA-confirmed resistant calls that were previously
  invisible.
* Established a real **sensitivity** measurement for this tool's fosfomycin
  logic in *K. pneumoniae* (as opposed to specificity only, which the CREC set
  already covered for MBL/fosA3-driven resistance in *E. coli*).

## Reproducing

Accessions are in `accessions.tsv` (used by `scripts/run_validation.sh`); full
source rows including PMID/strain/testing method are in `ground_truth.tsv`.
