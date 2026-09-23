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

## Results

| Sample | Lab phenotype (FOS) | Predicted FOS | fosA locus | Notes |
|---|---|---|---|---|
| KP_S_01–10 | Susceptible | **Susceptible** ×10 | fosAKP (intrinsic) | ✅ |
| KP_I_01–03 | Intermediate | Susceptible ×3 | fosAKP (intrinsic) | see below |
| KP_R_01–08 | **Resistant** | Susceptible ×8 | fosAKP (intrinsic) | see below |

**Specificity: 10/10.** Every genuinely Susceptible isolate is correctly
called Susceptible, with no false Resistant/Indeterminate call anywhere in
the set.

**Sensitivity: 0/11 (Resistant + Intermediate).** None of the 11 non-susceptible
isolates carry an acquired fosA-family enzyme, and none carries a curated
fosfomycin-resistance point mutation or a detectable loss-of-function change
in any fosfomycin transport/regulatory gene this tool tracks. Every one of the
21 genomes carries only the intact, intrinsic `fosAKP` — genotype gives no
basis to distinguish the 11 resistant isolates from the 10 susceptible ones.

This is a genuine, significant finding, and it traces to two separate causes:

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

### 2. A real, unfixed reference-database gap

The fosfomycin transport/regulatory genes this tool checks for loss of
function and curated point mutations (`uhpT`, `uhpA`, `uhpB`, `uhpC`, `glpT`,
`cyaA`, `ptsI`, `galU`, `murA`) are **all sourced from *E. coli* K-12** in the
bundled nucleotide database — none has a *K. pneumoniae*-specific reference.
Direct BLAST search confirms *K. pneumoniae*'s own orthologs are present,
full-length (97–100% coverage), but only **84–89% nucleotide identity** to the
*E. coli* reference — below this tool's 90% default detection threshold. The
genes are not absent; they are invisible to the detector for this species.

This means fosfomycin resistance mechanisms in these genes are **not actually
checked** for *K. pneumoniae* — the `fosAKP`-only genotype these 21 isolates
show is not "checked and clean," it is "unchecked." The `Susceptible` calls
above are correct read literally (no *scored* resistance mechanism was found),
but they understate what the tool could in principle detect if it had the
right reference. See
[`docs/METHODS.md`](../../docs/METHODS.md#9-known-limits) — this is now
tracked as a known, high-priority limitation, with the concrete fix scoped
(species-specific chromosomal references for these 9 genes, following the same
approach already used for `ompC`/`ompF` in *E. coli* and `ftsI`/`ompK36`/
`ompK35`/`envZ` in *K. pneumoniae*).

### Is 0/11 surprising?

Partly expected, partly not. Fosfomycin resistance in Enterobacterales that
is *not* explained by an acquired `fosA`-family enzyme is a well-documented,
partly-unexplained area in the literature — much of it driven by promoter or
regulatory changes (including IS-element insertions disrupting `uhpT`
expression) rather than simple coding-sequence loss of function, and some
fraction remains mechanistically unexplained even with full genome sequencing.
A tool limited to gene-content and coding-sequence variant calling was never
going to catch all of that. But 0/11 — not even one confirmed hit — is stronger
than "partly unexplained," and given the reference-database gap above, this
result cannot yet be read as "genuinely mechanism-negative." It should be
re-run once *K. pneumoniae*-specific transport-gene references exist.

## What this set changed

* Fixed a real intrinsic/acquired fosA mis-naming bug (above), verified on 4
  real genomes and locked in with 9 unit tests.
* Identified a significant, previously-undocumented reference-database gap:
  fosfomycin transport-gene detection is *E. coli*-only and silently
  ineffective for *K. pneumoniae* at the default identity threshold.
* Established, for the first time, a real **sensitivity** measurement for this
  tool's fosfomycin logic (as opposed to specificity only, which the CREC set
  already covered for MBL/fosA3-driven resistance).

## Reproducing

Accessions are in `accessions.tsv` (used by `scripts/run_validation.sh`); full
source rows including PMID/strain/testing method are in `ground_truth.tsv`.
