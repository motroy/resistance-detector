# ESKAPE fosfomycin GOLD — *K. pneumoniae*, round 2 (replication set)

24 more real *K. pneumoniae* clinical isolates with **measured fosfomycin
susceptibility phenotypes**, drawn from the same two source CSVs as
[`ESKAPE_fos_GOLD_Kpneumoniae`](../ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md)
(`ESKAPE_fosfomycin_GOLD_validation_set.csv`, `bvbrc_amr_LAB_CONFIRMED_with_accessions_FIXED.csv`),
but a disjoint set of accessions — none of these 24 were used in that
first round. Selected as a stratified sample (fixed seed) of 12 Resistant, 6
Intermediate, 6 Susceptible, deliberately weighted toward the
non-susceptible isolates: the first round found 0/11 sensitivity, and the
question here is whether that holds on an independent sample or was an
artefact of the first 11.

## Method

```bash
fos-cazavi batch -i <genomes>/ -o bioproject_tests/ESKAPE_fos_Kpneumoniae_round2 \
    --organism Klebsiella_pneumoniae
```

`ground_truth.tsv` carries the full source record for every isolate.
Reference data: AMRFinderPlus 2026-08-07.1 (unchanged from round 1).

## Results

| Sample group | Lab phenotype (FOS) | Predicted FOS | fosA locus |
|---|---|---|---|
| KP2_S_01–06 | Susceptible | **Susceptible** ×6 | fosAKP (intrinsic) |
| KP2_I_01–06 | Intermediate | Susceptible ×6 | fosAKP (intrinsic) |
| KP2_R_01–12 | **Resistant** | Susceptible ×12 | fosAKP (intrinsic) |

**Specificity: 6/6.** Every genuinely susceptible isolate is correctly
called Susceptible.

**Sensitivity: 0/18 (Resistant + Intermediate).** Identical to round 1's
0/11. Combined across both rounds: **0/29** non-susceptible *K. pneumoniae*
isolates from this data source have been called correctly by genotype
alone, out of 29 tested. Every one of the 24 genomes here carries only the
intact, intrinsic `fosAKP` and no curated fosfomycin resistance mutation —
same picture as round 1, and for the same documented reason: the
fosfomycin transport-gene references this tool checks (`uhpT`, `glpT`,
`cyaA`, `ptsI`, `galU`, etc.) are *E. coli*-only, and *K. pneumoniae*'s own
orthologs sit below the 90% identity detection threshold (see
[`ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md`](../ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md#2-a-real-unfixed-reference-database-gap)
and [`docs/METHODS.md`](../../docs/METHODS.md#9-known-limits) for the full
account). No new detection bug was found in this round — the intrinsic-naming
fix from round 1 held: all 24 genomes' single fosA locus resolved cleanly to
`fosAKP` with no ambiguous/acquired mis-calls.

**This round's contribution**: it turns a single-sample finding (n=11) into
a replicated one (n=29, across two independently drawn, non-overlapping
samples from two different source files). The 0/11 in round 1 was not
sampling noise — the reference-database gap is real and systematic, not
isolate-specific. This is the strongest evidence yet that the
species-specific transport-gene reference fix scoped in `docs/METHODS.md`
is the correct next investment, not a rare edge case.

### Ceftazidime-avibactam (not validated here, but observed)

This dataset has no measured CAZ/AVI phenotype (it is a fosfomycin-only
extract), so these are unvalidated genotype calls, included for reference
only. Six isolates carry a metallo-beta-lactamase — `blaNDM-5`
(KP2_I_02, KP2_I_04, KP2_I_05, KP2_I_06) or `blaVIM-1`/`blaVIM-19`
(KP2_R_09, KP2_R_10) — all at 100% identity, and are called CAZ/AVI
Resistant on that basis (avibactam does not inhibit metallo-enzymes). The
remaining 18 are called Susceptible. Several isolates also carry an
`ompK35` premature-stop loss-of-function (porin loss), which the tool
records but does not score toward CAZ/AVI resistance without an
accompanying beta-lactamase — consistent with the documented
porin-amplification rule.

## Reproducing

```bash
scripts/run_validation.sh -s ESKAPE_fos_Kpneumoniae_round2 -j 4
```

Accessions are in `accessions.tsv`; full source rows (measured phenotype,
testing method, standard, PMID, strain, source file) are in
`ground_truth.tsv`.
