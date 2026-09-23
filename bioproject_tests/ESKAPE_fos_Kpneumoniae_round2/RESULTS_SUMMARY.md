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

**Sensitivity: 0/18 (Resistant + Intermediate), unchanged by the
transport-gene reference fix.** This round originally replicated round 1's
0/11 exactly, at a time when the fosfomycin transport-gene references this
tool checks (`uhpT`, `glpT`, `cyaA`, `ptsI`, `galU`, etc.) were *E.
coli*-only and *K. pneumoniae*'s own orthologs sat below the 90% identity
detection threshold — a real, systematic detection gap, confirmed not to be
sampling noise by this very replication. That gap is now fixed: a
*K. pneumoniae*-specific reference was added for each of the 9 genes (see
[`ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md`](../ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md#2-a-real-reference-database-gap--now-fixed)),
and re-running this set with the fix applied still finds **zero** LOF in any
of the 9 genes across all 18 non-susceptible isolates here — the genes are no
longer invisible, they were checked and found intact. This round's 18
isolates' resistance is not explained by a coding-sequence change in these
specific genes; round 1's two `uhpB`-truncated isolates show that mechanism
is real and detectable when present, it simply is not what is driving
resistance in *this* round's isolates. No new detection bug was found in
this round either time — the intrinsic-naming fix from round 1 held: all 24
genomes' single fosA locus resolved cleanly to `fosAKP` with no
ambiguous/acquired mis-calls.

**This round's contribution**: it turned a single-sample finding (n=11) into
a replicated one (n=29 combined), which is what justified investing in the
reference fix rather than dismissing 0/11 as one unlucky sample. Now that the
fix is in, this round also shows the fix's limit honestly: closing a real
detection gap does not manufacture sensitivity that was never there to find —
combined sensitivity across both rounds is **2/29** (both from round 1).

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
