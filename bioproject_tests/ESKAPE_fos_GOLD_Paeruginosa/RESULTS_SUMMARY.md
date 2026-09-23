# ESKAPE fosfomycin GOLD — *P. aeruginosa* subset

15 real *P. aeruginosa* clinical isolates with measured fosfomycin phenotypes
(10 lab-Susceptible, 4 Intermediate, 1 Resistant), from the same gold-standard
ESKAPE fosfomycin AST source as the *K. pneumoniae* set beside this one — every
isolate with a usable phenotype and NCBI assembly accession for this species
across both source CSVs was included (there were only 24 in total; this is 15
of them, deduplicated).

**This set is why the tool's *P. aeruginosa* fosfomycin logic changed.**

## Method

```bash
fos-cazavi batch -i <genomes>/ -o bioproject_tests/ESKAPE_fos_GOLD_Paeruginosa \
    --organism Pseudomonas_aeruginosa
```

`ground_truth.tsv` carries the full source record for every isolate. Reference
data: AMRFinderPlus 2026-08-07.1.

## What running this set found

Before this set was run, the tool asserted `Resistant` for **every**
*P. aeruginosa* isolate's fosfomycin call, unconditionally, on the stated
grounds that the species is intrinsically resistant with no defined clinical
breakpoint. Run against these 15 real, MIC-tested isolates, that produced:

| Predicted (old logic) | Lab phenotype | Count |
|---|---|---|
| Resistant | Susceptible | **10** |
| Resistant | Intermediate | 4 |
| Resistant | Resistant | 1 |

**10 of 15 — every genuinely susceptible isolate in the set — contradicted.**
A literature check confirmed why: EUCAST publishes only an epidemiological
cut-off (ECOFF) for *Pseudomonas* spp. fosfomycin, explicitly *not* a clinical
breakpoint, citing insufficient outcome data; CLSI does not cover this
species/route at all (EUCAST, *Use of fosfomycin i.v. breakpoints*, May 2024).
Whatever breakpoint produced the lab phenotypes in this dataset is not a
validated clinical one — but asserting the opposite extreme (always resistant)
with equal confidence was just as unsupported, and this real data is what
exposed it: the "intrinsic resistance" claim was falsified on 10/15 isolates
on the very first real validation run against measured phenotypes.

## Fix

`fos_cazavi/phenotype.py` no longer asserts a species-level `Resistant` call
for fosfomycin in *P. aeruginosa*. The rule now:

* **No concrete mechanism found** (only the intrinsic chromosomal `fosA`, no
  acquired enzyme, no transport-gene loss of function): **Indeterminate**,
  stating that no validated clinical breakpoint exists for this species.
  Neither Susceptible nor Resistant is supportable from genotype alone.
* **A concrete mechanism found** (an acquired fosA-family enzyme beyond the
  intrinsic copy, or a transport-gene knockout): still **Resistant** — a real
  mechanism is real evidence regardless of the breakpoint question.

## Results after the fix

| Sample | Lab phenotype | Predicted FOS |
|---|---|---|
| PA_S_01–10 | Susceptible | Indeterminate ×10 |
| PA_I_01–04 | Intermediate | Indeterminate ×4 |
| PA_R_01 | **Resistant** | Indeterminate |

All 15 carry only the intrinsic `fosA_PA1129`, no acquired enzyme. `Indeterminate`
no longer contradicts the 14 non-resistant isolates (it asserts nothing), and
it is the honest answer for `PA_R_01` too: even the one genuinely resistant
isolate here shows no detectable mechanism this tool tracks, which is
consistent with the broader literature that most non-MBL *P. aeruginosa*
fosfomycin resistance is driven by efflux, OprD, or target-modification
mechanisms this tool does not assess — see
[`docs/METHODS.md`](../../docs/METHODS.md#9-known-limits).

This trades a confident, wrong answer for an honest, unhelpful one. That is
the right trade for a tool whose output feeds clinical or research decisions:
`Indeterminate` cannot itself mislead a reader the way a false `Resistant`
label can.

## Reproducing

Accessions are in `accessions.tsv`; full source rows in `ground_truth.tsv`.
