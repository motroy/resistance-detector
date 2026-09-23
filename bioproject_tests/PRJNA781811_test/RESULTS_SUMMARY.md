# PRJNA781811 — FOS-CAZAVI Resistance Detection Results

18 *Klebsiella pneumoniae* / *K. variicola* assemblies from BioProject PRJNA781811
(listed in `PRJNA781811.ncbi_datasets.tsv`; see
[Arena et al. 2022, Front. Microbiol. 13:983294](https://doi.org/10.3389/fmicb.2022.983294)),
run through the current pipeline and compared against the paper's own
**measured MICs** (Table 3, `Arena2022_Table3_measured_MICs.xlsx`), mapped
from the paper's `GMR###` strain IDs to NCBI accessions via each assembly's
BioSample "Sample name" attribute. This was previously a genotype-only
comparison; the measured-phenotype table turns it into a real accuracy
check for both drugs, like the CREC and ESKAPE-GOLD sets.

## Command used

```bash
datasets download genome accession <accession> --include genome
unzip ncbi_dataset.zip

fos-cazavi fos-cazavi-all \
    -a ncbi_dataset/data/<accession>/<accession>_*_genomic.fna \
    -o <accession> \
    -d fos_cazavi/data/example_database.fasta \
    --organism Klebsiella_pneumoniae
```

Reference data: AMRFinderPlus 2026-08-07.1 (see `fos_cazavi/data/DATA_VERSION.txt`).
Breakpoints applied to the paper's raw MICs: EUCAST Enterobacterales
fosfomycin (IV) S≤32/R>32 mg/L; ceftazidime-avibactam S≤8/R>8 mg/L (read
from the ceftazidime component of the paper's "ceftazidime/avibactam 4mg/L"
notation, e.g. `2/4` → ceftazidime MIC 2).

## Results

| Accession | GMR strain | Beta-lactamases | fosA-family | FOS measured | FOS predicted | CAZ/AVI measured | CAZ/AVI predicted |
|---|---|---|---|---|---|---|---|
| GCA_027151785.1 | GMR149 | SHV-like | fosAKP (intrinsic) + `uhpB` LOF | **Resistant** | **Resistant** ✅ | Susceptible | Susceptible ✅ |
| GCA_027151795.1 | GMR147 | SHV-like | fosAKP (intrinsic) | Susceptible | Susceptible ✅ | Susceptible | Susceptible ✅ |
| GCA_027151835.1 | GMR152 | KPC-3, SHV-like | fosAKP (intrinsic) | **Resistant** | Susceptible ❌ | Susceptible | Susceptible ✅ |
| GCA_027151845.1 | GMR150 | CTX-M-65, SHV-like | **fosA3** + intrinsic fosAKP | **Resistant** | **Resistant** ✅ | Susceptible | Susceptible ✅ |
| GCA_027151875.1 | GMR153 | KPC-3, CTX-M-15, CMY-like, SHV-like | fosAKP (intrinsic) | Susceptible | Susceptible ✅ | Susceptible | Susceptible ✅ |
| GCA_027151985.1 | GMR148 | SHV-like | fosAKP (intrinsic) | Susceptible | Susceptible ✅ | Susceptible | Susceptible ✅ |
| GCA_027151995.1 | GMR145 | CTX-M-15, SHV-like | fosAKP (intrinsic) | Susceptible | Susceptible ✅ | Susceptible | Susceptible ✅ |
| GCA_027152005.1 | GMR144 | CTX-M-15, SHV-like | fosAKP (intrinsic) | Susceptible | Susceptible ✅ | Susceptible | Susceptible ✅ |
| GCA_027152065.1 | GMR140 | KPC-3, SHV-like | fosAKP (intrinsic) | **Resistant** | Susceptible ❌ | **Resistant** | Susceptible ❌ |
| GCA_027152105.1 | GMR146 | CTX-M-15, SHV-like | fosAKP (intrinsic) | Susceptible | Susceptible ✅ | Susceptible | Susceptible ✅ |
| GCA_027152185.1 | GMR139 | SHV-like | fosAKP (intrinsic) | **Resistant** | Susceptible ❌ | Susceptible | Susceptible ✅ |
| GCA_027152205.1 | GMR142 | SHV-like | fosAKP (intrinsic) | **Resistant** | Susceptible ❌ | Susceptible | Susceptible ✅ |
| GCA_027152215.1 | GMR141 | SHV-like | fosAKP (intrinsic) | Susceptible | Susceptible ✅ | Susceptible | Susceptible ✅ |
| GCA_027152225.1 | GMR136 | KPC-2, CTX-M-15, SHV-like | fosAKP (intrinsic) | **Resistant** | Susceptible ❌ | Susceptible | Susceptible ✅ |
| GCA_027152245.1 | GMR135 | CTX-M-15, SHV-like | fosAKP (intrinsic) | Susceptible | Susceptible ✅ | Susceptible | Susceptible ✅ |
| GCA_027152405.1 | GMR134 | SHV-like | fosAKP (intrinsic) | **Resistant** | Susceptible ❌ | Susceptible | Susceptible ✅ |
| GCA_027152445.1 | GMR132 | KPC-3, SHV-like | fosAKP (intrinsic) + `glpT` LOF | **Resistant** | **Resistant** ✅ | Susceptible | Susceptible ✅ |
| GCA_027152495.1 (*K. variicola*) | GMR133 | SHV-like | fosA9, no intrinsic copy | Susceptible | Indeterminate — | Susceptible | Susceptible ✅ |

## Fosfomycin: 11/18 exact match (3/9 Resistant, 8/8 Susceptible + 1 Indeterminate)

**Specificity holds**: every genuinely susceptible isolate is called
Susceptible (or, for the one *K. variicola* genome, honestly `Indeterminate`
rather than a false Resistant — see below). **Sensitivity is 3/9**: the
`uhpB`/`glpT`/`fosA3` mechanisms this tool can detect account for 3 of the 9
resistant isolates; the other 6 carry no acquired fosA-family enzyme and no
loss of function in any of the 9 transport/regulatory genes this tool
tracks. This is not a new finding — it is the same, already-documented
pattern from the ESKAPE-fosfomycin GOLD sets (sensitivity 2/29 there), now
confirmed on a third, independent, real-MIC dataset: most fosfomycin
resistance in this data is not explained by a coding-sequence change this
tool can see, consistent with the literature on promoter/IS-element-driven
`uhpT` regulation. See
[`ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md`](../ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md).

**GCA_027152495.1 (*K. variicola*, GMR133)** is measured Susceptible and
called `Indeterminate` — not a miss in the sense of asserting the wrong
answer, but a case where this tool declines to assert Susceptible because
its single fosA-family hit cannot be resolved as intrinsic-vs-acquired for
a species with no matching intrinsic reference (see "Changes" below). Read
against ground truth, `Indeterminate` here happens to be the conservative
side of correct.

## Ceftazidime-avibactam: 17/18 exact match — one real miss, and it's instructive

**GCA_027152065.1 (GMR140)** is measured CAZ/AVI-resistant (ceftazidime MIC
>8) but called Susceptible: it carries `blaKPC-3` (H274Y only, no
Omega-loop/237-243/insertion-loop escape mutation) and `blaSHV-like`
(S235G, K236E), with `ompK36` carrying an in-frame insertion (not loss of
function) and `ompK35` truncated by a premature stop.

The striking part: **two other isolates in this same set —
GCA_027151835.1 (GMR152) and GCA_027152445.1 (GMR132) — have the identical
detectable genotype** (same `blaKPC-3` H274Y, same `blaSHV-like` S235G/K236E,
same `ompK36` insertion, same `ompK35` premature stop at residue 89) and are
both measured CAZ/AVI-**susceptible**. Genotype-only calling cannot
distinguish these three isolates from each other; whatever separates the
resistant one is not a coding-sequence difference this tool - or, on this
evidence, probably any gene-content-based caller - can see. This is a
concrete, real-world instance of the "Expression" limitation already
documented in `docs/METHODS.md` (section 9): gene copy number and promoter
strength affect real MICs and are not measured here.

This is also why `ompK35` truncation is deliberately *not* scored as
contributory evidence (`docs/METHODS.md`, section 5): it is present in 4 of
these 18 genomes (GMR152, GMR140, GMR132, and GMR150) but only one of those
four is actually CAZ/AVI-resistant. Scoring it would have flagged three
genuinely susceptible isolates as uncertain to catch one resistant one - a
net loss on this evidence, and this real data is a direct confirmation of
that design choice, not a reason to revisit it.

## Changes from earlier runs of this BioProject

These results have gone through four corrections/additions since first committed:

* **GCA_027152215.1** was originally reported as carrying `fosA5` at 96.19%
  identity and called fosfomycin-resistant. Allele assignment was fixed to
  choose by bitscore across all fosA references rather than from a truncated
  hit list, and the same locus matched intrinsic `fosAKP` better — susceptible.
* Once fixed, **GCA_027151985.1**, **GCA_027152185.1** and **GCA_027152495.1**
  still had their single fosA locus reported as an ambiguous acquired allele
  (`fosA10-like`, `fosA5-like`, `fosA9-like`) rather than the intrinsic gene,
  because the intrinsic reference lost the bitscore race by a hair at the
  identical span — a real hit-selection bug, not a labelling choice. All three
  were reported `Indeterminate`.
* That bug is now fixed (`BlastDetector._prefer_intrinsic_naming()`; see
  [`ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md`](../ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md)
  for the full account, found via real MIC-tested genomes). **GCA_027151985.1**
  and **GCA_027152185.1** now correctly resolve to intrinsic `fosAKP` —
  susceptible. **GCA_027152495.1** does not: it is *K. variicola*, a different
  species from the `fosAKP` reference strain, and its native copy sits outside
  the fix's identity tolerance — it remains genuinely ambiguous and stays
  `Indeterminate`.
* blaKPC alleles are now named from the observed protein changes against
  KPC-2, rather than from whichever KPC reference happened to win the BLAST
  hit.
* **GCA_027151785.1** and **GCA_027152445.1** report fosfomycin `Resistant`
  since the *K. pneumoniae*-specific transport-gene references were added
  (`uhpB` and `glpT` premature stops respectively, both GAMMA-confirmed) —
  see
  [`ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md`](../ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md#2-a-real-reference-database-gap--now-fixed).
  Both are now known to be genuinely correct: measured Resistant.
* **Beta-lactamase families widened** (`blaCTX-M`, `blaVEB`) after a Kleborate
  cross-check found gaps — GCA_027151845.1's `blaCTX-M-65` is now detected
  (previously invisible; see
  [`Kleborate_cross_check/RESULTS_SUMMARY.md`](../Kleborate_cross_check/RESULTS_SUMMARY.md)).
  No phenotype changes resulted for this set.
* **Measured phenotypes added** (this update): Table 3 from the source paper
  gives real fosfomycin and ceftazidime-avibactam MICs for every isolate,
  mapped to accessions via BioSample "Sample name". This turns the
  comparison from genotype-only into a real accuracy check; see the two
  sections above for the result.

See [../../docs/METHODS.md](../../docs/METHODS.md) for how allele assignment
and the phenotype rules work.

## Reproducing

`ground_truth.tsv` carries the measured phenotype, raw MIC, GMR strain ID and
ST for every isolate. Source: `Arena2022_Table3_measured_MICs.xlsx` (Table 3
of the paper, as supplied), mapped to accessions by matching each assembly's
BioSample "Sample name" attribute (fetched via the NCBI Datasets API) to the
paper's `GMR###` strain column.
