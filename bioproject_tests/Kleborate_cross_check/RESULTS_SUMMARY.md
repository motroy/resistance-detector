# Cross-check against Kleborate

An independent, third-party tool comparison, run on request: does
[Kleborate](https://github.com/klebgenomics/Kleborate) (v3.2.4) — the
community-standard genomic surveillance tool for the *Klebsiella
pneumoniae* species complex — agree with this tool's fosfomycin calls, and
does it even assess ceftazidime-avibactam-relevant genes at all?

Run across all 63 *K. pneumoniae* / *K. variicola* genomes in this repo's
validation set with either a measured phenotype or a real clinical source
(`ESKAPE_fos_GOLD_Kpneumoniae`, `ESKAPE_fos_Kpneumoniae_round2`,
`PRJNA781811_test`).

## What Kleborate actually checks

Kleborate's `klebsiella_pneumo_complex__amr` module BLASTs against a
curated subset of [CARD](https://card.mcmaster.ca/) (v3.2.9), plus
dedicated modules for OmpK35/36, QRDR (fluoroquinolone), SHV, colistin
(MgrB/PmrB) mutations. Relevant to this tool's two drugs:

* **Fosfomycin**: a `Fcyn_acquired` column covering acquired `fosA2/A3/A4/
  A7/A8`, `fosB`, `fosC/C2`, `fosD`, `fosK`, `fosL1`, `fosX` (23 alleles
  total). It does **not** track the plain `fosA`/`fosA5`/`fosA6`/`fosA9`-
  `fosA13` names this tool also carries, and — importantly — **it has no
  equivalent of `uhpT`/`uhpA`/`uhpB`/`uhpC`/`glpT`/`cyaA`/`ptsI`/`galU`/
  `murA` at all**: no chromosomal fosfomycin transport/regulatory gene
  tracking whatsoever, acquired or otherwise. So yes, this is a real
  independent check on the acquired-fosA side, and no, it cannot corroborate
  or contradict this tool's new `uhpB`/`glpT` loss-of-function findings
  (see `../ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md`) — that
  mechanism is simply outside what Kleborate looks for.
* **Ceftazidime-avibactam**: **yes**, Kleborate does search for the
  relevant genes, under `Bla_Carb_acquired` (carbapenemases — KPC, the
  OXA-48-like/OXA-23-like/other class D families, NDM, VIM, IMP, and more:
  620 CARD alleles) and `Bla_ESBL_acquired`/`Bla_ESBL_inhR_acquired`
  (CTX-M, SHV-ESBL variants, VEB, etc.: 540+23 alleles). It reports gene
  **presence**, the same way this tool's `Gene`/`Allele` columns do — it
  does not predict a CAZ/AVI phenotype itself (no Omega-loop/insertion-loop
  escape-mutation interpretation, no avibactam-inhibition reasoning). So it
  is a genuine, useful cross-check on this tool's gene-detection *input*,
  not a second phenotype predictor to compare final calls against.

## Setup (not a standard dependency of this repo)

Kleborate is **not** part of this project's environment — it was installed
only for this one-off comparison, in an isolated virtualenv, and is not
required to run or test this tool:

```bash
python3 -m venv kleborate_venv
kleborate_venv/bin/pip install "setuptools==68.2.2"  # kleborate's `mash` dependency's
                                                       # setup.py needs pre-modern distutils
kleborate_venv/bin/pip install kleborate
kleborate_venv/bin/pip install "kaptive==3.2.2"       # kleborate 3.2.4 pins no kaptive
                                                       # version; the latest (3.3.x) removed
                                                       # the `kaptive.database` module it imports
kleborate_venv/bin/pip install pandas                 # declared as a dependency but not
                                                       # actually installed by pip
apt-get install -y mash minimap2                      # external binaries kleborate shells out to
```

Kleborate 3.2.4's own CLI has a real bug when modules are selected with
`-m module_a,module_b` (not a preset): the header-filtering logic in
`__main__.py` only keeps the *first* listed module's output columns,
silently dropping the rest with no error. The module functions themselves
are correct — confirmed by calling `get_results()` directly. This
comparison was run that way (see `run_kleborate.py`, described below)
rather than through the CLI.

## Results

Full per-genome comparison: [`comparison.tsv`](comparison.tsv).

### Fosfomycin: acquired fosA-family gene presence

**62/63 concordant.** The one disagreement:

| Sample | Kleborate | This tool |
|---|---|---|
| GCA_027152495.1 (*K. variicola*) | No acquired gene (`Fcyn_acquired: -`) | `fosA-like`, reported **Indeterminate** (not Resistant) |

This is not a contradiction — it is independent corroboration. This tool
already flags this exact isolate as genuinely ambiguous: it is the sole
*K. variicola* genome in the set, its single fosA-family hit is not close
enough to the (*K. pneumoniae*-sourced) `fosAKP` reference to resolve as
intrinsic, and it is reported `Indeterminate` rather than asserted as
acquired (see `../PRJNA781811_test/RESULTS_SUMMARY.md`). Kleborate,
independently, also does not call it an acquired gene. Neither tool
confidently asserts resistance here, which is exactly the honest outcome.

### Ceftazidime-avibactam-relevant beta-lactamases

Comparing gene *families* found (not exact allele numbers — this tool only
does allele-level typing for blaKPC, by design; see
`docs/METHODS.md`), two real gaps turned up:

1. **`blaCTX-M-65`** (GCA_027151845.1) — undetected. This tool carried only
   one CTX-M reference (CTX-M-15, phylogenetic group 1); CTX-M-65 is group 9,
   ~80% nucleotide identity to CTX-M-15 — below even this tool's relaxed
   85% screening cutoff.
2. **`blaVEB-1`** (KP2_R_04, KP_R_06) — undetected. No VEB reference existed
   in the database at all.
3. **`blaOXA-23`-like** (KP2_R_03) — undetected. Already a known, documented
   exclusion (`docs/METHODS.md`, "Acquired class D oxacillinases other than
   the OXA-48-like group... are not in the database") — this cross-check
   found a concrete real-world instance of it, not a new gap.

**Fixed** (1 and 2): `blaCTX-M` and `blaVEB` were expanded from single
representative alleles to full families in
[`fos_cazavi/build_data.py`](../../fos_cazavi/build_data.py) (the same
approach already used for the metallo-beta-lactamase families). CTX-M
already had its family-level phenotype treatment in place
(`PERMEABILITY_AMPLIFIED_FAMILIES` keys on `gene_family()`, which strips
the allele number regardless of phylogenetic group), so this needed no
phenotype-logic change — re-running the full 112-genome validation suite
confirmed **zero phenotype changes anywhere** and 16 genomes across four
different validation sets (not just the K. pneumoniae ones) gained a real
or more precisely named beta-lactamase call. VEB was added for
detection/visibility only, not yet wired into the permeability-amplification
evidence — that would need the same kind of literature check already done
for the other families in that list, deliberately left as follow-up rather
than guessed.

**Not fixed** (3): OXA-23-like carbapenemases have different avibactam
susceptibility than the OXA-48-like family this tool already treats as
avibactam-inhibited (class D carbapenemases are not uniformly
avibactam-susceptible), so adding it requires the same literature
verification this tool's other phenotype rules were built on, not a
database-only change. Left as a documented, scoped follow-up.

After the fix: **62/63 concordant** on beta-lactamase gene families (the
one remaining gap is the already-documented OXA-23-like exclusion).

## Reproducing

```bash
python3 -m venv /tmp/kleborate_venv
/tmp/kleborate_venv/bin/pip install "setuptools==68.2.2" && \
  /tmp/kleborate_venv/bin/pip install kleborate && \
  /tmp/kleborate_venv/bin/pip install "kaptive==3.2.2" pandas
apt-get install -y mash minimap2

/tmp/kleborate_venv/bin/python3 run_kleborate.py   # writes /tmp/kleborate_direct_results.tsv
```

`run_kleborate.py` in this directory calls Kleborate's species and
K. pneumoniae-AMR module functions directly (working around the CLI bug
described above) across every genome in the three sets listed at the top.
