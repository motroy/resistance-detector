# CREC fosA3 — validation against measured MICs

Ten carbapenem-resistant *Escherichia coli* (CREC) clinical isolates from
Zhang *et al.*, *Spread of the fosfomycin resistance fosA3 gene via the IS26
mobile element between plasmids and the chromosome of carbapenem-resistant
Escherichia coli in China*, J Glob Antimicrob Resist 42 (2025) 80–87.

**This is the only set here with measured MICs for both of this tool's drugs.**
Everything else in `bioproject_tests/` compares genotype against a study's
reported genotype or its narrative phenotype; this one compares against broth
MICs in the paper's Table 1, with accessions from Table S4.

## Method

```bash
fos-cazavi fos-cazavi-all -a <isolate>.fna -o <isolate> \
    -d fos_cazavi/data/example_database.fasta \
    --organism Escherichia
```

Assemblies were resolved from the Table S4 WGS accessions (see
`accessions.tsv`). Reference data: AMRFinderPlus 2026-08-07.1.

## Results

MICs are from Table 1; FOS ≥256 and CZA >128 mg/L are resistant.

| Isolate | FOS MIC | Predicted FOS | CZA MIC | Predicted CAZ/AVI | Carbapenemase found |
|---|---|---|---|---|---|
| E2129 | >256 | **Resistant** ✅ | >128 | **Resistant** ✅ | blaNDM-5 |
| E1892 | 256 | **Resistant** ✅ | >128 | **Resistant** ✅ | blaNDM-5 |
| E1923 | >256 | **Resistant** ✅ | >128 | **Resistant** ✅ | blaNDM-9 |
| E1985 | 256 | **Resistant** ✅ | >128 | **Resistant** ✅ | blaNDM-5 |
| E2216 | >256 | **Resistant** ✅ | >128 | **Resistant** ✅ | blaNDM-5 |
| E2000 | 256 | **Resistant** ✅ | >128 | **Resistant** ✅ | blaNDM-5 |
| E2111 | >256 | **Resistant** ✅ | >128 | **Resistant** ✅ | blaNDM-5 |
| E2109 | >256 | **Resistant** ✅ | >128 | **Resistant** ✅ | blaNDM-5 |
| E2130 | >256 | **Resistant** ✅ | >128 | **Resistant** ✅ | blaNDM-5 |
| E2257 | >256 | **Resistant** ✅ | >128 | Indeterminate ⚠️ | none |

**Fosfomycin: 10/10 correct.** Every isolate carries `fosA3`, which the tool
calls as an acquired enzyme. *E. coli* has no intrinsic chromosomal fosA, so
there is none of the intrinsic-versus-acquired ambiguity that makes lone fosA
hits uncertain in *Klebsiella*.

**Ceftazidime-avibactam: 9/10 resistant, 1 indeterminate, 0 wrong.** Nine
isolates carry an NDM metallo-beta-lactamase, which avibactam does not inhibit.

## E2257, the one that is not a clean call

E2257 is ceftazidime-avibactam resistant (MIC >128) and carbapenem resistant
(IPM 16, MEM 64, ETP 128) but carries **no carbapenemase** — neither the paper's
Table S4 nor this tool finds one, and a relaxed BLAST search of the assembly
finds no carbapenemase fragment at all.

What it does carry is `blaCMY-2` (AmpC), `blaCTX-M-15`, and a **premature stop
at residue 257 of the 362-residue porin OmpF**. AmpC plus lost outer-membrane
permeability is a documented route to carbapenem and ceftazidime-avibactam
resistance with no carbapenemase involved.

Avibactam does inhibit CMY-2 and CTX-M, so the tool does not assert resistance.
But it no longer reports susceptible either: the call is **Indeterminate**, with
the porin loss and the enzymes it would amplify named in the evidence. Against a
measured MIC of >128 that is the honest answer — the mechanism is visible, its
sufficiency is not established.

Finding this changed the tool: *E. coli* `ompC` and `ompF` were not in the
reference database at all, so the porin defect was invisible. They are now
included, and the porin-amplification rule — previously gated on blaKPC alone —
covers any beta-lactamase whose activity reduced permeability can amplify
(KPC, CMY, CTX-M, SHV, PDC).

Porin status across the set, for context (4 of 10 carry an OmpF defect; the
other three also carry NDM, so it adds context rather than changing their call):

| Isolate | OmpC | OmpF |
|---|---|---|
| E1892 | intact | premature stop at 83 |
| E2129 | intact | premature stop at 83 |
| E2130 | intact | premature stop at 241 |
| E2257 | intact | premature stop at 257 |
| others | intact | intact |

## Limits

Ten isolates from one study and one country, all fosA3-positive and almost all
NDM-positive. It establishes that the fosfomycin side works against real MICs in
*E. coli* and that MBL-driven CAZ/AVI calls are right, but it contains no
fosfomycin-susceptible isolate and no CAZ/AVI-susceptible isolate, so it cannot
measure specificity. The recipient strains in Table 1 (J53, EC600 — FOS 1–2,
CZA 0.25–0.5 mg/L) would serve that purpose but have no deposited genomes.
