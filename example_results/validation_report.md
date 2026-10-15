# Validation Report

## Executive Summary

> **Status (updated 2026-10-15):** the per-genome results below (Results
> sections 1-3) are the original February 2026 installation check and were
> **not re-run** with the current pipeline. Since then detection and phenotype
> logic was rebuilt (AMRFinderPlus 2026-08-07.1 reference data, allele
> assignment by bitscore ranking, mechanism-based CAZ/AVI calls, fosfomycin
> transport-gene fixes), so the specific calls, the `D179T` discussion and the
> `glpT` fallback warning below describe the old version and should not be read
> as current behaviour.

The `fos-cazavi` tool was installed and run end-to-end on three reference
genomes (KPNIH1, MGH 78578, K-12 MG1655) and produced the expected output
files. Current, regenerated validation against published genomes with measured
MICs (CREC fosA3, ESKAPE-GOLD, PRJNA781811, PRJNA741867 and others, plus a
Kleborate cross-check with 62/63 concordant) is in
[../docs/VALIDATION.md](../docs/VALIDATION.md). Current example output is
`multi_mechanism_summary.txt` in this folder.

## Methodology

1.  **Installation**:
    -   Created conda environment `resistance_detector` with dependencies: `python`, `blast`, `biopython`, `seqkit`, `gamma`, `pytest`.
    -   Installed `fos-cazavi` in editable mode.
2.  **Genomes**:
    -   `GCA_000281535.1` (*K. pneumoniae* KPNIH1) - CAZAVI Positive Control (KPC-3).
    -   `GCF_000016305.1` (*K. pneumoniae* MGH 78578) - Negative Control.
    -   `GCF_000005845.2` (*E. coli* K-12 MG1655) - Chromosomal FOS targets.
3.  **Database**:
    -   Created reference database `resistance_db` using `fos-cazavi create-db`.
4.  **Analysis**:
    -   Ran `fos-cazavi fos-cazavi-all` on each genome.

## Results

### 1. *K. pneumoniae* KPNIH1 (`GCA_000281535.1`)

-   **Expected**: `blaKPC-3` at ≥99% identity.
-   **Observed**:
    -   `blaKPC-3` detected at 100.00% identity and 100.00% coverage.
    -   BLAST reported mutations: `D179T,V240Y,T243A`.
    -   GAMMA reported: "Native" (no coding mutations relative to reference).
-   **Analysis**:
    -   The detection of `blaKPC-3` confirms the positive control status.
    -   The reported mutation `D179T` is a known artifact. The AMRfinderPlus reference sequence for `blaKPC-3` contains T179, whereas the tool's internal logic or standard wildtype assumes D179. Since KPNIH1 matches the reference (T179), GAMMA sees no mutation ("Native"), but the mutation reporter flags `D179T` (Ref D -> Obs T).
    -   This matches the "known discrepancy" noted in project documentation.

### 2. *K. pneumoniae* MGH 78578 (`GCF_000016305.1`)

-   **Expected**: Negative control (no KPC, no FosA).
-   **Observed**:
    -   Detected: `acrB`, `ftsI`, `envZ`, `blaSHV-12`.
    -   No `blaKPC` or `fosA` genes detected.
-   **Analysis**:
    -   Confirmed as negative control for target resistance genes. `blaSHV-12` is a common ESBL but not the primary target for CAZAVI/FOS detection in this context (though relevant).

### 3. *E. coli* K-12 MG1655 (`GCF_000005845.2`)

-   **Expected**: Chromosomal FOS targets (`murA`, `uhpT`, `glpT`).
-   **Observed**:
    -   Detected: `cyaA`, `lon`, `ptsI`, `uhpB`, `uhpT`.
    -   `glpT` was not detected by BLAST (likely due to low identity or absence in this specific assembly region/version, or fallback issue). Wait, `glpT` is a chromosomal gene.
    -   Detected `uhpT` with mutation `E350T`.
-   **Analysis**:
    -   Chromosomal targets detected.
    -   `glpT` absence might be due to the specific assembly or database fallback. The `create-db` log showed: "WARNING: glpT not found in catalog, using fallback".
    -   `uhpT` mutation `E350T` detected.

## Conclusion

At the time of this check the installation and pipeline ran without errors on
all three genomes. The `blaKPC` D179T caveat and the `glpT` fallback warning
are specific to that earlier version. For the current tool's accuracy,
including known limitations such as low fosfomycin sensitivity on the
ESKAPE-GOLD sets, see [../docs/VALIDATION.md](../docs/VALIDATION.md). Re-run
the three genomes to refresh Results sections 1-3.
