# Validation Report

Re-run on 2026-10-15 with the current pipeline (AMRFinderPlus reference data
2026-08-07.1, bundled database).

## Executive Summary

`fos-cazavi fos-cazavi-all` ran without errors on all three reference genomes
and the calls match what is known about each strain:

| Genome | Fosfomycin | CAZ/AVI | Key finding |
|---|---|---|---|
| *K. pneumoniae* KPNIH1 | Susceptible | Susceptible | `blaKPC-3` found at 99.89% identity with only the defining `H274Y` change (no Omega-loop change, so no avibactam escape) |
| *K. pneumoniae* MGH 78578 | Susceptible | Susceptible | No `blaKPC`; chromosomal `blaSHV-like` only; no `fosA3` |
| *E. coli* K-12 MG1655 | Susceptible | Susceptible | All 12 chromosomal targets found, including `glpT`; no acquired genes |

The two issues recorded in the February version of this report no longer
occur: there is no spurious `D179T` call on KPNIH1, and `glpT` is detected.

## Methodology

1.  **Environment**: micromamba env with `blast`, `biopython`, `seqkit>=2.8`,
    `gamma`, `pytest` (as in `environment.yaml`); `pip install -e .`.
2.  **Genomes** (NCBI RefSeq/GenBank FTP):
    -   `GCA_000281535.2` (*K. pneumoniae* KPNIH1), the KPC-3 positive
        control. **Note:** `GCA_000281535.1`, used in the earlier report, has
        no sequence file on NCBI; `.2` is the same KPNIH1 strain.
    -   `GCF_000016305.1` (*K. pneumoniae* MGH 78578), negative control.
    -   `GCF_000005845.2` (*E. coli* K-12 MG1655), chromosomal targets.
3.  **Analysis**: `fos-cazavi fos-cazavi-all -a <fasta> -o <prefix>
    --organism <Klebsiella_pneumoniae|Escherichia> -t 4`, using the bundled
    reference data (`create-db` needs NCBI Entrez access and was not used).

## Results

### 1. *K. pneumoniae* KPNIH1 (`GCA_000281535.2`)

-   **Expected**: `blaKPC-3`, no avibactam-escape change.
-   **Observed**:
    -   `blaKPC-3` (allele assigned against `blaKPC-2`): 99.89% identity,
        100% coverage; only protein change `H274Y` (BLAST and GAMMA agree).
        `H274Y` is what makes a KPC a KPC-3 and is not a resistance marker.
    -   `blaSHV-like` (S235G, K236E) also present.
    -   Calls: fosfomycin Susceptible, CAZ/AVI Susceptible.
    -   `ompK35` has a premature stop at residue 89 (loss of function). It is
        reported but does not trigger a CAZ/AVI resistant call on its own.
-   **Caveat**: the many 50%-confidence `ompK35` changes (BLAST only, not
    confirmed by GAMMA) are downstream effects of the frameshift/truncation,
    not independent mutations.

### 2. *K. pneumoniae* MGH 78578 (`GCF_000016305.1`)

-   **Expected**: no KPC, no `fosA3`.
-   **Observed**: no `blaKPC` and no acquired fosfomycin enzyme. A `blaSHV-12`
    exact match and a `blaSHV-like` copy (S235G, K236E) were detected, plus the
    chromosomal targets. `ompK35` has a premature stop at residue 309.
    Both drugs called Susceptible.

### 3. *E. coli* K-12 MG1655 (`GCF_000005845.2`)

-   **Expected**: chromosomal fosfomycin/CAZ-AVI targets, nothing acquired.
-   **Observed**: `cyaA`, `lon`, `ptsI`, `uhpA/B/C/T`, `glpT`, `murA`, `ompC`,
    `ompF`, `galU` all found at 98.8-100% identity and 100% coverage with no
    protein changes reported. No acquired genes. Both drugs Susceptible.
-   **Caveat**: the amplicon table lists two `RT_rrsA` hits spanning hundreds
    of kb (coordinates only; the primer pair matches multiple rRNA operons).
    These are not used for any call.

## Conclusion

The tool installs from `environment.yaml`, runs end-to-end, and gives
biologically correct results on the three reference genomes. These are
reference strains, so this is a sanity check rather than an accuracy
measurement; for accuracy against measured MICs see
[../docs/VALIDATION.md](../docs/VALIDATION.md).
