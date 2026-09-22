# Output files

| File | Description |
|------|-------------|
| `*_results.tsv` | One row per detected gene copy: gene, assigned allele, identity, coverage, completeness, protein changes, loss of function, copy number |
| `*_summary.tsv` | One row per gene, with loci, allele(s), changes and the predicted phenotypes as trailing columns |
| `*_summary.json` | The same content as a structured document, including full per-locus detail and the reference data version |
| `*_summary.txt` | Human-readable summary, phenotypes first |
| `*_all_results.tsv` | Everything (genes, protein changes, GAMMA alignments, amplicons) in one long table |
| `*_genes.fasta` | The extracted gene sequences, with the assigned allele in the header |
| `*_unified_mutations.tsv` | Protein changes with their cross-caller confidence |
| `*_protein_mutations.tsv` | Raw GAMMA results |
| `*_amplicons.tsv` | Amplicon coordinates (`seqkit`) |
| `*_blast.txt`, `*_gamma.gamma`, `*_gamma.psl` | Raw tool output |
| `*_analysis.log` | Command, parameters and tool versions |

## Key columns in `*_results.tsv`

| Column | Meaning |
|---|---|
| `Gene` | The reference gene the locus was called against |
| `Allele` | The allele assigned from the observed protein changes (`blaKPC-31`, `novel blaKPC variant (...)`, or the gene name when no typing applies) |
| `Complete` | `no (contig boundary)` means the gene runs off the end of a contig and could not be fully assessed |
| `Reported_Mutations` | Changes this tool reports as resistance mutations for fosfomycin or ceftazidime-avibactam |
| `Other_Drug_Mutations` | Curated mutations found at this locus but described for a different drug (carbapenems, tigecycline, …), shown as context only |
| `All_Protein_Changes` | Every difference from the reference protein, including neutral ones |
| `Loss_Of_Function` | Premature stop, frameshift or truncation, with the residue numbers |

`Reported_Mutations` and `All_Protein_Changes` are deliberately separate. A
difference from a reference is not the same thing as a resistance mutation, and
the output never blurs the two.

## Predicted phenotypes

Each summary carries a call for both drugs, with the evidence that produced it.

| Call | Meaning |
|---|---|
| `Resistant` | A mechanism with established published evidence is present |
| `Indeterminate` | Something relevant was found but its effect is not established, or a target gene could not be assessed |
| `Susceptible` | The known mechanisms were looked for and not found |

**Fosfomycin** is Resistant for an acquired fosA-family enzyme (the fosA family,
fosB, fosC2, fosL1/L2), or loss of function in `uhpT`, `uhpA`, `uhpB`, `uhpC`,
`glpT`, `cyaA`, `ptsI` or `galU`, or a point mutation in those genes curated
*for fosfomycin* in the declared organism. The intrinsic chromosomal enzymes
(`fosAKP` in *K. pneumoniae*, `fosA_PA1129` in *P. aeruginosa*) are never scored
as resistance, and a lone fosA hit with no intrinsic copy recognised gives
Indeterminate rather than Resistant — see [METHODS.md](METHODS.md#6-fosfomycin).

**Ceftazidime-avibactam** is Resistant for a metallo-beta-lactamase (avibactam
does not inhibit those), or for a blaKPC carrying a documented escape variant or
an in-frame indel in the Omega loop. A blaKPC change in a hotspot that is not
documented gives Indeterminate, as do the contributory chromosomal changes —
curated avibactam-annotated mutations in `envZ`, `ftsI` or `ompK36`, and
`ompK36` loss of function alongside a blaKPC. OXA-48-like enzymes are inhibited
by avibactam and do not on their own produce a Resistant call.

Numbering for class A beta-lactamases is standardised Ambler numbering. See
[METHODS.md](METHODS.md) for the full rules, the evidence behind them, and the
method's limits.

## Confidence in `*_unified_mutations.tsv`

| Confidence | Meaning |
|---|---|
| 100% | The BLAST-based caller and GAMMA both report this change at this residue |
| 50% | Only one of the two reports it |

The `Reported_As_Resistance_Mutation` column says whether the change is one the
phenotype logic acts on, or simply a sequence difference.
