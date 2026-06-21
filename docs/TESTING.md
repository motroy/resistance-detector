# Testing

## Running Tests

Install dependencies and run the full test suite:

```bash
pip install "fos-cazavi[dev]"
python -m pytest tests/ -v
```

## Test Results

```
============================= test session starts ==============================
platform linux -- Python 3.11.14, pytest-9.0.2
collected 305 items

305 passed in 8.6s
```

**All 305 tests pass.** No external tools (BLAST, GAMMA, seqkit) are required to run the tests — all tool-dependent logic is tested via synthetic inputs.

## Test Coverage Summary

| Test Module | What Is Covered |
|---|---|
| `test_blast_parsing.py` | BLAST output parsing: identity/coverage filtering, gene name extraction from various ID formats, multi-hit handling, copy-number annotation, edge cases |
| `test_genome_creation.py` | Reference database loading (30 genes verified), `introduce_mutation()` utility, synthetic contig construction for all FOS and CAZAVI genome scenarios |
| `test_gamma_parsing.py` | GAMMA Codon_Changes field parsing for KPC (D179Y, V240G, T243M), OXA-48 (P68A, Y211S), CMY-178 (N70T), porins (OmpK35/36), AcrB, gene name normalization, and fosA variant names (fosA3/4/5/7/11) |
| `test_mutation_detection.py` | `detect_mutations()` for all gene families: wildtype no-call verification plus every documented resistance mutation across FOS and CAZAVI pathways |
| `test_dual_detection.py` | Dual-method (GAMMA + SeqKit amplicon) confidence scoring: 100% when both methods agree, 50% when only one detects a mutation |
| `test_cli_summary.py` | Copy-number aggregation and machine-readable JSON/TSV summary generation |
| `test_phenotype.py` | Genotype-to-phenotype prediction rules for FOS and CAZ/AVI (`fos_cazavi/phenotype.py`) |

## Genes and Mutations Covered by Tests

**Fosfomycin (FOS) — Plasmidic:**

| Gene | Mutations Tested |
|---|---|
| fosA3, fosA4, fosA5, fosA7, fosA11 | K90E, H119Q (+ wildtype) |
| fosAKP | I91V (+ wildtype) |

**Fosfomycin (FOS) — Chromosomal:**

| Gene | Mutations Tested |
|---|---|
| murA | D369N, L370I |
| uhpB | G469R, H350Y, H350Q |
| uhpC | F384L |
| uhpA | D54N, R139C, R139H |
| uhpT | G55D, W198\*, E258\*, W350\* |
| glpT | E44\*, W88\*, G90D, W234\*, R362C, R362H |
| cyaA | G463D, G463\* |
| ptsI | H191Y, H191Q |
| galU | R282V |
| lon | Q558\* |

**Ceftazidime-Avibactam (CAZAVI) — Plasmidic:**

| Gene | Mutations Tested |
|---|---|
| blaKPC-2, -3, -31, -190 | D179Y, V240G, T243M (incl. double mutants) |
| blaOXA-48 | P68A, Y211S (incl. double mutant) |
| blaCMY-178 | N70T |
| blaSHV-12 | G238S, E240K |

**Ceftazidime-Avibactam (CAZAVI) — Chromosomal:**

| Gene | Mutations Tested |
|---|---|
| ompK36 | G134D, G134\*, D135\*, G213D |
| ompK35 | G134D, D135\*, D181G, D181\* |
| acrB | G617D, G617N, F626L, A628T, A628V |
| mexR | W69G, W69\*, A75V |
| nalD | Q153\*, L174R |
| ftsI | A333V, A333T, Y350C, Y350S, S357N |
| envZ | G244S, G244D, T324I, T324A |

See [VALIDATION.md](VALIDATION.md) for end-to-end validation against real, published genomes.
