# Testing

```bash
pip install "fos-cazavi[dev]"
python3 -m pytest tests/ -q
```

Tests that need BLAST+ skip themselves automatically when it is not installed.

## What is tested

### `tests/test_variants.py` — the variant caller, directly

* Ambler numbering: positions 58 and 253 are absent, anchors are correct, and
  the mapping round-trips.
* Substitutions, in-frame deletions (including adjacent ones merged into a
  range) and insertions are each reported as what they are.
* A deletion does not shift the calls that follow it — the failure mode that
  fabricates downstream substitutions.
* Premature stops, frameshifts and truncations are detected, and a small
  in-frame deletion is *not* mistaken for a knockout.
* Gene spans are recovered correctly from partial alignments, from the reverse
  strand, and are flagged when they run off a contig.

### `tests/test_pipeline.py` — the real pipeline on synthetic genomes

Genes are taken from the bundled reference data, mutated at a known position,
embedded in a contig, and run through BLAST and the caller. Checks include:

* wild-type KPC-2 is typed as KPC-2 and called susceptible;
* H274Y is typed as KPC-3 and is *not* treated as a resistance marker;
* D179Y is typed as KPC-33 and called resistant;
* an Omega-loop deletion is called resistant on the mechanism;
* an undocumented hotspot change gives Indeterminate;
* a gene on the reverse strand gives the same call;
* NDM-1 gives a resistant CAZ/AVI call even alongside wild-type KPC;
* acquired fosA3 gives fosfomycin resistance, intrinsic fosAKP does not;
* a nonsense mutation in `uhpT` gives fosfomycin resistance, and the same gene
  cut by a contig boundary does not;
* curated chromosomal mutations are not reported without `--organism`;
* the summary files are written and parse correctly.

### `tests/test_scenarios.py` — ground-truth scenarios

`create_test_genomes.py` writes 14 synthetic genomes together with the
phenotype each is built to produce (`expected_results.tsv`). The test
regenerates them, runs the pipeline over each, and asserts the declared
expectation for both drugs. Regenerate them yourself with:

```bash
python3 create_test_genomes.py test_data
```

## Validation on real genomes

See [VALIDATION.md](VALIDATION.md).
