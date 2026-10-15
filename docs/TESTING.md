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
  strand, and are flagged when they run off a contig - with the missing part
  padded so the observed part stays in frame, on either strand.
* A residue that was not observed (`X`) is never reported as a change.

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
* the summary files are written and parse correctly;
* a blaKPC with a premature stop or frameshift, a disrupted NDM and a disrupted
  fosA3 give Indeterminate rather than Resistant, and a disrupted NDM does not
  hide an intact one;
* a contig cut mid-codon keeps the frame (KPC-33 is still typed exactly), a
  gene a few bases short of a contig start is not called a premature stop,
  and a KPC split across two contigs is reported as partial rather than
  missed.

### `tests/test_inputs.py` — inputs, installation and failure behaviour

* Empty, non-FASTA, FASTQ, protein, all-`N` and duplicate-contig inputs are
  rejected with a clear message and exit code 2, writing no results.
* Gzip, Windows line endings, lower case and IUPAC codes are accepted and
  normalised.
* The BLAST index: an existing fresh index is used; otherwise one is built in
  the cache without writing next to the database; an edited database or an
  index older than its FASTA is never used.
* File names containing spaces and shell metacharacters (`;`) work, and are
  not executed by GAMMA's shell command.
* One bad assembly does not stop a batch; clashing sample names are refused;
  gzipped assemblies are found.

### `tests/test_batch.py` — batch mode and the combined tables

Serial and parallel runs give byte-identical tables; the combined columns are
complete; OXA-48-like enzymes are listed as carbapenemases.

### `tests/test_scenarios.py` — ground-truth scenarios

`create_test_genomes.py` writes 14 synthetic genomes together with the
phenotype each is built to produce (`expected_results.tsv`). The test
regenerates them, runs the pipeline over each, and asserts the declared
expectation for both drugs. Regenerate them yourself with:

```bash
python3 create_test_genomes.py test_data
```

## What CI checks

Every push and pull request (`.github/workflows/ci.yml`):

* unit tests on Python 3.10-3.12 without the bioinformatics tools;
* the full suite in the conda environment with BLAST+, GAMMA and seqkit;
* the **installed wheel**: built, installed non-editable, made read-only,
  and run end-to-end from an unrelated directory on a file name containing a
  space, with an empty input checked to fail;
* `flake8`, failing on any pyflakes finding.

## Validation on real genomes

See [VALIDATION.md](VALIDATION.md). `scripts/run_validation.sh` re-downloads
all 112 validation assemblies and re-runs them; the committed results under
`bioproject_tests/` are what it produces, so `git diff` after a run shows
exactly what a code change did to real-genome calls.
