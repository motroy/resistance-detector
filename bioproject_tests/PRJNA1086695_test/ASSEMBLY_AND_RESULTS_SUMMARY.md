# PRJNA1086695 — myloasm Assembly + FOS-CAZAVI Resistance Detection

Two isolates from PRJNA1086695 (see [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11590670/)) were assembled with
[myloasm](https://myloasm-docs.github.io/) and run through the `fos-cazavi` resistance detector.

## 1. Assembly (myloasm)

myloasm v0.5.1 (built from source, `bluenote-1577/myloasm`) was run on each `.fastq.gz` file:

```
myloasm SRR28296939.fastq.gz -o SRR28296939_myloasm -t 4 \
    --singleton-coverage-threshold 0 --secondary-coverage-threshold 0
myloasm SRR28296940.fastq.gz -o SRR28296940_myloasm -t 4 \
    --singleton-coverage-threshold 0 --secondary-coverage-threshold 0
```

**Important caveat about the input data and a myloasm bug encountered:**

- The "reads" in these `.fastq.gz` files are not raw long reads — each record is a full pre-assembled
  contig (SPAdes-style `NODE_x_length_y_cov_z` headers, lengths up to ~365 kb) packaged into FASTQ
  records with uniform quality strings. There is essentially no redundant read coverage of any region
  (almost no overlaps between records), which is exactly the opposite of what myloasm (an
  overlap/consensus long-read assembler) expects.
- With default settings myloasm filtered out nearly all "reads" as low-coverage singletons and then
  **crashed** (`capacity overflow` panic in `skani::triangle::triangle_return`) while trying to
  dereplicate an empty contig set.
- Setting `--singleton-coverage-threshold 0 --secondary-coverage-threshold 0` let myloasm retain the
  unitigs, but the final skani-based dereplication step (`assembly_primary.fa`) is still broken on this
  input: for SRR28296939 it collapsed 119 contigs (5.88 Mb) down to a single 14.7 kb contig, and for
  SRR28296940 it crashed with the same panic after polishing completed.
- This looks like a genuine myloasm/skani edge-case bug (not something tunable via CLI flags — disabling
  dereplication thresholds did not change the outcome) triggered by assembling a near-complete set of
  non-overlapping, already-assembled contigs rather than real noisy long reads.
- **Workaround used:** the pre-dereplication, pre-polish unitig set
  (`<run>/3-mapping/final_contigs_nopolish.fa`) was used as the final assembly for downstream analysis,
  since it is the last stage at which all contig content survives intact. Full myloasm run directories
  (logs, graphs, intermediate stages) are kept in this folder for inspection.

| Isolate | Assembly used | Contigs | Total bases | N50-ish (largest) |
|---|---|---|---|---|
| SRR28296939 | `SRR28296939_myloasm_assembly.fasta.gz` | 119 | 5,880,009 bp | 365,202 bp |
| SRR28296940 | `SRR28296940_myloasm_assembly.fasta.gz` | 122 | 5,886,378 bp | 365,202 bp |

Both assemblies are ~5.9 Mb, consistent with a single *Klebsiella*-sized bacterial genome.

## 2. Resistance Detection (fos-cazavi)

Run with the bundled reference database and the current pipeline:

```bash
fos-cazavi fos-cazavi-all \
    -a <sample>_myloasm_assembly.fasta \
    -o <sample> \
    -d fos_cazavi/data/example_database.fasta \
    --organism Klebsiella_pneumoniae
```

Reference data: AMRFinderPlus 2026-08-07.1.

### Results

| Sample | Beta-lactamases | blaKPC changes (Ambler) | Predicted FOS | Predicted CAZ/AVI |
|---|---|---|---|---|
| SRR28296939 | **blaKPC-179**, (intrinsic fosAKP) | A133T; insS@180 | Susceptible | **Resistant** |
| SRR28296940 | blaKPC-2, (intrinsic fosAKP) | none | Susceptible | Susceptible |

SRR28296939 carries a KPC whose change set matches NCBI allele **KPC-179**
exactly: a single-residue insertion immediately after the Omega-loop residue
D179, plus A133T. NCBI curates KPC-179 as an *inhibitor-resistant
extended-spectrum class A beta-lactamase*, so the isolate is predicted
ceftazidime-avibactam-resistant.

The earlier committed run of these same assemblies reported a plain wild-type
blaKPC and a susceptible call, because insertions were not represented in the
variant caller and allele assignment did not exist. See
[../../docs/METHODS.md](../../docs/METHODS.md).

Neither isolate carries an acquired fosA-family enzyme or a loss-of-function
change in the fosfomycin uptake genes, so both are predicted
fosfomycin-susceptible.
