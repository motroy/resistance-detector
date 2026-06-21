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
| SRR28296939 | `SRR28296939_myloasm_assembly.fasta` | 119 | 5,880,009 bp | 365,202 bp |
| SRR28296940 | `SRR28296940_myloasm_assembly.fasta` | 122 | 5,886,378 bp | 365,202 bp |

Both assemblies are ~5.9 Mb, consistent with a single *Klebsiella*-sized bacterial genome.

## 2. Resistance Detection (fos-cazavi)

Run with the bundled reference database (`fos_cazavi/data/example_database.fasta`, BLAST-indexed) and
bundled primers/genes:

```
fos-cazavi fos-cazavi-all -a SRR28296939_myloasm_assembly.fasta -d resistance_db -o SRR28296939
fos-cazavi fos-cazavi-all -a SRR28296940_myloasm_assembly.fasta -d resistance_db -o SRR28296940
```

### Key findings

Both isolates carry **blaKPC-2** (carbapenemase, CAZAVI-relevant) plus the same set of chromosomal
porin/efflux/PBP genes (ompK35/ompK36/acrB/envZ/ftsI) and **fosAKP** (chromosomal fosfomycin gene, I91V).

| Gene | SRR28296939 | SRR28296940 |
|---|---|---|
| blaKPC-2 | 99.55% id / 100.34% cov — **V240G, T242G** | 100.00% id / 100.00% cov — wild-type |
| fosAKP | 98.81% id / 100.00% cov — I91V | 98.81% id / 100.00% cov — I91V |
| ompK36 | 92.22% id / 102.45% cov — wild-type | 92.22% id / 102.45% cov — wild-type |
| ompK35 | 99.35% id / 100.00% cov — G134A, D135V, D181V | same |
| acrB | 99.84% id / 100.00% cov — G617A, F626A, A628T/V | same |
| envZ | 99.93% id / 100.00% cov — G244S/D, T324Q | same |
| ftsI | 99.32% id / 100.00% cov — A333P, Y350L, S357Q | same |

**Notable difference between isolates:** SRR28296939's blaKPC-2 carries the V240G/T242G mutations
(confirmed independently by both GAMMA and the SeqKit amplicon method, see
`*_unified_mutations.tsv`), which are associated with reduced ceftazidime-avibactam susceptibility,
while SRR28296940's blaKPC-2 is wild-type at those positions. This is consistent with the two isolates
being related but phenotypically/genotypically distinct KPC-producing strains, matching the premise of
PRJNA1086695 (paired isolates from the same study/outbreak).

## Files in this folder

- `SRR28296939.fastq.gz`, `SRR28296940.fastq.gz` — original input "read" files
- `SRR28296939_myloasm/`, `SRR28296940_myloasm/` — full myloasm run directories (logs, graphs, stages)
- `SRR28296939_myloasm_assembly.fasta`, `SRR28296940_myloasm_assembly.fasta` — final assemblies used for
  resistance detection (pre-dereplication unitigs, see caveat above)
- `SRR28296939_*`, `SRR28296940_*` — fos-cazavi outputs (`_summary.txt`, `_results.tsv`, `_all_results.tsv`,
  `_unified_mutations.tsv`, `_genes.fasta`, `_blast.txt`, `_gamma.gamma`/`.psl`, `_amplicons.tsv`,
  `_seqkit_primers.tsv`, `_analysis.log`)
