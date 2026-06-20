# Comparison: fos-cazavi results vs. Hernández-García et al. 2022 (JCM 60:e02245-21)

"Impact of Ceftazidime-Avibactam Treatment in the Emergence of Novel KPC Variants
in the ST307-*Klebsiella pneumoniae* High-Risk Clone..."

## 1. Setup

- Tool installed from this repo (`pip install -e .`) with system dependencies:
  BLAST+ 2.12.0 (apt), seqkit v2.8.2, GAMMA.py (rastanton/GAMMA) + UCSC `blat` v39x1.
- Reference DB: bundled `fos_cazavi/data/example_database.fasta` (AMRfinderPlus-derived
  FOS/CAZAVI gene set), built into a BLAST nucleotide DB.
- Input: `test_data/PRJNA741867.ncbi_dataset.zip` → 6 *K. pneumoniae* ST307 assemblies
  (BioProject PRJNA741867, NCBI GenBank).
- Command per genome: `fos-cazavi fos-cazavi-all -a <assembly> -d resistance_db -o <prefix>`

## 2. Mapping assemblies to the paper's isolates

The zip contains no isolate names matching the paper directly, but BioSample
`collection_date`, `isolation_source`, and `host_disease` attributes match the
paper's Fig. 1 timeline and Table 2 exactly:

| Accession        | Strain | Collected   | Source            | → Paper isolate | Paper-reported KPC |
|------------------|--------|-------------|--------------------|------------------|---------------------|
| GCA_022423605.1  | S16    | 2020-03-02  | Broncho aspirate   | **A-1**          | KPC-3 (CAZ/AVI-S)   |
| GCA_022423565.1  | 1-A3   | 2020-03-25  | Rectal             | **A-2**          | KPC-46 (L168P)      |
| GCA_022423575.1  | 1-G9   | 2020-06-02  | Rectal             | **B-1**          | KPC-3 (CAZ/AVI-S)   |
| GCA_022423525.1  | 2-C6   | 2020-06-16  | Rectal             | **B-2**          | KPC-66 (E167_L168del) |
| GCA_022423535.1  | 1-I6   | 2020-06-09  | Rectal             | **C-1**          | KPC-3 (CAZ/AVI-S)   |
| GCA_022423505.1  | HC     | 2020-07-13  | Blood              | **C-2**          | KPC-92 (E167_N169delinsD) |

(dates/sources for each isolate cross-checked against Table 2 / Fig. 1 of the paper —
every accession matches exactly one isolate, with no ambiguity.)

## 3. Gene-level detection (BLAST) — all isolates

All 6 assemblies gave 8 gene hits passing thresholds (≥90% id / ≥80% cov), identical
gene content across isolates: `blaKPC-3`, `blaSHV-12`, `fosAKP`, `acrB`, `ftsI`,
`ompK35`, `ompK36`, `envZ`. This matches the paper's statement that "all KPC
K. pneumoniae isolates ... carried identical antimicrobial resistance genes" —
the GAMMA alignments for blaSHV-12 (Y3F,Q31L,S234G,K235E) and fosAKP (Q130P) are
**identical across all 6 genomes**, confirming clonal identity (ST307) independent
of the WGS-based MLST/SNP typing reported in the paper.

blaKPC-3 itself is detected at ~99.3–100% nucleotide identity / 100% coverage in
every isolate, correctly identifying the presence of a KPC-3-family carbapenemase
gene in all 6 genomes — consistent with the paper (all isolates are blaKPC-3-derived).

The tool's bundled database does not include blaCTX-M-15, blaTEM-1, blaOXA-1/-9, or
the SHV-28 allele reported in the paper, because `fos-cazavi` is scoped only to
FOS/CAZAVI-relevant genes (it is not a general AMR-resistome tool like ABRicate,
which the paper used). This is expected and not a discrepancy.

## 4. Mutation detection — the key comparison

### 4a. GAMMA raw protein-level alignment (ground truth within this pipeline)

| Isolate (mapped) | GAMMA Match_Type | GAMMA Codon_Changes (raw) | Paper-reported variant | Concordance |
|---|---|---|---|---|
| A-1 (S16)      | Native | (none)                              | KPC-3  (wild-type)         | ✅ exact match |
| A-2 (1-A3)     | Mutant | **L168P**                           | KPC-46 (L168P)              | ✅ **exact match** |
| B-1 (1-G9)     | Native | (none)                              | KPC-3  (wild-type)         | ✅ exact match |
| B-2 (2-C6)     | Indel  | 6 bp deletion at nt 496 + L166W     | KPC-66 (E167_L168del, 6 bp/2-codon deletion) | ✅ deletion size/location matches; GAMMA's automatic indel re-annotation (L166W) differs from the paper's manual HGVS notation, but both pinpoint the same 6-bp in-frame deletion in the Ω(X)-loop |
| C-1 (1-I6)     | Native | (none)                              | KPC-3  (wild-type)         | ✅ exact match |
| C-2 (HC)       | Indel  | 6 bp deletion at nt 501 + E167D,N169D | KPC-92 (E167_N169delinsD)  | ✅ deletion-insertion size/location matches; GAMMA reports it as a 6-bp deletion plus 2 flanking substitutions, the paper as a delins — same underlying 9-bp event in the same codons |

All six calls are in **complete qualitative agreement** with the paper: the three
pre-CAZ/AVI isolates (A-1/B-1/C-1) are correctly called wild-type blaKPC-3, and the
three post-CAZ/AVI isolates (A-2/B-2/C-2) all carry in-frame substitutions/deletions
clustered at codons 166–169 — precisely the X-loop/Ω-loop region (Arg164–Asp179)
the paper identifies as the mutational hotspot for blaKPC-46/-66/-92. The A-2 call
(L168P) is an **exact, unambiguous match** to KPC-46.

### 4b. What the pipeline's high-level summary/report actually shows

Despite GAMMA's raw output being correct, **none of these mutations is surfaced**
in `*_summary.txt`'s "MUTATION DETECTION SUMMARY (Dual-Method)" section ("No
resistance mutations detected by GAMMA or seqkit" for every isolate) or in
`*_protein_mutations.tsv`. Root cause:

1. `fos_cazavi/data/example_database_mutations.tsv` only defines blaKPC positions
   **179 (D179Y/N), 240 (V240G), 243 (T243M)**. It has no entries for codons
   166–169, so `_filter_known_mutations()` in `mutations.py` discards the real
   L168P / E167D / N169D substitutions as "not a known resistance position."
2. `_parse_gamma_codon_changes()` uses a regex (`^[A-Z*]\d+[A-Z*]$`) that only
   accepts simple single-codon substitution strings and explicitly drops any
   GAMMA description containing indel/deletion text — so the B-2 and C-2 raw
   calls ("6 bp Deletion at ...") are filtered out entirely before they even
   reach the known-mutation lookup.

**Net effect: as configured out of the box, `fos-cazavi` would not flag any of
the three CAZ/AVI-resistant isolates in this dataset as carrying a clinically
significant KPC mutation**, even though the underlying GAMMA alignment (which the
tool already runs) contains the correct information. This mirrors the paper's own
central finding — that novel ESBL-like blaKPC variants are easily "missed" — but
here the miss happens at the bioinformatic curation/parsing layer rather than the
phenotypic-test layer the paper studied.

### 4c. A second, unrelated issue: spurious BLAST-based mutation calls

The `*_results.tsv` / summary's BLAST-derived "Mutations:" column for blaKPC-3
reports *different, inconsistent* amino-acid changes at positions 179/240/243 for
every isolate (e.g. `D179T,V240Y,T243A` for A-1; `D179H,V240C,T243A` for C-1),
**including the three isolates GAMMA correctly calls wild-type with zero coding
changes**. This is a translation-frame artifact in `acquired.py`'s
`extract_hit_sequence()` / `detect_mutations()` (naive `Seq(...).translate()` of
the raw BLAST-aligned subsequence, without GAMMA's proper codon-frame handling) —
not a real biological signal. The GAMMA-based calls (section 4a) should be treated
as authoritative; the BLAST-based mutation column in the summary is unreliable for
blaKPC in this run and should not be used to compare against the paper.

## 5. Bottom line

| Aspect | Result |
|---|---|
| Correct isolate→gene assignment (blaKPC-3 present in all 6) | ✅ matches paper |
| Correct wild-type vs. mutant calls (GAMMA raw) | ✅ 6/6 matches paper |
| Correct mutation identity/location (GAMMA raw) | ✅ L168P exact match (KPC-46); 6-bp X-loop deletions co-localize with KPC-66/KPC-92 |
| Mutations surfaced in tool's summary/report | ❌ 0/3 — known-mutation DB and indel-discarding parser miss all three novel variants |
| BLAST-based "Mutations" column for blaKPC-3 | ❌ unreliable/spurious for this gene in this run |
| Clonal consistency (identical SHV/fosA/acrB/etc. alleles across isolates) | ✅ matches paper's WGS-based clonality finding |

**Conclusion:** the underlying gene/variant detection machinery in `fos-cazavi`
(BLAST for gene presence, GAMMA for protein-level alignment) successfully
reproduces the paper's core genomic findings when the raw GAMMA output is
inspected directly. However, the bundled `example_database_mutations.tsv`
curation (covering only D179Y/N, V240G, T243M) does not yet include the X-loop
positions 166–169 needed to automatically flag the novel blaKPC-46/-66/-92
variants reported in this study, so the tool's headline summary currently
under-reports CAZ/AVI resistance risk for isolates carrying these specific
variants. Expanding the mutation database to include indel-tolerant entries for
codons 164–179 (the full Ω/X-loop) would let `fos-cazavi` correctly flag this
exact scenario end-to-end.
