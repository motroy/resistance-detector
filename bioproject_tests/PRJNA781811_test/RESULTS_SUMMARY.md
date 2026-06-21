# PRJNA781811 — FOS-CAZAVI Resistance Detection Results

18 *Klebsiella pneumoniae* / *K. variicola* assemblies from BioProject PRJNA781811 (listed in
`PRJNA781811.ncbi_datasets.tsv`; see
[Arena et al. 2022, Front. Microbiol. 13:983294](https://doi.org/10.3389/fmicb.2022.983294), whose
Data Availability statement deposits its sequenced hypermucoviscous *K. pneumoniae*/*K. variicola*
bacteremia isolates under this BioProject) were downloaded with the NCBI `datasets` CLI and run
through `fos-cazavi fos-cazavi-all` using the bundled reference database
(`fos_cazavi/data/example_database.fasta` / `example_database_deduplicated.fasta` /
`example_database_mutations.tsv`) and default primers/genes.

## Command used

```bash
datasets download genome accession <18 GCA accessions, space-separated> \
    --include genome --filename PRJNA781811.ncbi_dataset.zip

unzip PRJNA781811.ncbi_dataset.zip   # produces ncbi_dataset/data/<accession>/<accession>_*_genomic.fna

fos-cazavi fos-cazavi-all \
    -a ncbi_dataset/data/<accession>/<accession>_*_genomic.fna \
    -o <accession> \
    -d fos_cazavi/data/example_database.fasta \
    --genes fos_cazavi/data/example_database_deduplicated.fasta \
    --mutations fos_cazavi/data/example_database_mutations.tsv
```

(Required external tools — BLAST+ 2.17.0, GAMMA, SeqKit, BLAT, plus GAMMA's `unidecode` Python
dependency — were installed locally via `install_deps.sh` to run the pipeline; note the script's
hardcoded BLAST+ download URL (`ncbi-blast-2.16.0+`) currently 404s because NCBI rotated `LATEST` to
2.17.0 — 2.17.0 was substituted manually. None of these tools/binaries are part of this repo or the
committed results.)

## Genomes analyzed

| Accession | Organism | Strain/Isolate |
|---|---|---|
| GCA_027151845.1 | *K. pneumoniae* | 57HMV004 |
| GCA_027152185.1 | *K. pneumoniae* | 41HMV001-R44 |
| GCA_027151985.1 | *K. pneumoniae* | 57HMV001 |
| GCA_027151795.1 | *K. pneumoniae* | 59HMV |
| GCA_027152215.1 | *K. pneumoniae* | 49HMV003 |
| GCA_027152405.1 | *K. pneumoniae* | 21HMV001 |
| GCA_027151785.1 | *K. pneumoniae* | 57HMV003 |
| GCA_027152065.1 | *K. pneumoniae* | 49HMV002 |
| GCA_027151835.1 | *K. pneumoniae* | 29HMV004 |
| GCA_027152205.1 | *K. pneumoniae* | HMV007LS |
| GCA_027152245.1 | *K. pneumoniae* | 43HMV001 |
| GCA_027152005.1 | *K. pneumoniae* | 46HMV001-A |
| GCA_027152105.1 | *K. pneumoniae* | 46HMV001-C |
| GCA_027151875.1 | *K. pneumoniae* | 29HMV007 |
| GCA_027151995.1 | *K. pneumoniae* | 46HMV001-B |
| GCA_027152225.1 | *K. pneumoniae* | 30HMV001 |
| GCA_027152445.1 | *K. pneumoniae* | 17HMV002 |
| GCA_027152495.1 | *K. variicola* | 03HMV002 |

## Predicted phenotype summary

All 18 isolates were predicted **Susceptible** to ceftazidime-avibactam (no blaKPC Omega-loop/X-loop
resistance marker detected in any genome — wild-type or absent blaKPC). Three isolates were predicted
**Resistant** to fosfomycin via an acquired plasmid-borne `fosA` variant; the remaining 15 carry only the
chromosomal, non-resistance-conferring `fosAKP`:

| Accession | FOS phenotype | Acquired fosA gene | Identity |
|---|---|---|---|
| GCA_027151845.1 | **Resistant** | fosA3 | 100.00% |
| GCA_027152185.1 | **Resistant** | fosA5 | 96.67% |
| GCA_027152215.1 | **Resistant** | fosA5 | 96.19% |
| (other 15 genomes) | Susceptible | — (chromosomal fosAKP only) | — |

Full per-genome predicted phenotype, gene, and mutation calls are in `PRJNA781811_summary_table.tsv`
and the per-accession `*_summary.txt`/`*_summary.json` files.

## Other resistance genes detected

Beyond the FOS-specific markers, the chromosomal background genes (acrB, envZ, ompK35, ompK36, ftsI,
blaSHV-12) were detected with their expected baseline mutation sets in all 18 genomes (consistent with
the chromosomal-background findings in the PRJNA595047/PRJNA741867 datasets in this repo). Five isolates
additionally carry an acquired carbapenemase:

| Accession | Carbapenemase | Identity |
|---|---|---|
| GCA_027151835.1 | blaKPC-3 | 100.00% |
| GCA_027151875.1 | blaKPC-3, blaCMY-178 | 100.00%, 98.87% |
| GCA_027152065.1 | blaKPC-3 | 100.00% |
| GCA_027152225.1 | blaKPC-2 | 100.00% |
| GCA_027152445.1 | blaKPC-3 | 100.00% |

No mutations were detected in any of the blaKPC alleles found here (all `Match_Type=Native` in GAMMA),
so the CAZ/AVI-susceptible phenotype call is consistent: wild-type blaKPC remains avibactam-inhibitable.

## Comparison with the source paper (Arena et al. 2022)

The paper enrolled 19 confirmed hypermucoviscous (HMV, "string test"-positive) bloodstream isolates
from a 2016–17 nationwide Italian surveillance (43 laboratories, 1,502 *K. pneumoniae* bacteremia
episodes screened): 18 *K. pneumoniae* + 1 *K. variicola*. The genomes were typed with Kleborate/
PathogenWatch and phenotyped by broth microdilution (cephalosporins, CAZ/AVI, colistin, fosfomycin,
etc.) plus *Galleria mellonella*/murine virulence models. Only 18 of the 19 assemblies are present in
`PRJNA781811.ncbi_datasets.tsv` (the 18 analyzed here) — one isolate from the paper's cohort does not
appear to have a corresponding public assembly accession in this BioProject.

**Acquired carbapenemases — count matches.** The paper reports 5/19 isolates (across the two major
ST307 and ST512 clones) carrying an acquired blaKPC carbapenemase: 2 isolates with blaKPC-3/blaKPC-2
among the six ST307 strains, and blaKPC-3 in all three ST512 strains. `fos-cazavi` independently
detected blaKPC in exactly 5/18 genomes here (4× blaKPC-3, 1× blaKPC-2, one genome also carrying
blaCMY-178) — consistent in count and allele mix with the paper's own genotyping.

**Ceftazidime-avibactam — a genuine discrepancy.** The paper found all 19 isolates phenotypically
resistant to plain ceftazidime, but only **one** isolate (an ST512 strain, "GMR140") was phenotypically
resistant to ceftazidime-avibactam (94.7% susceptible overall, MIC90 2/4 µg/mL). `fos-cazavi` predicted
**all 18** genomes here Susceptible to CAZ/AVI, because every detected blaKPC allele was GAMMA
`Match_Type=Native` (no Omega-loop/X-loop mutation). This means the tool likely missed the one
genuinely CAZ/AVI-resistant ST512 isolate: CAZ/AVI resistance in blaKPC-producing *K. pneumoniae* can
also arise from mechanisms this pipeline doesn't screen for (porin loss/ompK mutations, blaKPC gene
amplification/copy-number increase, efflux), not just canonical Omega-loop/X-loop substitutions — a
known limitation of mutation-marker-only CAZ/AVI prediction, worth flagging for future improvement.

**Fosfomycin — partial agreement, similar genotype/phenotype gap to the original study.** The paper
phenotypically found 47.4% (9/19) of isolates fosfomycin-resistant, but explicitly identified an
acquired `fosA3` gene in only **one** of them (the ST11 isolate). The other resistant isolates (all
three ST512 strains, plus single ST29/ST35/new-ST37 isolates) had no acquired fosfomycin-resistance
gene detected by the paper's own Kleborate-based pipeline either — i.e., the source study itself
reports a similar genotype-phenotype gap for fosfomycin, most likely from un-screened chromosomal
mutations or efflux-mediated resistance. `fos-cazavi` here detected acquired `fosA` in 3/18 genomes
(1× fosA3 at 100% identity — consistent with the paper's ST11 fosA3 call; 2× fosA5 at ~96% identity,
not narratively reported in the paper's main text but plausibly present in its supplementary Kleborate
output). Both studies therefore likely under-call true fosfomycin resistance when relying only on
acquired-gene detection.

## Files in this folder

- `PRJNA781811.ncbi_dataset.zip` — original NCBI `datasets` download (unzip to get
  `ncbi_dataset/data/<accession>/<accession>_*_genomic.fna` assemblies plus
  `data_summary.tsv`/`assembly_data_report.jsonl` metadata; not kept unzipped here to avoid duplicating
  ~100 MB of genome data already in the zip)
- `<accession>_*` — fos-cazavi outputs per genome: `_summary.txt` (human-readable report),
  `_summary.json`/`_summary.tsv` (machine-readable summary), `_results.tsv`/`_all_results.tsv` (BLAST
  gene detection), `_unified_mutations.tsv` (cross-method mutation confidence, where present),
  `_protein_mutations.tsv` (GAMMA mutation calls), `_genes.fasta` (extracted gene sequences),
  `_blast.txt` (raw BLAST output), `_gamma.gamma`/`_gamma.psl` (raw GAMMA output),
  `_amplicons.tsv`/`_seqkit_primers.tsv` (SeqKit amplicon detection, where amplicons were found),
  `_analysis.log` (tool versions/parameters log)
- `PRJNA781811_summary_table.tsv` — one-row-per-genome rollup of predicted phenotypes and detected
  genes/mutations across all 18 assemblies
