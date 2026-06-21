# PRJNA781811 — FOS-CAZAVI Resistance Detection Results

18 *Klebsiella pneumoniae* / *K. variicola* assemblies from BioProject PRJNA781811 (listed in
`PRJNA781811.ncbi_datasets.tsv`; see
[paper](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2026.1694693/full))
were downloaded with the NCBI `datasets` CLI and run through `fos-cazavi fos-cazavi-all` using the
bundled reference database (`fos_cazavi/data/example_database.fasta` /
`example_database_deduplicated.fasta` / `example_database_mutations.tsv`) and default primers/genes.

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
