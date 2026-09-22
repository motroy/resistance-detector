#!/bin/bash
# Batch analysis of multiple assemblies for FOS-CAZAVI resistance.
#
# This is a thin wrapper around `fos-cazavi batch`, which does the work:
# it runs the assemblies in parallel worker processes and writes one combined
# table for the whole collection.

set -euo pipefail

DATABASE=""
INPUT_DIR=""
OUTPUT_DIR="resistance_results"
ORGANISM=""
MIN_ID=90
MIN_COV=80
JOBS=""
THREADS=1

show_help() {
    cat << EOF
Usage: $(basename "$0") -i INPUT_DIR [OPTIONS]

Batch analysis of multiple assemblies for FOS-CAZAVI resistance.

Required arguments:
  -i    Input directory containing assemblies (*.fasta, *.fa, *.fna)

Optional arguments:
  -d    Resistance gene database (FASTA) [default: the bundled database]
  -o    Output directory (default: resistance_results)
  -s    Organism: Escherichia, Klebsiella_pneumoniae or Pseudomonas_aeruginosa.
        Required for curated chromosomal point mutations.
  -j    Assemblies to analyse in parallel (default: cores - 1)
  -t    Threads per assembly for blastn/seqkit (default: 1)
        Total CPU use is roughly -j x -t; keep the product <= your core count.
  --min_id   Minimum percent identity (default: 90)
  --min_cov  Minimum percent coverage (default: 80)
  -h    Show this help

Outputs, in the output directory:
  <sample>_*                    per-sample results, as for a single run
  batch_combined_summary.tsv    one row per sample, with both phenotype calls
  batch_combined_genes.tsv      one row per detected gene copy, all samples

Example:
  $(basename "$0") -i assemblies/ -o results/ -s Klebsiella_pneumoniae -j 8
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d) DATABASE="$2"; shift 2 ;;
        -i) INPUT_DIR="$2"; shift 2 ;;
        -o) OUTPUT_DIR="$2"; shift 2 ;;
        -s) ORGANISM="$2"; shift 2 ;;
        -j) JOBS="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        --min_id) MIN_ID="$2"; shift 2 ;;
        --min_cov) MIN_COV="$2"; shift 2 ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Unknown option: $1" >&2; show_help; exit 1 ;;
    esac
done

if [[ -z "$INPUT_DIR" ]]; then
    echo "ERROR: -i INPUT_DIR is required" >&2
    show_help
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: input directory not found: $INPUT_DIR" >&2
    exit 1
fi

args=(batch -i "$INPUT_DIR" -o "$OUTPUT_DIR" -t "$THREADS"
      --min_id "$MIN_ID" --min_cov "$MIN_COV")
[[ -n "$DATABASE" ]] && args+=(-d "$DATABASE")
[[ -n "$ORGANISM" ]] && args+=(--organism "$ORGANISM")
[[ -n "$JOBS" ]] && args+=(-j "$JOBS")

if command -v fos-cazavi >/dev/null 2>&1; then
    exec fos-cazavi "${args[@]}"
else
    exec python3 -m fos_cazavi.cli "${args[@]}"
fi
