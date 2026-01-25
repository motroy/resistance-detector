#!/bin/bash
# Batch analysis script for multiple assemblies

set -e

# Default parameters
DATABASE=""
INPUT_DIR=""
OUTPUT_DIR="resistance_results"
MIN_ID=90
MIN_COV=80
THREADS=1

# Help message
show_help() {
    cat << EOF
Usage: $(basename $0) -d DATABASE -i INPUT_DIR [OPTIONS]

Batch analysis of multiple assemblies for FOS-CAZAVI resistance

Required arguments:
  -d    Path to resistance gene database (FASTA)
  -i    Input directory containing assembly files (*.fasta, *.fa, *.fna)

Optional arguments:
  -o    Output directory (default: resistance_results)
  --min_id    Minimum identity threshold (default: 90)
  --min_cov   Minimum coverage threshold (default: 80)
  -t    Number of parallel jobs (default: 1)
  -h    Show this help message

Examples:
  $(basename $0) -d resistance_db.fasta -i assemblies/ -o results
  $(basename $0) -d resistance_db.fasta -i assemblies/ -t 4 --min_id 95
