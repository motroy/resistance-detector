#!/bin/bash
# The script that gets recorded to produce demo/validation-run.cast (and .gif).
#
# It runs the real validation across all six genome sets, logged through
# `dochist` so the session doubles as the provenance record. Re-record with:
#
#   scripts/record_demo.sh --record
#
# Requires: asciinema, agg (https://github.com/asciinema/agg), dochist.

set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ "${1:-}" == "--record" ]]; then
    mkdir -p demo
    asciinema rec --overwrite --idle-time-limit 2 --cols 100 --rows 32 \
        -c "bash scripts/record_demo.sh" demo/validation-run.cast
    agg --speed 2 --font-size 15 --theme asciinema \
        demo/validation-run.cast demo/validation-run.gif
    ls -lh demo/validation-run.cast demo/validation-run.gif
    exit 0
fi

BOLD=$'\033[1m'; GREEN=$'\033[1;32m'; DIM=$'\033[2m'; OFF=$'\033[0m'

say() { printf '%s\n' "${DIM}# $*${OFF}"; sleep 0.6; }
run() { printf '%s$%s %s\n' "$GREEN" "$OFF" "$*"; sleep 0.4; eval "$@"; echo; sleep 0.6; }

printf '%sfos-cazavi — full validation run, logged for provenance%s\n\n' "$BOLD" "$OFF"

say "Start a dochist session: every command and every artifact it produces is recorded."
run "dochist init fos-cazavi-validation -d 'Validation across six published genome sets'"

say "FAIR metadata for the session."
run "dochist meta set license MIT && dochist meta set author motroy"

say "Snapshot the software environment."
run "dochist env snapshot"

say "Run each validation set as its own recorded command: 52 assemblies, 3 species."
for entry in \
    "PRJNA741867|Klebsiella_pneumoniae|PRJNA741867_test_results" \
    "PRJNA595047|Klebsiella_pneumoniae|PRJNA595047_test" \
    "PRJNA781811|Klebsiella_pneumoniae|PRJNA781811_test" \
    "PRJNA1086695|Klebsiella_pneumoniae|PRJNA1086695_test" \
    "Paeruginosa_ML_subset|Pseudomonas_aeruginosa|Paeruginosa_ML_subset" \
    "CREC_fosA3_China|Escherichia|CREC_fosA3_China" ; do
    IFS='|' read -r set_name organism out_dir <<< "$entry"
    run "dochist run -- 'fos-cazavi batch -i .validation_genomes/${set_name} -o bioproject_tests/${out_dir} --organism ${organism} -j 6'"
done

say "Combine every set into one table."
run "dochist run -- 'fos-cazavi combine -i bioproject_tests/*_test bioproject_tests/*_test_results bioproject_tests/Paeruginosa_ML_subset bioproject_tests/CREC_fosA3_China -o bioproject_tests/all_validation'"

say "What did the session record?"
run "dochist log"

say "The headline table: one row per sample, across every set."
run "cut -f1,3,4,6 bioproject_tests/all_validation_combined_summary.tsv | head -14"

say "Generate the FAIR provenance document and a curated rerun script."
run "dochist report --output docs/PROVENANCE.md && dochist extract --output scripts/rerun_validation.sh"

run "ls -lh docs/PROVENANCE.md scripts/rerun_validation.sh"

printf '%sDone.%s\n' "$BOLD" "$OFF"
sleep 1
