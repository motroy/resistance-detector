#!/bin/bash
# Run every validation set through `fos-cazavi batch` and combine the results.
#
# Assemblies are fetched from NCBI on first run and cached, so a clean checkout
# reproduces the whole thing with one command. Each set is a directory of
# assemblies analysed with the organism it belongs to; results land in
# bioproject_tests/<set>/ next to the per-set write-ups.
#
# Usage:
#   scripts/run_validation.sh [-j JOBS] [-g GENOME_CACHE] [-s SET]
#
#   -j  assemblies analysed in parallel (default: cores - 1)
#   -g  where to cache downloaded assemblies (default: .validation_genomes)
#   -s  run only this set (default: all)

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GENOME_CACHE="${REPO}/.validation_genomes"
JOBS=""
ONLY_SET=""

while getopts "j:g:s:h" opt; do
    case "$opt" in
        j) JOBS="$OPTARG" ;;
        g) GENOME_CACHE="$OPTARG" ;;
        s) ONLY_SET="$OPTARG" ;;
        h) sed -n '2,15p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) exit 1 ;;
    esac
done

# set  |  organism  |  output directory
SETS=(
    "PRJNA741867|Klebsiella_pneumoniae|PRJNA741867_test_results"
    "PRJNA595047|Klebsiella_pneumoniae|PRJNA595047_test"
    "PRJNA781811|Klebsiella_pneumoniae|PRJNA781811_test"
    "PRJNA1086695|Klebsiella_pneumoniae|PRJNA1086695_test"
    "Paeruginosa_ML_subset|Pseudomonas_aeruginosa|Paeruginosa_ML_subset"
    "CREC_fosA3_China|Escherichia|CREC_fosA3_China"
    "ESKAPE_fos_GOLD_Kpneumoniae|Klebsiella_pneumoniae|ESKAPE_fos_GOLD_Kpneumoniae"
    "ESKAPE_fos_GOLD_Paeruginosa|Pseudomonas_aeruginosa|ESKAPE_fos_GOLD_Paeruginosa"
    "ESKAPE_fos_Kpneumoniae_round2|Klebsiella_pneumoniae|ESKAPE_fos_Kpneumoniae_round2"
)

fos_cazavi() {
    if command -v fos-cazavi >/dev/null 2>&1; then
        fos-cazavi "$@"
    else
        (cd "$REPO" && python3 -m fos_cazavi.cli "$@")
    fi
}

fetch_set() {
    # Assembly accessions per set live in bioproject_tests/<dir>/accessions.tsv
    # (two columns: sample name, accession). Sets predating that file list their
    # accessions in their RESULTS_SUMMARY.md and are cached by hand.
    local set_name="$1" out_dir="$2" target="$3"
    local manifest="${REPO}/bioproject_tests/${out_dir}/accessions.tsv"

    # Some sets ship their assemblies in the repository (built from reads, so
    # there is no assembly accession to fetch). Unpack those instead.
    mkdir -p "$target"
    local packaged found=0
    for packaged in "${REPO}/bioproject_tests/${out_dir}"/*_assembly.fasta.gz; do
        [[ -e "$packaged" ]] || continue
        found=1
        local sample
        sample="$(basename "$packaged" _myloasm_assembly.fasta.gz)"
        [[ -f "${target}/${sample}.fna" ]] || {
            echo "  unpacking ${sample} from the repository"
            gunzip -c "$packaged" > "${target}/${sample}.fna"
        }
    done
    [[ "$found" == 1 ]] && return 0

    if [[ ! -f "$manifest" ]]; then
        echo "  no accessions.tsv for ${set_name}; expecting cached assemblies" >&2
        return 0
    fi

    mkdir -p "$target"
    while IFS=$'\t' read -r sample accession; do
        [[ -z "${sample:-}" || "$sample" == \#* ]] && continue
        [[ -f "${target}/${sample}.fna" ]] && continue
        echo "  fetching ${sample} (${accession})"
        if [[ "$accession" == GC[AF]_* ]]; then
            curl -sL -o "${target}/${sample}.zip" \
                "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${accession}/download?include_annotation_type=GENOME_FASTA"
            unzip -qo "${target}/${sample}.zip" -d "${target}/${sample}_tmp"
            find "${target}/${sample}_tmp" -name '*.fna' -exec mv {} "${target}/${sample}.fna" \;
            rm -rf "${target}/${sample}_tmp" "${target}/${sample}.zip"
        else
            # One or more nucleotide accessions, comma separated
            : > "${target}/${sample}.fna"
            for part in ${accession//,/ }; do
                curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${part}&rettype=fasta&retmode=text" \
                    >> "${target}/${sample}.fna"
                sleep 0.5
            done
        fi
    done < <(tail -n +1 "$manifest")
}

echo "fos-cazavi validation run"
echo "  repository:    ${REPO}"
echo "  genome cache:  ${GENOME_CACHE}"

for entry in "${SETS[@]}"; do
    IFS='|' read -r set_name organism out_dir <<< "$entry"
    [[ -n "$ONLY_SET" && "$ONLY_SET" != "$set_name" ]] && continue

    genomes="${GENOME_CACHE}/${set_name}"
    results="${REPO}/bioproject_tests/${out_dir}"

    echo
    echo "=== ${set_name} (${organism}) ==="
    fetch_set "$set_name" "$out_dir" "$genomes"

    if ! compgen -G "${genomes}/*.f*a" >/dev/null; then
        echo "  no assemblies in ${genomes}; skipping" >&2
        continue
    fi

    args=(batch -i "$genomes" -o "$results" --organism "$organism")
    [[ -n "$JOBS" ]] && args+=(-j "$JOBS")
    fos_cazavi "${args[@]}"
done

echo
echo "=== combining every set ==="
mapfile -t result_dirs < <(for entry in "${SETS[@]}"; do
    IFS='|' read -r _ _ out_dir <<< "$entry"
    echo "${REPO}/bioproject_tests/${out_dir}"
done)
fos_cazavi combine -i "${result_dirs[@]}" -o "${REPO}/bioproject_tests/all_validation"

echo
echo "Done. Master tables:"
echo "  bioproject_tests/all_validation_combined_summary.tsv"
echo "  bioproject_tests/all_validation_combined_genes.tsv"
