"""Batch-run Kleborate's species + K. pneumoniae AMR modules directly,
bypassing a CLI aggregation bug in kleborate 3.2.4's `-m` custom module
list handling (the module functions themselves work fine when called
directly; only __main__.py's header-filtering for multi-module `-m` runs
drops everything but the first module's columns).

Run with the venv's own interpreter (see RESULTS_SUMMARY.md for setup):
    /tmp/kleborate_venv/bin/python3 run_kleborate.py
"""
import argparse
import csv
import glob
import sys

from kleborate.modules.enterobacterales__species import enterobacterales__species as species_mod
from kleborate.modules.klebsiella_pneumo_complex__amr import klebsiella_pneumo_complex__amr as amr_mod

parser = argparse.ArgumentParser()
species_mod.add_cli_options(parser)
amr_mod.add_cli_options(parser)
args = parser.parse_args([])

genome_globs = [
    '.validation_genomes/ESKAPE_fos_GOLD_Kpneumoniae/*.fna',
    '.validation_genomes/ESKAPE_fos_Kpneumoniae_round2/*.fna',
    '.validation_genomes/PRJNA781811/*.fna',
]
assemblies = sorted(p for g in genome_globs for p in glob.glob(g))
print(f"{len(assemblies)} assemblies to process", file=sys.stderr)

rows = []
for i, assembly in enumerate(assemblies, 1):
    strain = assembly.split('/')[-1].rsplit('.', 1)[0]
    print(f"[{i}/{len(assemblies)}] {strain}", file=sys.stderr)
    results = {'strain': strain}
    sp = species_mod.get_results(assembly, None, args, results)
    results.update(sp)
    amr = amr_mod.get_results(assembly, None, args, results)
    row = {'strain': strain, 'species': sp.get('species', '')}
    row.update(amr)
    rows.append(row)

fieldnames = ['strain', 'species'] + [k for k in rows[0] if k not in ('strain', 'species')]
with open('/tmp/kleborate_direct_results.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    w.writeheader()
    for row in rows:
        w.writerow({k: (';'.join(v) if isinstance(v, list) else v) for k, v in row.items()})

print(f"Wrote {len(rows)} rows to /tmp/kleborate_direct_results.tsv", file=sys.stderr)
