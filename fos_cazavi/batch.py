"""Running many assemblies, and combining their results into one table.

Two levels of parallelism are available and they trade off against each other:

* ``--jobs`` runs that many assemblies at once, each in its own process.
* ``--threads`` is handed to the external tools (``blastn``, ``seqkit``) within
  one assembly.

For a collection, ``--jobs`` is almost always the better lever: one assembly
takes a couple of seconds and BLAST scales poorly across threads on a query that
small. Total CPU use is roughly ``jobs x threads``, so keep the product at or
below the number of cores.
"""

import contextlib
import csv
import json
import os
import subprocess
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .references import gene_family, is_mbl

# One row per sample.
COMBINED_COLUMNS = [
    'Sample',
    'Organism',
    'Predicted_Phenotype_Fosfomycin',
    'Predicted_Phenotype_Ceftazidime_Avibactam',
    'Acquired_Genes',
    'Carbapenemases',
    'Fosfomycin_Enzymes',
    'Loss_Of_Function',
    'Reported_Mutations',
    'Genes_Detected',
    'Incomplete_Genes',
    'Fosfomycin_Evidence',
    'Ceftazidime_Avibactam_Evidence',
    'Reference_Data',
]

# One row per detected gene copy, across all samples.
GENE_COLUMNS = [
    'Sample', 'Organism', 'Gene', 'Allele', 'Copy_Number', 'Contig',
    'Start', 'End', 'Identity%', 'Coverage%', 'Complete',
    'Reported_Mutations', 'Loss_Of_Function',
]

ASSEMBLY_SUFFIXES = ('.fasta', '.fa', '.fna', '.fas', '.fsa')


def find_assemblies(inputs):
    """Expand files and directories into a sorted list of assembly paths."""
    assemblies = []
    for item in inputs:
        path = Path(item)
        if path.is_dir():
            assemblies.extend(
                sorted(child for child in path.iterdir()
                       if child.suffix.lower() in ASSEMBLY_SUFFIXES))
        elif path.exists():
            assemblies.append(path)
        else:
            print(f"WARNING: no such file or directory: {item}", file=sys.stderr)
    return assemblies


def prepare_blast_database(database):
    """Build the BLAST index once, before any worker starts.

    Workers would otherwise race to create it and corrupt each other's index.
    """
    index_files = [f"{database}.{extension}" for extension in ('nhr', 'nin', 'nsq')]
    if all(Path(path).exists() for path in index_files):
        return
    print(f"Creating BLAST database from {database}...")
    subprocess.run(['makeblastdb', '-in', str(database), '-dbtype', 'nucl'],
                   check=True, capture_output=True)


def _sample_name(assembly):
    name = Path(assembly).name
    for suffix in ASSEMBLY_SUFFIXES:
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return name


def run_one(task):
    """Analyse a single assembly. Runs in a worker process.

    Unless the run is verbose, everything the pipeline prints is captured into
    ``<prefix>_run.log`` instead of the terminal: fifty samples x twenty lines
    of tool chatter buries the one line per sample that actually matters.
    """
    assembly = task['assembly']
    prefix = task['prefix']
    verbose = task.get('verbose', False)

    log_handle = None
    try:
        if verbose:
            redirect = contextlib.nullcontext()
        else:
            log_handle = open(f"{prefix}_run.log", 'w')
            redirect = contextlib.redirect_stdout(log_handle)

        with redirect:
            return _analyse(task, assembly, prefix)
    except Exception:                                    # noqa: BLE001
        return {'sample': _sample_name(assembly), 'prefix': prefix,
                'error': traceback.format_exc(), 'fos': '-', 'cazavi': '-'}
    finally:
        if log_handle is not None:
            log_handle.close()


def _analyse(task, assembly, prefix):
    """The actual pipeline for one assembly, imported per worker process."""
    import argparse

    from .acquired import run_acquired_detection
    from .cli import write_summary
    from .mutations import run_mutation_detection
    from .phenotype import predict_phenotypes
    from .utils import setup_logger

    setup_logger(prefix, argparse.Namespace(**task),
                 console=task.get('verbose', False))

    blast_results = run_acquired_detection(
        assembly, task['database'], prefix, task['min_id'], task['min_cov'],
        task['mutations'], task['organism'], task['threads'])

    gamma_results, amplicon_results, _, unified_results = run_mutation_detection(
        assembly, prefix, task['genes'], task['primers'],
        blast_results=blast_results, mutation_db_file=task['mutations'],
        organism=task['organism'], threads=task['threads'])

    write_summary(prefix, assembly, blast_results, gamma_results,
                  amplicon_results, None, unified_results,
                  organism=task['organism'])

    phenotypes = predict_phenotypes(blast_results, unified_results, task['organism'])
    return {
        'sample': _sample_name(assembly),
        'prefix': prefix,
        'error': None,
        'fos': phenotypes['fosfomycin']['phenotype'],
        'cazavi': phenotypes['ceftazidime_avibactam']['phenotype'],
    }


def summarise(summary_json):
    """Turn one sample's JSON summary into a combined-table row."""
    with open(summary_json) as handle:
        summary = json.load(handle)

    phenotypes = summary.get('predicted_phenotypes', {})
    fos = phenotypes.get('fosfomycin', {})
    cazavi = phenotypes.get('ceftazidime_avibactam', {})

    acquired, carbapenemases, fos_enzymes = [], [], []
    loss_of_function, reported, incomplete = [], [], []
    gene_copies = 0

    for entry in summary.get('genes', []):
        for locus in entry['loci']:
            gene_copies += 1
            allele = locus.get('allele') or entry['gene']
            if entry['gene'].startswith('bla'):
                acquired.append(allele)
                if is_mbl(entry['gene']) or gene_family(entry['gene']) == 'blaKPC':
                    carbapenemases.append(allele)
            elif entry['gene'].lower().startswith('fos'):
                acquired.append(allele)
                fos_enzymes.append(allele)
            if locus.get('loss_of_function'):
                loss_of_function.append(f"{entry['gene']}: {locus['loss_of_function']}")
            reported.extend(locus.get('reported_mutations') or [])
            if not locus.get('complete', True):
                incomplete.append(entry['gene'])

    def joined(values):
        return ';'.join(sorted(set(values))) if values else '-'

    return {
        'Sample': summary.get('sample', Path(summary_json).name),
        'Organism': summary.get('organism') or '-',
        'Predicted_Phenotype_Fosfomycin': fos.get('phenotype', '-'),
        'Predicted_Phenotype_Ceftazidime_Avibactam': cazavi.get('phenotype', '-'),
        'Acquired_Genes': joined(acquired),
        'Carbapenemases': joined(carbapenemases),
        'Fosfomycin_Enzymes': joined(fos_enzymes),
        'Loss_Of_Function': joined(loss_of_function),
        'Reported_Mutations': joined(reported),
        'Genes_Detected': str(gene_copies),
        'Incomplete_Genes': joined(incomplete),
        'Fosfomycin_Evidence': '; '.join(fos.get('evidence', [])) or '-',
        'Ceftazidime_Avibactam_Evidence': '; '.join(cazavi.get('evidence', [])) or '-',
        'Reference_Data': (summary.get('reference_data') or '-').replace('\n', ' '),
    }


def gene_rows(summary_json):
    """One row per detected gene copy, for the long-format table."""
    with open(summary_json) as handle:
        summary = json.load(handle)

    rows = []
    for entry in summary.get('genes', []):
        for locus in entry['loci']:
            rows.append({
                'Sample': summary.get('sample', ''),
                'Organism': summary.get('organism') or '-',
                'Gene': entry['gene'],
                'Allele': locus.get('allele') or entry['gene'],
                'Copy_Number': str(entry['copy_number']),
                'Contig': locus.get('contig', '-'),
                'Start': str(locus.get('start', '-')),
                'End': str(locus.get('end', '-')),
                'Identity%': str(locus.get('identity', '-')),
                'Coverage%': str(locus.get('coverage', '-')),
                'Complete': 'yes' if locus.get('complete', True) else 'no',
                'Reported_Mutations': ';'.join(locus.get('reported_mutations') or []) or '-',
                'Loss_Of_Function': locus.get('loss_of_function') or '-',
            })
    return rows


def find_summaries(inputs):
    """Every ``*_summary.json`` under the given files or directories."""
    found = []
    for item in inputs:
        path = Path(item)
        if path.is_dir():
            found.extend(sorted(path.glob('*_summary.json')))
        elif path.name.endswith('_summary.json'):
            found.append(path)
    return found


def write_combined(summaries, output_prefix):
    """Write the combined per-sample and per-gene tables."""
    rows, genes = [], []
    for summary_json in summaries:
        try:
            rows.append(summarise(summary_json))
            genes.extend(gene_rows(summary_json))
        except Exception as error:                        # noqa: BLE001
            print(f"WARNING: could not read {summary_json}: {error}", file=sys.stderr)

    rows.sort(key=lambda row: row['Sample'])
    genes.sort(key=lambda row: (row['Sample'], row['Gene']))

    combined = f"{output_prefix}_combined_summary.tsv"
    with open(combined, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=COMBINED_COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} sample rows to {combined}")

    per_gene = f"{output_prefix}_combined_genes.tsv"
    with open(per_gene, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=GENE_COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(genes)
    print(f"Wrote {len(genes)} gene rows to {per_gene}")

    return combined, per_gene


def run_batch(assemblies, output_dir, database, genes, primers, mutations,
              min_id, min_cov, organism, jobs, threads, verbose=False):
    """Analyse every assembly, then write the combined tables."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    prepare_blast_database(database)

    tasks = []
    for assembly in assemblies:
        tasks.append({
            'assembly': str(assembly),
            'prefix': str(output_dir / _sample_name(assembly)),
            'database': str(database),
            'genes': str(genes),
            'primers': str(primers),
            'mutations': mutations,
            'min_id': min_id,
            'min_cov': min_cov,
            'organism': organism,
            'threads': threads,
            'verbose': verbose,
        })

    print(f"Analysing {len(tasks)} assemblies with {jobs} parallel job(s), "
          f"{threads} thread(s) each")

    total = len(tasks)
    results = []

    def report(result):
        done = len(results)
        if result['error']:
            print(f"  [{done:>3}/{total}] {result['sample']:<28} FAILED", flush=True)
        else:
            print(f"  [{done:>3}/{total}] {result['sample']:<28} "
                  f"FOS={result['fos']:<14} CAZ/AVI={result['cazavi']}", flush=True)

    if jobs <= 1:
        for task in tasks:
            results.append(run_one(task))
            report(results[-1])
    else:
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            futures = [pool.submit(run_one, task) for task in tasks]
            for future in as_completed(futures):
                results.append(future.result())
                report(results[-1])

    failures = [r for r in results if r['error']]
    for failure in failures:
        print(f"ERROR: {failure['sample']} failed:\n{failure['error']}",
              file=sys.stderr)

    succeeded = [r for r in results if not r['error']]
    print(f"Completed {len(succeeded)}/{total} assemblies")

    summaries = [Path(f"{r['prefix']}_summary.json") for r in succeeded]
    summaries = [path for path in summaries if path.exists()]
    write_combined(summaries, str(output_dir / 'batch'))

    return 1 if failures else 0


def default_jobs():
    """A sensible default: use the machine, but leave a core free."""
    return max(1, (os.cpu_count() or 2) - 1)
