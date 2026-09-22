import argparse
import json
import sys
from collections import Counter
from pathlib import Path

from .acquired import run_acquired_detection
from .db import create_db
from .mutations import run_mutation_detection
from .phenotype import predict_phenotypes
from .references import (
    DEFAULT_DATABASE, SUPPORTED_ORGANISMS, data_version,
)
from .utils import log_tool_versions, setup_logger

_DATA_DIR = Path(__file__).parent / 'data'
_DEFAULT_GENES = str(DEFAULT_DATABASE)
_DEFAULT_PRIMERS = str(_DATA_DIR / 'primers.tsv')


def write_combined_tsv(output_prefix, assembly, blast_results, gamma_results,
                       amplicon_results, unified_results):
    output_file = f"{output_prefix}_all_results.tsv"
    print(f"Writing combined summary to {output_file}...")

    blast_copy_counts = Counter(r['gene'] for r in blast_results) if blast_results else Counter()
    gamma_copy_counts = Counter(r['protein'] for r in gamma_results) if gamma_results else Counter()

    with open(output_file, 'w') as handle:
        handle.write('\t'.join([
            'Sample', 'Type', 'Gene', 'Allele', 'Method', 'Result_Details',
            'Identity_Confidence', 'Contig', 'Start', 'End', 'Copy_Number', 'Notes',
        ]) + '\n')

        sample = Path(assembly).name

        for result in blast_results or []:
            notes = []
            if result['changes']:
                notes.append(f"protein changes: {';'.join(result['changes'])}")
            if result['lof_description']:
                notes.append(result['lof_description'])
            if not result['complete']:
                notes.append('gene incomplete (contig boundary)')
            handle.write('\t'.join([
                sample,
                'Acquired gene' if result['acquired'] else 'Chromosomal gene',
                result['gene'],
                result['allele'],
                'BLAST',
                f"Coverage: {result['coverage']}%",
                f"{result['identity']}%",
                result['contig'],
                str(result['start']),
                str(result['end']),
                str(blast_copy_counts[result['gene']]),
                '; '.join(notes) if notes else 'no protein changes vs reference',
            ]) + '\n')

        for item in unified_results or []:
            handle.write('\t'.join([
                sample,
                'Protein change',
                item['gene'],
                '-',
                '+'.join(item['methods']),
                item['mutation'],
                f"{item['confidence']}%",
                item['contig'],
                '-', '-', '-',
                'reported as resistance mutation' if item['reported']
                else 'sequence difference, not a curated resistance mutation',
            ]) + '\n')

        for result in gamma_results or []:
            if result['mutations']:
                continue
            handle.write('\t'.join([
                sample, 'Gene alignment', result['protein'], '-', 'GAMMA',
                'No codon changes detected', f"{result['identity']:.2f}%",
                result['contig'], str(result['contig_start']), str(result['contig_end']),
                str(gamma_copy_counts[result['protein']]),
                f"Coverage: {result['coverage']:.2f}%",
            ]) + '\n')

        for amplicon in amplicon_results or []:
            handle.write('\t'.join([
                sample, 'Amplicon', amplicon['pair_id'], '-', 'SeqKit/amplicon',
                f"Length: {amplicon['length']}bp", '-', amplicon['contig'],
                str(amplicon['start']), str(amplicon['end']), '-',
                'Genes in region: ' + ('; '.join(amplicon['mutations_found'])
                                       if amplicon['mutations_found'] else 'none'),
            ]) + '\n')


def _aggregate_genes(blast_results):
    by_gene = {}
    for result in blast_results or []:
        entry = by_gene.setdefault(result['gene'], {
            'gene': result['gene'],
            'copy_number': 0,
            'loci': [],
        })
        entry['copy_number'] += 1
        entry['loci'].append({
            'contig': result['contig'],
            'start': result['start'],
            'end': result['end'],
            'allele': result['allele'],
            'identity': result['identity'],
            'coverage': result['coverage'],
            'complete': result['complete'],
            'changes': result['changes'],
            'reported_mutations': result['reported_mutations'],
            'loss_of_function': result['lof_description'],
        })
    return list(by_gene.values())


def write_machine_summary(output_prefix, assembly, blast_results, gamma_results,
                          amplicon_results, unified_results, organism=None):
    sample = Path(assembly).name
    gene_entries = _aggregate_genes(blast_results)
    phenotypes = predict_phenotypes(blast_results, unified_results, organism)

    summary = {
        'sample': sample,
        'assembly': assembly,
        'organism': organism,
        'reference_data': data_version(),
        'total_genes_detected': len(blast_results) if blast_results else 0,
        'unique_genes': len(gene_entries),
        'genes': gene_entries,
        'protein_changes': unified_results or [],
        'gene_alignments': gamma_results or [],
        'amplicons': amplicon_results or [],
        'predicted_phenotypes': phenotypes,
    }

    json_file = f"{output_prefix}_summary.json"
    print(f"Writing machine-readable JSON summary to {json_file}...")
    with open(json_file, 'w') as handle:
        json.dump(summary, handle, indent=2, default=str)

    tsv_file = f"{output_prefix}_summary.tsv"
    print(f"Writing machine-readable TSV summary to {tsv_file}...")
    fos = phenotypes['fosfomycin']
    cazavi = phenotypes['ceftazidime_avibactam']
    phenotype_columns = [
        fos['phenotype'], '; '.join(fos['evidence']),
        cazavi['phenotype'], '; '.join(cazavi['evidence']),
        phenotypes['disclaimer'],
    ]

    with open(tsv_file, 'w') as handle:
        handle.write('\t'.join([
            'Sample', 'Organism', 'Gene', 'Allele', 'Copy_Number', 'Loci',
            'Max_Identity%', 'Max_Coverage%', 'Protein_Changes',
            'Reported_Mutations', 'Loss_Of_Function',
            'Predicted_Phenotype_Fosfomycin', 'Fosfomycin_Evidence',
            'Predicted_Phenotype_Ceftazidime_Avibactam', 'Ceftazidime_Avibactam_Evidence',
            'Phenotype_Disclaimer',
        ]) + '\n')

        if not gene_entries:
            handle.write('\t'.join([sample, organism or '-', '-', '-', '0',
                                    '-', '-', '-', '-', '-', '-'] + phenotype_columns) + '\n')

        for entry in gene_entries:
            loci = ';'.join(f"{l['contig']}:{l['start']}-{l['end']}" for l in entry['loci'])
            alleles = ';'.join(sorted({l['allele'] for l in entry['loci']}))
            changes = sorted({change for l in entry['loci'] for change in l['changes']})
            reported = sorted({m for l in entry['loci'] for m in l['reported_mutations']})
            lof = ';'.join(sorted({l['loss_of_function'] for l in entry['loci']
                                   if l['loss_of_function']}))
            handle.write('\t'.join([
                sample,
                organism or '-',
                entry['gene'],
                alleles,
                str(entry['copy_number']),
                loci,
                f"{max(float(l['identity']) for l in entry['loci']):.2f}",
                f"{max(float(l['coverage']) for l in entry['loci']):.2f}",
                ';'.join(changes) if changes else '-',
                ';'.join(reported) if reported else '-',
                lof or '-',
            ] + phenotype_columns) + '\n')


def write_summary(output_prefix, assembly, blast_results, gamma_results,
                  amplicon_results, seqkit_mut_results=None, unified_results=None,
                  organism=None):
    write_combined_tsv(output_prefix, assembly, blast_results, gamma_results,
                       amplicon_results, unified_results)
    write_machine_summary(output_prefix, assembly, blast_results, gamma_results,
                          amplicon_results, unified_results, organism)

    summary_file = f"{output_prefix}_summary.txt"
    print(f"Writing summary to {summary_file}...")

    with open(summary_file, 'w') as handle:
        handle.write('=' * 70 + '\n')
        handle.write('FOS-CAZAVI Resistance Detection Summary\n')
        handle.write('=' * 70 + '\n\n')
        handle.write(f"Assembly: {assembly}\n")
        organism_note = organism or ('not specified - chromosomal point mutations '
                                     'were not evaluated')
        handle.write(f"Organism: {organism_note}\n")
        handle.write(f"Reference data: {data_version()}\n\n")

        phenotypes = predict_phenotypes(blast_results, unified_results, organism)
        handle.write('PREDICTED PHENOTYPES (genotype-based):\n')
        handle.write('-' * 50 + '\n')
        for drug, key in (('Fosfomycin (FOS)', 'fosfomycin'),
                          ('Ceftazidime-Avibactam (CAZ/AVI)', 'ceftazidime_avibactam')):
            prediction = phenotypes[key]
            handle.write(f"  {drug}: {prediction['phenotype']}\n")
            for evidence in prediction['evidence']:
                handle.write(f"    - {evidence}\n")
        handle.write(f"  Note: {phenotypes['disclaimer']}\n\n")

        if blast_results is not None:
            acquired = [r for r in blast_results if r['acquired']]
            chromosomal = [r for r in blast_results if not r['acquired']]

            handle.write(f"Total genes detected: {len(blast_results)}\n")
            handle.write('Method: BLAST+ with codon-aware variant calling\n\n')

            if acquired:
                handle.write('ACQUIRED RESISTANCE GENES:\n')
                handle.write('-' * 50 + '\n')
                for result in acquired:
                    handle.write(f"  {result['allele']} (gene {result['gene']}, "
                                 f"copy number {result['copy_number']}): "
                                 f"{result['identity']}% identity, "
                                 f"{result['coverage']}% coverage\n")
                    if result['changes']:
                        handle.write(f"    Protein changes vs {result['gene']}: "
                                     f"{', '.join(result['changes'])}\n")
                    if not result['complete']:
                        handle.write('    WARNING: gene runs off a contig boundary\n')
                handle.write('\n')

            if chromosomal:
                handle.write('CHROMOSOMAL TARGET GENES:\n')
                handle.write('-' * 50 + '\n')
                for result in chromosomal:
                    handle.write(f"  {result['gene']} (copy number "
                                 f"{result['copy_number']}): {result['identity']}% identity, "
                                 f"{result['coverage']}% coverage\n")
                    if result['reported_mutations']:
                        handle.write('    Curated resistance mutations: '
                                     f"{', '.join(result['reported_mutations'])}\n")
                    if result['lof_description']:
                        handle.write(f"    Loss of function: {result['lof_description']}\n")
                    if not result['complete']:
                        handle.write('    WARNING: gene runs off a contig boundary\n')
                handle.write('\n')

        if unified_results:
            confirmed = [item for item in unified_results if item['confidence'] == 100]
            handle.write('PROTEIN CHANGES (BLAST caller cross-checked with GAMMA):\n')
            handle.write(f"  {len(confirmed)} of {len(unified_results)} changes were "
                         f"reported by both callers\n")
            handle.write('-' * 70 + '\n')
            by_gene = {}
            for item in unified_results:
                by_gene.setdefault(item['gene'], []).append(item)
            for gene in sorted(by_gene):
                handle.write(f"\n  Gene: {gene}\n")
                for item in by_gene[gene]:
                    flag = '' if item['reported'] else '  (not a curated resistance mutation)'
                    handle.write(f"    {item['mutation']:<20} "
                                 f"confidence {item['confidence']:>3}%  "
                                 f"[{'+'.join(item['methods'])}]{flag}\n")
            handle.write('\n')

        if amplicon_results:
            handle.write('\nDETECTED AMPLICONS (coordinates only):\n')
            handle.write('-' * 50 + '\n')
            for amplicon in amplicon_results:
                handle.write(f"  {amplicon['pair_id']}: {amplicon['contig']}:"
                             f"{amplicon['start']}-{amplicon['end']} "
                             f"({amplicon['length']} bp)\n")


def handle_create_db(args):
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    setup_logger(str(output_dir / args.prefix), args)
    create_db(args.email, args.output, args.prefix)


def handle_acquired(args):
    logger = setup_logger(args.output, args)
    log_tool_versions(logger)
    blast_results = run_acquired_detection(
        args.assembly, args.database, args.output, args.min_id, args.min_cov,
        args.mutations, getattr(args, 'organism', None))
    write_summary(args.output, args.assembly, blast_results, None, None,
                  organism=getattr(args, 'organism', None))


def handle_mutations(args):
    logger = setup_logger(args.output, args)
    log_tool_versions(logger)
    gamma_results, amplicon_results, _, unified_results = run_mutation_detection(
        args.assembly, args.output, args.genes, args.primers,
        blast_results=None, mutation_db_file=getattr(args, 'mutations', None),
        organism=getattr(args, 'organism', None))
    write_summary(args.output, args.assembly, None, gamma_results, amplicon_results,
                  None, unified_results, organism=getattr(args, 'organism', None))


def handle_all(args):
    logger = setup_logger(args.output, args)
    log_tool_versions(logger)

    if args.genes == _DEFAULT_GENES and Path(args.database).exists():
        args.genes = args.database

    if not args.organism:
        print("NOTE: --organism was not given, so curated chromosomal point "
              "mutations will not be evaluated. Acquired genes and blaKPC "
              "variants are organism independent and are still reported.",
              file=sys.stderr)

    blast_results = run_acquired_detection(
        args.assembly, args.database, args.output, args.min_id, args.min_cov,
        args.mutations, args.organism)

    gamma_results, amplicon_results, _, unified_results = run_mutation_detection(
        args.assembly, args.output, args.genes, args.primers,
        blast_results=blast_results, mutation_db_file=args.mutations,
        organism=args.organism)

    write_summary(args.output, args.assembly, blast_results, gamma_results,
                  amplicon_results, None, unified_results, organism=args.organism)


def main():
    parser = argparse.ArgumentParser(description='FOS-CAZAVI Resistance Detector CLI')
    subparsers = parser.add_subparsers(dest='command', required=True)

    parser_db = subparsers.add_parser('create-db', help='Create reference database')
    parser_db.add_argument('-e', '--email', required=True, help='Email for NCBI Entrez')
    parser_db.add_argument('-o', '--output', default='.', help='Output directory [default: .]')
    parser_db.add_argument('-p', '--prefix', default='resistance_db',
                           help='Output filename prefix [default: resistance_db]')
    parser_db.set_defaults(func=handle_create_db)

    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument('-a', '--assembly', required=True, help='Input assembly file (FASTA)')
    parent.add_argument('-o', '--output', required=True, help='Output prefix')
    parent.add_argument('--organism', choices=SUPPORTED_ORGANISMS,
                        help='Organism of the sample. Required for curated chromosomal '
                             'point mutations, which are only meaningful against a '
                             'species-matched reference.')

    parser_acq = subparsers.add_parser('fos-cazavi-acquired', parents=[parent],
                                       help='Detect resistance genes and call their variants')
    parser_acq.add_argument('-d', '--database', default=_DEFAULT_GENES,
                            help='Resistance gene database (FASTA) [default: bundled]')
    parser_acq.add_argument('--mutations', help='Point mutation definitions file (TSV)')
    parser_acq.add_argument('--min_id', type=float, default=90.0,
                            help='Minimum percent identity (default: 90)')
    parser_acq.add_argument('--min_cov', type=float, default=80.0,
                            help='Minimum percent coverage (default: 80)')
    parser_acq.set_defaults(func=handle_acquired)

    parser_mut = subparsers.add_parser('fos-cazavi-mutations', parents=[parent],
                                       help='Run the GAMMA/seqkit analyses on their own')
    parser_mut.add_argument('--genes', default=_DEFAULT_GENES,
                            help='Nucleotide CDS database for GAMMA (FASTA) [default: bundled]')
    parser_mut.add_argument('--primers', default=_DEFAULT_PRIMERS,
                            help='Primer definitions file (TSV) [default: bundled]')
    parser_mut.add_argument('--mutations', help='Point mutation definitions file (TSV)')
    parser_mut.set_defaults(func=handle_mutations)

    parser_all = subparsers.add_parser('fos-cazavi-all', parents=[parent],
                                       help='Run the full detection pipeline')
    parser_all.add_argument('-d', '--database', default=_DEFAULT_GENES,
                            help='Resistance gene database (FASTA) [default: bundled]')
    parser_all.add_argument('--mutations', help='Point mutation definitions file (TSV)')
    parser_all.add_argument('--primers', default=_DEFAULT_PRIMERS,
                            help='Primer definitions file (TSV) [default: bundled]')
    parser_all.add_argument('--genes', default=_DEFAULT_GENES,
                            help='Nucleotide CDS database for GAMMA (FASTA) [default: bundled]')
    parser_all.add_argument('--min_id', type=float, default=90.0,
                            help='Minimum percent identity (default: 90)')
    parser_all.add_argument('--min_cov', type=float, default=80.0,
                            help='Minimum percent coverage (default: 80)')
    parser_all.set_defaults(func=handle_all)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
