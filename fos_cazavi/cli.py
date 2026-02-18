import argparse
import sys
from pathlib import Path
from .db import create_db
from .acquired import run_acquired_detection
from .mutations import run_mutation_detection
from .utils import setup_logger, log_tool_versions

_DATA_DIR = Path(__file__).parent / 'data'
_DEFAULT_PROTEINS = str(_DATA_DIR / 'cazavi_proteins.fasta')
_DEFAULT_PRIMERS = str(_DATA_DIR / 'primers.tsv')


def write_combined_tsv(output_prefix, assembly, blast_results, miniprot_results, amplicon_results, unified_results):
    """Write a combined TSV summary of all results"""
    output_file = f"{output_prefix}_all_results.tsv"
    print(f"Writing combined summary to {output_file}...")

    with open(output_file, 'w') as f:
        # Header
        f.write('\t'.join([
            'Sample', 'Type', 'Gene_Protein', 'Method', 'Result_Details',
            'Identity_Confidence', 'Contig', 'Start', 'End', 'Notes'
        ]) + '\n')

        sample = Path(assembly).name

        # 1. BLAST Results (Acquired Genes)
        if blast_results:
            for r in blast_results:
                f.write('\t'.join([
                    sample,
                    'Acquired Gene',
                    r['gene'],
                    'BLAST',
                    f"Coverage: {r['coverage']}%",
                    f"{r['identity']}%",
                    r['contig'],
                    str(r['start']),
                    str(r['end']),
                    f"Mutations: {r['mutations']}"
                ]) + '\n')

        # 2. Unified Mutation Results (Miniprot + SeqKit)
        if unified_results:
            for r in unified_results:
                mp = r['miniprot_detail']
                # Prefer miniprot location if available
                contig = mp['contig'] if mp else '-'
                start = str(mp['contig_start']) if mp else '-'
                end = str(mp['contig_end']) if mp else '-'

                f.write('\t'.join([
                    sample,
                    'Mutation',
                    r['gene'],
                    '+'.join(r['methods']),
                    r['mutation'],
                    f"{r['confidence']}%",
                    contig,
                    start,
                    end,
                    'Dual-method consensus'
                ]) + '\n')

        # 3. Miniprot Results (Wildtype/No mutation detected)
        if miniprot_results:
            for r in miniprot_results:
                if not r['mutations']:
                    f.write('\t'.join([
                        sample,
                        'Protein Alignment',
                        r['protein'],
                        'Miniprot',
                        'Wildtype / No mutations detected',
                        f"{r['identity']:.2f}%",
                        r['contig'],
                        str(r['contig_start']),
                        str(r['contig_end']),
                        f"Coverage: {r['coverage']:.2f}%"
                    ]) + '\n')

        # 4. Amplicon Results (Locations)
        if amplicon_results:
            for r in amplicon_results:
                muts = '; '.join(r['mutations_found']) if r['mutations_found'] else 'None'
                f.write('\t'.join([
                    sample,
                    'Amplicon',
                    r['pair_id'],
                    'SeqKit/Bed',
                    f"Length: {r['length']}bp",
                    '-',
                    r['contig'],
                    str(r['start']),
                    str(r['end']),
                    f"Genes/Mutations in region: {muts}"
                ]) + '\n')


def write_summary(output_prefix, assembly, blast_results, miniprot_results,
                  amplicon_results, seqkit_mut_results=None, unified_results=None):
    """Write a summary of detected resistance mechanisms"""

    # Also write the combined TSV
    write_combined_tsv(output_prefix, assembly, blast_results, miniprot_results, amplicon_results, unified_results)

    summary_file = f"{output_prefix}_summary.txt"

    print(f"Writing summary to {summary_file}...")

    with open(summary_file, 'w') as f:
        f.write("=" * 70 + '\n')
        f.write("FOS-CAZAVI Resistance Detection Summary\n")
        f.write("=" * 70 + '\n\n')

        f.write(f"Assembly: {assembly}\n")

        if blast_results is not None:
            f.write(f"Total genes detected: {len(blast_results)}\n")
            f.write("Method: BLAST+\n\n")

            # Group by drug class
            fos_genes = [r for r in blast_results if r['gene'].startswith('fos')]
            kpc_genes = [r for r in blast_results if 'KPC' in r['gene'].upper()]
            oxa_genes = [r for r in blast_results if 'OXA' in r['gene'].upper()]
            other_genes = [r for r in blast_results
                          if r not in fos_genes + kpc_genes + oxa_genes]

            if fos_genes:
                f.write("FOSFOMYCIN RESISTANCE GENES:\n")
                f.write("-" * 50 + '\n')
                for gene in fos_genes:
                    f.write(f"  {gene['gene']}: {gene['identity']}% identity, "
                           f"{gene['coverage']}% coverage\n")
                    if gene['mutations'] != '-':
                        f.write(f"    Mutations: {gene['mutations']}\n")
                f.write('\n')

            if kpc_genes:
                f.write("CEFTAZIDIME-AVIBACTAM RESISTANCE (KPC):\n")
                f.write("-" * 50 + '\n')
                for gene in kpc_genes:
                    f.write(f"  {gene['gene']}: {gene['identity']}% identity, "
                           f"{gene['coverage']}% coverage\n")
                    if gene['mutations'] != '-':
                        f.write(f"    Mutations: {gene['mutations']}\n")
                f.write('\n')

            if oxa_genes:
                f.write("CEFTAZIDIME-AVIBACTAM RESISTANCE (OXA):\n")
                f.write("-" * 50 + '\n')
                for gene in oxa_genes:
                    f.write(f"  {gene['gene']}: {gene['identity']}% identity, "
                           f"{gene['coverage']}% coverage\n")
                    if gene['mutations'] != '-':
                        f.write(f"    Mutations: {gene['mutations']}\n")
                f.write('\n')

            if other_genes:
                f.write("OTHER RESISTANCE GENES:\n")
                f.write("-" * 50 + '\n')
                for gene in other_genes:
                    f.write(f"  {gene['gene']}: {gene['identity']}% identity, "
                           f"{gene['coverage']}% coverage\n")
                    if gene['mutations'] != '-':
                        f.write(f"    Mutations: {gene['mutations']}\n")
                f.write('\n')

            if not blast_results:
                f.write("No resistance genes detected\n")

        # ── Unified dual-method mutation detection ────────────────────────────
        if unified_results:
            f.write("\n")
            f.write("MUTATION DETECTION SUMMARY (Dual-Method):\n")
            f.write("Methods: Miniprot (protein alignment) + SeqKit (amplicon)\n")
            f.write("Confidence: 100% = detected by both methods; "
                    "50% = detected by one method\n")
            f.write("-" * 70 + '\n')

            # Group by gene for readability
            from collections import defaultdict
            by_gene = defaultdict(list)
            for r in unified_results:
                by_gene[r['gene']].append(r)

            for gene in sorted(by_gene):
                f.write(f"\n  Gene: {gene}\n")
                for r in by_gene[gene]:
                    methods_str = '+'.join(r['methods'])
                    f.write(f"    {r['mutation']:<15} "
                            f"Confidence: {r['confidence']:>3}%  "
                            f"[{methods_str}]\n")
            f.write('\n')
        elif miniprot_results is not None or (seqkit_mut_results is not None):
            # At least one method ran but found no cross-confirmed mutations
            f.write("\n")
            f.write("MUTATION DETECTION SUMMARY (Dual-Method):\n")
            f.write("-" * 70 + '\n')
            f.write("  No resistance mutations detected by miniprot or seqkit.\n\n")

        if amplicon_results:
            f.write("\n")
            f.write("DETECTED AMPLICONS:\n")
            f.write("Method: Seqkit/Amplicon\n")
            f.write("-" * 50 + '\n')
            for amp in amplicon_results:
                f.write(f"  Pair: {amp['pair_id']}\n")
                f.write(f"    Location: {amp['contig']}:{amp['start']}-{amp['end']} ({amp['length']} bp)\n")
                f.write(f"    Primers: {amp['f_primer']} -> {amp['r_primer']}\n")
                if amp['mutations_found']:
                    f.write("    Mutations/Genes in amplicon:\n")
                    for m in amp['mutations_found']:
                        f.write(f"      - {m}\n")
                else:
                    f.write("    No resistance genes/mutations detected in amplicon\n")

        if miniprot_results:
            f.write("\n")
            f.write("PROTEIN MUTATION ANALYSIS (Miniprot):\n")
            f.write("Method: Miniprot (protein-to-genome alignment)\n")
            f.write("-" * 50 + '\n')
            f.write("Reference: doi.org/10.3389/fcimb.2025.1645042\n\n")

            # Group by protein type
            porins = [r for r in miniprot_results if 'Omp' in r['protein']]
            efflux = [r for r in miniprot_results if 'acr' in r['protein'].lower()]

            if porins:
                f.write("Outer Membrane Porins:\n")
                for r in porins:
                    f.write(f"  {r['protein']}: {r['identity']:.1f}% identity, {r['coverage']:.1f}% coverage\n")
                    f.write(f"    Location: {r['contig']}:{r['contig_start']}-{r['contig_end']}\n")
                    if r['mutations']:
                        f.write(f"    Mutations: {', '.join(r['mutations'])}\n")
                    else:
                        f.write("    No mutations detected\n")

            if efflux:
                f.write("\nEfflux Pump Components:\n")
                for r in efflux:
                    f.write(f"  {r['protein']}: {r['identity']:.1f}% identity, {r['coverage']:.1f}% coverage\n")
                    f.write(f"    Location: {r['contig']}:{r['contig_start']}-{r['contig_end']}\n")
                    if r['mutations']:
                        f.write(f"    Mutations: {', '.join(r['mutations'])}\n")
                    else:
                        f.write("    No mutations detected\n")


def handle_create_db(args):
    # Setup logging
    logger = setup_logger(args.output, args)
    create_db(args.email, args.output)


def handle_acquired(args):
    logger = setup_logger(args.output, args)
    log_tool_versions(logger)

    blast_results = run_acquired_detection(
        args.assembly,
        args.database,
        args.output,
        args.min_id,
        args.min_cov,
        args.mutations
    )
    # Write summary only for acquired part
    write_summary(args.output, args.assembly, blast_results, None, None)


def handle_mutations(args):
    logger = setup_logger(args.output, args)
    log_tool_versions(logger)

    # Run miniprot + seqkit mutation detection with cross-method confidence scoring
    miniprot_results, amplicon_results, seqkit_mut_results, unified_results = \
        run_mutation_detection(
            args.assembly,
            args.output,
            args.proteins,
            args.primers,
            blast_results=None
        )
    write_summary(
        args.output, args.assembly,
        None, miniprot_results, amplicon_results,
        seqkit_mut_results, unified_results
    )


def handle_all(args):
    logger = setup_logger(args.output, args)
    log_tool_versions(logger)

    # Auto-detect protein database if using default and custom BLAST DB provided
    if args.proteins == _DEFAULT_PROTEINS:
        db_path = Path(args.database)
        # Check for _proteins.fasta sibling
        if db_path.suffix == '.fasta':
            prot_path = db_path.parent / (db_path.stem + '_proteins.fasta')
        else:
            prot_path = Path(str(db_path) + '_proteins.fasta')

        if prot_path.exists():
            print(f"Auto-detected protein database: {prot_path}")
            args.proteins = str(prot_path)

    # Run BLAST
    blast_results = run_acquired_detection(
        args.assembly,
        args.database,
        args.output,
        args.min_id,
        args.min_cov,
        args.mutations
    )

    # Run Mutations (Miniprot + SeqKit) with BLAST results for cross-ref
    miniprot_results, amplicon_results, seqkit_mut_results, unified_results = \
        run_mutation_detection(
            args.assembly,
            args.output,
            args.proteins,
            args.primers,
            blast_results=blast_results
        )

    write_summary(
        args.output, args.assembly,
        blast_results, miniprot_results, amplicon_results,
        seqkit_mut_results, unified_results
    )


def main():
    parser = argparse.ArgumentParser(description='FOS-CAZAVI Resistance Detector CLI')
    subparsers = parser.add_subparsers(dest='command', required=True)

    # create-db
    parser_db = subparsers.add_parser('create-db', help='Create reference database')
    parser_db.add_argument('-e', '--email', required=True, help='Email for NCBI Entrez')
    parser_db.add_argument('-o', '--output', default='resistance_db', help='Output prefix')
    parser_db.set_defaults(func=handle_create_db)

    # Common arguments for detection commands
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-a', '--assembly', required=True, help='Input assembly file (FASTA)')
    parent_parser.add_argument('-o', '--output', required=True, help='Output prefix')

    # acquired
    parser_acq = subparsers.add_parser('fos-cazavi-acquired', parents=[parent_parser], help='Detect acquired resistance genes')
    parser_acq.add_argument('-d', '--database', required=True, help='Resistance gene database (FASTA)')
    parser_acq.add_argument('--mutations', help='Mutation definitions file (TSV)')
    parser_acq.add_argument('--min_id', type=float, default=90.0, help='Minimum percent identity (default: 90)')
    parser_acq.add_argument('--min_cov', type=float, default=80.0, help='Minimum percent coverage (default: 80)')
    parser_acq.set_defaults(func=handle_acquired)

    # mutations
    parser_mut = subparsers.add_parser('fos-cazavi-mutations', parents=[parent_parser], help='Detect resistance mutations (Miniprot/SeqKit dual-method)')
    parser_mut.add_argument('--proteins', default=_DEFAULT_PROTEINS, help=f'Protein sequences for miniprot (FASTA) [default: bundled]')
    parser_mut.add_argument('--primers', default=_DEFAULT_PRIMERS, help=f'Primers definitions file (TSV) [default: bundled]')
    parser_mut.set_defaults(func=handle_mutations)

    # all
    parser_all = subparsers.add_parser('fos-cazavi-all', parents=[parent_parser], help='Run full detection pipeline')
    parser_all.add_argument('-d', '--database', required=True, help='Resistance gene database (FASTA)')
    parser_all.add_argument('--mutations', help='Mutation definitions file (TSV)')
    parser_all.add_argument('--primers', default=_DEFAULT_PRIMERS, help=f'Primers definitions file (TSV) [default: bundled]')
    parser_all.add_argument('--proteins', default=_DEFAULT_PROTEINS, help=f'Protein sequences for miniprot (FASTA) [default: bundled]')
    parser_all.add_argument('--min_id', type=float, default=90.0, help='Minimum percent identity (default: 90)')
    parser_all.add_argument('--min_cov', type=float, default=80.0, help='Minimum percent coverage (default: 80)')
    parser_all.set_defaults(func=handle_all)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
