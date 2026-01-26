import argparse
import sys
from pathlib import Path
from .db import create_db
from .acquired import run_acquired_detection
from .mutations import run_mutation_detection

def write_summary(output_prefix, assembly, blast_results, miniprot_results, amplicon_results):
    """Write a summary of detected resistance mechanisms"""
    summary_file = f"{output_prefix}_summary.txt"

    print(f"Writing summary to {summary_file}...")

    with open(summary_file, 'w') as f:
        f.write("=" * 70 + '\n')
        f.write("FOS-CAZAVI Resistance Detection Summary\n")
        f.write("=" * 70 + '\n\n')

        f.write(f"Assembly: {assembly}\n")

        if blast_results is not None:
            f.write(f"Total genes detected: {len(blast_results)}\n\n")

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

        if amplicon_results:
            f.write("\n")
            f.write("DETECTED AMPLICONS:\n")
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
            f.write("PROTEIN MUTATION ANALYSIS (miniprot):\n")
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
    create_db(args.email, args.output)

def handle_acquired(args):
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
    # This command runs miniprot and amplicons.
    # It does not run BLAST, so amplicon cross-referencing is limited.
    miniprot_results, amplicon_results = run_mutation_detection(
        args.assembly,
        args.output,
        args.proteins,
        args.primers,
        blast_results=None
    )
    write_summary(args.output, args.assembly, None, miniprot_results, amplicon_results)

def handle_all(args):
    # Run BLAST
    blast_results = run_acquired_detection(
        args.assembly,
        args.database,
        args.output,
        args.min_id,
        args.min_cov,
        args.mutations
    )

    # Run Mutations (Miniprot + Amplicons) with BLAST results for cross-ref
    miniprot_results, amplicon_results = run_mutation_detection(
        args.assembly,
        args.output,
        args.proteins,
        args.primers,
        blast_results=blast_results
    )

    write_summary(args.output, args.assembly, blast_results, miniprot_results, amplicon_results)

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
    parser_mut = subparsers.add_parser('fos-cazavi-mutations', parents=[parent_parser], help='Detect resistance mutations (Miniprot/Amplicons)')
    parser_mut.add_argument('--proteins', help='Protein sequences for miniprot (FASTA)')
    parser_mut.add_argument('--primers', help='Primers definitions file (TSV)')
    parser_mut.set_defaults(func=handle_mutations)

    # all
    parser_all = subparsers.add_parser('fos-cazavi-all', parents=[parent_parser], help='Run full detection pipeline')
    parser_all.add_argument('-d', '--database', required=True, help='Resistance gene database (FASTA)')
    parser_all.add_argument('--mutations', help='Mutation definitions file (TSV)')
    parser_all.add_argument('--primers', help='Primers definitions file (TSV)')
    parser_all.add_argument('--proteins', help='Protein sequences for miniprot (FASTA)')
    parser_all.add_argument('--min_id', type=float, default=90.0, help='Minimum percent identity (default: 90)')
    parser_all.add_argument('--min_cov', type=float, default=80.0, help='Minimum percent coverage (default: 80)')
    parser_all.set_defaults(func=handle_all)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
