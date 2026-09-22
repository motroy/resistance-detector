"""Orthogonal mutation detection with GAMMA, plus amplicon mapping with seqkit.

GAMMA aligns the assembly against the same nucleotide reference database at the
codon level and reports its own list of changes.  Running it alongside the
BLAST-based caller gives a second, independently implemented opinion: a change
reported by both tools is far less likely to be an artefact of either one.

Two things are deliberately *not* done here:

* A change is never invented from the presence of an amplicon.  Several of the
  bundled primers are laboratory mutagenesis primers; a PCR product from them
  says nothing about the genotype of a clinical isolate, so amplicons are
  reported as coordinates only.
* GAMMA's positions are converted to the same numbering the rest of the tool
  uses before the two callers are compared, so agreement means agreement about
  the same residue.
"""

import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from .references import gene_family
from .utils import load_primers
from .variants import sequential_to_ambler

_SUBSTITUTION = re.compile(r'^([A-Z*])(\d+)([A-Z*])$')
_INDEL = re.compile(r'^(\d+)\s*bp\s+(Deletion|Insertion|Duplication)\s+at\s+(\d+)$',
                    re.IGNORECASE)

# Families whose positions are reported in Ambler numbering elsewhere in the
# tool, so GAMMA's sequential positions must be converted before comparison.
_AMBLER_FAMILIES = ('blaKPC', 'blaSHV', 'blaOXA', 'blaCTX-M')


class MutationDetector:
    def __init__(self, assembly, output_prefix, genes_file=None, primers_file=None,
                 mutation_db_file=None, organism=None, threads=1):
        self.assembly = assembly
        self.output_prefix = output_prefix
        self.genes_file = genes_file
        self.primers_file = primers_file
        self.organism = organism
        self.threads = max(1, int(threads or 1))
        self.gamma_results = []
        self.amplicon_results = []
        self.unified_results = []
        self.primers = load_primers(primers_file) if primers_file else {}

    # ------------------------------------------------------------------ GAMMA

    def run_gamma(self):
        if not self.genes_file or not Path(self.genes_file).exists():
            return

        print("Running GAMMA for independent codon-level mutation analysis...")
        output_prefix = f"{self.output_prefix}_gamma"
        output_file = f"{output_prefix}.gamma"

        try:
            subprocess.run(['GAMMA.py', self.assembly, self.genes_file, output_prefix],
                           capture_output=True, text=True, check=True)
        except FileNotFoundError:
            print("WARNING: GAMMA not found, skipping the second-opinion analysis")
            return
        except subprocess.CalledProcessError as error:
            print(f"WARNING: GAMMA failed: {error.stderr}", file=sys.stderr)
            return

        if not Path(output_file).exists():
            print("WARNING: GAMMA output file not found")
            return

        with open(output_file) as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                gene = row['Gene'].rstrip('‡')
                start, stop = int(row['Start']), int(row['Stop'])
                changes = self._parse_gamma_description(row.get('Description', ''), gene)
                self.gamma_results.append({
                    'protein': gene,
                    'contig': row['Contig'],
                    'contig_start': min(start, stop),
                    'contig_end': max(start, stop),
                    'identity': float(row['Codon_Percent']) * 100,
                    'coverage': float(row['Percent_Length']) * 100,
                    'match_type': row['Match_Type'],
                    'mutations': changes,
                })

        print(f"Found {len(self.gamma_results)} gene alignments")

    @staticmethod
    def _parse_gamma_description(description, gene_name):
        """Turn GAMMA's Description field into comparable change labels.

        GAMMA writes substitutions as ``D179Y`` in sequential numbering and
        indels as ``6 bp Deletion at 496`` (a nucleotide offset).  Substitutions
        are renumbered for class A beta-lactamases; indels are kept as a
        descriptive label because the two callers describe them differently.
        """
        if not description:
            return []

        ambler = gene_family(gene_name) in _AMBLER_FAMILIES
        changes = []
        for token in description.split(','):
            token = token.strip()
            if not token or token in ('0', '-'):
                continue

            substitution = _SUBSTITUTION.match(token)
            if substitution:
                reference_aa, position, variant_aa = substitution.groups()
                position = int(position)
                if ambler:
                    position = sequential_to_ambler(position)
                changes.append(f"{reference_aa}{position}{variant_aa}")
                continue

            indel = _INDEL.match(token)
            if indel:
                length, kind, offset = indel.groups()
                codon = (int(offset) - 1) // 3 + 1
                if ambler:
                    codon = sequential_to_ambler(codon)
                changes.append(f"{int(length)}bp{kind.lower()}@codon{codon}")
                continue

            changes.append(token)
        return changes

    # --------------------------------------------------------------- amplicons

    def detect_amplicons(self):
        """Map primer-pair amplicon coordinates with seqkit (reporting only)."""
        if not self.primers:
            return

        print("Running SeqKit for amplicon coordinate mapping...")
        pairs = defaultdict(dict)
        for name, info in self.primers.items():
            pair_id = info.get('pair_id')
            if not pair_id or pair_id == '-':
                continue
            if name.endswith('_F') or 'Fwd' in name or '-F' in name or '_F' in name:
                pairs[pair_id]['F'] = info['seq']
                pairs[pair_id]['F_name'] = name
            elif name.endswith('_R') or 'Rev' in name or '-R' in name or '_R' in name:
                pairs[pair_id]['R'] = info['seq']
                pairs[pair_id]['R_name'] = name

        complete_pairs = {pair_id: pair for pair_id, pair in pairs.items()
                          if 'F' in pair and 'R' in pair}
        if not complete_pairs:
            print("No complete primer pairs identified for amplicon detection")
            return

        primer_file = f"{self.output_prefix}_seqkit_primers.tsv"
        with open(primer_file, 'w') as handle:
            for pair_id, pair in complete_pairs.items():
                handle.write(f"{pair_id}\t{pair['F']}\t{pair['R']}\n")

        try:
            completed = subprocess.run(
                ['seqkit', 'amplicon', '-j', str(self.threads),
                 '-p', primer_file, self.assembly, '--bed'],
                capture_output=True, text=True, check=True)
        except FileNotFoundError:
            print("WARNING: seqkit not found, skipping amplicon mapping")
            return
        except subprocess.CalledProcessError as error:
            print(f"ERROR running seqkit: {error.stderr}", file=sys.stderr)
            return

        for line in completed.stdout.strip().split('\n'):
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 6:
                continue
            contig, start, end, pair_id = fields[0], int(fields[1]), int(fields[2]), fields[3]
            self.amplicon_results.append({
                'pair_id': pair_id,
                'contig': contig,
                'start': start,
                'end': end,
                'length': end - start,
                'f_primer': complete_pairs.get(pair_id, {}).get('F_name', '?'),
                'r_primer': complete_pairs.get(pair_id, {}).get('R_name', '?'),
                'mutations_found': [],
            })

        print(f"  SeqKit: mapped coordinates for {len(self.amplicon_results)} amplicons")

    def analyze_amplicons(self, blast_results=None):
        """Note which detected genes fall inside each mapped amplicon."""
        if not self.amplicon_results or not blast_results:
            return

        for amplicon in self.amplicon_results:
            for result in blast_results:
                if result['contig'] != amplicon['contig']:
                    continue
                start = min(result['start'], result['end']) - 1
                end = max(result['start'], result['end'])
                if max(amplicon['start'], start) < min(amplicon['end'], end):
                    changes = ';'.join(result['changes']) if result['changes'] else 'no changes'
                    amplicon['mutations_found'].append(f"{result['allele']}: {changes}")

    # ------------------------------------------------------------------ merging

    def merge_detection_results(self, blast_results=None):
        """Cross-reference the BLAST caller and GAMMA.

        Confidence is 100% when both callers report the same change at the same
        residue of the same gene copy, and 50% when only one does.
        """
        gamma_by_gene = defaultdict(list)
        for result in self.gamma_results:
            gamma_by_gene[gene_family(result['protein'])].append(result)

        unified = []
        seen = set()

        for result in blast_results or []:
            family = gene_family(result['gene'])
            gamma_changes = {change
                             for gamma in gamma_by_gene.get(family, [])
                             for change in gamma['mutations']}
            gamma_detail = gamma_by_gene.get(family, [None])[0]

            for change in result['changes']:
                key = (family, change, result['contig'], result['start'])
                if key in seen:
                    continue
                seen.add(key)
                confirmed = change in gamma_changes
                unified.append({
                    'gene': result['gene'],
                    'family': family,
                    'mutation': change,
                    'confidence': 100 if confirmed else 50,
                    'methods': ['blast', 'gamma'] if confirmed else ['blast'],
                    'reported': change in result['reported_mutations'] or result['acquired'],
                    'gamma_detail': gamma_detail if confirmed else None,
                    'contig': result['contig'],
                })

        # Changes GAMMA found in genes the BLAST caller did not report at all.
        blast_families = {gene_family(result['gene']) for result in blast_results or []}
        for family, gamma_results in gamma_by_gene.items():
            if family in blast_families:
                continue
            for gamma in gamma_results:
                for change in gamma['mutations']:
                    unified.append({
                        'gene': gamma['protein'],
                        'family': family,
                        'mutation': change,
                        'confidence': 50,
                        'methods': ['gamma'],
                        'reported': False,
                        'gamma_detail': gamma,
                        'contig': gamma['contig'],
                    })

        unified.sort(key=lambda item: (item['gene'], item['mutation']))
        self.unified_results = unified
        return unified

    # ------------------------------------------------------------------ reports

    def write_unified_report(self):
        report_file = f"{self.output_prefix}_unified_mutations.tsv"
        if not self.unified_results:
            if Path(report_file).exists():
                Path(report_file).unlink()
            return

        print(f"Writing unified mutation report to {report_file}...")
        with open(report_file, 'w') as handle:
            handle.write('\t'.join([
                'Gene', 'Mutation', 'Confidence(%)', 'Methods', 'Contig',
                'Reported_As_Resistance_Mutation',
            ]) + '\n')
            for item in self.unified_results:
                handle.write('\t'.join([
                    item['gene'],
                    item['mutation'],
                    str(item['confidence']),
                    '+'.join(item['methods']),
                    item['contig'],
                    'yes' if item['reported'] else 'no',
                ]) + '\n')

    def write_gamma_report(self):
        if not self.gamma_results:
            return
        report_file = f"{self.output_prefix}_protein_mutations.tsv"
        print(f"Writing GAMMA results to {report_file}...")
        with open(report_file, 'w') as handle:
            handle.write('\t'.join(['Gene', 'Contig', 'Start', 'End', 'Identity',
                                    'Coverage', 'Match_Type', 'Changes', 'Method']) + '\n')
            for result in self.gamma_results:
                handle.write('\t'.join([
                    result['protein'], result['contig'], str(result['contig_start']),
                    str(result['contig_end']), f"{result['identity']:.2f}",
                    f"{result['coverage']:.2f}", result['match_type'],
                    ';'.join(result['mutations']) if result['mutations'] else '-',
                    'GAMMA',
                ]) + '\n')

    def write_amplicon_report(self):
        if not self.amplicon_results:
            return
        report_file = f"{self.output_prefix}_amplicons.tsv"
        print(f"Writing amplicon results to {report_file}...")
        with open(report_file, 'w') as handle:
            handle.write('\t'.join(['Pair_ID', 'Contig', 'Start', 'End', 'Length',
                                    'Genes_In_Region', 'Method']) + '\n')
            for amplicon in self.amplicon_results:
                handle.write('\t'.join([
                    amplicon['pair_id'], amplicon['contig'], str(amplicon['start']),
                    str(amplicon['end']), str(amplicon['length']),
                    ';'.join(amplicon['mutations_found']) if amplicon['mutations_found'] else '-',
                    'SeqKit/amplicon',
                ]) + '\n')

    def run(self, blast_results=None):
        self.run_gamma()
        unified = self.merge_detection_results(blast_results)
        self.detect_amplicons()
        self.analyze_amplicons(blast_results)
        self.write_gamma_report()
        self.write_amplicon_report()
        self.write_unified_report()
        return self.gamma_results, self.amplicon_results, [], unified


def run_mutation_detection(assembly, output, genes, primers, blast_results=None,
                           mutation_db_file=None, organism=None, threads=1):
    detector = MutationDetector(assembly, output, genes, primers, mutation_db_file,
                                organism, threads)
    return detector.run(blast_results)
