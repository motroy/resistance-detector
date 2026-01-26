import sys
import subprocess
from pathlib import Path
from collections import defaultdict
import re
from .utils import load_primers

class MutationDetector:
    def __init__(self, assembly, output_prefix, proteins_file=None, primers_file=None):
        self.assembly = assembly
        self.output_prefix = output_prefix
        self.proteins_file = proteins_file
        self.primers_file = primers_file
        self.miniprot_results = []
        self.amplicon_results = []

        if primers_file:
            self.primers = load_primers(primers_file)
        else:
            self.primers = {}

    def run_miniprot(self):
        """
        Run miniprot to detect mutations in CAZAVI resistance proteins.
        """
        if not self.proteins_file or not Path(self.proteins_file).exists():
            return

        print(f"Running miniprot for protein mutation analysis...")

        miniprot_output = f"{self.output_prefix}_miniprot.paf"

        cmd = [
            'miniprot',
            '-c',  # Output cs tag
            '--outs=0.95',  # High identity threshold
            '-t', '1',
            self.assembly,
            self.proteins_file
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            with open(miniprot_output, 'w') as f:
                f.write(result.stdout)

            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue

                fields = line.split('\t')
                if len(fields) < 12:
                    continue

                protein_name = fields[0]
                protein_len = int(fields[1])
                protein_start = int(fields[2])
                protein_end = int(fields[3])
                contig = fields[5]
                contig_start = int(fields[7])
                contig_end = int(fields[8])
                matches = int(fields[9])
                alignment_len = int(fields[10])

                identity = (matches / alignment_len * 100) if alignment_len > 0 else 0
                coverage = ((protein_end - protein_start) / protein_len * 100) if protein_len > 0 else 0

                mutations = []
                cs_tag = None
                for field in fields[12:]:
                    if field.startswith('cs:Z:'):
                        cs_tag = field[5:]
                        break

                if cs_tag:
                    mutations = self.parse_miniprot_cs(cs_tag, protein_name)

                self.miniprot_results.append({
                    'protein': protein_name,
                    'contig': contig,
                    'contig_start': contig_start,
                    'contig_end': contig_end,
                    'protein_start': protein_start,
                    'protein_end': protein_end,
                    'identity': identity,
                    'coverage': coverage,
                    'mutations': mutations,
                    'cs_tag': cs_tag
                })

            print(f"Found {len(self.miniprot_results)} protein alignments")
            for r in self.miniprot_results:
                if r['mutations']:
                    print(f"  {r['protein']}: {len(r['mutations'])} mutations detected")

        except FileNotFoundError:
            print("WARNING: miniprot not found, skipping protein analysis")
        except subprocess.CalledProcessError as e:
            print(f"WARNING: miniprot failed: {e.stderr}", file=sys.stderr)

    def parse_miniprot_cs(self, cs_tag, protein_name):
        """Parse miniprot cs tag to extract mutations."""
        mutations = []
        position = 1  # 1-indexed amino acid position
        pattern = r'(:[0-9]+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[a-z]+|~[a-z]+[0-9]+[a-z]+|[a-z]+[A-Z\*])'

        for match in re.finditer(pattern, cs_tag):
            op = match.group(1)

            if op.startswith(':'):
                num_identical = int(op[1:])
                position += num_identical

            elif op.startswith('*'):
                ref_aa = op[1].upper()
                query_aa = op[2].upper()
                mutations.append(f"{ref_aa}{position}{query_aa}")
                position += 1

            elif op.startswith('+'):
                inserted = op[1:]
                mutations.append(f"ins{position}_{inserted}")

            elif op.startswith('-'):
                deleted_nt = op[1:]
                deleted_aa_count = len(deleted_nt) // 3
                if deleted_aa_count > 0:
                    mutations.append(f"del{position}_{deleted_aa_count}aa")
                position += deleted_aa_count

            elif op.startswith('~'):
                pass

            elif op[:-1].islower() and op[-1].isupper():
                query_aa = op[-1]
                mutations.append(f"?{position}{query_aa}")
                position += 1

            elif op[-1] == '*':
                mutations.append(f"?{position}*")
                position += 1

        return mutations

    def detect_amplicons(self):
        """Detect amplicons using seqkit amplicon"""
        if not self.primers:
            return

        print("Detecting amplicons with seqkit...")

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

        if not pairs:
            print("No primer pairs identified for amplicon detection")
            return

        seqkit_primer_file = f"{self.output_prefix}_seqkit_primers.tsv"
        with open(seqkit_primer_file, 'w') as f:
            for pair_id, p in pairs.items():
                if 'F' in p and 'R' in p:
                    f.write(f"{pair_id}\t{p['F']}\t{p['R']}\n")

        amplicon_output = f"{self.output_prefix}_amplicons.fasta"
        cmd = [
            'seqkit', 'amplicon',
            '-p', seqkit_primer_file,
            self.assembly,
            '--bed'
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            for line in result.stdout.strip().split('\n'):
                if not line: continue
                fields = line.split('\t')
                if len(fields) < 6: continue

                contig = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                pair_id = fields[3]
                strand = fields[5]

                # Store f_primer and r_primer names if needed, but pairs dict has them
                f_primer = pairs[pair_id].get('F_name', '?')
                r_primer = pairs[pair_id].get('R_name', '?')

                self.amplicon_results.append({
                    'pair_id': pair_id,
                    'contig': contig,
                    'start': start,
                    'end': end,
                    'length': end - start,
                    'f_primer': f_primer,
                    'r_primer': r_primer,
                    'mutations_found': []
                })

            print(f"Found {len(self.amplicon_results)} amplicons")

        except subprocess.CalledProcessError as e:
            print(f"ERROR running seqkit: {e.stderr}", file=sys.stderr)
        except Exception as e:
            print(f"ERROR processing amplicons: {e}", file=sys.stderr)

    def analyze_amplicons(self, blast_results=None):
        """
        Check if detected genes fall within amplicons.
        blast_results: list of dictionaries from BlastDetector
        """
        if not self.amplicon_results:
            return

        if not blast_results:
            # Cannot cross-reference without blast results
            # But maybe we can just verify the amplicon exists?
            return

        print("Checking for mutations within amplicons...")

        for amp in self.amplicon_results:
            amp_contig = amp['contig']
            amp_start = amp['start']
            amp_end = amp['end']

            for res in blast_results:
                if res['contig'] != amp_contig:
                    continue

                res_start = res['start']
                res_end = res['end']

                start = min(res_start, res_end) - 1
                end = max(res_start, res_end)

                if max(amp_start, start) < min(amp_end, end):
                    mut_str = res['mutations']
                    if mut_str != '-':
                        amp['mutations_found'].append(f"{res['gene']}: {mut_str}")
                    else:
                        amp['mutations_found'].append(f"{res['gene']}: (wildtype)")

    def write_miniprot_report(self):
        """Write miniprot mutation detection results"""
        if not self.miniprot_results:
            return

        report_file = f"{self.output_prefix}_protein_mutations.tsv"
        print(f"Writing protein mutation results to {report_file}...")

        with open(report_file, 'w') as f:
            f.write('\t'.join(['Protein', 'Contig', 'Start', 'End', 'Identity',
                              'Coverage', 'Mutations', 'CS_Tag']) + '\n')

            for r in self.miniprot_results:
                mutations_str = ';'.join(r['mutations']) if r['mutations'] else '-'
                f.write('\t'.join([
                    r['protein'],
                    r['contig'],
                    str(r['contig_start']),
                    str(r['contig_end']),
                    f"{r['identity']:.2f}",
                    f"{r['coverage']:.2f}",
                    mutations_str,
                    r['cs_tag'] or '-'
                ]) + '\n')

    def write_amplicon_report(self):
        """Write amplicon detection results to TSV file"""
        if not self.amplicon_results:
            return

        report_file = f"{self.output_prefix}_amplicons.tsv"

        print(f"Writing amplicon results to {report_file}...")

        with open(report_file, 'w') as f:
            f.write('\t'.join(['Pair_ID', 'Contig', 'Start', 'End', 'Length',
                              'Mutations_Found']) + '\n')

            for amp in self.amplicon_results:
                mut_str = ';'.join(amp['mutations_found']) if amp['mutations_found'] else '-'
                f.write('\t'.join([
                    amp['pair_id'],
                    amp['contig'],
                    str(amp['start']),
                    str(amp['end']),
                    str(amp['length']),
                    mut_str
                ]) + '\n')

    def run(self, blast_results=None):
        self.run_miniprot()
        self.detect_amplicons()
        self.analyze_amplicons(blast_results)

        self.write_miniprot_report()
        self.write_amplicon_report()
        return self.miniprot_results, self.amplicon_results

def run_mutation_detection(assembly, output, proteins, primers, blast_results=None):
    detector = MutationDetector(assembly, output, proteins, primers)
    return detector.run(blast_results)
