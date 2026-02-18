import sys
import os
import subprocess
from pathlib import Path
from collections import defaultdict
import re
from .utils import load_primers, detect_mutations as _detect_mutations, KNOWN_MUTATIONS


class MutationDetector:
    def __init__(self, assembly, output_prefix, proteins_file=None, primers_file=None):
        self.assembly = assembly
        self.output_prefix = output_prefix
        self.proteins_file = proteins_file
        self.primers_file = primers_file
        self.miniprot_results = []
        self.amplicon_results = []
        self.seqkit_mut_results = []
        self.unified_results = []

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
            #'--no-cs',  # Output cs tag
            '--outs=0.95',  # High identity threshold
            '--outc=0.1', # Output an alignment only if FLOAT fraction of the query protein is aligned [defaul=0.1]
            '--trans', # Output translated protein sequences on ‘##STA’ lines.
            '-t', '1',
            '-j', '0',
            '-C', '0',
            '-S',
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

    @staticmethod
    def _normalize_gene_name(gene_name):
        """
        Normalize gene name for cross-method comparison.

        Handles both primer-style names ('uhpB', 'blaKPC') and miniprot
        protein FASTA ID style ('acrB_WP_002892069.1', 'blaKPC-3').
        Returns a lowercase base gene name for comparison.
        """
        # Strip accession suffix (e.g., acrB_WP_002892069.1 -> acrB)
        base = gene_name.split('_')[0]
        # Strip numeric variant suffix after hyphen (e.g., blaKPC-3 -> blaKPC)
        base = base.split('-')[0]
        # For bla genes, strip trailing digits from the class part (e.g., blaOXA -> blaOXA)
        if base.lower().startswith('bla'):
            base = 'bla' + base[3:].rstrip('0123456789')
        return base.lower()

    def _extract_pair_id_from_header(self, header, valid_pairs):
        """
        Extract pair_id from a seqkit amplicon FASTA header.

        Seqkit amplicon headers embed the primer pair name in the header
        (e.g., '>contig1_uhpB_ver:100-500'). Scans known pair IDs.
        """
        for pair_id in valid_pairs:
            if pair_id in header:
                return pair_id
        return None

    def _build_primer_pairs(self):
        """
        Group primers into forward/reverse pairs from self.primers.
        Returns {pair_id: {'F': seq, 'R': seq, 'gene': str, 'mutation': str}}.
        """
        pairs = defaultdict(dict)
        for name, info in self.primers.items():
            pair_id = info.get('pair_id')
            if not pair_id or pair_id == '-':
                continue

            gene = info.get('gene', '') or ''
            mutation_anno = info.get('mutation')

            if pair_id not in pairs:
                pairs[pair_id] = {
                    'gene': gene,
                    'mutation': mutation_anno,
                }

            # Determine forward/reverse using the same heuristics as detect_amplicons
            if name.endswith('_F') or 'Fwd' in name or '-F' in name or '_F' in name:
                pairs[pair_id]['F'] = info['seq']
                pairs[pair_id]['F_name'] = name
            elif name.endswith('_R') or 'Rev' in name or '-R' in name or '_R' in name:
                pairs[pair_id]['R'] = info['seq']
                pairs[pair_id]['R_name'] = name
        return pairs

    def detect_seqkit_mutations(self):
        """
        Detect resistance mutations using the seqkit amplicon method.

        Runs seqkit amplicon to extract amplicon sequences, then:
        - For mutation-annotated primer pairs: amplicon presence confirms
          the annotated mutation (primer is mutation-specific).
        - For gene-verification primer pairs: translates the amplicon in
          all 3 frames and calls detect_mutations() to find mutations.

        Returns a list of dicts:
            {gene, pair_id, mutations (list), method='seqkit'}
        """
        if not self.primers:
            return []

        print("Running SeqKit for targeted mutation detection (amplicon extraction)...")

        pairs = self._build_primer_pairs()

        # Keep only pairs that have both primers and a gene annotation
        _skip_purposes = ('deletion', 'cloning', 'pBAD', 'pDS', 'RT-qPCR', 'Quantification')
        valid_pairs = {}
        for pid, p in pairs.items():
            if 'F' not in p or 'R' not in p:
                continue
            gene = p.get('gene', '') or ''
            if not gene or gene in ('-', ''):
                continue
            # Skip cloning/RT-qPCR pairs that are not for mutation/gene verification
            purpose_hint = p.get('mutation', '') or ''
            if any(skip in purpose_hint for skip in _skip_purposes):
                continue
            valid_pairs[pid] = p

        if not valid_pairs:
            print("  No valid primer pairs for seqkit mutation detection")
            return []

        seqkit_primer_file = f"{self.output_prefix}_seqkit_mut_primers.tsv"
        with open(seqkit_primer_file, 'w') as f:
            for pair_id, p in valid_pairs.items():
                f.write(f"{pair_id}\t{p['F']}\t{p['R']}\n")

        seqkit_mut_results = []

        try:
            result = subprocess.run(
                ['seqkit', 'amplicon', '-p', seqkit_primer_file, self.assembly],
                capture_output=True, text=True
            )

            if not result.stdout.strip():
                print("  SeqKit: no amplicons found")
                return []

            # Parse FASTA output
            current_header = None
            current_seq_parts = []
            amplicons = []

            for line in result.stdout.strip().split('\n'):
                if line.startswith('>'):
                    if current_header is not None and current_seq_parts:
                        amplicons.append((current_header, ''.join(current_seq_parts)))
                    current_header = line[1:]
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line.strip())

            if current_header is not None and current_seq_parts:
                amplicons.append((current_header, ''.join(current_seq_parts)))

            print(f"  SeqKit (Targeted): Extracted {len(amplicons)} amplicons for mutation analysis")

            for header, seq in amplicons:
                pair_id = self._extract_pair_id_from_header(header, valid_pairs)
                if not pair_id:
                    continue

                p = valid_pairs[pair_id]
                gene = p.get('gene', '')
                mutation_anno = p.get('mutation')

                if not gene or gene in ('-', ''):
                    continue

                mutations_found = []

                if mutation_anno and mutation_anno.strip() and mutation_anno.strip() not in ('-', ''):
                    # Mutation-specific primer: amplicon presence = mutation confirmed.
                    # Annotation format: "<gene> <MutName>" e.g. "uhpB G469R"
                    parts = mutation_anno.strip().split()
                    # Exclude non-mutation annotations
                    skip_words = {'deletion', 'cloning', 'verification', 'expression', 'expression(reverse)'}
                    if len(parts) >= 2 and not any(w.lower() in skip_words for w in parts):
                        mut_name = parts[-1]  # e.g. "G469R"
                        mutations_found = [mut_name]
                else:
                    # Gene-level verification primer: translate & check for mutations.
                    # Try all 3 reading frames to account for different amplicon starts.
                    seen = set()
                    for frame in range(3):
                        frame_seq = seq[frame:]
                        frame_muts = _detect_mutations(gene, frame_seq, KNOWN_MUTATIONS)
                        for m in frame_muts:
                            if m not in seen:
                                seen.add(m)
                                mutations_found.append(m)

                if mutations_found:
                    seqkit_mut_results.append({
                        'gene': gene,
                        'pair_id': pair_id,
                        'mutations': mutations_found,
                        'method': 'seqkit',
                    })

        except FileNotFoundError:
            print("WARNING: seqkit not found, skipping seqkit mutation analysis")
        except subprocess.CalledProcessError as e:
            print(f"WARNING: seqkit amplicon failed: {e.stderr}", file=sys.stderr)
        finally:
            try:
                os.remove(seqkit_primer_file)
            except Exception:
                pass

        self.seqkit_mut_results = seqkit_mut_results
        gene_count = len({r['gene'] for r in seqkit_mut_results})
        if gene_count > 0:
            print(f"  SeqKit (Targeted): Detected mutations in {gene_count} gene(s) across {len(seqkit_mut_results)} amplicon(s)")
        else:
            print("  SeqKit (Targeted): No specific mutations detected in extracted amplicons")
        return seqkit_mut_results

    def merge_detection_results(self, seqkit_mut_results):
        """
        Merge miniprot and seqkit mutation results with confidence scores.

        Confidence scoring:
          - 100% : mutation found by BOTH miniprot AND seqkit
          -  50% : mutation found by only ONE method

        Returns a list of dicts (sorted by gene then mutation):
        {
            'gene'            : normalized gene name (str),
            'mutation'        : mutation string e.g. 'D179Y' (str),
            'confidence'      : 50 or 100 (int),
            'methods'         : list of method names, e.g. ['miniprot', 'seqkit'],
            'miniprot_detail' : matching miniprot result dict, or None,
            'seqkit_detail'   : matching seqkit result dict, or None,
        }
        """
        # Index miniprot results by normalized (gene, mutation)
        miniprot_by_key = {}
        for r in self.miniprot_results:
            gene_norm = self._normalize_gene_name(r['protein'])
            for mut in r.get('mutations', []):
                key = (gene_norm, mut)
                if key not in miniprot_by_key:
                    miniprot_by_key[key] = r

        # Index seqkit results by normalized (gene, mutation)
        seqkit_by_key = {}
        for r in seqkit_mut_results:
            gene_norm = self._normalize_gene_name(r['gene'])
            for mut in r.get('mutations', []):
                key = (gene_norm, mut)
                if key not in seqkit_by_key:
                    seqkit_by_key[key] = r

        all_keys = set(miniprot_by_key.keys()) | set(seqkit_by_key.keys())

        unified = []
        for (gene, mutation) in sorted(all_keys):
            in_miniprot = (gene, mutation) in miniprot_by_key
            in_seqkit = (gene, mutation) in seqkit_by_key

            methods = []
            if in_miniprot:
                methods.append('miniprot')
            if in_seqkit:
                methods.append('seqkit')

            confidence = 100 if (in_miniprot and in_seqkit) else 50

            unified.append({
                'gene': gene,
                'mutation': mutation,
                'confidence': confidence,
                'methods': methods,
                'miniprot_detail': miniprot_by_key.get((gene, mutation)),
                'seqkit_detail': seqkit_by_key.get((gene, mutation)),
            })

        self.unified_results = unified
        return unified

    def write_unified_report(self):
        """Write the unified dual-method mutation report (TSV)."""
        if not self.unified_results:
            return

        report_file = f"{self.output_prefix}_unified_mutations.tsv"
        print(f"Writing unified mutation report to {report_file}...")

        with open(report_file, 'w') as f:
            f.write('\t'.join([
                'Gene', 'Mutation', 'Confidence(%)', 'Methods',
                'Miniprot_Contig', 'Miniprot_Identity(%)', 'Miniprot_Coverage(%)',
                'SeqKit_PairID',
            ]) + '\n')

            for r in self.unified_results:
                mp = r['miniprot_detail']
                sk = r['seqkit_detail']

                mp_contig = mp['contig'] if mp else '-'
                mp_identity = f"{mp['identity']:.2f}" if mp else '-'
                mp_coverage = f"{mp['coverage']:.2f}" if mp else '-'
                sk_pair = sk['pair_id'] if sk else '-'

                f.write('\t'.join([
                    r['gene'],
                    r['mutation'],
                    str(r['confidence']),
                    '+'.join(r['methods']),
                    mp_contig,
                    mp_identity,
                    mp_coverage,
                    sk_pair,
                ]) + '\n')

    def detect_amplicons(self):
        """Detect amplicon coordinates using seqkit amplicon (BED output)"""
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

            print(f"  SeqKit (Mapping): Mapped coordinates for {len(self.amplicon_results)} amplicons")

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
                              'Coverage', 'Mutations', 'CS_Tag', 'Method']) + '\n')

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
                    r['cs_tag'] or '-',
                    'Miniprot'
                ]) + '\n')

    def write_amplicon_report(self):
        """Write amplicon detection results to TSV file"""
        if not self.amplicon_results:
            return

        report_file = f"{self.output_prefix}_amplicons.tsv"

        print(f"Writing amplicon results to {report_file}...")

        with open(report_file, 'w') as f:
            f.write('\t'.join(['Pair_ID', 'Contig', 'Start', 'End', 'Length',
                              'Mutations_Found', 'Method']) + '\n')

            for amp in self.amplicon_results:
                mut_str = ';'.join(amp['mutations_found']) if amp['mutations_found'] else '-'
                f.write('\t'.join([
                    amp['pair_id'],
                    amp['contig'],
                    str(amp['start']),
                    str(amp['end']),
                    str(amp['length']),
                    mut_str,
                    'Seqkit/Amplicon'
                ]) + '\n')

    def run(self, blast_results=None):
        """
        Run the full mutation detection pipeline:
          1. miniprot  – protein-to-genome alignment, CS-tag mutation parsing
          2. seqkit    – amplicon extraction + sequence-level mutation detection
          3. merge     – cross-reference both methods, assign confidence scores
          4. seqkit BED – amplicon location detection (cross-ref with BLAST)

        Returns:
            (miniprot_results, amplicon_results, seqkit_mut_results, unified_results)
        """
        self.run_miniprot()
        seqkit_mut_results = self.detect_seqkit_mutations()
        unified_results = self.merge_detection_results(seqkit_mut_results)

        self.detect_amplicons()
        self.analyze_amplicons(blast_results)

        self.write_miniprot_report()
        self.write_amplicon_report()
        self.write_unified_report()

        return self.miniprot_results, self.amplicon_results, seqkit_mut_results, unified_results


def run_mutation_detection(assembly, output, proteins, primers, blast_results=None):
    detector = MutationDetector(assembly, output, proteins, primers)
    return detector.run(blast_results)
