import sys
import os
import csv
import subprocess
from pathlib import Path
from collections import defaultdict
import re
from .utils import load_primers, load_mutation_db, detect_mutations as _detect_mutations, KNOWN_MUTATIONS


class MutationDetector:
    def __init__(self, assembly, output_prefix, genes_file=None, primers_file=None, mutation_db_file=None):
        self.assembly = assembly
        self.output_prefix = output_prefix
        self.genes_file = genes_file
        self.primers_file = primers_file
        self.gamma_results = []
        self.amplicon_results = []
        self.seqkit_mut_results = []
        self.unified_results = []

        if primers_file:
            self.primers = load_primers(primers_file)
        else:
            self.primers = {}

        # Load mutation database for filtering GAMMA results
        # Normalize keys to lowercase so lookups via _normalize_gene_name succeed
        raw_db = load_mutation_db(mutation_db_file) if mutation_db_file else KNOWN_MUTATIONS
        self.mutation_db = {
            self._normalize_gene_name(k): v
            for k, v in raw_db.items()
        }

    def run_gamma(self):
        """
        Run GAMMA to detect mutations in resistance genes.

        GAMMA performs protein-level alignment using nucleotide CDS sequences
        and reports mutations in the Codon_Changes column of the .gamma output.
        """
        if not self.genes_file or not Path(self.genes_file).exists():
            return

        print("Running GAMMA for resistance gene mutation analysis...")

        gamma_output_prefix = f"{self.output_prefix}_gamma"
        gamma_output_file = f"{gamma_output_prefix}.gamma"

        cmd = ['GAMMA.py', self.assembly, self.genes_file, gamma_output_prefix]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            if not Path(gamma_output_file).exists():
                print("WARNING: GAMMA output file not found")
                return

            with open(gamma_output_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    # Strip ambiguous match marker (‡) if present
                    gene_name = row['Gene'].rstrip('\u2021')
                    contig = row['Contig']
                    start = int(row['Start'])
                    stop = int(row['Stop'])
                    match_type = row['Match_Type']
                    codon_changes = row.get('Codon_Changes', '0')
                    codon_percent = float(row['Codon_Percent'])
                    percent_length = float(row['Percent_Length'])

                    identity = codon_percent * 100
                    coverage = percent_length * 100

                    raw_mutations = self._parse_gamma_codon_changes(codon_changes)
                    mutations = self._filter_known_mutations(raw_mutations, gene_name)

                    self.gamma_results.append({
                        'protein': gene_name,
                        'contig': contig,
                        'contig_start': min(start, stop),
                        'contig_end': max(start, stop),
                        'identity': identity,
                        'coverage': coverage,
                        'mutations': mutations,
                        'match_type': match_type,
                    })

            print(f"Found {len(self.gamma_results)} gene alignments")
            for r in self.gamma_results:
                if r['mutations']:
                    print(f"  {r['protein']}: {len(r['mutations'])} mutations detected")

        except FileNotFoundError:
            print("WARNING: GAMMA not found, skipping gene mutation analysis")
        except subprocess.CalledProcessError as e:
            print(f"WARNING: GAMMA failed: {e.stderr}", file=sys.stderr)

    @staticmethod
    def _parse_gamma_codon_changes(codon_changes):
        """Parse GAMMA's Codon_Changes field into a list of mutation strings.

        GAMMA reports non-degenerate codon differences as comma-separated
        substitution strings like 'D179Y,V240G' or '0' for wildtype matches.
        Only standard amino acid substitution strings (e.g. 'D179Y') are returned;
        indels, frameshifts, and truncation annotations are skipped.
        """
        if not codon_changes or codon_changes.strip() in ('0', '-', ''):
            return []
        mutations = []
        for mut in codon_changes.split(','):
            mut = mut.strip()
            if re.match(r'^[A-Z*]\d+[A-Z*]$', mut):
                mutations.append(mut)
        return mutations

    def _filter_known_mutations(self, raw_mutations, gene_name):
        """
        Filter raw mutations to only known resistance mutations from mutation_db.

        raw_mutations: list of strings from _parse_gamma_codon_changes, e.g. ['D369N', 'G463D']
        gene_name:     gene FASTA header (e.g. 'acrB_WP_002892069.1' or 'blaKPC-3')

        Returns a list of mutation name strings (from mutation_db) for confirmed
        resistance mutations only.  Unknown positional variants or non-resistance
        positions are dropped.
        """
        gene_norm = self._normalize_gene_name(gene_name)

        if gene_norm not in self.mutation_db:
            return []

        known = self.mutation_db[gene_norm]
        filtered = []

        for mut_str in raw_mutations:
            # Only process standard substitution strings like "D369N", "G463*"
            m = re.match(r'^([A-Z*])(\d+)([A-Z*])$', mut_str)
            if not m:
                continue  # skip frameshifts, indels, ? entries

            ref, pos, var = m.group(1), int(m.group(2)), m.group(3)

            if pos not in known:
                continue  # not a known resistance position

            known_info = known[pos]
            if var in known_info.get('variants', []):
                # Known resistance variant at a known position
                filtered.append(known_info.get('name', mut_str))

        return filtered

    @staticmethod
    def _normalize_gene_name(gene_name):
        """
        Normalize gene name for cross-method comparison.

        Handles both primer-style names ('uhpB', 'blaKPC') and GAMMA
        gene FASTA ID style ('acrB_WP_002892069.1', 'blaKPC-3').
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
        Merge GAMMA and seqkit mutation results with confidence scores.

        Confidence scoring:
          - 100% : mutation found by BOTH GAMMA AND seqkit
          -  50% : mutation found by only ONE method

        Returns a list of dicts (sorted by gene then mutation):
        {
            'gene'         : normalized gene name (str),
            'mutation'     : mutation string e.g. 'D179Y' (str),
            'confidence'   : 50 or 100 (int),
            'methods'      : list of method names, e.g. ['gamma', 'seqkit'],
            'gamma_detail' : matching GAMMA result dict, or None,
            'seqkit_detail': matching seqkit result dict, or None,
        }
        """
        # Index GAMMA results by normalized (gene, mutation)
        gamma_by_key = {}
        for r in self.gamma_results:
            gene_norm = self._normalize_gene_name(r['protein'])
            for mut in r.get('mutations', []):
                key = (gene_norm, mut)
                if key not in gamma_by_key:
                    gamma_by_key[key] = r

        # Index seqkit results by normalized (gene, mutation)
        seqkit_by_key = {}
        for r in seqkit_mut_results:
            gene_norm = self._normalize_gene_name(r['gene'])
            for mut in r.get('mutations', []):
                key = (gene_norm, mut)
                if key not in seqkit_by_key:
                    seqkit_by_key[key] = r

        all_keys = set(gamma_by_key.keys()) | set(seqkit_by_key.keys())

        unified = []
        for (gene, mutation) in sorted(all_keys):
            in_gamma = (gene, mutation) in gamma_by_key
            in_seqkit = (gene, mutation) in seqkit_by_key

            methods = []
            if in_gamma:
                methods.append('gamma')
            if in_seqkit:
                methods.append('seqkit')

            confidence = 100 if (in_gamma and in_seqkit) else 50

            unified.append({
                'gene': gene,
                'mutation': mutation,
                'confidence': confidence,
                'methods': methods,
                'gamma_detail': gamma_by_key.get((gene, mutation)),
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
                'GAMMA_Contig', 'GAMMA_Identity(%)', 'GAMMA_Coverage(%)',
                'SeqKit_PairID',
            ]) + '\n')

            for r in self.unified_results:
                gm = r['gamma_detail']
                sk = r['seqkit_detail']

                gm_contig = gm['contig'] if gm else '-'
                gm_identity = f"{gm['identity']:.2f}" if gm else '-'
                gm_coverage = f"{gm['coverage']:.2f}" if gm else '-'
                sk_pair = sk['pair_id'] if sk else '-'

                f.write('\t'.join([
                    r['gene'],
                    r['mutation'],
                    str(r['confidence']),
                    '+'.join(r['methods']),
                    gm_contig,
                    gm_identity,
                    gm_coverage,
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

    def write_gamma_report(self):
        """Write GAMMA mutation detection results."""
        if not self.gamma_results:
            return

        report_file = f"{self.output_prefix}_protein_mutations.tsv"
        print(f"Writing GAMMA mutation results to {report_file}...")

        with open(report_file, 'w') as f:
            f.write('\t'.join(['Gene', 'Contig', 'Start', 'End', 'Identity',
                              'Coverage', 'Match_Type', 'Mutations', 'Method']) + '\n')

            for r in self.gamma_results:
                mutations_str = ';'.join(r['mutations']) if r['mutations'] else '-'
                f.write('\t'.join([
                    r['protein'],
                    r['contig'],
                    str(r['contig_start']),
                    str(r['contig_end']),
                    f"{r['identity']:.2f}",
                    f"{r['coverage']:.2f}",
                    r['match_type'],
                    mutations_str,
                    'GAMMA'
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
          1. GAMMA    – protein-level gene alignment, Codon_Changes mutation parsing
          2. seqkit   – amplicon extraction + sequence-level mutation detection
          3. merge    – cross-reference both methods, assign confidence scores
          4. seqkit BED – amplicon location detection (cross-ref with BLAST)

        Returns:
            (gamma_results, amplicon_results, seqkit_mut_results, unified_results)
        """
        self.run_gamma()
        seqkit_mut_results = self.detect_seqkit_mutations()
        unified_results = self.merge_detection_results(seqkit_mut_results)

        self.detect_amplicons()
        self.analyze_amplicons(blast_results)

        self.write_gamma_report()
        self.write_amplicon_report()
        self.write_unified_report()

        return self.gamma_results, self.amplicon_results, seqkit_mut_results, unified_results


def run_mutation_detection(assembly, output, genes, primers, blast_results=None, mutation_db_file=None):
    detector = MutationDetector(assembly, output, genes, primers, mutation_db_file)
    return detector.run(blast_results)
