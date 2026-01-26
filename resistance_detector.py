#!/usr/bin/env python3
"""
FOS-CAZAVI Resistance Detector
Detects fosfomycin and ceftazidime-avibactam resistance genes and mutations
from WGS assemblies or reads.

Author: Motroy
Similar to Kleborate approach
"""

import argparse
import subprocess
import sys
import os
import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import tempfile


class ResistanceDetector:
    """Main class for detecting resistance genes and mutations"""
    
    # Default known resistance mutations (fallback)
    KNOWN_MUTATIONS = {
        # ... (Same as before, kept as fallback/initialization)
        # But for brevity in this update, I will keep them but allow overwrite.
    }
    
    def __init__(self, assembly, database, output_prefix, min_identity=90, min_coverage=80, mutation_db=None, primers_file=None, proteins_file=None):
        self.assembly = assembly
        self.database = database
        self.output_prefix = output_prefix
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.mutation_db = mutation_db
        self.primers_file = primers_file
        self.proteins_file = proteins_file
        self.results = []
        self.detected_genes = []
        self.primer_results = []
        self.amplicon_results = []
        self.miniprot_results = []
        self.primers = {}
        
        # Initialize default mutations (fallback)
        self.init_default_mutations()
        
        # Load external mutations if provided
        if self.mutation_db:
            self.load_mutation_db(self.mutation_db)
        else:
            # Try to auto-discover
            default_mut_file = Path(str(self.database).replace('.fasta', '') + '_mutations.tsv')
            if default_mut_file.exists():
                print(f"Auto-detected mutation file: {default_mut_file}")
                self.load_mutation_db(default_mut_file)

        # Load primers if provided
        if self.primers_file:
            self.load_primers(self.primers_file)
        else:
            pass
    
    def load_primers(self, primers_file):
        """Load primers from TSV file"""
        if not primers_file or not Path(primers_file).exists():
            return

        print(f"Loading primers from {primers_file}...")

        try:
            self.primers = {}
            with open(primers_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    # Support different column naming conventions
                    name = row.get('Primer', row.get('name'))
                    seq = row.get('Nucleotide_sequence', row.get('seq', row.get('sequence')))
                    purpose = row.get('Purpose', row.get('purpose', ''))
                    mutation = row.get('Mutation', row.get('mutation', ''))
                    gene = row.get('Gene', row.get('gene', ''))
                    pair_id = row.get('Pair_ID', row.get('pair_id', ''))

                    if name and seq:
                        self.primers[name] = {
                            'seq': seq.strip(),
                            'purpose': purpose.strip(),
                            'mutation': mutation.strip() if mutation and mutation.strip() != '-' else None,
                            'gene': gene.strip() if gene and gene.strip() != '-' else None,
                            'pair_id': pair_id.strip() if pair_id and pair_id.strip() != '-' else None
                        }

            print(f"Loaded {len(self.primers)} primers")

        except Exception as e:
            print(f"ERROR loading primers file: {e}", file=sys.stderr)

    def init_default_mutations(self):
        """Initialize hardcoded mutations as fallback"""
        self.KNOWN_MUTATIONS = {
            'blaKPC': {
                179: {'ref': 'D', 'variants': ['Y', 'N'], 'name': 'D179Y/N'},
                240: {'ref': 'V', 'variants': ['G'], 'name': 'V240G'},
                243: {'ref': 'T', 'variants': ['M'], 'name': 'T243M'}
            },
            'blaOXA-48': {
                68: {'ref': 'P', 'variants': ['A'], 'name': 'P68A'},
                211: {'ref': 'Y', 'variants': ['S'], 'name': 'Y211S'}
            },
            'murA': {
                369: {'ref': 'D', 'variants': ['N'], 'name': 'D369N'},
                370: {'ref': 'L', 'variants': ['I'], 'name': 'L370I'}
            },
            'uhpT': {
                55: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G55D/*'},
                198: {'ref': 'W', 'variants': ['*', 'R'], 'name': 'W198*/R'},
                258: {'ref': 'E', 'variants': ['*', 'K'], 'name': 'E258*/K'},
                350: {'ref': 'W', 'variants': ['*', 'R'], 'name': 'W350*/R'}
            },
            'glpT': {
                44: {'ref': 'E', 'variants': ['*', 'K'], 'name': 'E44*/K'},
                88: {'ref': 'W', 'variants': ['*', 'R', 'G'], 'name': 'W88*/R/G'},
                90: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G90D/*'},
                234: {'ref': 'W', 'variants': ['*', 'R'], 'name': 'W234*/R'},
                362: {'ref': 'R', 'variants': ['C', 'H', '*'], 'name': 'R362C/H/*'}
            },
            'uhpA': {
                54: {'ref': 'D', 'variants': ['N', 'A'], 'name': 'D54N/A'},
                139: {'ref': 'R', 'variants': ['C', 'H'], 'name': 'R139C/H'}
            },
            'uhpB': {
                469: {'ref': 'G', 'variants': ['R'], 'name': 'G469R'},
                350: {'ref': 'H', 'variants': ['Y', 'Q'], 'name': 'H350Y/Q'}
            },
            'uhpC': {
                384: {'ref': 'F', 'variants': ['L'], 'name': 'F384L'}
            },
            'cyaA': {
                463: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G463D/*'}
            },
            'ptsI': {
                191: {'ref': 'H', 'variants': ['Y', 'Q'], 'name': 'H191Y/Q'}
            },
            'galU': {
                282: {'ref': 'R', 'variants': ['V'], 'name': 'R282V'}
            },
            'lon': {
                558: {'ref': 'Q', 'variants': ['*'], 'name': 'Q558*'}
            },
            'fosAKP': {
                91: {'ref': 'I', 'variants': ['V'], 'name': 'I91V'}
            },
            'fosA': {
                90: {'ref': 'K', 'variants': ['E', 'Q'], 'name': 'K90E/Q'},
                119: {'ref': 'H', 'variants': ['Q', 'R'], 'name': 'H119Q/R'}
            },
            'acrB': {
                617: {'ref': 'G', 'variants': ['D', 'N'], 'name': 'G617D/N'},
                626: {'ref': 'F', 'variants': ['L'], 'name': 'F626L'},
                628: {'ref': 'A', 'variants': ['T', 'V'], 'name': 'A628T/V'}
            },
            'ompK36': {
                134: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G134D/*'},
                135: {'ref': 'D', 'variants': ['*', 'N'], 'name': 'D135*/N'},
                213: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G213D/*'}
            },
            'ftsI': {
                333: {'ref': 'A', 'variants': ['V', 'T'], 'name': 'A333V/T'},
                350: {'ref': 'Y', 'variants': ['C', 'S'], 'name': 'Y350C/S'},
                357: {'ref': 'S', 'variants': ['N'], 'name': 'S357N'}
            },
            'envZ': {
                244: {'ref': 'G', 'variants': ['S', 'D'], 'name': 'G244S/D'},
                324: {'ref': 'T', 'variants': ['I', 'A'], 'name': 'T324I/A'}
            }
        }

    def load_mutation_db(self, mutation_file):
        """Load mutations from TSV file, replacing defaults"""
        if not mutation_file or not Path(mutation_file).exists():
            return
            
        print(f"Loading mutations from {mutation_file}...")
        
        # Clear defaults to avoid mixing incompatible coordinate systems
        self.KNOWN_MUTATIONS = defaultdict(dict)
        
        try:
            with open(mutation_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    gene = row['Gene']
                    # Handle base gene name logic
                    # The file should contain base names, but let's be safe
                    
                    pos = int(row['Position'])
                    ref = row['Ref']
                    var = row['Variant']
                    
                    # Handle variants (could be 'Y/N' or just 'Y')
                    vars_list = var.split('/')
                    
                    if gene not in self.KNOWN_MUTATIONS:
                        self.KNOWN_MUTATIONS[gene] = {}
                    
                    if pos not in self.KNOWN_MUTATIONS[gene]:
                        self.KNOWN_MUTATIONS[gene][pos] = {
                            'ref': ref,
                            'variants': set()
                        }
                    
                    # Warn on conflict but proceed
                    if self.KNOWN_MUTATIONS[gene][pos]['ref'] != ref:
                        # This might happen if multiple entries define same pos with different ref? 
                        # Unlikely in valid DB, but possible if mixing sources.
                        pass
                    
                    self.KNOWN_MUTATIONS[gene][pos]['variants'].update(vars_list)
            
            # Post-process to format names
            for gene in self.KNOWN_MUTATIONS:
                for pos in self.KNOWN_MUTATIONS[gene]:
                    entry = self.KNOWN_MUTATIONS[gene][pos]
                    variants_str = '/'.join(sorted(entry['variants']))
                    entry['name'] = f"{entry['ref']}{pos}{variants_str}"
                    entry['variants'] = list(entry['variants'])
                    
            print(f"Loaded mutation definitions for {len(self.KNOWN_MUTATIONS)} genes")
            
        except Exception as e:
            print(f"ERROR loading mutation file: {e}", file=sys.stderr)

    def check_dependencies(self):
        """Check if required tools are installed"""
        required = ['blastn', 'makeblastdb']
        missing = []
        
        for tool in required:
            try:
                subprocess.run([tool, '-version'], 
                             capture_output=True, 
                             check=True)
            except (subprocess.CalledProcessError, FileNotFoundError):
                missing.append(tool)
        
        if missing:
            print(f"ERROR: Missing required tools: {', '.join(missing)}", 
                  file=sys.stderr)
            print("Please install BLAST+ toolkit", file=sys.stderr)
            sys.exit(1)
    
    def prepare_database(self):
        """Check if BLAST database exists, create if needed"""
        db_files = [f"{self.database}.{ext}" for ext in ['nhr', 'nin', 'nsq']]
        
        if not all(Path(f).exists() for f in db_files):
            print(f"Creating BLAST database from {self.database}...")
            cmd = ['makeblastdb', '-in', self.database, 
                   '-dbtype', 'nucl', '-parse_seqids']
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                print("Database created successfully")
            except subprocess.CalledProcessError as e:
                print(f"ERROR creating database: {e.stderr.decode()}", 
                      file=sys.stderr)
                sys.exit(1)
    
    def run_blast(self):
        """Run BLAST search for resistance genes"""
        print(f"Running BLAST search (min_id={self.min_identity}%, min_cov={self.min_coverage}%)...")
        
        blast_output = f"{self.output_prefix}_blast.txt"
        
        cmd = [
            'blastn',
            '-query', self.assembly,
            '-db', self.database,
            '-outfmt', '6 qseqid sseqid pident length qstart qend qlen sstart send slen',
            '-evalue', '1e-20',
            '-max_target_seqs', '5'
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, 
                                  text=True, check=True)
            
            with open(blast_output, 'w') as f:
                f.write(result.stdout)
            
            return self.parse_blast_output(result.stdout)
            
        except subprocess.CalledProcessError as e:
            print(f"ERROR running BLAST: {e.stderr}", file=sys.stderr)
            sys.exit(1)
    
    def parse_blast_output(self, blast_output):
        """Parse BLAST results and filter by identity and coverage"""
        hits = []
        
        for line in blast_output.strip().split('\n'):
            if not line:
                continue
            
            fields = line.split('\t')
            query_id = fields[0]
            subject_id = fields[1]
            pident = float(fields[2])
            length = int(fields[3])
            qstart = int(fields[4])
            qend = int(fields[5])
            qlen = int(fields[6])
            sstart = int(fields[7])
            send = int(fields[8])
            slen = int(fields[9])
            
            # Calculate coverage based on subject (reference gene)
            coverage = (length / slen) * 100
            
            if pident >= self.min_identity and coverage >= self.min_coverage:
                hits.append({
                    'query_id': query_id,
                    'subject_id': subject_id,
                    'gene': self.extract_gene_name(subject_id),
                    'identity': pident,
                    'coverage': coverage,
                    'qstart': qstart,
                    'qend': qend,
                    'qlen': qlen,
                    'sstart': sstart,
                    'send': send,
                    'slen': slen
                })
        
        print(f"Found {len(hits)} gene hits passing thresholds")
        return hits
    
    def extract_gene_name(self, subject_id):
        """Extract gene name from BLAST subject ID"""
        # Handle various formats: >fosA3_reference, >gb|XXX|fosA3, etc.
        # Updated to handle formats like >fosA3 ... or >WP_...|...|fosA3|...
        # If headers are from our new create script, they might look like:
        # >murA reference sequence (AMRfinderPlus:...)
        # or >fosA3 reference sequence (AMRfinderPlus)
        
        # Strategy: 
        # 1. Try splitting by pipe if present
        if '|' in subject_id:
            parts = subject_id.split('|')
            # Look for gene name in parts?
            # AMR_CDS.fa format: >Prot|Nuc|1|1|Gene|...
            if len(parts) > 4:
                return parts[4]
        
        # 2. Try splitting by space (SeqIO often keeps just ID part)
        # 3. Try splitting by underscore if standard format
        
        # If ID is just "murA"
        if '_' not in subject_id and '|' not in subject_id:
            return subject_id
            
        # Fallback to taking first part before _ or space
        # But some genes have underscores e.g. fosA_3? No usually fosA3.
        # But subject_id might be "fosA3_reference".
        return subject_id.split('|')[-1].split('_')[0]
    
    def extract_hit_sequence(self, hit):
        """Extract the sequence of a BLAST hit from the assembly"""
        query_id = hit['query_id']
        qstart = hit['qstart']
        qend = hit['qend']
        
        # Load assembly and find the contig
        for record in SeqIO.parse(self.assembly, 'fasta'):
            if record.id == query_id:
                # Extract sequence (accounting for 1-based coordinates)
                if qstart < qend:
                    seq = record.seq[qstart-1:qend]
                else:
                    # Reverse strand
                    seq = record.seq[qend-1:qstart].reverse_complement()
                
                return str(seq)
        
        return None
    
    def detect_mutations(self, gene_name, sequence):
        """Detect known mutations in a gene sequence"""
        mutations_found = []
        
        # Get base gene name (remove variant numbers)
        base_gene = gene_name.split('-')[0]
        if base_gene.startswith('bla'):
            base_gene = 'bla' + base_gene[3:].rstrip('0123456789')
        
        # Try exact match first
        if gene_name in self.KNOWN_MUTATIONS:
            mutation_dict = self.KNOWN_MUTATIONS[gene_name]
        elif base_gene in self.KNOWN_MUTATIONS:
            mutation_dict = self.KNOWN_MUTATIONS[base_gene]
        else:
            return mutations_found
        
        # Translate to protein
        try:
            # Ensure sequence length is multiple of 3
            if len(sequence) % 3 != 0:
                sequence = sequence[:-(len(sequence) % 3)]
            
            protein = str(Seq(sequence).translate())
            
            # Check each known mutation position
            for pos, mut_info in mutation_dict.items():
                if pos <= len(protein):
                    observed_aa = protein[pos-1]
                    ref_aa = mut_info['ref']
                    
                    if observed_aa in mut_info['variants']:
                        mutations_found.append(mut_info['name'])
                    elif observed_aa != ref_aa:
                        # Novel mutation at known position
                        mutations_found.append(f"{ref_aa}{pos}{observed_aa}")
        
        except Exception as e:
            print(f"Warning: Could not translate {gene_name}: {e}", 
                  file=sys.stderr)
        
        return mutations_found
    
    def detect_amplicons(self):
        """Detect amplicons using seqkit amplicon"""
        if not self.primers:
            return

        print("Detecting amplicons with seqkit...")

        # 1. Generate seqkit-compatible primer file
        # Format: Name <tab> Fwd <tab> Rev
        # We need to group primers by Pair_ID from our primers dict

        # Group primers by Pair_ID
        pairs = defaultdict(dict)
        for name, info in self.primers.items():
            pair_id = info.get('pair_id')
            if not pair_id or pair_id == '-':
                continue

            # Identify if F or R based on name suffix or explicit logic
            # Using simple heuristic: ends with _F or _R, or similar
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
                    # Seqkit expects: Name Fwd Rev
                    f.write(f"{pair_id}\t{p['F']}\t{p['R']}\n")

        # 2. Run seqkit amplicon
        amplicon_output = f"{self.output_prefix}_amplicons.fasta"
        # We use --bed to get coordinates easily, or just Parse FASTA headers if they contain info?
        # seqkit amplicon outputs FASTA.
        # But we need coordinates to check for internal mutations (BLAST hits).
        # seqkit amplicon --bed outputs BED format.

        cmd = [
            'seqkit', 'amplicon',
            '-p', seqkit_primer_file,
            self.assembly,
            '--bed'
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            # Parse BED output
            # BED columns: chrom, start(0-based), end(0-based), name, score, strand, amplicon_seq?
            # Wait, seqkit amplicon --bed:
            # col 1: chrom
            # col 2: start (0-based)
            # col 3: end (0-based)
            # col 4: primer name (pair_id)
            # col 5: score
            # col 6: strand
            # col 7: amplicon sequence (if --bed is used with recent versions? checks docs: "output in BED6+1 format with amplicon as the 7th column")

            for line in result.stdout.strip().split('\n'):
                if not line: continue
                fields = line.split('\t')
                if len(fields) < 6: continue

                contig = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                pair_id = fields[3]
                strand = fields[5]
                # seq = fields[6] if len(fields) > 6 else ""

                self.amplicon_results.append({
                    'pair_id': pair_id,
                    'contig': contig,
                    'start': start, # 0-based
                    'end': end,     # 0-based
                    'length': end - start,
                    'mutations_found': []
                })

            print(f"Found {len(self.amplicon_results)} amplicons")

        except subprocess.CalledProcessError as e:
            print(f"ERROR running seqkit: {e.stderr}", file=sys.stderr)
        except Exception as e:
            print(f"ERROR processing amplicons: {e}", file=sys.stderr)

    def write_amplicon_report(self):
        """Write amplicon detection results to TSV file"""
        report_file = f"{self.output_prefix}_amplicons.tsv"

        print(f"Writing amplicon results to {report_file}...")

        with open(report_file, 'w') as f:
            # Header
            f.write('\t'.join(['Pair_ID', 'Contig', 'Start', 'End', 'Length',
                              'Mutations_Found']) + '\n')

            # Results
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

    def run_miniprot(self):
        """
        Run miniprot to detect mutations in CAZAVI resistance proteins.
        Uses protein-to-genome alignment with cs tag parsing for mutation detection.
        Reference: doi.org/10.3389/fcimb.2025.1645042
        """
        if not self.proteins_file or not Path(self.proteins_file).exists():
            return

        print(f"Running miniprot for protein mutation analysis...")

        # Run miniprot with cs tag output
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

            # Write raw output
            with open(miniprot_output, 'w') as f:
                f.write(result.stdout)

            # Parse PAF output with cs tags
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

                # Calculate identity
                identity = (matches / alignment_len * 100) if alignment_len > 0 else 0
                coverage = ((protein_end - protein_start) / protein_len * 100) if protein_len > 0 else 0

                # Parse cs tag for mutations
                mutations = []
                cs_tag = None
                for field in fields[12:]:
                    if field.startswith('cs:Z:'):
                        cs_tag = field[5:]
                        break

                if cs_tag:
                    mutations = self.parse_miniprot_cs(cs_tag, protein_name)

                # Store result
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

            # Summarize mutations found
            for r in self.miniprot_results:
                if r['mutations']:
                    print(f"  {r['protein']}: {len(r['mutations'])} mutations detected")

        except FileNotFoundError:
            print("WARNING: miniprot not found, skipping protein analysis")
        except subprocess.CalledProcessError as e:
            print(f"WARNING: miniprot failed: {e.stderr}", file=sys.stderr)

    def parse_miniprot_cs(self, cs_tag, protein_name):
        """
        Parse miniprot cs tag to extract mutations.

        cs tag format:
        - ":[0-9]+" represents identical amino acids
        - "[acgtn]+[A-Z*]" represents substitution (ref codons to query aa)
        - "+[A-Z]+" represents insertion to reference
        - "-[acgtn]+" represents deletion from reference
        - "~[acgtn]+[0-9]+[acgtn]+" represents intron
        - "*[A-Z][A-Z]" represents amino acid substitution
        """
        import re
        mutations = []
        position = 1  # 1-indexed amino acid position

        # Split cs tag into operations
        # Pattern matches: :num, *XY (aa sub), +AA (ins), -nt (del), ~intron
        pattern = r'(:[0-9]+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[a-z]+|~[a-z]+[0-9]+[a-z]+|[a-z]+[A-Z\*])'

        for match in re.finditer(pattern, cs_tag):
            op = match.group(1)

            if op.startswith(':'):
                # Identical residues
                num_identical = int(op[1:])
                position += num_identical

            elif op.startswith('*'):
                # Amino acid substitution *XY means ref X -> query Y
                ref_aa = op[1].upper()
                query_aa = op[2].upper()
                mutations.append(f"{ref_aa}{position}{query_aa}")
                position += 1

            elif op.startswith('+'):
                # Insertion in query
                inserted = op[1:]
                mutations.append(f"ins{position}_{inserted}")
                # Insertions don't advance reference position

            elif op.startswith('-'):
                # Deletion from reference
                deleted_nt = op[1:]
                deleted_aa_count = len(deleted_nt) // 3
                if deleted_aa_count > 0:
                    mutations.append(f"del{position}_{deleted_aa_count}aa")
                position += deleted_aa_count

            elif op.startswith('~'):
                # Intron - skip
                pass

            elif op[:-1].islower() and op[-1].isupper():
                # Codon substitution: ref codons -> query aa
                # e.g., "gat" + "Y" means GAT -> TAT (D -> Y)
                ref_codons = op[:-1]
                query_aa = op[-1]
                # This represents a substitution at this position
                mutations.append(f"?{position}{query_aa}")
                position += 1

            elif op[-1] == '*':
                # Stop codon
                mutations.append(f"?{position}*")
                position += 1

        return mutations

    def write_miniprot_report(self):
        """Write miniprot mutation detection results"""
        if not self.miniprot_results:
            return

        report_file = f"{self.output_prefix}_protein_mutations.tsv"
        print(f"Writing protein mutation results to {report_file}...")

        with open(report_file, 'w') as f:
            # Header
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

    def analyze_hits(self, hits):
        """Analyze BLAST hits and detect mutations"""
        print("Analyzing hits and detecting mutations...")
        
        for hit in hits:
            sequence = self.extract_hit_sequence(hit)
            
            if sequence:
                mutations = self.detect_mutations(hit['gene'], sequence)
                
                result = {
                    'contig': hit['query_id'],
                    'gene': hit['gene'],
                    'identity': f"{hit['identity']:.2f}",
                    'coverage': f"{hit['coverage']:.2f}",
                    'mutations': ','.join(mutations) if mutations else '-',
                    'sequence': sequence,
                    'start': hit['qstart'],
                    'end': hit['qend']
                }
                
                self.results.append(result)
                self.detected_genes.append(SeqRecord(
                    Seq(sequence),
                    id=f"{hit['query_id']}_{hit['gene']}",
                    description=f"identity={hit['identity']:.2f}% coverage={hit['coverage']:.2f}% mutations={result['mutations']}"
                ))

        # Analyze amplicons for overlaps with detected genes
        self.analyze_amplicons()

    def analyze_amplicons(self):
        """Check if detected genes fall within amplicons"""
        if not self.amplicon_results or not self.results:
            return

        print("Checking for mutations within amplicons...")

        for amp in self.amplicon_results:
            amp_contig = amp['contig']
            amp_start = amp['start']
            amp_end = amp['end']

            for res in self.results:
                if res['contig'] != amp_contig:
                    continue

                # Convert BLAST 1-based coords to 0-based range
                res_start = res['start']
                res_end = res['end']

                start = min(res_start, res_end) - 1
                end = max(res_start, res_end)

                # Check overlap
                # Overlap if max(starts) < min(ends)
                if max(amp_start, start) < min(amp_end, end):
                    # Overlap found
                    # Report gene and mutations
                    mut_str = res['mutations']
                    if mut_str != '-':
                        amp['mutations_found'].append(f"{res['gene']}: {mut_str}")
                    else:
                        amp['mutations_found'].append(f"{res['gene']}: (wildtype)")

    def write_report(self):
        """Write results to TSV file"""
        report_file = f"{self.output_prefix}_results.tsv"
        
        print(f"Writing results to {report_file}...")
        
        with open(report_file, 'w') as f:
            # Header
            f.write('\t'.join(['Contig', 'Gene', 'Identity%', 'Coverage%', 
                              'Mutations']) + '\n')
            
            # Results
            for result in self.results:
                f.write('\t'.join([
                    result['contig'],
                    result['gene'],
                    result['identity'],
                    result['coverage'],
                    result['mutations']
                ]) + '\n')
        
        print(f"Detected {len(self.results)} resistance genes")
    
    def write_sequences(self):
        """Write detected gene sequences to FASTA file"""
        fasta_file = f"{self.output_prefix}_genes.fasta"
        
        if self.detected_genes:
            print(f"Writing {len(self.detected_genes)} gene sequences to {fasta_file}...")
            SeqIO.write(self.detected_genes, fasta_file, 'fasta')
        else:
            print("No genes detected to write")
    
    def write_summary(self):
        """Write a summary of detected resistance mechanisms"""
        summary_file = f"{self.output_prefix}_summary.txt"
        
        print(f"Writing summary to {summary_file}...")
        
        # Group by drug class
        fos_genes = [r for r in self.results if r['gene'].startswith('fos')]
        kpc_genes = [r for r in self.results if 'KPC' in r['gene'].upper()]
        oxa_genes = [r for r in self.results if 'OXA' in r['gene'].upper()]
        other_genes = [r for r in self.results 
                      if r not in fos_genes + kpc_genes + oxa_genes]
        
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + '\n')
            f.write("FOS-CAZAVI Resistance Detection Summary\n")
            f.write("=" * 70 + '\n\n')
            
            f.write(f"Assembly: {self.assembly}\n")
            f.write(f"Total genes detected: {len(self.results)}\n\n")
            
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
            
            if not self.results:
                f.write("No resistance genes detected\n")

            if self.amplicon_results:
                f.write("\n")
                f.write("DETECTED AMPLICONS:\n")
                f.write("-" * 50 + '\n')
                for amp in self.amplicon_results:
                    f.write(f"  Pair: {amp['pair_id']}\n")
                    f.write(f"    Location: {amp['contig']}:{amp['start']}-{amp['end']} ({amp['length']} bp)\n")
                    f.write(f"    Primers: {amp['f_primer']} -> {amp['r_primer']}\n")
                    if amp['mutations_found']:
                        f.write("    Mutations/Genes in amplicon:\n")
                        for m in amp['mutations_found']:
                            f.write(f"      - {m}\n")
                    else:
                        f.write("    No resistance genes/mutations detected in amplicon\n")

            if self.miniprot_results:
                f.write("\n")
                f.write("PROTEIN MUTATION ANALYSIS (miniprot):\n")
                f.write("-" * 50 + '\n')
                f.write("Reference: doi.org/10.3389/fcimb.2025.1645042\n\n")

                # Group by protein type
                porins = [r for r in self.miniprot_results if 'Omp' in r['protein']]
                efflux = [r for r in self.miniprot_results if 'acr' in r['protein'].lower()]

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

    def run(self):
        """Run the complete detection pipeline"""
        print("=" * 70)
        print("FOS-CAZAVI Resistance Detector")
        print("=" * 70)
        
        self.check_dependencies()
        self.prepare_database()
        
        # Run BLAST first to get acquired genes
        hits = self.run_blast()
        if hits:
            self.analyze_hits(hits)

        # Then run amplicon detection
        self.detect_amplicons()
        if self.amplicon_results:
            self.analyze_amplicons()

        # Run miniprot for protein mutation analysis
        self.run_miniprot()

        # Write reports
        self.write_report()
        self.write_sequences()
        self.write_summary()
        if self.amplicon_results:
            self.write_amplicon_report()
        if self.miniprot_results:
            self.write_miniprot_report()

        print("=" * 70)
        print("Analysis complete!")
        print(f"Results written to {self.output_prefix}_*")
        print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description='Detect fosfomycin and ceftazidime-avibactam resistance genes and mutations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -a assembly.fasta -d resistance_db.fasta -o sample1
  %(prog)s -a assembly.fasta -d resistance_db.fasta -o sample1 --min_id 95 --min_cov 90
  %(prog)s -a assembly.fasta -d resistance_db.fasta -o sample1 --mutations resistance_db_mutations.tsv

Output files:
  <prefix>_results.tsv     - Tab-delimited results
  <prefix>_genes.fasta     - FASTA file with detected gene sequences
  <prefix>_summary.txt     - Human-readable summary
  <prefix>_blast.txt       - Raw BLAST output

Author: Motroy
        """
    )
    
    parser.add_argument('-a', '--assembly', required=True,
                       help='Input assembly file (FASTA format)')
    parser.add_argument('-d', '--database', required=True,
                       help='Resistance gene database (FASTA format)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output prefix for result files')
    parser.add_argument('--mutations',
                       help='Mutation definitions file (TSV)')
    parser.add_argument('--primers',
                       help='Primers definitions file (TSV)')
    parser.add_argument('--proteins',
                       help='Protein sequences for miniprot mutation analysis (FASTA)')
    parser.add_argument('--min_id', type=float, default=90.0,
                       help='Minimum percent identity (default: 90)')
    parser.add_argument('--min_cov', type=float, default=80.0,
                       help='Minimum percent coverage (default: 80)')
    parser.add_argument('-v', '--version', action='version',
                       version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.assembly).exists():
        print(f"ERROR: Assembly file not found: {args.assembly}", 
              file=sys.stderr)
        sys.exit(1)
    
    if not Path(args.database).exists():
        print(f"ERROR: Database file not found: {args.database}", 
              file=sys.stderr)
        sys.exit(1)
    
    # Run detector
    detector = ResistanceDetector(
        args.assembly,
        args.database,
        args.output,
        args.min_id,
        args.min_cov,
        args.mutations,
        args.primers,
        args.proteins
    )
    detector.run()


if __name__ == '__main__':
    main()
