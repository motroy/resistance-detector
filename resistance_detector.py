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
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import tempfile


class ResistanceDetector:
    """Main class for detecting resistance genes and mutations"""
    
    # Known resistance mutations
    # Based on literature: fosfomycin resistance mutations and carbapenem resistance mutations
    KNOWN_MUTATIONS = {
        # Carbapenem resistance - KPC mutations (ceftazidime-avibactam resistance)
        'blaKPC': {
            179: {'ref': 'D', 'variants': ['Y', 'N'], 'name': 'D179Y/N'},
            240: {'ref': 'V', 'variants': ['G'], 'name': 'V240G'},
            243: {'ref': 'T', 'variants': ['M'], 'name': 'T243M'}
        },
        # OXA-48 mutations
        'blaOXA-48': {
            68: {'ref': 'P', 'variants': ['A'], 'name': 'P68A'},
            211: {'ref': 'Y', 'variants': ['S'], 'name': 'Y211S'}
        },
        # Fosfomycin target enzyme - MurA mutations (D369N, L370I confer resistance)
        'murA': {
            369: {'ref': 'D', 'variants': ['N'], 'name': 'D369N'},
            370: {'ref': 'L', 'variants': ['I'], 'name': 'L370I'}
        },
        # Fosfomycin transporter - UhpT mutations (loss of function confers resistance)
        'uhpT': {
            # Key positions where mutations affect transporter function
            55: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G55D/*'},
            198: {'ref': 'W', 'variants': ['*', 'R'], 'name': 'W198*/R'},
            258: {'ref': 'E', 'variants': ['*', 'K'], 'name': 'E258*/K'},
            350: {'ref': 'W', 'variants': ['*', 'R'], 'name': 'W350*/R'}
        },
        # Fosfomycin transporter - GlpT mutations (loss of function confers resistance)
        'glpT': {
            # Key positions where mutations affect transporter function
            44: {'ref': 'E', 'variants': ['*', 'K'], 'name': 'E44*/K'},
            88: {'ref': 'W', 'variants': ['*', 'R', 'G'], 'name': 'W88*/R/G'},
            90: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G90D/*'},
            234: {'ref': 'W', 'variants': ['*', 'R'], 'name': 'W234*/R'},
            362: {'ref': 'R', 'variants': ['C', 'H', '*'], 'name': 'R362C/H/*'}
        },
        # UhpA mutations (transcriptional activator)
        'uhpA': {
            54: {'ref': 'D', 'variants': ['N', 'A'], 'name': 'D54N/A'},
            139: {'ref': 'R', 'variants': ['C', 'H'], 'name': 'R139C/H'}
        },
        # UhpB mutations (sensor kinase)
        'uhpB': {
            469: {'ref': 'G', 'variants': ['R'], 'name': 'G469R'},
            350: {'ref': 'H', 'variants': ['Y', 'Q'], 'name': 'H350Y/Q'}
        },
        # UhpC mutations (membrane sensor)
        'uhpC': {
            384: {'ref': 'F', 'variants': ['L'], 'name': 'F384L'}
        },
        # CyaA mutations (adenylate cyclase - affects uhpT expression)
        'cyaA': {
            463: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G463D/*'}
        },
        # PtsI mutations (phosphoenolpyruvate:sugar phosphotransferase)
        'ptsI': {
            191: {'ref': 'H', 'variants': ['Y', 'Q'], 'name': 'H191Y/Q'}
        },
        # Other fosfomycin resistance genes
        'galU': {
            282: {'ref': 'R', 'variants': ['V'], 'name': 'R282V'}
        },
        'lon': {
            558: {'ref': 'Q', 'variants': ['*'], 'name': 'Q558*'}
        },
        # FosA variants mutations
        'fosAKP': {
            91: {'ref': 'I', 'variants': ['V'], 'name': 'I91V'}
        },
        'fosA': {
            90: {'ref': 'K', 'variants': ['E', 'Q'], 'name': 'K90E/Q'},
            119: {'ref': 'H', 'variants': ['Q', 'R'], 'name': 'H119Q/R'}
        }
    }
    
    def __init__(self, assembly, database, output_prefix, min_identity=90, min_coverage=80):
        self.assembly = assembly
        self.database = database
        self.output_prefix = output_prefix
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.results = []
        self.detected_genes = []
        
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
        gene = subject_id.split('|')[-1].split('_')[0]
        return gene
    
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
        
        if base_gene not in self.KNOWN_MUTATIONS:
            return mutations_found
        
        # Translate to protein
        try:
            # Ensure sequence length is multiple of 3
            if len(sequence) % 3 != 0:
                sequence = sequence[:-(len(sequence) % 3)]
            
            protein = str(Seq(sequence).translate())
            
            # Check each known mutation position
            for pos, mut_info in self.KNOWN_MUTATIONS[base_gene].items():
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
                    'sequence': sequence
                }
                
                self.results.append(result)
                self.detected_genes.append(SeqRecord(
                    Seq(sequence),
                    id=f"{hit['query_id']}_{hit['gene']}",
                    description=f"identity={hit['identity']:.2f}% coverage={hit['coverage']:.2f}% mutations={result['mutations']}"
                ))
    
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
    
    def run(self):
        """Run the complete detection pipeline"""
        print("=" * 70)
        print("FOS-CAZAVI Resistance Detector")
        print("=" * 70)
        
        self.check_dependencies()
        self.prepare_database()
        hits = self.run_blast()
        
        if hits:
            self.analyze_hits(hits)
            self.write_report()
            self.write_sequences()
            self.write_summary()
        else:
            print("No resistance genes detected")
            # Still create empty output files
            with open(f"{self.output_prefix}_results.tsv", 'w') as f:
                f.write('\t'.join(['Contig', 'Gene', 'Identity%', 
                                  'Coverage%', 'Mutations']) + '\n')
            with open(f"{self.output_prefix}_summary.txt", 'w') as f:
                f.write("No resistance genes detected\n")
        
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
        args.min_cov
    )
    detector.run()


if __name__ == '__main__':
    main()
