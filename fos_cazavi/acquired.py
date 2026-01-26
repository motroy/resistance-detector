import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from .utils import detect_mutations, check_dependencies, load_mutation_db, KNOWN_MUTATIONS

class BlastDetector:
    def __init__(self, assembly, database, output_prefix, min_identity=90, min_coverage=80, mutation_db_file=None):
        self.assembly = assembly
        self.database = database
        self.output_prefix = output_prefix
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.results = []
        self.detected_genes = []

        # Load mutations
        if mutation_db_file:
            self.mutation_db = load_mutation_db(mutation_db_file)
        else:
            # Try to auto-discover
            default_mut_file = Path(str(self.database).replace('.fasta', '') + '_mutations.tsv')
            if default_mut_file.exists():
                print(f"Auto-detected mutation file: {default_mut_file}")
                self.mutation_db = load_mutation_db(default_mut_file)
            else:
                self.mutation_db = KNOWN_MUTATIONS

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
        if not check_dependencies(['blastn', 'makeblastdb']):
            sys.exit(1)

        self.prepare_database()

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
        if '|' in subject_id:
            parts = subject_id.split('|')
            if len(parts) > 4:
                return parts[4]

        if '_' not in subject_id and '|' not in subject_id:
            return subject_id

        return subject_id.split('|')[-1].split('_')[0]

    def extract_hit_sequence(self, hit):
        """Extract the sequence of a BLAST hit from the assembly"""
        query_id = hit['query_id']
        qstart = hit['qstart']
        qend = hit['qend']

        for record in SeqIO.parse(self.assembly, 'fasta'):
            if record.id == query_id:
                if qstart < qend:
                    seq = record.seq[qstart-1:qend]
                else:
                    seq = record.seq[qend-1:qstart].reverse_complement()
                return str(seq)
        return None

    def analyze_hits(self, hits):
        """Analyze BLAST hits and detect mutations"""
        print("Analyzing hits and detecting mutations...")

        for hit in hits:
            sequence = self.extract_hit_sequence(hit)

            if sequence:
                mutations = detect_mutations(hit['gene'], sequence, self.mutation_db)

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

    def write_report(self):
        """Write results to TSV file"""
        report_file = f"{self.output_prefix}_results.tsv"

        print(f"Writing results to {report_file}...")

        with open(report_file, 'w') as f:
            f.write('\t'.join(['Contig', 'Gene', 'Identity%', 'Coverage%',
                              'Mutations']) + '\n')

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

    def run(self):
        hits = self.run_blast()
        if hits:
            self.analyze_hits(hits)
        self.write_report()
        self.write_sequences()
        return self.results

def run_acquired_detection(assembly, database, output, min_id, min_cov, mutation_db=None):
    detector = BlastDetector(assembly, database, output, min_id, min_cov, mutation_db)
    return detector.run()
