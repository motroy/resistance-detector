import sys
import subprocess
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from .utils import detect_mutations_aligned, check_dependencies, load_mutation_db, KNOWN_MUTATIONS

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
            '-outfmt', '6 qseqid sseqid pident length qstart qend qlen sstart send slen qseq sseq',
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
            qseq = fields[10] if len(fields) > 10 else ''
            sseq = fields[11] if len(fields) > 11 else ''

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
                    'slen': slen,
                    'qseq': qseq,
                    'sseq': sseq
                })

        print(f"Found {len(hits)} gene hits passing thresholds")
        return self.filter_redundant_hits(hits)

    def filter_redundant_hits(self, hits):
        """Filter hits to keep only the best match per genomic location"""
        if not hits:
            return []

        # Sort by coverage (desc), identity (desc), then length (desc)
        # We use alignment length on query as a proxy for match quality if others are tied
        hits.sort(key=lambda x: (x['coverage'], x['identity'], abs(x['qend'] - x['qstart'])), reverse=True)

        kept_hits = []
        for hit in hits:
            # Normalize coordinates
            start = min(hit['qstart'], hit['qend'])
            end = max(hit['qstart'], hit['qend'])

            # Check overlap with kept hits
            is_redundant = False
            for kept in kept_hits:
                if hit['query_id'] != kept['query_id']:
                    continue

                k_start = min(kept['qstart'], kept['qend'])
                k_end = max(kept['qstart'], kept['qend'])

                # Calculate overlap
                overlap_start = max(start, k_start)
                overlap_end = min(end, k_end)
                overlap_len = max(0, overlap_end - overlap_start + 1)

                if overlap_len > 0:
                    # Calculate overlap percentage relative to the NEW hit
                    # If the new hit significantly overlaps with an existing better hit, discard it.
                    hit_len = end - start + 1
                    overlap_pct = (overlap_len / hit_len) * 100
                    if overlap_pct > 50:
                        is_redundant = True
                        break

            if not is_redundant:
                kept_hits.append(hit)

        print(f"Filtered to {len(kept_hits)} hits after redundancy check")
        return kept_hits

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
        """Extract the sequence of a BLAST hit, oriented to match the subject
        (reference gene) sense.

        BLAST reports qstart/qend and sstart/send independently; either can be
        descending depending on which strand the hit lies on. The extracted
        sequence must be reverse-complemented whenever the query and subject
        strands differ (not merely whenever qstart > qend), otherwise genes
        aligned to the subject's minus strand get translated on the wrong
        strand entirely.
        """
        query_id = hit['query_id']
        qstart = hit['qstart']
        qend = hit['qend']
        sstart = hit['sstart']
        send = hit['send']

        for record in SeqIO.parse(self.assembly, 'fasta'):
            if record.id == query_id:
                q_lo, q_hi = min(qstart, qend), max(qstart, qend)
                seq = record.seq[q_lo-1:q_hi]
                query_plus = qstart <= qend
                subject_plus = sstart <= send
                if query_plus != subject_plus:
                    seq = seq.reverse_complement()
                return str(seq)
        return None

    def analyze_hits(self, hits):
        """Analyze BLAST hits and detect mutations"""
        print("Analyzing hits and detecting mutations...")

        for hit in hits:
            sequence = self.extract_hit_sequence(hit)

            if sequence:
                if hit.get('qseq') and hit.get('sseq'):
                    mutations = detect_mutations_aligned(
                        hit['gene'], hit['qseq'], hit['sseq'],
                        hit['sstart'], hit['send'], self.mutation_db
                    )
                else:
                    mutations = []

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

        # Annotate each result with how many distinct loci carry the same gene
        copy_counts = Counter(r['gene'] for r in self.results)
        for result in self.results:
            result['copy_number'] = copy_counts[result['gene']]

    def write_report(self):
        """Write results to TSV file"""
        report_file = f"{self.output_prefix}_results.tsv"

        print(f"Writing results to {report_file}...")

        with open(report_file, 'w') as f:
            f.write('\t'.join(['Contig', 'Gene', 'Identity%', 'Coverage%',
                              'Mutations', 'Method', 'Copy_Number']) + '\n')

            for result in self.results:
                f.write('\t'.join([
                    result['contig'],
                    result['gene'],
                    result['identity'],
                    result['coverage'],
                    result['mutations'],
                    'BLAST',
                    str(result['copy_number'])
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
