"""
Create reference database for FOS-CAZAVI resistance genes
Uses AMRfinderPlus data to source sequences and mutation definitions.

Database preparation for GAMMA:
  After building the nucleotide FASTA, run GAMMA_DB_Maker.py to validate
  reading frames, remove nonstandard bases, and deduplicate sequences:

    python3 GAMMA_DB_Maker.py <output_prefix>.fasta

  The formatted database (<output_prefix>_Formatted.fasta) is then ready
  for use with GAMMA.  See https://github.com/rastanton/GAMMA_DB_Maker
"""

import csv
import subprocess
import urllib.request
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
import time

# URLs for AMRfinderPlus data
AMRFINDER_URLS = {
    'catalog': 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt',
    'sequences': 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS.fa',
    'mutations': 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt-mutation.tsv'
}

# Genes of interest
# Chromosomal genes requiring specific reference sequences (wildtype)
CHROMOSOMAL_GENES = {
    # Escherichia coli genes (Fosfomycin resistance)
    'murA': {'taxgroup': 'Escherichia'},
    'uhpT': {'taxgroup': 'Escherichia'},
    'glpT': {'taxgroup': 'Escherichia'},
    'uhpA': {'taxgroup': 'Escherichia'},
    'uhpB': {'taxgroup': 'Escherichia'},
    'uhpC': {'taxgroup': 'Escherichia'},
    'cyaA': {'taxgroup': 'Escherichia'},
    'ptsI': {'taxgroup': 'Escherichia'},
    'galU': {'taxgroup': 'Escherichia'},
    'lon': {'taxgroup': 'Escherichia'},
    
    # Klebsiella pneumoniae genes (Ceftazidime-Avibactam resistance)
    'acrB': {'taxgroup': 'Klebsiella_pneumoniae'},
    'ompK36': {'taxgroup': 'Klebsiella_pneumoniae'},
    'ftsI': {'taxgroup': 'Klebsiella_pneumoniae'},
    'envZ': {'taxgroup': 'Klebsiella_pneumoniae'},
    'ompK35': {'taxgroup': 'Klebsiella_pneumoniae'},

    # Pseudomonas aeruginosa
    'mexR': {'taxgroup': 'Pseudomonas_aeruginosa'},
    'nalD': {'taxgroup': 'Pseudomonas_aeruginosa'},
}

# Acquired genes (fetched from AMR_CDS.fa by family/name)
ACQUIRED_GENES = [
    'fosA3', 'fosA4', 'fosA5', 'fosA7', 'fosA11',
    'blaKPC-2', 'blaKPC-3', 'blaOXA-48',
    'blaKPC-31', 'blaKPC-190', 'blaCMY-178', 'blaSHV-12'
]

# Fallback for genes missing in AMRfinderPlus catalog
# Gene -> (Accession, Start, Stop, Strand)
FALLBACK_GENES = {
    # E. coli K-12 MG1655 (U00096.3)
    'glpT': ('U00096.3', 2351016, 2352374, '-'),
    'uhpB': ('U00096.3', 3848634, 3850136, '-'),
    'uhpC': ('U00096.3', 3847305, 3848624, '-'),
    'galU': ('U00096.3', 1291457, 1292365, '+'),

    # K. pneumoniae MGH 78578 (CP000647.1)
    'acrB': ('CP000647.1', 484760, 487906, '-'),
}

# Manual genes to fetch directly (not in catalog logic)
MANUAL_GENES = {
    # K. pneumoniae chromosomal fosA (FosA6), used as fosAKP
    'fosAKP': ('KU254579.1', 59422, 59841, '+')
}

# Manual mutations (for acquired genes or those missing in AMRfinderPlus)
# Format: Gene -> Position -> {ref, variants, name}
MANUAL_MUTATIONS = {
    'blaKPC': [
        {'pos': 179, 'ref': 'D', 'variants': ['Y', 'N'], 'name': 'D179Y/N'},
        {'pos': 240, 'ref': 'V', 'variants': ['G'], 'name': 'V240G'},
        {'pos': 243, 'ref': 'T', 'variants': ['M'], 'name': 'T243M'}
    ],
    'blaOXA-48': [
        {'pos': 68, 'ref': 'P', 'variants': ['A'], 'name': 'P68A'},
        {'pos': 211, 'ref': 'Y', 'variants': ['S'], 'name': 'Y211S'}
    ],
}


class DatabaseBuilder:
    def __init__(self, email, output_prefix):
        self.email = email
        self.output_fasta = f"{output_prefix}.fasta"
        self.output_proteins = f"{output_prefix}_proteins.fasta"
        self.output_mutations = f"{output_prefix}_mutations.tsv"
        Entrez.email = email
        self.sequences = []
        self.protein_sequences = []
        self.mutations = [] # List of dicts
        self.download_dir = Path("amrfinder_data")
        self.download_dir.mkdir(exist_ok=True)
        
        self.catalog_file = self.download_dir / "ReferenceGeneCatalog.txt"
        self.sequences_file = self.download_dir / "AMR_CDS.fa"
        self.mutations_file = self.download_dir / "AMRProt-mutation.tsv"

    def download_files(self):
        """Download AMRfinderPlus files if not present"""
        print("Checking AMRfinderPlus data files...")
        for name, url in AMRFINDER_URLS.items():
            filepath = self.download_dir / Path(url).name
            if not filepath.exists():
                print(f"Downloading {name} from {url}...")
                urllib.request.urlretrieve(url, filepath)
            else:
                print(f"Using existing {name}: {filepath}")

    def load_catalog(self):
        """Parse ReferenceGeneCatalog.txt"""
        print("Parsing ReferenceGeneCatalog...")
        self.catalog_data = []
        with open(self.catalog_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                self.catalog_data.append(row)

    def load_mutation_defs(self):
        """Parse AMRProt-mutation.tsv"""
        print("Parsing mutation definitions...")
        self.mutation_defs = []
        with open(self.mutations_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                self.mutation_defs.append(row)

    def fetch_chromosomal_genes(self):
        """Fetch chromosomal reference sequences based on catalog"""
        print("Processing chromosomal genes...")
        
        for gene, info in CHROMOSOMAL_GENES.items():
            print(f"Looking up {gene} ({info['taxgroup']})...")
            
            candidates = [
                row for row in self.catalog_data
                if row['gene_family'] == gene 
                and info['taxgroup'] in row['whitelisted_taxa']
                and row['type'] == 'AMR' 
                and (row['subtype'] == 'POINT' or row['subtype'] == 'POINT_DISRUPT')
            ]
            
            if not candidates:
                candidates = [
                    row for row in self.catalog_data
                    if gene in row['allele']
                    and info['taxgroup'] in row['whitelisted_taxa']
                    and row['type'] == 'AMR' 
                    and (row['subtype'] == 'POINT' or row['subtype'] == 'POINT_DISRUPT')
                ]

            if candidates:
                ref_entry = candidates[0]
                
                if ref_entry['refseq_nucleotide_accession']:
                    accession = ref_entry['refseq_nucleotide_accession']
                    start = int(ref_entry['refseq_start'])
                    stop = int(ref_entry['refseq_stop'])
                    strand = ref_entry['refseq_strand'] # + or -
                elif ref_entry['genbank_nucleotide_accession']:
                    print(f"  WARNING: RefSeq missing for {gene}, using GenBank")
                    accession = ref_entry['genbank_nucleotide_accession']
                    start = int(ref_entry['genbank_start'])
                    stop = int(ref_entry['genbank_stop'])
                    strand = ref_entry['genbank_strand']
                else:
                    print(f"  ERROR: No accession found for {gene}")
                    continue
                
                print(f"  Found reference: {accession} ({start}-{stop}, {strand})")
                
                record = self.fetch_sequence_from_ncbi(accession, start, stop, strand, gene)
                if record:
                    self.sequences.append(record)
                    self.extract_mutations_for_gene(gene, info['taxgroup'])

                    # Translate and add to protein sequences
                    try:
                        prot_seq = record.seq.translate(table=11, cds=True)
                        prot_record = SeqRecord(prot_seq, id=gene, description=f"{gene} reference protein (translated)")
                        self.protein_sequences.append(prot_record)
                    except Exception as e:
                        print(f"  WARNING: Could not translate {gene}: {e}")

            elif gene in FALLBACK_GENES:
                print(f"  WARNING: {gene} not found in catalog, using fallback")
                acc, start, stop, strand = FALLBACK_GENES[gene]
                print(f"  Using fallback: {acc} ({start}-{stop}, {strand})")
                record = self.fetch_sequence_from_ncbi(acc, start, stop, strand, gene)
                if record:
                    self.sequences.append(record)
                    # We can't extract mutations easily without catalog protein accessions,
                    # unless we add them manually.

                    # Translate and add to protein sequences
                    try:
                        prot_seq = record.seq.translate(table=11, cds=True)
                        prot_record = SeqRecord(prot_seq, id=gene, description=f"{gene} reference protein (translated)")
                        self.protein_sequences.append(prot_record)
                    except Exception as e:
                        print(f"  WARNING: Could not translate {gene}: {e}")
            else:
                print(f"  WARNING: {gene} not found in AMRfinderPlus catalog for {info['taxgroup']}")

    def fetch_manual_genes(self):
        """Fetch manually defined genes"""
        print("Processing manual genes...")
        for gene, (acc, start, stop, strand) in MANUAL_GENES.items():
            print(f"Fetching manual gene {gene} from {acc}...")
            record = self.fetch_sequence_from_ncbi(acc, start, stop, strand, gene)
            if record:
                self.sequences.append(record)

                # Translate and add to protein sequences
                try:
                    prot_seq = record.seq.translate(table=11, cds=True)
                    prot_record = SeqRecord(prot_seq, id=gene, description=f"{gene} reference protein (translated)")
                    self.protein_sequences.append(prot_record)
                except Exception as e:
                    print(f"  WARNING: Could not translate {gene}: {e}")

    def fetch_sequence_from_ncbi(self, accession, start, stop, strand, gene_name):
        """Fetch sequence from NCBI"""
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=accession,
                rettype="fasta",
                retmode="text",
                seq_start=start,
                seq_stop=stop
            )
            record = SeqIO.read(handle, "fasta")
            handle.close()

            if strand == '-':
                record.seq = record.seq.reverse_complement()
            
            record.id = gene_name
            record.description = f"{gene_name} reference sequence (AMRfinderPlus:{accession}:{start}-{stop})"
            
            if not record.seq.startswith("ATG") and not record.seq.startswith("GTG") and not record.seq.startswith("TTG"):
                 print(f"  WARNING: Sequence for {gene_name} does not start with ATG/GTG/TTG: {record.seq[:10]}...")

            return record
        except Exception as e:
            print(f"  FAILED to fetch {gene_name}: {e}")
            return None

    def extract_mutations_for_gene(self, gene_name, taxgroup):
        """Extract mutations from parsed mutation file"""
        prot_accessions = set()
        for row in self.catalog_data:
            if gene_name in row['allele'] and taxgroup in row['whitelisted_taxa']:
                 if row['refseq_protein_accession']:
                     prot_accessions.add(row['refseq_protein_accession'])

        count = 0
        for row in self.mutation_defs:
            if row['accession_version'] in prot_accessions:
                pos = row['mutation_position']
                symbol = row['standard_mutation_symbol']
                try:
                    mutation_part = symbol.split('_')[-1] # L370I
                    ref = mutation_part[0]
                    var = mutation_part[-1]
                    num = ''.join(filter(str.isdigit, mutation_part))
                    if num != pos:
                        pass
                    
                    self.mutations.append({
                        'Gene': gene_name,
                        'Position': pos,
                        'Ref': ref,
                        'Variant': var,
                        'Name': symbol
                    })
                    count += 1
                except:
                    print(f"  Skipping parse of symbol {symbol}")
        
        print(f"  Added {count} mutation definitions for {gene_name}")

    def fetch_acquired_genes(self):
        """Fetch acquired genes from AMR_CDS.fa"""
        print("Processing acquired genes...")
        
        found_genes = set()
        
        with open(self.sequences_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                header_parts = record.description.split('|')
                if len(header_parts) > 4:
                    gene_symbol = header_parts[4]
                    
                    if gene_symbol in ACQUIRED_GENES:
                        if gene_symbol not in found_genes:
                            record.id = gene_symbol
                            record.description = f"{gene_symbol} reference sequence (AMRfinderPlus)"
                            self.sequences.append(record)
                            found_genes.add(gene_symbol)
                            print(f"  Found {gene_symbol}")

                            # Translate and add to protein sequences
                            try:
                                prot_seq = record.seq.translate(table=11, cds=True)
                                prot_record = SeqRecord(prot_seq, id=gene_symbol, description=f"{gene_symbol} reference protein (translated)")
                                self.protein_sequences.append(prot_record)
                            except Exception as e:
                                print(f"  WARNING: Could not translate {gene_symbol}: {e}")
        
        for gene in ACQUIRED_GENES:
            if gene not in found_genes:
                print(f"  WARNING: Acquired gene {gene} not found in AMR_CDS.fa")

    def add_manual_mutations(self):
        """Add manual/hardcoded mutations for genes not covered by AMRfinderPlus point mutations"""
        print("Adding manual mutation definitions...")
        for gene, mutations in MANUAL_MUTATIONS.items():
            for m in mutations:
                self.mutations.append({
                    'Gene': gene,
                    'Position': m['pos'],
                    'Ref': m['ref'],
                    'Variant': '/'.join(m['variants']),
                    'Name': m['name']
                })

    def prepare_gamma_db(self):
        """Run GAMMA_DB_Maker on the nucleotide FASTA to prepare it for GAMMA.

        GAMMA_DB_Maker validates reading frames, removes sequences with
        nonstandard bases, verifies start/stop codons, and deduplicates
        entries.  The formatted output (<prefix>_Formatted.fasta) is the
        database to pass to GAMMA via the --genes argument.

        Requires GAMMA_DB_Maker.py to be on PATH or in the current directory.
        See https://github.com/rastanton/GAMMA_DB_Maker
        """
        import shutil
        gamma_db_maker = shutil.which('GAMMA_DB_Maker.py') or 'GAMMA_DB_Maker.py'

        print(f"Running GAMMA_DB_Maker on {self.output_fasta}...")
        try:
            subprocess.run(
                ['python3', gamma_db_maker, self.output_fasta],
                check=True
            )
            formatted = Path(self.output_fasta).stem + '_Formatted.fasta'
            if Path(formatted).exists():
                print(f"GAMMA-formatted database written to {formatted}")
                print(f"Use '--genes {formatted}' when running GAMMA mutation detection.")
            else:
                print("WARNING: GAMMA_DB_Maker ran but formatted output not found; "
                      f"using {self.output_fasta} directly.")
        except FileNotFoundError:
            print("WARNING: GAMMA_DB_Maker.py not found. Skipping database formatting. "
                  "Download it from https://github.com/rastanton/GAMMA_DB_Maker and run:\n"
                  f"  python3 GAMMA_DB_Maker.py {self.output_fasta}")
        except subprocess.CalledProcessError as e:
            print(f"WARNING: GAMMA_DB_Maker failed: {e}")

    def build(self):
        self.download_files()
        self.load_catalog()
        self.load_mutation_defs()

        self.fetch_chromosomal_genes()
        self.fetch_manual_genes()
        self.fetch_acquired_genes()
        self.add_manual_mutations()

        print(f"Writing sequences to {self.output_fasta}...")
        SeqIO.write(self.sequences, self.output_fasta, "fasta")

        if self.protein_sequences:
            print(f"Writing protein sequences to {self.output_proteins}...")
            SeqIO.write(self.protein_sequences, self.output_proteins, "fasta")

        print(f"Writing mutations to {self.output_mutations}...")
        with open(self.output_mutations, 'w') as f:
            writer = csv.DictWriter(f, fieldnames=['Gene', 'Position', 'Ref', 'Variant', 'Name'], delimiter='\t')
            writer.writeheader()
            writer.writerows(self.mutations)

        self.prepare_gamma_db()

        print("Done!")

def create_db(email, output_prefix):
    builder = DatabaseBuilder(email, output_prefix)
    builder.build()
