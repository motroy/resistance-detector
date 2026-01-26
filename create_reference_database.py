#!/usr/bin/env python3
"""
Create reference database for FOS-CAZAVI resistance genes
Uses AMRfinderPlus data to source sequences and mutation definitions.
"""

import argparse
import sys
import os
import csv
import urllib.request
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
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
    'envZ': {'taxgroup': 'Klebsiella_pneumoniae'}
}

# Acquired genes (fetched from AMR_CDS.fa by family/name)
ACQUIRED_GENES = [
    'fosA3', 'fosA4', 'fosA5', 'fosA7',
    'blaKPC-2', 'blaKPC-3', 'blaOXA-48'
]

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
    # Fallback for chromosomal if not found in AMRfinderPlus
    # (But we prefer AMRfinderPlus if available)
}


class DatabaseBuilder:
    def __init__(self, email, output_prefix):
        self.email = email
        self.output_fasta = f"{output_prefix}.fasta"
        self.output_mutations = f"{output_prefix}_mutations.tsv"
        Entrez.email = email
        self.sequences = []
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
            
            # Find in catalog
            # We look for entries with matching gene symbol and taxgroup
            # Ideally type=AMR, subtype=POINT
            candidates = [
                row for row in self.catalog_data
                if row['gene_family'] == gene 
                and info['taxgroup'] in row['whitelisted_taxa']
                and row['type'] == 'AMR' 
                and row['subtype'] == 'POINT'
            ]
            
            if not candidates:
                # Try relaxed search (gene symbol in product name or allele?)
                # Some genes might be listed under specific alleles e.g. murA_C115D
                candidates = [
                    row for row in self.catalog_data
                    if gene in row['allele']
                    and info['taxgroup'] in row['whitelisted_taxa']
                    and row['type'] == 'AMR' 
                    and row['subtype'] == 'POINT'
                ]

            if candidates:
                # Pick the first one as reference source
                # Usually all mutations for a gene in a taxgroup point to the same reference accession
                ref_entry = candidates[0]
                
                # Try RefSeq first, then GenBank
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
                
                # Fetch sequence
                record = self.fetch_sequence_from_ncbi(accession, start, stop, strand, gene)
                if record:
                    self.sequences.append(record)
                    
                    # Add mutations from catalog/mutation file
                    self.extract_mutations_for_gene(gene, info['taxgroup'])
            else:
                print(f"  WARNING: {gene} not found in AMRfinderPlus catalog for {info['taxgroup']}")
                # Here we could fall back to hardcoded if we wanted, but user said 'use them'
                # and implied sticking to AMRfinderPlus.

    def fetch_sequence_from_ncbi(self, accession, start, stop, strand, gene_name):
        """Fetch sequence from NCBI"""
        try:
            # Entrez uses 1-based coordinates
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

            # Reverse complement if needed
            if strand == '-':
                record.seq = record.seq.reverse_complement()
            
            record.id = gene_name
            record.description = f"{gene_name} reference sequence (AMRfinderPlus:{accession}:{start}-{stop})"
            
            # Check for valid start codon (simple check)
            if not record.seq.startswith("ATG") and not record.seq.startswith("GTG") and not record.seq.startswith("TTG"):
                 print(f"  WARNING: Sequence for {gene_name} does not start with ATG/GTG/TTG: {record.seq[:10]}...")

            return record
        except Exception as e:
            print(f"  FAILED to fetch {gene_name}: {e}")
            return None

    def extract_mutations_for_gene(self, gene_name, taxgroup):
        """Extract mutations from parsed mutation file"""
        # Find all mutations for this gene in the mutation file
        # The mutation file uses 'accession_version' (protein) to link.
        # We need to map gene_name -> protein_accession(s) first.
        
        # Get protein accessions from catalog
        prot_accessions = set()
        for row in self.catalog_data:
            if gene_name in row['allele'] and taxgroup in row['whitelisted_taxa']:
                 if row['refseq_protein_accession']:
                     prot_accessions.add(row['refseq_protein_accession'])

        count = 0
        for row in self.mutation_defs:
            if row['accession_version'] in prot_accessions:
                # Parse mutation: e.g. murA_L370I -> pos=370, ref=L, var=I
                # 'mutation_position': '370'
                # 'standard_mutation_symbol': 'murA_L370I'
                # 'mutated_protein_name': ...
                
                pos = row['mutation_position']
                symbol = row['standard_mutation_symbol']
                # Parse symbol usually: Gene_RefPosVar
                # But sometimes complex.
                # Let's trust the position column.
                # But we need ref and var.
                # Try to parse symbol: murA_L370I
                try:
                    mutation_part = symbol.split('_')[-1] # L370I
                    ref = mutation_part[0]
                    var = mutation_part[-1]
                    # verify position matches
                    num = ''.join(filter(str.isdigit, mutation_part))
                    if num != pos:
                        # Sometimes symbol is different format?
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
        
        # Build index of AMR_CDS.fa
        # Format: >ProtAcc|NucAcc|...|GeneSymbol|...
        
        found_genes = set()
        
        with open(self.sequences_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                header_parts = record.description.split('|')
                if len(header_parts) > 4:
                    gene_symbol = header_parts[4]
                    
                    # Check if this gene is in our list
                    if gene_symbol in ACQUIRED_GENES:
                        # Avoid duplicates
                        if gene_symbol not in found_genes:
                            record.id = gene_symbol
                            record.description = f"{gene_symbol} reference sequence (AMRfinderPlus)"
                            self.sequences.append(record)
                            found_genes.add(gene_symbol)
                            print(f"  Found {gene_symbol}")
        
        # Check missing
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
                    'Variant': '/'.join(m['variants']), # Support multiple variants in string?
                    'Name': m['name']
                })
                # Note: resistance_detector.py needs to parse 'Variant' carefully if it contains '/'
                # Currently it expects single char check usually, but the new logic will iterate variants.

    def build(self):
        self.download_files()
        self.load_catalog()
        self.load_mutation_defs()
        
        self.fetch_chromosomal_genes()
        self.fetch_acquired_genes()
        self.add_manual_mutations()
        
        # Write outputs
        print(f"Writing sequences to {self.output_fasta}...")
        SeqIO.write(self.sequences, self.output_fasta, "fasta")
        
        print(f"Writing mutations to {self.output_mutations}...")
        with open(self.output_mutations, 'w') as f:
            writer = csv.DictWriter(f, fieldnames=['Gene', 'Position', 'Ref', 'Variant', 'Name'], delimiter='\t')
            writer.writeheader()
            writer.writerows(self.mutations)
            
        print("Done!")

def main():
    parser = argparse.ArgumentParser(description='Create reference database using AMRfinderPlus data')
    parser.add_argument('-e', '--email', required=True, help='Email for NCBI Entrez')
    parser.add_argument('-o', '--output', default='resistance_db', help='Output prefix')
    args = parser.parse_args()
    
    builder = DatabaseBuilder(args.email, args.output)
    builder.build()

if __name__ == '__main__':
    main()
