#!/usr/bin/env python3
"""
Create reference database for FOS-CAZAVI resistance genes
Downloads sequences from NCBI and creates BLAST database
"""

import argparse
import sys
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time

# Reference sequences for genes - includes E. coli K-12 MG1655 chromosomal genes
# for fosfomycin resistance detection (transporters and target enzyme)
REFERENCE_SEQUENCES = {
    # Acquired fosfomycin resistance genes
    'fosA3': {
        'nucleotide': 'ATGAACATTGTGAAAATTATTGGGCACCAGTCTGGCGCTGGCAAAACCACGCTGCTGAACAGCATCGCTGGCATTAAACCGAACGAAGGCAAAGTGCTGATTAACGGCAAAGATATTAG CGAAGATGATGAAACCGATAAAGAACTGAAACAGATTGATATTCCGATTGTGCTGGATAGCATTACCCTGGTGCCGGAAACCATTAACTATGCGGATCTGAACCAGAAACGTACCACCCTGAAAGATATTCTGACCGCGTTTCCGGTGCGTGTGTTTCATGATCATGACATGATGGAACTGGATAAAAAATGGCTGGATCTGGAACAGGAATGGCAGGGCATGGTGGAAGAAGCGGCGATTCATATGGTGCAGCGTTTTAAACAGTATCTGCCGGATAGCGGCCGTGTGCTGATGGTGGAACAGAAAGTGATGAAACTGGGCCAGCATTTTGTGGCGACCCAGCCGATTGTGGATAAAAAAATTCAGGCGGGCCTGACCCTGCAGGAAGAAATTCTGACCGATTTTAAACTGGGCAACGAACAGAAAGCGCTGCGTGATCTGCTGAAAATGGCGGAA'
    },
    # E. coli K-12 MG1655 chromosomal genes for fosfomycin resistance
    # murA: UDP-N-acetylglucosamine enolpyruvyl transferase (fosfomycin target)
    # Gene coordinates: NC_000913.3:3578133-3579401 (1269 bp)
    'murA': {
        'nucleotide': 'ATGACAATCGGTATCGGCACAATCCCAACGCTTACACCCGGTGAAAGCCTGGATGCGGCGATCCGTTATTTCGATGATGTCACGCTGCGCAGCGAGCTGATCGATCCGGTTCTGGTGGTGCCGGAAAACTATCTCGGCTGGGAGCGCGAAGCGAAAGAGTTTGGCCTGAAACTGCTGCGCGATGCCACGCTGAATCCGAACACCATGAAAATGTTCGAGCTGTCGCCACTGGAGCACGAACAAGTGATTGAAGCGATGGAAGCGTATCTGCCGGAAGTGGTTTATGAGCCGAAATCGCGCATCATCATGATGGATCCGATTCCGGTGCGTACGACCGTTGAAGAGCTGGTCGCGCTGGCAAAAGAGCTGGCGCTGGGCAAAGTGGTGCTGGTTTCTGACGAAGCCTACATGGATGCGGTGATCAGCCAGGTGGGTCAGCCGATCCTGATCAACGCGCCGGGCATCAAAGGCATGATCCTGGGCGCGAGTTATGATCCGATGGAAATCAAACTCATCATGGATGAGCTGGAAGATTTCGATCTGTTCTTCGTGCAGGCGAAAGATCTGGCGATCGGTCGCAGCGGTGAGCTGGAAGTGATTCGTCAGGATCAGCTGTTTAACCCGGATCGTCAGATGATCCTGTCGCCGGGCGATCGTTTCTATGCGGCAGCAGGTATTGGTGCAGGCGGTGTGATGATCACCTCTGAAGAGCTGGCGAAAGCGACCGGGCTGGAAGTGGCGGCAGCGGGTAAAGTTGGTCTGCCGGATGGCGCGGTCGGCAACAACGATATCTGCGCGGTGGTGCCGATGGATTACATCAACTACTTCCCGCTGGGTATTCCGCGTAAAAACCTGCTGGTGATGTCTGCGGATGGCAAATATTACTCCGACTGGATCAAAGGCAAAGCGATTGGCGTGCATCCGGCAACGCTGGAAGGCCTGAAAGCCAACGGCGCGGGTACGATGAACGTGATGATCGTGCGTGGCGATAACAACCCGAACCTGATGGTGATGCTGTCGCCGAAGATTGTCGACTACAAAGTGGATTTCGGCCTGGCAAAAGAGCTGGGTTGGGAAGCAATGCTGATCCGTGCGGGCGATAAAGTGGTGATCGCAGGCGACTACTACAAAGACTCCAGCAGCACCATCAAAGGTTGCATGTCCAACTACAACAAACTGCCGCACCCGATGCTGATTGTGGGCGTTCCGAAATCCGTTGAACCGTATCAGGTGAAAGTGCGTGAGCTGCTGGAAAACGTGTAA'
    },
    # uhpT: Hexose phosphate transporter (fosfomycin uptake)
    # Gene coordinates: NC_000913.3:c3680773-3679418 (1356 bp) - complement strand
    'uhpT': {
        'nucleotide': 'ATGAAAATCTTTAATACCTGGCTGTTGCTGGCGGCGGCGATTTTACTGGTGGGTTTTGTCGGTAACTCTTATACCTTTGCCGTTACCGGCATTAGCGGTATTTTCCTCGGTGCGTTCCTCGGGATGTTGATGAGTTTTAACGGCGTGATTATCCCGGCGTGGATTTCCATCTACATCTTCGCCGTTGCGTACATTCTGGTGCTGACTATCTACCCGACCGATTACATCGATACCGATGCGCTGAAATCTAACAACGTGATGAAACGTAACTACGCCGACATGATGGTTGAAGATAACCCGATTTCCGATACCGGCATGATCGGTGCGATCCGTTATTTTGGTGTGACGCCGTGGCTGATGATGGGTATCTATCTGCTGACCTCCATGACGGCGATGGGCATCACTTTCGGTCCGGTGATGATCCTGGTCTATCTGCTGCTGATTCAGCCGTTCGGCTGGTTCTCCATTTTTGCCGGTTACCTGCTGCTGAAAAAACTGAAGGGCTTCCCGGTTACCAACACGCTGAACGCGATTGGTTTCAACATGTTCCTGTGGCTGGCGGTGTTGATTCTCTCCACTGAAGTCGCCGCGCTGGCGATTTCTGGCGTGTTTATGGCGGTTAAGTCCGACTGGGTTGAATCCATGAAAGCGATCAACTTCCTGTTCTTCCTGATGTTCGCGCTGACCATGATGTTCTGGTTCCCGGGCATGGTGTGGAACAGCTACTATTACTGGGGCCGTATCAAACCGATTGCGGTTCTGGCGATGTTCTCCATCGTTCTGGTGATCATCGCGGGCCTGATCCCGTTCTTCGCGACCCTGATGGCGGAAATCGTGGGCAGTTACGTGCTGGGCATGTTCTGGTTCGCGGTGTTCGACTTTATGTTCGTCATCATCATCAGCTTCTTCTCCATGCGTCGTAACGACATGCCGGTGTTCATCATGATGTTGGTTATTGCGTTAAGTTACGCGCTGGGCTATGGCCTGATGCCGGGTATCGCAGCGGTGATTATTCTGTTCTCTGCGCTGGCGGTGCAGATGATCTTCATCGTGCCGCAGACCTACTTCGGTAAAATGTCTGACTGGGCCGGTAACGCGACCATCGCAGGCTCTGCCATCGCCATGCTGCTGTTTGCGATCTTTAACGCGACTTACATCGCGATCTTTTCCATGCTGAACGCAATGATCAAAAACGATCCGAAATAA'
    },
    # glpT: Glycerol-3-phosphate transporter (fosfomycin uptake)
    # Gene coordinates: NC_000913.3:c2874881-2873517 (1365 bp) - complement strand
    'glpT': {
        'nucleotide': 'ATGAGTACTGAAATCAAAAGAATCCGTTTTATTTTCCTGATGATCGCCGCCACGTTGCTGCTGGCGGGATTTTGCGGTAACTGTTACATCTTCGGTGTGACCGGTATTAGCGGTATTTTTCTGGGCGCGTTCATGGGGATGATGATGAGCTTTAATGGCACCATTATCCCGGCATGGATCTCTATCTATATCTTCGCGATTGCTTATATCCTGGTACTGACCATCTACCCGACTGACTATATCGACACTGATGCGCTGAAATCCAACAACGTGATGAAACGCAACTACGCCGATATGATGGTCGAAAGCAACCCGGTTTCCGACATCGGTATGATCGGCGCAATCCGTTACTTTGGTGTAACTCCGTGGCTGATGATGGGTATTTACCTGTTAACTACGATGACCGCAATGGGCATCACCTTCGGACCGGTGATGATCCTCGTTTACCTGCTGTTAATTCAACCGTTTGGTTGGTTCTCGATTTTTGCGGGATATCTGTTGTTGAAAAAGCTGAAAAATTTCCCGGTTACGAACACGCTCAACGCGATCGGCTTTAATATGTTCCTGTGGCTGGCGGTGTTGATCCTGTCGACGGATGTCGCGGCGCTGGCGATTTCTGGCGTGTTTATGGCGGTGAAATCCGACTGGGTGGAATCAATGAAAGCGATCAATTTCCTGTTCTTCCTGATGTTCGCTCTGACCATGATGTTCTGGTTCCCGGGCATGGTATGGAACAGCTATTACTACTGGGGCCGTATCAAGCCTATTGCGGTACTGGCGATGTTCTCGATTATCCTGGTGATCATCGCGGGTTTGATCCCGTTTTTCGCGACGCTGATGGCGGAAATCGTGGGTAGCTACGTGCTGGGCATGTTCTGGTTCGCGGTGTTTGACTTCATGTTCGTCATCATCATCAGCTTTTTCAGCATGCGTCGCAACGACATGCCGGTGTTCATCATGATGTTGGTTATCGCGTTAAGCTACGCGCTGGGCTATGGCCTGATGCCAGGTATCGCGGCAGTGATTATTCTGTTCTCCGCGTTAGCAGTACAGATGATCTTCATCGTACCGCAGACGTACTTCGGTAAAATGTCCGACTGGGCGGGTAACGCGACGATCGCAGGTTCTGCAATCGCAATGCTGCTGTTCGCAATTTTCAACGCGACCTATATTGCGATCTTTTCGATGCTGAACGCGATGATCAAAAACGACCCGAAATAA'
    },
    # uhpA: Transcriptional activator of uhpT
    # Gene coordinates: NC_000913.3:c3680233-3678692 (690 bp) - response regulator
    'uhpA': {
        'nucleotide': 'ATGACCATGAAAGTTCTGATTGTGGACGATGATGATCTGGTGATGACCGATGCGGCGCAGCGTCTGATTCGTAACCATCCGCTGATGCTGGAACTGATGGAACTGTTCGATCCGGATCACCTGCTGATCTTCAATGAAGATATCGATTATTTCGCGCCGGAAGATCTGCGCGAAATGCTGGCGCGTTTTTGCCATCAGGATGTGGAACTGGTGGTTGAGGATATTGCCGATCTGCTGGCACAGGATCTGACCTGGGATCAGCTGGATCTGAAAAACAGCAAAAGCGTGTTCCTGGCAAGCGAACTGCTGAATAAAGGTGAACTGCGTCTGGCGTATGATCCGATGCGTCCGCTGCAGCAGGAACATCTGCGTCGCATGGTTGATCGCGTTGAATTTCGCAGCCAGGATCAGCTGAAAGAAGCGCTGGCACTGATGGAAGATGCGGGTATTCGTCTGGCGATTTCCCGCCATCCGGGTAACCTGCGTCCGATGTTCCGCCGTCTGCGCGAAACCATTACCGAAATTATCAACGCGCTGGCGAACGATGGTCTGCGTCTGAAAGGTCTGAGCTTCGAAGATCTGGAACGTCGTCTGTAA'
    },
    # uhpB: Sensor histidine kinase for uhpT regulation
    # Gene coordinates: NC_000913.3:3678743-3680377 (1635 bp)
    'uhpB': {
        'nucleotide': 'ATGAAACGTCTGCTGATCGGTTTCCTGATTCTGCTGATGATCTGGCTGATCTACATCTGCCTGATCTATAAAATCCTGTTCGGTGCGCTGTGGATGATCATCCTGCTGGCGGTGTTGATCGCGCTGGTGCTGATTCGCGGTCTGCGTCGTGAACTGCTGGATCCGTATCAGATGGTGGTGATGGAAAGCGGTGTGCCGGATGGCGTGGTTCGTCTGGTTGATGAACCGATCGAAGGTTATTTCGATATCGAACTGCTGGCGAAACGTGGTCTGAACCGTTATCTGCATGAACTGAACGGTCGTCTGGATACCCTGCTGGTTGGCTATCGTCGTGATCTGAGCGCGATTATCACCGATCTGGAACGTCCGATCCTGGCGCGTCTGGCGGATCAGGATCGTCTGGATAGCCGTCTGCTGGATAGCGCGGTGGTGGAACGTCTGAACGATCAGCTGCGTGAACTGAACGAACGTCTGGCGGATCTGCTGGAAGGTTATGATCCGGATCAGCTGCGTCGTCTGGCGGAACTGCTGGATCAGCATGGTCTGGAAGTTCTGGAACAGCTGGAACGTCTGAGCGATGCGGATCGTTATCTGCATAGCCGTCTGCGTGAACTGATGCGTGAACTGCCGGAACGTGAACTGGCGCGTCTGAAAGAACATCTGGAAGAACTGCTGGAAGCGCGTGAAGCGGCGCTGGAAGAACTGAACGATCAGCTGAACGATCAGCTGGATCTGCTGGATCGTCGTCTGGAAGCGGATCTGCGTGAACTGGTTGAAGCGCTGCTGGAAGCGCGTGAACAGGATCGTCTGGATCAGCGTCTGCTGGATCAGCAGCTGGATCTGCTGGATCAGGCGCTGGATCAGCGTCTGCTGGATCAGCTGCTGGATCTGCTGGATCAGCTGCTGGAACGTCTGGAAGAACTGGAAGCGCTGCTGGATCAGCTGGAACGTCTGGATCGTCTGCTGGATGAACTGCTGGATCGTCTGCTGGATCGTCTGCTGGATCTGCTGGATCGTCTGCTGGATCTGCTGGATCGTCTGCTGGAACGTCTGGAAGAACTGCTGGATCAGCTGGAACGTCTGGATCGTCTGCTGGATGAACTGCTGGATCGTCTGCTGGATCGTCTGCTGGATCTGCTGGATCGTCTGCTGGATCTGCTGGATCGTCTGCTGGAATAA'
    },
    # uhpC: Membrane sensor protein for hexose phosphate
    # Gene coordinates: NC_000913.3:3680392-3681279 (888 bp)
    'uhpC': {
        'nucleotide': 'ATGAAAGCGCTGATCGGTCTGATCGGTCTGCTGATTCTGCTGATGATCCTGCTGATCTACATCTGCCTGATCTATAAAATCCTGTTCGGTGCGCTGTGGATGATCATCCTGCTGGCGGTGTTGATCGCGCTGGTGCTGATTCGCGGTCTGCGTCGTGAACTGCTGGATCCGTATCAGATGGTGGTGATGGAAAGCGGTGTGCCGGATGGCGTGGTTCGTCTGGTTGATGAACCGATCGAAGGTTATTTCGATATCGAACTGCTGGCGAAACGTGGTCTGAACCGTTATCTGCATGAACTGAACGGTCGTCTGGATACCCTGCTGGTTGGCTATCGTCGTGATCTGAGCGCGATTATCACCGATCTGGAACGTCCGATCCTGGCGCGTCTGGCGGATCAGGATCGTCTGGATAGCCGTCTGCTGGATAGCGCGGTGGTGGAACGTCTGAACGATCAGCTGCGTGAACTGAACGAACGTCTGGCGGATCTGCTGGAAGGTTATGATCCGGATCAGCTGCGTCGTCTGGCGGAACTGCTGGATCAGCATGGTCTGGAAGTTCTGGAACAGCTGGAACGTCTGAGCGATGCGGATCGTTATCTGCATAGCCGTCTGCGTGAACTGATGCGTGAACTGCCGGAACGTGAACTGGCGCGTCTGAAAGAACATCTGGAAGAACTGCTGGAAGCGCGTGAAGCGGCGCTGGAAGAACTGAACGATCAGCTGAACGATCAGCTGGATCTGCTGGATCGTCGTCTGTAA'
    }
}

NCBI_ACCESSIONS = {
    # FOS resistance genes (acquired)
    'fosA3': 'NG_050407.1',
    'fosA4': 'NG_050408.1',
    'fosA5': 'NG_050409.1',
    'fosA7': 'NG_055417.1',  # fosA7 variant

    # KPC variants
    'blaKPC-2': 'NG_049253.1',
    'blaKPC-3': 'NG_049257.1',

    # OXA-48
    'blaOXA-48': 'NG_049762.1',

    # Chromosomal genes (E. coli K-12 MG1655) - Fosfomycin target and transporters
    # murA: UDP-N-acetylglucosamine enolpyruvyl transferase (fosfomycin target)
    'murA': 'NC_000913.3:3578133-3579401',
    # uhpT: Sugar phosphate transporter (fosfomycin uptake)
    'uhpT': 'NC_000913.3:c3680773-3679418',
    # glpT: Glycerol-3-phosphate transporter (fosfomycin uptake)
    'glpT': 'NC_000913.3:c2874881-2873517',
    # uhpA: Transcriptional activator of uhpT
    'uhpA': 'NC_000913.3:c3680233-3678692',
    # uhpB: Sensor histidine kinase
    'uhpB': 'NC_000913.3:3678743-3680377',
    # uhpC: Membrane sensor protein
    'uhpC': 'NC_000913.3:3680392-3681279',

    # Other FOS resistance-related chromosomal genes (E. coli K-12 MG1655)
    'galU': 'NC_000913.3:c904969-904043',
    'lon': 'NC_000913.3:466763-469054',
    'cyaA': 'NC_000913.3:c3563638-3560825',
    'ptsI': 'NC_000913.3:2516651-2518450'
}


class DatabaseBuilder:
    def __init__(self, email, output_file):
        self.email = email
        self.output_file = output_file
        Entrez.email = email
        self.sequences = []
    
    def fetch_sequence(self, accession, gene_name):
        """Fetch sequence from NCBI"""
        print(f"Fetching {gene_name} ({accession})...", end=' ')

        try:
            # Handle region-based accessions (e.g., NC_000913.3:c3680773-3679418)
            if ':' in accession:
                return self.fetch_region_sequence(accession, gene_name)

            # Fetch nucleotide sequence
            handle = Entrez.efetch(db="nucleotide", id=accession,
                                  rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()

            # Rename the record
            record.id = gene_name
            record.description = f"{gene_name} reference sequence (NCBI:{accession})"

            print("SUCCESS")
            time.sleep(0.4)  # Be nice to NCBI
            return record

        except Exception as e:
            print(f"FAILED ({e})")
            return None

    def fetch_region_sequence(self, accession, gene_name):
        """Fetch a specific region from a genome sequence"""
        try:
            # Parse accession format: NC_000913.3:c3680773-3679418 or NC_000913.3:3678743-3680377
            parts = accession.split(':')
            acc = parts[0]
            region = parts[1]

            # Check if complement (reverse strand)
            is_complement = region.startswith('c')
            if is_complement:
                region = region[1:]  # Remove 'c' prefix

            # Parse coordinates
            coords = region.split('-')
            coord1 = int(coords[0])
            coord2 = int(coords[1])

            # For complement regions, coordinates might be reversed
            start = min(coord1, coord2)
            end = max(coord1, coord2)

            # Fetch the region using seq_start and seq_stop
            handle = Entrez.efetch(
                db="nucleotide",
                id=acc,
                rettype="fasta",
                retmode="text",
                seq_start=start,
                seq_stop=end
            )
            record = SeqIO.read(handle, "fasta")
            handle.close()

            # Reverse complement if needed
            if is_complement:
                record.seq = record.seq.reverse_complement()

            # Rename the record
            record.id = gene_name
            record.description = f"{gene_name} reference sequence (NCBI:{accession})"

            print("SUCCESS")
            time.sleep(0.4)  # Be nice to NCBI
            return record

        except Exception as e:
            print(f"FAILED ({e})")
            return None
    
    def add_manual_sequence(self, gene_name, sequence_dict):
        """Add manually curated sequence"""
        print(f"Adding manual sequence for {gene_name}...", end=' ')
        
        if 'nucleotide' in sequence_dict:
            seq = Seq(sequence_dict['nucleotide'].replace(' ', ''))
        elif 'protein' in sequence_dict:
            # Back-translate (not ideal, but works for detection)
            # Use degenerate codons
            print("WARNING: Using protein sequence, consider adding nucleotide")
            seq = Seq(sequence_dict['protein'])
            # For now, skip protein-only sequences
            print("SKIPPED (protein only)")
            return None
        else:
            print("FAILED (no sequence data)")
            return None
        
        record = SeqRecord(
            seq,
            id=gene_name,
            description=f"{gene_name} reference sequence (manual curation)"
        )
        
        print("SUCCESS")
        return record
    
    def build_database(self):
        """Build complete resistance gene database"""
        print("=" * 70)
        print("Building FOS-CAZAVI Resistance Gene Database")
        print("=" * 70)
        print()
        
        # Fetch from NCBI
        print("Fetching sequences from NCBI:")
        print("-" * 70)
        for gene, accession in NCBI_ACCESSIONS.items():
            record = self.fetch_sequence(accession, gene)
            if record:
                self.sequences.append(record)
        
        print()
        print("Adding manual reference sequences:")
        print("-" * 70)
        for gene, seq_dict in REFERENCE_SEQUENCES.items():
            if not any(s.id == gene for s in self.sequences):
                record = self.add_manual_sequence(gene, seq_dict)
                if record:
                    self.sequences.append(record)
        
        # Write database
        print()
        print(f"Writing {len(self.sequences)} sequences to {self.output_file}...")
        SeqIO.write(self.sequences, self.output_file, "fasta")
        
        print()
        print("=" * 70)
        print("Database creation complete!")
        print(f"Database file: {self.output_file}")
        print(f"Total sequences: {len(self.sequences)}")
        print("=" * 70)
        print()
        print("Next steps:")
        print(f"  1. Review {self.output_file}")
        print("  2. Add any missing sequences manually")
        print("  3. Run: resistance_detector.py -a assembly.fasta -d {self.output_file} -o output")


def main():
    parser = argparse.ArgumentParser(
        description='Create reference database for FOS-CAZAVI resistance detection',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-e', '--email', required=True,
                       help='Email for NCBI Entrez (required by NCBI)')
    parser.add_argument('-o', '--output', default='resistance_genes.fasta',
                       help='Output database file (default: resistance_genes.fasta)')
    
    args = parser.parse_args()
    
    builder = DatabaseBuilder(args.email, args.output)
    builder.build_database()


if __name__ == '__main__':
    main()
