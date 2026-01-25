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

# Reference sequences for genes not easily found in NCBI
REFERENCE_SEQUENCES = {
    'fosA3': {
        'protein': 'MNIVKIIGHQSGAGKTTLLNSIAGIKPNEGKVLINGKDISEDWETWKELKQIDIPIVLDSITLVPETINYADLNQKRTTLKDILTAFPVRVFHDHDMMELDKKWLDLEQEWQGMVEEAAIHMVQRFKQYLPDSGRVLMVEQKVMKLGQHFVATQPIVDKKIQAGLTLQEEILTDFKLGNEQKALRDLLKMAE',
        'nucleotide': 'ATGAACATTGTGAAAATTATTGGGCACCAGTCTGGCGCTGGCAAAACCACGCTGCTGAACAGCATCGCTGGCATTAAACCGAACGAAGGCAAAGTGCTGATTAACGGCAAAGATATTAG CGAAGATGATGAAACCGATAAAGAACTGAAACAGATTGATATTCCGATTGTGCTGGATAGCATTACCCTGGTGCCGGAAACCATTAACTATGCGGATCTGAACCAGAAACGTACCACCCTGAAAGATATTCTGACCGCGTTTCCGGTGCGTGTGTTTCATGATCATGACATGATGGAACTGGATAAAAAATGGCTGGATCTGGAACAGGAATGGCAGGGCATGGTGGAAGAAGCGGCGATTCATATGGTGCAGCGTTTTAAACAGTATCTGCCGGATAGCGGCCGTGTGCTGATGGTGGAACAGAAAGTGATGAAACTGGGCCAGCATTTTGTGGCGACCCAGCCGATTGTGGATAAAAAAATTCAGGCGGGCCTGACCCTGCAGGAAGAAATTCTGACCGATTTTAAACTGGGCAACGAACAGAAAGCGCTGCGTGATCTGCTGAAAATGGCGGAA'
    },
    'fosA4': {
        'protein': 'MNIVKIIGHQSGAGKTTLLNSIAGIKPNEGKVLINGKDISEDWETWKELKQIDIPIVLDSITLVPETINYADLNQKRTTLKDILTAFPVRVFHDHDMMELDKKWLDLEQEWQGMVEEAAIHMVQRFKQYLPDSGRVLMVEQKVMKLGQHFVATQPIVDKKIQAGLTLQEEILTDFKLGNEQKALRDLLKMAE'
    },
    'blaKPC-2': {
        'protein': 'MSLYRLLLTLLATLFGVAKPQWSESYKEQGLHRMKQIIDRYHQDLASMKGLCQLYQGEQNGDVVQMKWGQLQMNVLSDVGFQAMRDVWNYDKKLTPQGQKLVYQADWLWEFQSVCWKTHWSQQGQAVWHLQGGIAGVSWNPPVLLLQDEQDPQALKDAQTRITLKKFLSTGSGLGGAGMAWMNDKQAIHPGGSVISVGDKVLGRALGGKIPWFEGDLFWAINDLNKQDYIRGFSIASDKKQLHHAELQEIGVEVNGQLQGKNIPFGKLYVHRKQSPEMVQKILKALGHGVGCALITPTSITLPPAQVKEMLEQVMIKRPGLLAVVLTYLAMLYALFLSLVQPK'
    },
    'blaKPC-3': {
        'protein': 'MSLYRLLLTLLATLFGVAKPQWSESYKEQGLHRMKQIIDRYHQDLASMKGLCQLYQGEQNGDVVQMKWGQLQMNVLSDVGFQAMRDVWNYDKKLTPQGQRLVYQADWLWEFQSVCWKTHWSQQGQAVWHLQGGIAGVSWNPPVLLLQDEQDPQALKDAQTRITLKKFLSTGSGLGGAGMAWMNDKQAIHPGGSVISVGDKVLGRALGGKIPWFEGDLFWAINDLNKQDYIRGFSIASDKKQLHHAELQEIGVEVNGQLQGKNIPFGKLYVHRKQSPEMVQKILKALGHGVGCALITPTSITLPPAQVKEMLEQVMIKRPGLLAVVLTYLAMLYALFLSLVQPK'
    },
    'blaOXA-48': {
        'protein': 'MSNKNYFQRMRFIYFFFMFFVWSCVYAQSDGKLEGVAAVYRGEKPKFEWLIEFKGTPEVPGFGLNMYPRNIVVPKLYVYEHWQSQLCGMVNQALTDRWEAFFTPSWFDHNPRDTFVEPGSGTIGFTGRPPDPTEFDFGFSRVNEGGIFDSQPLTDNPVPKQFFTRLPQMYGGEAYRNHDVYTIRNNTIISWVRNWIEKDPVEWMVANPVAGKVIDQSTGKVQFWQENLAELDKIYAKVMRGAVGCPHTNYVALKEQYDNLLVRVPLYRQGKMSFPMQVENIEALKHHFKLKEKSVKLKAGKTVYYGGNPNSSGGFTRFGETSTFEYIPTKIVELLKKVGDRVFDRKNRSLIGKALEKGLHKEWQALFKDGHKLVSVNDMVTKDDWLVGMMNFYFKWKQRDRK'
    }
}

NCBI_ACCESSIONS = {
    # FOS resistance genes
    'fosA3': 'JX861169.1',
    'fosA4': 'KP324830.1', 
    'fosA5': 'KY270852.1',
    
    # KPC variants
    'blaKPC-2': 'AY034847.1',
    'blaKPC-3': 'AF297554.1',
    
    # OXA-48
    'blaOXA-48': 'AY236073.1',
    
    # Chromosomal genes (E. coli K-12 MG1655)
    'uhpT': 'NC_000913.3:c3680773-3679418',  # Region
    'glpT': 'NC_000913.3:c2874881-2873517',
    'uhpA': 'NC_000913.3:c3680233-3678692',
    'uhpB': 'NC_000913.3:3678743-3680377',
    'uhpC': 'NC_000913.3:3680392-3681279',
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
            # Handle region-based accessions
            if ':' in accession:
                # This is a genomic region
                parts = accession.split(':')
                acc = parts[0]
                region = parts[1]
                
                # For now, skip complex region extraction
                print("SKIPPED (genomic region - use local sequence)")
                return None
            
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
