#!/usr/bin/env python3
"""Test script to verify mutation detection functionality"""

import sys
sys.path.insert(0, '/home/user/resistance-detector')

from Bio.Seq import Seq

# Test the ResistanceDetector mutation detection
class TestMutationDetection:
    # Known resistance mutations (from resistance_detector.py)
    KNOWN_MUTATIONS = {
        'murA': {
            369: {'ref': 'D', 'variants': ['N'], 'name': 'D369N'},
            370: {'ref': 'L', 'variants': ['I'], 'name': 'L370I'}
        },
        'uhpT': {
            55: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G55D/*'},
            198: {'ref': 'W', 'variants': ['*', 'R'], 'name': 'W198*/R'},
        },
        'glpT': {
            44: {'ref': 'E', 'variants': ['*', 'K'], 'name': 'E44*/K'},
            88: {'ref': 'W', 'variants': ['*', 'R', 'G'], 'name': 'W88*/R/G'},
        }
    }
    
    def detect_mutations(self, gene_name, sequence):
        """Detect known mutations in a gene sequence"""
        mutations_found = []
        
        base_gene = gene_name.split('-')[0]
        if base_gene.startswith('bla'):
            base_gene = 'bla' + base_gene[3:].rstrip('0123456789')
        
        if base_gene not in self.KNOWN_MUTATIONS:
            return mutations_found
        
        try:
            if len(sequence) % 3 != 0:
                sequence = sequence[:-(len(sequence) % 3)]
            
            protein = str(Seq(sequence).translate())
            
            for pos, mut_info in self.KNOWN_MUTATIONS[base_gene].items():
                if pos <= len(protein):
                    observed_aa = protein[pos-1]
                    ref_aa = mut_info['ref']
                    
                    if observed_aa in mut_info['variants']:
                        mutations_found.append(mut_info['name'])
                    elif observed_aa != ref_aa:
                        mutations_found.append(f"{ref_aa}{pos}{observed_aa}")
        
        except Exception as e:
            print(f"Warning: Could not translate {gene_name}: {e}")
        
        return mutations_found

# Run tests
if __name__ == '__main__':
    detector = TestMutationDetection()
    
    # Test murA with D369N mutation
    # Create a sequence that would produce D at position 369, then mutate it to N
    print("Testing murA mutation detection...")
    
    # Example: Create a murA-like sequence with D369N mutation
    # Position 369 = codon at (369-1)*3 = 1104
    # D = GAT/GAC, N = AAT/AAC
    
    # Read murA sequence and modify it
    with open('/home/user/resistance-detector/test/full_resistance_db.fasta') as f:
        content = f.read()
        # Extract murA sequence
        for entry in content.split('>')[1:]:
            if entry.startswith('murA'):
                lines = entry.strip().split('\n')
                murA_seq = ''.join(lines[1:])
                break
    
    print(f"murA sequence length: {len(murA_seq)} bp")
    
    # Translate to check current amino acids at positions 369 and 370
    protein = str(Seq(murA_seq).translate())
    print(f"Protein length: {len(protein)} aa")
    
    if len(protein) >= 370:
        print(f"Position 369: {protein[368]}")
        print(f"Position 370: {protein[369]}")
    
    # Test mutation detection on the reference sequence
    mutations = detector.detect_mutations('murA', murA_seq)
    print(f"Mutations detected in reference: {mutations}")
    
    # Create a mutant sequence with D369N
    # D (GAT) -> N (AAT)
    mutant_seq = murA_seq[:1104] + 'A' + murA_seq[1105:]  # Change first base of codon 369
    mutant_protein = str(Seq(mutant_seq).translate())
    print(f"\nMutant position 369: {mutant_protein[368]}")
    
    mutations = detector.detect_mutations('murA', mutant_seq)
    print(f"Mutations detected in mutant: {mutations}")
    
    print("\n" + "="*50)
    print("SUMMARY: Mutation detection test completed")
    print("="*50)
