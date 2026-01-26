import sys
import os
import unittest
from pathlib import Path

# Add parent dir to path
sys.path.append(str(Path(__file__).parent.parent))

from resistance_detector import ResistanceDetector

class TestAmplicons(unittest.TestCase):
    def test_amplicon_detection(self):
        # Create dummy assembly with primers forming an amplicon
        assembly = "test/dummy_amp_assembly.fasta"

        # uhpB_F -> uhpB_R (uhpB_ver)
        # Sequence: uhpB_F + 100bp + uhpB_R (rc)
        seq_f = "ACTGGGCGTCAGTAACGACG" # uhpB_F
        seq_r = "ATGGCGCATCGGCAGGCGCT" # uhpB_R
        from Bio.Seq import Seq
        seq_r_rc = str(Seq(seq_r).reverse_complement())

        middle = "A" * 100
        full_seq = seq_f + middle + seq_r_rc

        with open(assembly, 'w') as f:
            f.write(f'>contig1\nNNNN{full_seq}NNNN\n')

        database = "example_database.fasta"
        output = "test/test_amp_output"
        primers_file = "primers.tsv"

        detector = ResistanceDetector(assembly, database, output, primers_file=primers_file)

        # 1. Detect amplicons (uses seqkit)
        detector.detect_amplicons()

        self.assertEqual(len(detector.amplicon_results), 1)
        amp = detector.amplicon_results[0]
        self.assertEqual(amp['pair_id'], 'uhpB_ver')
        self.assertEqual(amp['contig'], 'contig1')

        # Seqkit output might be slightly different in coordinates/logic
        # For seqkit amplicon -p primers.tsv assembly.fasta --bed
        # It finds the amplicon.
        # Check start/end.
        # F primer start is 4. Length is 20.
        # R primer start RC is 4 + 20 + 100. Length is 20.
        # End is 4 + 20 + 100 + 20 = 144.
        # Start should be 4.

        self.assertEqual(amp['start'], 4)
        self.assertEqual(amp['end'], 144)

        # 2. Mock BLAST results (gene within amplicon)
        # Gene is roughly in the middle
        gene_start = 4 + len(seq_f) + 10
        gene_end = gene_start + 50

        mock_result = {
            'contig': 'contig1',
            'gene': 'uhpB',
            'identity': '100.00',
            'coverage': '100.00',
            'mutations': 'G469R',
            'sequence': 'AAA',
            'start': gene_start + 1, # 1-based
            'end': gene_end + 1      # 1-based
        }
        detector.results = [mock_result]

        # 3. Analyze overlap
        detector.analyze_amplicons()

        # 4. Verify mutation reported in amplicon
        self.assertEqual(len(amp['mutations_found']), 1)
        self.assertIn('uhpB: G469R', amp['mutations_found'][0])

        # Clean up
        os.remove(assembly)
        if os.path.exists(f"{output}_seqkit_primers.tsv"):
            os.remove(f"{output}_seqkit_primers.tsv")

if __name__ == '__main__':
    unittest.main()
