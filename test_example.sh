#!/bin/bash
# Test script demonstrating the resistance detector

set -e

echo "================================================================================"
echo "FOS-CAZAVI Resistance Detector - Example/Test"
echo "================================================================================"
echo ""

# Check if tools exist
if [[ ! -f "resistance_detector.py" ]]; then
    echo "ERROR: resistance_detector.py not found in current directory"
    exit 1
fi

echo "Creating test assembly..."
# Create a mock assembly with resistance genes
cat > test_assembly.fasta << 'TESTASM'
>contig1_with_fosA3
ATGAACATTGTGAAAATTATTGGGCACCAGTCTGGCGCTGGCAAAACCACGCTGCTGAACAGCATCGCTGGCATTAAACCGAACGAAGGCAAAGTGCTGATTAACGGCAAAGATATTAG
CGAAGATGATGAAACCGATAAAGAACTGAAACAGATTGATATTCCGATTGTGCTGGATAGCATTACCCTGGTGCCGGAAACCATTAACTATGCGGATCTGAACCAGAAACGTACCACCCTG
AAAGATATTCTGACCGCGTTTCCGGTGCGTGTGTTTCATGATCATGACATGATGGAACTGGATAAAAAATGGCTGGATCTGGAACAGGAATGGCAGGGCATGGTGGAAGAAGCGGCGATT
CATATGGTGCAGCGTTTTAAACAGTATCTGCCGGATAGCGGCCGTGTGCTGATGGTGGAACAGAAAGTGATGAAACTGGGCCAGCATTTTGTGGCGACCCAGCCGATTGTGGATAAAAAAATT
CAGGCGGGCCTGACCCTGCAGGAAGAAATTCTGACCGATTTTAAACTGGGCAACGAACAGAAAGCGCTGCGTGATCTGCTGAAAATGGCGGAATAAGGCGATTAAGCGCGTTAATCGATG
>contig2_with_blaKPC-3_D179Y
ATGTCACTGTATCGCCTTCTCCTTATTGCTATTGCTGTCGGCACTGCGCCGCCGCAATGGAGCGAGTTATAAAAGAACAAGGCCTGCATCGGATGAAGCAAATAATCGATCGGTATCATCAGGATCTGGCCTCGATGAAGGGGCTTTGCCAACTGTACCAGGGGGAAAACAACGGTGATGTTGTCCAGATATGGAAATGGGGGGCTCAACAAATGAATGTGTTATCTGATGTGGGGTTTCAGGCCATGCGTGATGTCTGGAATTATGATAAAAAACTCACGCCGCCCGGGCAACGACTGGTATACCAGGCAGATTGGCTCTGGGTTTTTGAACAATCTTGCTGGAAAACACACTGGAGCCAGCAAGGACAGGCAGTCGTTTGGCATCTGCAAGGAGGAATTGCCGGTGTATCGTGGAATCCCCGTGTTGTTGTTGCAGGATGAGCAAGACCCACAGGCTTTGAAAGATGCACAGACCCGTATCACCCTGAAGAAATTCCTGTCAACAGGCTCAGGGCTGGGTGGAGCCGGTATGGCATGGATGAATGACAAGCAAGCAATTCATCCGGGGGGATCGGTTATCTCGGTTGGAGACAAAGTGCTGGGTAGAGCCTTGGGAGGTAAGATTCCGTGGTTCGAAGGGGACTTGTGGTGGGCCATTAATGACTTGAATAAGCAAGATTATATCCGTGGATTCTCAATAGCATCGGACAAGAAGCAGTTGCACCATGCAGAGTTACAGGAGATCGGTGTGGAGGTAAATGGCCAATTGCAGGGAAAAAATATCCCGTTCGGGAAATTATACGTACACCGGAAGCAATCTCCGGAAATGGTTCAGAAGATTCTGAAAGCATTGGGTGCTCATGGAGTGGGGTGCGCATTGATTACACCGACTTCCAGCATCACTTTGCCGCCGGCACAGGTGAAGGAAATGCTGGAACAAGTCATGATTAAGCGTCCAGGTCTTTTGGCGGTTGTGCTGACTTATCTGGCGATGCTGTATGCGCTGTTCCTGTCTTTAGTGCAGCCAAAATAG
>contig3_random_sequence
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
CGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
TESTASM

echo "Test assembly created: test_assembly.fasta"
echo ""

echo "Running resistance detector..."
python3 resistance_detector.py \
    -a test_assembly.fasta \
    -d example_database.fasta \
    -o test_output

echo ""
echo "================================================================================"
echo "Test Results"
echo "================================================================================"
echo ""

echo "1. Results table (test_output_results.tsv):"
echo "-------------------------------------------"
cat test_output_results.tsv
echo ""

echo "2. Summary (test_output_summary.txt):"
echo "--------------------------------------"
cat test_output_summary.txt
echo ""

echo "3. Detected gene sequences (test_output_genes.fasta):"
echo "------------------------------------------------------"
if [[ -f "test_output_genes.fasta" ]]; then
    head -20 test_output_genes.fasta
    echo "..."
else
    echo "No gene sequences file generated"
fi
echo ""

echo "================================================================================"
echo "Test complete!"
echo "================================================================================"
echo ""
echo "Generated files:"
ls -lh test_output_* test_assembly.fasta
echo ""
echo "You can now use the tool with your own assemblies:"
echo "  python3 resistance_detector.py -a your_assembly.fasta -d database.fasta -o output"
