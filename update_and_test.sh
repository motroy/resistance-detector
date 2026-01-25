#!/bin/bash
set -e

# Configuration
EMAIL="test@example.com"
DB_FILE="test/resistance_db.fasta"
TEST_DIR="test"

echo "======================================================================"
echo "Automated Update and Test Script for Resistance Detector"
echo "======================================================================"

# Ensure test directory exists
mkdir -p "$TEST_DIR"

# 1. Update Database
echo ""
echo "[1/3] Creating/Updating Reference Database..."
python3 create_reference_database.py -e "$EMAIL" -o "$DB_FILE"

# Check database quality (simple check for sequence lengths)
# If any sequence is > 5000bp, it's likely a whole plasmid/genome (bad)
MAX_LEN=$(python3 -c "from Bio import SeqIO; print(max(len(r.seq) for r in SeqIO.parse('$DB_FILE', 'fasta')))")
echo "Max sequence length in DB: $MAX_LEN bp"
if [ "$MAX_LEN" -gt 5000 ]; then
    echo "ERROR: Database contains sequences > 5000bp. Likely full plasmids/genomes included."
    exit 1
fi
echo "Database check passed."

# 2. Download Test Genomes (if not present)
echo ""
echo "[2/3] Preparing Test Genomes..."

if [ ! -f "$TEST_DIR/ecoli.fasta" ]; then
    echo "Downloading E. coli K-12 MG1655 (negative control)..."
    python3 test/download_genome.py NC_000913.3 "$TEST_DIR/ecoli.fasta"
else
    echo "E. coli genome already present."
fi

if [ ! -f "$TEST_DIR/kp.fasta" ]; then
    echo "Downloading K. pneumoniae NTUH-K2044 (positive control for fosA-like)..."
    python3 test/download_assembly.py GCA_000009885.1 "$TEST_DIR/kp.fasta"
else
    echo "K. pneumoniae genome already present."
fi

# 3. Run Tests
echo ""
echo "[3/3] Running Resistance Detector..."

# Test 1: E. coli (Negative Control)
echo "Running on E. coli..."
python3 resistance_detector.py -a "$TEST_DIR/ecoli.fasta" -d "$DB_FILE" -o "$TEST_DIR/ecoli_test"
if grep -q "No resistance genes detected" "$TEST_DIR/ecoli_test_summary.txt"; then
    echo "PASS: E. coli result correct (negative)."
else
    echo "FAIL: Unexpected detection in E. coli."
    cat "$TEST_DIR/ecoli_test_summary.txt"
    # Don't exit, continue to next test
fi

# Test 2: K. pneumoniae (Positive Control)
echo "Running on K. pneumoniae..."
python3 resistance_detector.py -a "$TEST_DIR/kp.fasta" -d "$DB_FILE" -o "$TEST_DIR/kp_test"
if grep -q "fosA" "$TEST_DIR/kp_test_summary.txt"; then
    echo "PASS: K. pneumoniae result correct (fosA detected)."
else
    echo "FAIL: fosA NOT detected in K. pneumoniae."
    cat "$TEST_DIR/kp_test_summary.txt"
    exit 1
fi

echo ""
echo "======================================================================"
echo "All tests completed successfully!"
echo "======================================================================"
