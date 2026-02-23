
"""
Tests for redundant hit filtering in BlastDetector.
"""
import sys, os
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from fos_cazavi.acquired import BlastDetector

def make_detector():
    d = BlastDetector.__new__(BlastDetector)
    d.assembly = "fake.fasta"
    d.database = "fake_db.fasta"
    d.output_prefix = "fake_out"
    d.min_identity = 90
    d.min_coverage = 80
    d.results = []
    d.detected_genes = []
    d.mutation_db = {}
    return d

def blast_line(query, subject, pident, length, qs, qe, qlen, ss, se, slen):
    return f"{query}\t{subject}\t{pident}\t{length}\t{qs}\t{qe}\t{qlen}\t{ss}\t{se}\t{slen}"

def test_redundant_hits_same_location():
    d = make_detector()
    # Three hits at the same location (100-1100 on contig1)
    # Hit 1: 95% identity, 100% coverage
    line1 = blast_line("contig1", "geneA_var1", 95.0, 1000, 100, 1100, 5000, 1, 1000, 1000)
    # Hit 2: 100% identity, 100% coverage (The best one)
    line2 = blast_line("contig1", "geneA_var2", 100.0, 1000, 100, 1100, 5000, 1, 1000, 1000)
    # Hit 3: 98% identity, 100% coverage
    line3 = blast_line("contig1", "geneA_var3", 98.0, 1000, 100, 1100, 5000, 1, 1000, 1000)

    output = "\n".join([line1, line2, line3])

    hits = d.parse_blast_output(output)

    # Currently, this should return 3 hits.
    # The requirement is to return only the top match (Hit 2).
    # So we assert that len(hits) == 1 and hits[0]['subject_id'] == 'geneA_var2'

    assert len(hits) == 1
    assert hits[0]['subject_id'] == "geneA_var2"

def test_redundant_hits_overlapping_location():
    d = make_detector()
    # Hit 1: 100-1100 (Best match for this region)
    line1 = blast_line("contig1", "geneA", 100.0, 1000, 100, 1100, 5000, 1, 1000, 1000)

    # Hit 2: 105-1095 (Included in Hit 1, lower coverage/identity potentially)
    line2 = blast_line("contig1", "geneA_frag", 99.0, 990, 105, 1095, 5000, 1, 990, 990)

    output = "\n".join([line1, line2])
    hits = d.parse_blast_output(output)

    assert len(hits) == 1
    assert hits[0]['subject_id'] == "geneA"

def test_distinct_locations_kept():
    d = make_detector()
    # Hit 1: 100-1100
    line1 = blast_line("contig1", "geneA", 100.0, 1000, 100, 1100, 5000, 1, 1000, 1000)

    # Hit 2: 2000-3000 (Distinct location)
    line2 = blast_line("contig1", "geneB", 100.0, 1000, 2000, 3000, 5000, 1, 1000, 1000)

    output = "\n".join([line1, line2])
    hits = d.parse_blast_output(output)

    assert len(hits) == 2
