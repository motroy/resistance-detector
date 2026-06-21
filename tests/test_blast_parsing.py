"""
Unit tests for BlastDetector.parse_blast_output() and extract_gene_name().
No external tools required – BLAST output is provided as synthetic strings.
"""
import sys, os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from fos_cazavi.acquired import BlastDetector


# ── helper ─────────────────────────────────────────────────────────────────────

def make_detector(min_identity=90.0, min_coverage=80.0):
    """Return a BlastDetector wired to a fake assembly/db (not used in parsing)."""
    d = BlastDetector.__new__(BlastDetector)
    d.assembly = "fake.fasta"
    d.database = "fake_db.fasta"
    d.output_prefix = "fake_out"
    d.min_identity = min_identity
    d.min_coverage = min_coverage
    d.results = []
    d.detected_genes = []
    d.mutation_db = {}
    return d


def blast_line(query, subject, pident, length, qs, qe, qlen, ss, se, slen):
    return f"{query}\t{subject}\t{pident}\t{length}\t{qs}\t{qe}\t{qlen}\t{ss}\t{se}\t{slen}"


# ══════════════════════════════════════════════════════════════════════════════

class TestParseBlastOutput:

    def test_perfect_hit_passes(self):
        d = make_detector()
        line = blast_line("contig1", "fosA3_reference", 100.0, 576, 1, 576, 576, 1, 576, 576)
        hits = d.parse_blast_output(line)
        assert len(hits) == 1
        assert hits[0]['identity'] == 100.0
        assert abs(hits[0]['coverage'] - 100.0) < 0.1

    def test_high_identity_high_coverage_passes(self):
        d = make_detector(min_identity=90, min_coverage=80)
        # 92% identity, 85% coverage (length=870, slen=1024)
        line = blast_line("c1", "blaKPC-3_reference", 92.0, 870, 1, 870, 900, 1, 870, 1024)
        hits = d.parse_blast_output(line)
        assert len(hits) == 1

    def test_low_identity_filtered(self):
        d = make_detector(min_identity=90)
        line = blast_line("c1", "blaKPC-3_reference", 80.0, 1035, 1, 1035, 1035, 1, 1035, 1035)
        hits = d.parse_blast_output(line)
        assert hits == []

    def test_low_coverage_filtered(self):
        d = make_detector(min_coverage=80)
        # 50% coverage: length=100, slen=200
        line = blast_line("c1", "fosA3_reference", 99.0, 100, 1, 100, 576, 1, 100, 200)
        hits = d.parse_blast_output(line)
        assert hits == []

    def test_multiple_hits_filtered_correctly(self):
        d = make_detector(min_identity=90, min_coverage=80)
        lines = "\n".join([
            blast_line("c1", "fosA3_reference",  100.0, 576,  1,  576, 576, 1, 576, 576),   # pass
            blast_line("c2", "blaKPC-3_reference", 85.0, 1035, 1, 1035,1035, 1,1035,1035),  # fail id
            blast_line("c3", "blaOXA-48_reference", 95.0, 200, 1,  200, 200, 1, 200,1194),  # fail cov
        ])
        hits = d.parse_blast_output(lines)
        assert len(hits) == 1
        assert hits[0]['query_id'] == "c1"

    def test_empty_output(self):
        d = make_detector()
        hits = d.parse_blast_output("")
        assert hits == []

    def test_contig_start_end_stored(self):
        d = make_detector()
        line = blast_line("chr1", "uhpB_reference", 99.0, 951, 100, 1050, 2000, 1, 951, 951)
        hits = d.parse_blast_output(line)
        assert hits[0]['qstart'] == 100
        assert hits[0]['qend'] == 1050

    def test_gene_name_extracted(self):
        d = make_detector()
        line = blast_line("c1", "blaKPC-2_reference", 100.0, 1035, 1, 1035, 1035, 1, 1035, 1035)
        hits = d.parse_blast_output(line)
        assert hits[0]['gene'] != ""   # gene name extracted

    @pytest.mark.parametrize("subject,expected_prefix", [
        ("fosA3_reference",    "fosA3"),
        ("blaKPC-3_reference", "blaKPC"),
        ("blaOXA-48_reference","blaOXA"),
        ("uhpB_reference",     "uhpB"),
        ("blaCMY-178_reference","blaCMY"),
        ("blaSHV-12_reference", "blaSHV"),
        ("ompK35_reference",    "ompK35"),
        ("mexR_reference",      "mexR"),
        ("nalD_reference",      "nalD"),
    ])
    def test_gene_name_extraction(self, subject, expected_prefix):
        d = make_detector()
        line = blast_line("c1", subject, 100.0, 1035, 1, 1035, 1035, 1, 1035, 1035)
        hits = d.parse_blast_output(line)
        assert hits[0]['gene'].startswith(expected_prefix)

    def test_redundant_hits_same_location(self):
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

        # Should return only the top match (Hit 2)
        assert len(hits) == 1
        assert hits[0]['subject_id'] == "geneA_var2"

    def test_redundant_hits_overlapping_location(self):
        d = make_detector()
        # Hit 1: 100-1100 (Best match for this region)
        line1 = blast_line("contig1", "geneA", 100.0, 1000, 100, 1100, 5000, 1, 1000, 1000)

        # Hit 2: 105-1095 (Included in Hit 1, lower coverage/identity potentially)
        line2 = blast_line("contig1", "geneA_frag", 99.0, 990, 105, 1095, 5000, 1, 990, 990)

        output = "\n".join([line1, line2])
        hits = d.parse_blast_output(output)

        assert len(hits) == 1
        assert hits[0]['subject_id'] == "geneA"

    def test_distinct_locations_kept(self):
        d = make_detector()
        # Hit 1: 100-1100
        line1 = blast_line("contig1", "geneA", 100.0, 1000, 100, 1100, 5000, 1, 1000, 1000)

        # Hit 2: 2000-3000 (Distinct location)
        line2 = blast_line("contig1", "geneB", 100.0, 1000, 2000, 3000, 5000, 1, 1000, 1000)

        output = "\n".join([line1, line2])
        hits = d.parse_blast_output(output)

        assert len(hits) == 2


class TestCopyNumber:

    def test_multiple_loci_same_gene_get_copy_number(self):
        d = make_detector()
        d.extract_hit_sequence = lambda hit: "ATGC"
        # Same gene detected at two distinct loci on contig1
        line1 = blast_line("contig1", "fosA3_reference", 100.0, 576, 1, 576, 576, 1, 576, 576)
        line2 = blast_line("contig1", "fosA3_reference", 95.0, 576, 2000, 2576, 576, 1, 576, 576)
        hits = d.parse_blast_output("\n".join([line1, line2]))
        d.analyze_hits(hits)

        assert len(d.results) == 2
        assert all(r['copy_number'] == 2 for r in d.results)

    def test_distinct_genes_get_copy_number_one(self):
        d = make_detector()
        d.extract_hit_sequence = lambda hit: "ATGC"
        line1 = blast_line("contig1", "geneA", 100.0, 1000, 1, 1000, 1000, 1, 1000, 1000)
        line2 = blast_line("contig1", "geneB", 100.0, 1000, 2000, 3000, 1000, 1, 1000, 1000)
        hits = d.parse_blast_output("\n".join([line1, line2]))
        d.analyze_hits(hits)

        assert len(d.results) == 2
        assert all(r['copy_number'] == 1 for r in d.results)


class TestExtractGeneName:

    def setup_method(self):
        self.d = make_detector()

    def test_simple_id(self):
        assert self.d.extract_gene_name("fosA3") == "fosA3"

    def test_underscore_id(self):
        result = self.d.extract_gene_name("fosA3_reference")
        assert result.startswith("fosA3")

    def test_pipe_delimited_genbank(self):
        sid = "gb|CP001234|1|1|blaKPC-3"
        result = self.d.extract_gene_name(sid)
        assert "blaKPC" in result
