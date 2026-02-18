"""
Unit tests for MutationDetector._parse_gamma_codon_changes() and
GAMMA output parsing integration.

Tests the GAMMA .gamma file Codon_Changes field parsing for all CAZAVI
resistance mutations. No external tools are required.
"""
import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from fos_cazavi.mutations import MutationDetector


def make_detector():
    from fos_cazavi.utils import KNOWN_MUTATIONS
    d = MutationDetector.__new__(MutationDetector)
    d.assembly = "fake.fasta"
    d.output_prefix = "fake"
    d.genes_file = None
    d.primers_file = None
    d.gamma_results = []
    d.amplicon_results = []
    d.primers = {}
    d.mutation_db = {
        MutationDetector._normalize_gene_name(k): v
        for k, v in KNOWN_MUTATIONS.items()
    }
    return d


# ── Helpers ────────────────────────────────────────────────────────────────────

def parse(codon_changes):
    """Parse a Codon_Changes string through _parse_gamma_codon_changes."""
    return MutationDetector._parse_gamma_codon_changes(codon_changes)


def filter_known(d, raw_mutations, gene_name):
    """Filter raw mutations through _filter_known_mutations."""
    return d._filter_known_mutations(raw_mutations, gene_name)


# ══════════════════════════════════════════════════════════════════════════════
# _parse_gamma_codon_changes – field parsing
# ══════════════════════════════════════════════════════════════════════════════

class TestParseGammaCodonChanges:

    def test_wildtype_zero_returns_empty(self):
        assert parse('0') == []

    def test_empty_string_returns_empty(self):
        assert parse('') == []

    def test_dash_returns_empty(self):
        assert parse('-') == []

    def test_single_substitution(self):
        assert parse('D179Y') == ['D179Y']

    def test_two_substitutions(self):
        result = parse('D179Y,T243M')
        assert 'D179Y' in result
        assert 'T243M' in result
        assert len(result) == 2

    def test_three_substitutions(self):
        result = parse('G148N,G540S,D749E')
        assert result == ['G148N', 'G540S', 'D749E']

    def test_stop_codon_variant(self):
        assert parse('Q558*') == ['Q558*']

    def test_whitespace_around_entries_stripped(self):
        result = parse(' D179Y , T243M ')
        assert 'D179Y' in result
        assert 'T243M' in result

    def test_non_substitution_entries_skipped(self):
        # Lowercase or indel strings are not standard substitution format
        result = parse('ins100,D179Y,del5aa')
        assert result == ['D179Y']

    def test_wildtype_space_returns_empty(self):
        assert parse('  0  ') == []


# ══════════════════════════════════════════════════════════════════════════════
# KPC mutations – parse + filter
# ══════════════════════════════════════════════════════════════════════════════

class TestKpcGamma:
    """blaKPC GAMMA Codon_Changes parsing for D179Y, V240G, T243M."""

    def test_D179Y_detected(self):
        d = make_detector()
        raw = parse('D179Y')
        muts = filter_known(d, raw, 'blaKPC-3')
        assert 'D179Y/N' in muts

    def test_V240G_detected(self):
        d = make_detector()
        raw = parse('V240G')
        muts = filter_known(d, raw, 'blaKPC-3')
        assert 'V240G' in muts

    def test_T243M_detected(self):
        d = make_detector()
        raw = parse('T243M')
        muts = filter_known(d, raw, 'blaKPC-3')
        assert 'T243M' in muts

    def test_D179Y_T243M_double(self):
        d = make_detector()
        raw = parse('D179Y,T243M')
        muts = filter_known(d, raw, 'blaKPC-3')
        assert 'D179Y/N' in muts
        assert 'T243M' in muts

    def test_wildtype_no_mutations(self):
        d = make_detector()
        raw = parse('0')
        muts = filter_known(d, raw, 'blaKPC-3')
        assert muts == []

    def test_unknown_position_filtered_out(self):
        d = make_detector()
        raw = parse('A100T')
        muts = filter_known(d, raw, 'blaKPC-3')
        assert muts == []


# ══════════════════════════════════════════════════════════════════════════════
# OXA-48 mutations
# ══════════════════════════════════════════════════════════════════════════════

class TestOxa48Gamma:
    """blaOXA-48 GAMMA parsing for P68A and Y211S."""

    def test_P68A_detected(self):
        d = make_detector()
        raw = parse('P68A')
        muts = filter_known(d, raw, 'blaOXA-48')
        assert 'P68A' in muts

    def test_Y211S_detected(self):
        d = make_detector()
        raw = parse('Y211S')
        muts = filter_known(d, raw, 'blaOXA-48')
        assert 'Y211S' in muts

    def test_P68A_Y211S_double(self):
        d = make_detector()
        raw = parse('P68A,Y211S')
        muts = filter_known(d, raw, 'blaOXA-48')
        assert 'P68A' in muts
        assert 'Y211S' in muts


# ══════════════════════════════════════════════════════════════════════════════
# CMY-178 – N70T
# ══════════════════════════════════════════════════════════════════════════════

class TestCmy178Gamma:

    def test_N70T_detected(self):
        d = make_detector()
        raw = parse('N70T')
        muts = filter_known(d, raw, 'blaCMY-178')
        assert 'N70T' in muts


# ══════════════════════════════════════════════════════════════════════════════
# OmpK36 / OmpK35 porin mutations
# ══════════════════════════════════════════════════════════════════════════════

class TestOmpK36Gamma:

    def test_G134D(self):
        d = make_detector()
        raw = parse('G134D')
        muts = filter_known(d, raw, 'OmpK36')
        assert 'G134D/*' in muts

    def test_G213D(self):
        d = make_detector()
        raw = parse('G213D')
        muts = filter_known(d, raw, 'OmpK36')
        assert 'G213D/*' in muts


class TestOmpK35Gamma:

    def test_G134D(self):
        d = make_detector()
        raw = parse('G134D')
        muts = filter_known(d, raw, 'OmpK35')
        assert 'G134D/*' in muts

    def test_D181G(self):
        d = make_detector()
        raw = parse('D181G')
        muts = filter_known(d, raw, 'OmpK35')
        assert 'D181G/*' in muts


# ══════════════════════════════════════════════════════════════════════════════
# AcrB efflux mutations
# ══════════════════════════════════════════════════════════════════════════════

class TestAcrBGamma:

    def test_G617D(self):
        d = make_detector()
        raw = parse('G617D')
        muts = filter_known(d, raw, 'acrB')
        assert 'G617D/N' in muts

    def test_F626L(self):
        d = make_detector()
        raw = parse('F626L')
        muts = filter_known(d, raw, 'acrB')
        assert 'F626L' in muts

    def test_A628T(self):
        d = make_detector()
        raw = parse('A628T')
        muts = filter_known(d, raw, 'acrB')
        assert 'A628T/V' in muts


# ══════════════════════════════════════════════════════════════════════════════
# _filter_known_mutations – gene name normalization
# ══════════════════════════════════════════════════════════════════════════════

class TestFilterKnownMutationsGeneNorm:

    def test_accession_suffix_stripped(self):
        """acrB_WP_002892069.1 normalizes to acrB for lookup."""
        d = make_detector()
        raw = parse('G617D')
        muts = filter_known(d, raw, 'acrB_WP_002892069.1')
        assert 'G617D/N' in muts

    def test_variant_suffix_stripped(self):
        """blaKPC-3 normalizes to blaKPC for lookup."""
        d = make_detector()
        raw = parse('D179Y')
        muts = filter_known(d, raw, 'blaKPC-3')
        assert 'D179Y/N' in muts

    def test_unknown_gene_returns_empty(self):
        d = make_detector()
        raw = parse('A100T')
        muts = filter_known(d, raw, 'unknownGene')
        assert muts == []

    def test_ambiguous_marker_stripped(self):
        """Gene names with ‡ (GAMMA ambiguous marker) are handled."""
        d = make_detector()
        raw = parse('D179Y')
        # Simulate gene name with ambiguous marker stripped (done in run_gamma)
        gene = 'blaKPC-3\u2021'.rstrip('\u2021')
        muts = filter_known(d, raw, gene)
        assert 'D179Y/N' in muts


# ══════════════════════════════════════════════════════════════════════════════
# GAMMA result dict integration
# ══════════════════════════════════════════════════════════════════════════════

class TestGammaResultIntegration:
    """
    Simulate what run_gamma() does when processing a .gamma file row,
    without requiring GAMMA to be installed.
    """

    def _make_gamma_row(self, gene, codon_changes, codon_percent=1.0, percent_length=1.0,
                        contig='contig1', start=1000, stop=2000, match_type='mutant'):
        """Simulate a parsed GAMMA output row."""
        return {
            'Gene': gene,
            'Contig': contig,
            'Start': str(start),
            'Stop': str(stop),
            'Match_Type': match_type,
            'Codon_Changes': codon_changes,
            'Codon_Percent': str(codon_percent),
            'Percent_Length': str(percent_length),
        }

    def _process_row(self, d, row):
        """Process a GAMMA row the same way run_gamma() does."""
        gene_name = row['Gene'].rstrip('\u2021')
        start = int(row['Start'])
        stop = int(row['Stop'])
        codon_changes = row.get('Codon_Changes', '0')
        codon_percent = float(row['Codon_Percent'])
        percent_length = float(row['Percent_Length'])

        raw_mutations = MutationDetector._parse_gamma_codon_changes(codon_changes)
        mutations = d._filter_known_mutations(raw_mutations, gene_name)

        return {
            'protein': gene_name,
            'contig': row['Contig'],
            'contig_start': min(start, stop),
            'contig_end': max(start, stop),
            'identity': codon_percent * 100,
            'coverage': percent_length * 100,
            'mutations': mutations,
            'match_type': row['Match_Type'],
        }

    def test_kpc_mutant_row(self):
        d = make_detector()
        row = self._make_gamma_row('blaKPC-3', 'D179Y', match_type='mutant')
        result = self._process_row(d, row)
        assert 'D179Y/N' in result['mutations']
        assert result['match_type'] == 'mutant'
        assert result['identity'] == 100.0

    def test_wildtype_row_no_mutations(self):
        d = make_detector()
        row = self._make_gamma_row('blaKPC-3', '0', match_type='native')
        result = self._process_row(d, row)
        assert result['mutations'] == []

    def test_identity_coverage_calculation(self):
        d = make_detector()
        row = self._make_gamma_row('blaKPC-3', '0', codon_percent=0.985, percent_length=0.95)
        result = self._process_row(d, row)
        assert abs(result['identity'] - 98.5) < 0.001
        assert abs(result['coverage'] - 95.0) < 0.001

    def test_reverse_strand_coordinates_normalized(self):
        """When Stop < Start (reverse strand), contig_start/end are still min/max."""
        d = make_detector()
        row = self._make_gamma_row('blaKPC-3', '0', start=2000, stop=1000)
        result = self._process_row(d, row)
        assert result['contig_start'] == 1000
        assert result['contig_end'] == 2000

    def test_ambiguous_gene_marker_stripped(self):
        d = make_detector()
        row = self._make_gamma_row('blaKPC-3\u2021', 'D179Y')
        result = self._process_row(d, row)
        assert result['protein'] == 'blaKPC-3'
        assert 'D179Y/N' in result['mutations']
