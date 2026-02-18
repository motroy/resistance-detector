"""
Unit tests for MutationDetector.parse_miniprot_cs().
Tests the CS tag parsing for all CAZAVI protein mutations detectable via miniprot.
No external tools are required.
"""
import sys, os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from fos_cazavi.mutations import MutationDetector


def make_detector():
    from fos_cazavi.utils import KNOWN_MUTATIONS
    d = MutationDetector.__new__(MutationDetector)
    d.assembly = "fake.fasta"
    d.output_prefix = "fake"
    d.proteins_file = None
    d.primers_file = None
    d.miniprot_results = []
    d.amplicon_results = []
    d.primers = {}
    d.mutation_db = KNOWN_MUTATIONS
    return d


# ── CS tag builders ────────────────────────────────────────────────────────────

def identical(n):
    """CS operation for n identical residues."""
    return f":{n}"


def substitution(ref_aa, query_aa):
    """CS substitution operation (lowercase ref, uppercase query)."""
    return f"*{ref_aa.lower()}{query_aa.upper()}"


def build_cs(ops):
    """Concatenate a list of CS operations into a CS string."""
    return ''.join(ops)


# ══════════════════════════════════════════════════════════════════════════════
# KPC mutations
# ══════════════════════════════════════════════════════════════════════════════

class TestKpcMiniprotCS:
    """blaKPC cs-tag parsing for D179Y, V240G, T243M."""

    def test_D179Y_detected(self):
        d = make_detector()
        # 178 identical residues, then D->Y, then the rest
        cs = build_cs([identical(178), substitution('D', 'Y'), identical(167)])
        muts = d.parse_miniprot_cs(cs, 'blaKPC-3')
        assert 'D179Y' in muts

    def test_V240G_detected(self):
        d = make_detector()
        cs = build_cs([identical(239), substitution('V', 'G'), identical(106)])
        muts = d.parse_miniprot_cs(cs, 'blaKPC-3')
        assert 'V240G' in muts

    def test_T243M_detected(self):
        d = make_detector()
        cs = build_cs([identical(242), substitution('T', 'M'), identical(103)])
        muts = d.parse_miniprot_cs(cs, 'blaKPC-3')
        assert 'T243M' in muts

    def test_D179Y_T243M_double(self):
        d = make_detector()
        cs = build_cs([
            identical(178), substitution('D', 'Y'),  # pos 179
            identical(63),  substitution('T', 'M'),  # pos 243
        ])
        muts = d.parse_miniprot_cs(cs, 'blaKPC-3')
        assert 'D179Y' in muts
        assert 'T243M' in muts

    def test_wildtype_no_substitutions(self):
        d = make_detector()
        cs = build_cs([identical(345)])
        muts = d.parse_miniprot_cs(cs, 'blaKPC-3')
        assert muts == []


# ══════════════════════════════════════════════════════════════════════════════
# OXA-48 mutations
# ══════════════════════════════════════════════════════════════════════════════

class TestOxa48MiniprotCS:
    """blaOXA-48 cs-tag parsing for P68A and Y211S."""

    def test_P68A_detected(self):
        d = make_detector()
        cs = build_cs([identical(67), substitution('P', 'A'), identical(330)])
        muts = d.parse_miniprot_cs(cs, 'blaOXA-48')
        assert 'P68A' in muts

    def test_Y211S_detected(self):
        d = make_detector()
        cs = build_cs([identical(210), substitution('Y', 'S'), identical(187)])
        muts = d.parse_miniprot_cs(cs, 'blaOXA-48')
        assert 'Y211S' in muts

    def test_P68A_Y211S_double(self):
        d = make_detector()
        cs = build_cs([
            identical(67),  substitution('P', 'A'),  # pos 68
            identical(142), substitution('Y', 'S'),  # pos 211
        ])
        muts = d.parse_miniprot_cs(cs, 'blaOXA-48')
        assert 'P68A' in muts
        assert 'Y211S' in muts


# ══════════════════════════════════════════════════════════════════════════════
# CMY-178 – N70T
# ══════════════════════════════════════════════════════════════════════════════

class TestCmy178MiniprotCS:

    def test_N70T_detected(self):
        d = make_detector()
        cs = build_cs([identical(69), substitution('N', 'T'), identical(312)])
        muts = d.parse_miniprot_cs(cs, 'blaCMY-178')
        assert 'N70T' in muts


# ══════════════════════════════════════════════════════════════════════════════
# OmpK36 / OmpK35 porin mutations
# ══════════════════════════════════════════════════════════════════════════════

class TestOmpK36MiniprotCS:

    def test_G134D(self):
        d = make_detector()
        cs = build_cs([identical(133), substitution('G', 'D'), identical(295)])
        muts = d.parse_miniprot_cs(cs, 'OmpK36')
        assert 'G134D' in muts

    def test_G213D(self):
        d = make_detector()
        cs = build_cs([identical(212), substitution('G', 'D'), identical(217)])
        muts = d.parse_miniprot_cs(cs, 'OmpK36')
        assert 'G213D' in muts


class TestOmpK35MiniprotCS:

    def test_G134D(self):
        d = make_detector()
        cs = build_cs([identical(133), substitution('G', 'D'), identical(205)])
        muts = d.parse_miniprot_cs(cs, 'OmpK35')
        assert 'G134D' in muts

    def test_D181G(self):
        d = make_detector()
        cs = build_cs([identical(180), substitution('D', 'G'), identical(158)])
        muts = d.parse_miniprot_cs(cs, 'OmpK35')
        assert 'D181G' in muts


# ══════════════════════════════════════════════════════════════════════════════
# AcrB efflux mutations
# ══════════════════════════════════════════════════════════════════════════════

class TestAcrBMiniprotCS:

    def test_G617D(self):
        d = make_detector()
        cs = build_cs([identical(616), substitution('G', 'D'), identical(433)])
        muts = d.parse_miniprot_cs(cs, 'acrB')
        assert 'G617D' in muts

    def test_F626L(self):
        d = make_detector()
        cs = build_cs([identical(625), substitution('F', 'L'), identical(424)])
        muts = d.parse_miniprot_cs(cs, 'acrB')
        assert 'F626L' in muts

    def test_A628T(self):
        d = make_detector()
        cs = build_cs([identical(627), substitution('A', 'T'), identical(422)])
        muts = d.parse_miniprot_cs(cs, 'acrB')
        assert 'A628T' in muts


# ══════════════════════════════════════════════════════════════════════════════
# CS tag edge cases
# ══════════════════════════════════════════════════════════════════════════════

class TestCsTagEdgeCases:

    def test_empty_cs_tag(self):
        d = make_detector()
        muts = d.parse_miniprot_cs("", "blaKPC-3")
        assert muts == []

    def test_only_identical_residues(self):
        d = make_detector()
        muts = d.parse_miniprot_cs(":200", "fosA3")
        assert muts == []

    def test_position_accumulates_correctly(self):
        """Multiple identical blocks should not shift position counting."""
        d = make_detector()
        # 10 identical, 5 identical, D->Y at position 16
        cs = build_cs([identical(10), identical(5), substitution('D', 'Y'), identical(100)])
        muts = d.parse_miniprot_cs(cs, 'protein')
        assert 'D16Y' in muts

    def test_two_substitutions_positions(self):
        d = make_detector()
        cs = build_cs([
            identical(4), substitution('A', 'V'),   # pos 5
            identical(9), substitution('G', 'D'),   # pos 15
        ])
        muts = d.parse_miniprot_cs(cs, 'protein')
        assert 'A5V' in muts
        assert 'G15D' in muts
