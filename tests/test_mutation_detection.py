"""
Unit tests for detect_mutations() covering all CAVI/FOS genes listed in the gist:
https://gist.github.com/motroy/c10bfb64aa1f8a5a86bb2d3986c8c724#gistcomment-5951482

Each test:
  1. Takes the reference synthetic gene sequence from conftest.py
  2. Introduces the exact codon change for the target mutation
  3. Calls detect_mutations() and asserts the mutation is reported
  4. Optionally verifies that the wildtype sequence is NOT flagged
"""
import pytest
from fos_cazavi.utils import detect_mutations, KNOWN_MUTATIONS


def mutate(dna, codon_pos, new_codon):
    """Replace the codon at *codon_pos* (1-indexed) in *dna*."""
    start = (codon_pos - 1) * 3
    return dna[:start] + new_codon + dna[start + 3:]


# ── helpers ────────────────────────────────────────────────────────────────────

def check_detected(gene_name, dna, expected_mutation_name):
    """Assert that the expected mutation is reported.

    detect_mutations() returns combined names such as 'H350Y/Q' for any
    variant at position 350.  This helper therefore:
      1. First tries an exact substring match  ('H350Y' in 'H350Y/Q' = True)
      2. If that fails, checks that the position number AND the specific
         variant amino acid both appear in one of the returned names –
         covering cases like 'H350Q' vs 'H350Y/Q' where the variant
         character follows a '/' separator.
    """
    import re
    found = detect_mutations(gene_name, dna, KNOWN_MUTATIONS)
    # 1. Direct substring match
    if any(expected_mutation_name in m for m in found):
        return
    # 2. Position + variant char match (handles combined names)
    m_exp = re.match(r'^([A-Z\*?]?)(\d+)([A-Z\*])$', expected_mutation_name)
    if m_exp:
        pos = m_exp.group(2)
        var = m_exp.group(3)
        for found_name in found:
            if pos in found_name:
                idx = found_name.find(pos)
                suffix = found_name[idx + len(pos):]   # e.g. 'Y/Q' or 'S/D'
                variants_in_name = re.findall(r'[A-Z\*]', suffix)
                if var in variants_in_name:
                    return
    assert False, (
        f"Expected '{expected_mutation_name}' for {gene_name}, got: {found}"
    )


def check_not_detected(gene_name, dna, unexpected_name):
    """Assert that *unexpected_name* does NOT appear (wildtype check)."""
    found = detect_mutations(gene_name, dna, KNOWN_MUTATIONS)
    assert not any(unexpected_name in m for m in found), (
        f"Did not expect '{unexpected_name}' for wildtype {gene_name}, got: {found}"
    )


def check_no_mutations(gene_name, dna):
    """Assert that no resistance mutations are reported for a wildtype sequence."""
    found = detect_mutations(gene_name, dna, KNOWN_MUTATIONS)
    assert found == [], (
        f"Expected no mutations for wildtype {gene_name}, got: {found}"
    )


# ══════════════════════════════════════════════════════════════════════════════
# FOS – Plasmidic genes
# ══════════════════════════════════════════════════════════════════════════════

class TestFosPlasmidic:
    """fosA family – plasmid-borne glutathione S-transferases.

    Each test uses the actual variant gene name (fosA3, fosA4, …) as the
    gene_name argument so that the normalization path fosA3 -> fosA is
    exercised, matching what both BLAST and GAMMA return at runtime.
    """

    # fosA3 / fosA4 / fosA5 / fosA7 / fosA11 share the same mutation positions
    @pytest.mark.parametrize("gene", ["fosA3", "fosA4", "fosA5", "fosA7", "fosA11"])
    def test_wildtype_no_mutations(self, gene_seqs, gene):
        # gene_name is the variant name so that normalization (fosA3 -> fosA) is tested
        check_no_mutations(gene, gene_seqs[gene])

    @pytest.mark.parametrize("gene", ["fosA3", "fosA4", "fosA5", "fosA7", "fosA11"])
    def test_K90E(self, gene_seqs, gene):
        """fosA K90E – loss-of-function variant"""
        dna = mutate(gene_seqs[gene], 90, 'GAA')   # K(AAA) -> E(GAA)
        check_detected(gene, dna, 'K90E')

    @pytest.mark.parametrize("gene", ["fosA3", "fosA4", "fosA5", "fosA7", "fosA11"])
    def test_H119Q(self, gene_seqs, gene):
        """fosA H119Q"""
        dna = mutate(gene_seqs[gene], 119, 'CAA')  # H(CAT) -> Q(CAA)
        check_detected(gene, dna, 'H119Q')


# ══════════════════════════════════════════════════════════════════════════════
# FOS – Chromosomal genes
# ══════════════════════════════════════════════════════════════════════════════

class TestFosAKP:
    """fosAKP – chromosomal fosfomycin resistance in K. pneumoniae."""

    def test_wildtype(self, gene_seqs):
        check_no_mutations('fosAKP', gene_seqs['fosAKP'])

    def test_I91V(self, gene_seqs):
        """fosAKP I91V"""
        dna = mutate(gene_seqs['fosAKP'], 91, 'GTT')  # I(ATT) -> V(GTT)
        check_detected('fosAKP', dna, 'I91V')


class TestMurA:
    """murA – UDP-N-acetylglucosamine enolpyruvyl transferase (FOS target)."""

    def test_wildtype(self, gene_seqs):
        check_no_mutations('murA', gene_seqs['murA'])

    def test_D369N(self, gene_seqs):
        dna = mutate(gene_seqs['murA'], 369, 'AAT')  # D(GAT) -> N(AAT)
        check_detected('murA', dna, 'D369N')

    def test_L370I(self, gene_seqs):
        dna = mutate(gene_seqs['murA'], 370, 'ATT')  # L(CTG) -> I(ATT)
        check_detected('murA', dna, 'L370I')


class TestUhpB:
    """uhpB – sensor histidine kinase of the Uhp two-component system."""

    def test_wildtype(self, gene_seqs):
        check_no_mutations('uhpB', gene_seqs['uhpB'])

    def test_G469R(self, gene_seqs):
        dna = mutate(gene_seqs['uhpB'], 469, 'CGT')  # G(GGT) -> R(CGT)
        check_detected('uhpB', dna, 'G469R')

    def test_H350Y(self, gene_seqs):
        dna = mutate(gene_seqs['uhpB'], 350, 'TAT')  # H(CAT) -> Y(TAT)
        check_detected('uhpB', dna, 'H350Y')

    def test_H350Q(self, gene_seqs):
        dna = mutate(gene_seqs['uhpB'], 350, 'CAA')  # H(CAT) -> Q(CAA)
        check_detected('uhpB', dna, 'H350Q')


class TestUhpC:
    """uhpC – sugar phosphate permease (UhpT regulation)."""

    def test_wildtype(self, gene_seqs):
        check_no_mutations('uhpC', gene_seqs['uhpC'])

    def test_F384L(self, gene_seqs):
        dna = mutate(gene_seqs['uhpC'], 384, 'CTG')  # F(TTT) -> L(CTG)
        check_detected('uhpC', dna, 'F384L')


class TestUhpA:
    """uhpA – transcriptional activator of uhpT."""

    def test_wildtype(self, gene_seqs):
        check_no_mutations('uhpA', gene_seqs['uhpA'])

    def test_D54N(self, gene_seqs):
        dna = mutate(gene_seqs['uhpA'], 54, 'AAT')   # D(GAT) -> N(AAT)
        check_detected('uhpA', dna, 'D54N')

    def test_R139C(self, gene_seqs):
        dna = mutate(gene_seqs['uhpA'], 139, 'TGT')  # R(CGT) -> C(TGT)
        check_detected('uhpA', dna, 'R139C')

    def test_R139H(self, gene_seqs):
        dna = mutate(gene_seqs['uhpA'], 139, 'CAT')  # R(CGT) -> H(CAT)
        check_detected('uhpA', dna, 'R139H')


class TestUhpT:
    """uhpT – hexose phosphate transport protein (main FOS uptake route)."""

    def test_wildtype(self, gene_seqs):
        check_no_mutations('uhpT', gene_seqs['uhpT'])

    def test_G55D(self, gene_seqs):
        dna = mutate(gene_seqs['uhpT'], 55, 'GAT')   # G -> D
        check_detected('uhpT', dna, 'G55D')

    def test_W198_stop(self, gene_seqs):
        dna = mutate(gene_seqs['uhpT'], 198, 'TAA')  # W -> *
        check_detected('uhpT', dna, 'W198')

    def test_E258_stop(self, gene_seqs):
        dna = mutate(gene_seqs['uhpT'], 258, 'TAA')  # E -> *
        check_detected('uhpT', dna, 'E258')

    def test_W350_stop(self, gene_seqs):
        dna = mutate(gene_seqs['uhpT'], 350, 'TAA')  # W -> *
        check_detected('uhpT', dna, 'W350')


class TestGlpT:
    """glpT – glycerol-3-phosphate transporter (alternative FOS uptake)."""

    def test_wildtype(self, gene_seqs):
        check_no_mutations('glpT', gene_seqs['glpT'])

    def test_E44_stop(self, gene_seqs):
        dna = mutate(gene_seqs['glpT'], 44, 'TAA')   # E -> *
        check_detected('glpT', dna, 'E44')

    def test_W88_stop(self, gene_seqs):
        dna = mutate(gene_seqs['glpT'], 88, 'TAA')   # W -> *
        check_detected('glpT', dna, 'W88')

    def test_G90D(self, gene_seqs):
        dna = mutate(gene_seqs['glpT'], 90, 'GAT')   # G -> D
        check_detected('glpT', dna, 'G90D')

    def test_W234_stop(self, gene_seqs):
        dna = mutate(gene_seqs['glpT'], 234, 'TAA')  # W -> *
        check_detected('glpT', dna, 'W234')

    def test_R362C(self, gene_seqs):
        dna = mutate(gene_seqs['glpT'], 362, 'TGT')  # R -> C
        check_detected('glpT', dna, 'R362C')

    def test_R362H(self, gene_seqs):
        dna = mutate(gene_seqs['glpT'], 362, 'CAT')  # R -> H
        check_detected('glpT', dna, 'R362H')


class TestCyaA:
    """cyaA – adenylate cyclase (cAMP signalling, indirectly controls uhpT)."""

    def test_wildtype(self, gene_seqs):
        check_no_mutations('cyaA', gene_seqs['cyaA'])

    def test_G463D(self, gene_seqs):
        dna = mutate(gene_seqs['cyaA'], 463, 'GAT')  # G -> D
        check_detected('cyaA', dna, 'G463D')

    def test_G463_stop(self, gene_seqs):
        dna = mutate(gene_seqs['cyaA'], 463, 'TAA')  # G -> *
        check_detected('cyaA', dna, 'G463')


class TestPtsI:
    """ptsI – enzyme I of the PTS (phosphoenolpyruvate-protein kinase)."""

    def test_wildtype(self, gene_seqs):
        check_no_mutations('ptsI', gene_seqs['ptsI'])

    def test_H191Y(self, gene_seqs):
        dna = mutate(gene_seqs['ptsI'], 191, 'TAT')  # H -> Y
        check_detected('ptsI', dna, 'H191Y')

    def test_H191Q(self, gene_seqs):
        dna = mutate(gene_seqs['ptsI'], 191, 'CAA')  # H -> Q
        check_detected('ptsI', dna, 'H191Q')


class TestGalU:
    """galU – UTP-glucose-1-phosphate uridylyltransferase."""

    def test_wildtype(self, gene_seqs):
        check_no_mutations('galU', gene_seqs['galU'])

    def test_R282V(self, gene_seqs):
        dna = mutate(gene_seqs['galU'], 282, 'GTT')  # R(CGT) -> V(GTT)
        check_detected('galU', dna, 'R282V')


class TestLon:
    """lon – ATP-dependent serine protease."""

    def test_wildtype(self, gene_seqs):
        check_no_mutations('lon', gene_seqs['lon'])

    def test_Q558_stop(self, gene_seqs):
        dna = mutate(gene_seqs['lon'], 558, 'TAA')   # Q -> * (stop)
        check_detected('lon', dna, 'Q558')


# ══════════════════════════════════════════════════════════════════════════════
# CAZAVI – Plasmidic beta-lactamases
# ══════════════════════════════════════════════════════════════════════════════

class TestBlaKPC:
    """blaKPC variants – class A serine carbapenemases."""

    @pytest.mark.parametrize("gene", ["blaKPC-2", "blaKPC-3", "blaKPC-31", "blaKPC-190"])
    def test_wildtype_no_mutations(self, gene_seqs, gene):
        check_no_mutations('blaKPC', gene_seqs[gene])

    @pytest.mark.parametrize("gene", ["blaKPC-2", "blaKPC-3", "blaKPC-31"])
    def test_D179Y(self, gene_seqs, gene):
        """D179Y – key CAZAVI resistance mutation in KPC."""
        dna = mutate(gene_seqs[gene], 179, 'TAT')   # D(GAT) -> Y(TAT)
        check_detected('blaKPC', dna, 'D179Y')

    def test_blaKPC3_V240G(self, gene_seqs):
        """blaKPC-3 V240G"""
        dna = mutate(gene_seqs['blaKPC-3'], 240, 'GGT')  # V(GTT) -> G(GGT)
        check_detected('blaKPC', dna, 'V240G')

    def test_blaKPC3_D179Y_T243M_double(self, gene_seqs):
        """blaKPC-3 D179Y + T243M double mutant"""
        dna = mutate(gene_seqs['blaKPC-3'], 179, 'TAT')
        dna = mutate(dna, 243, 'ATG')               # T(ACT) -> M(ATG)
        muts = detect_mutations('blaKPC', dna, KNOWN_MUTATIONS)
        assert any('D179Y' in m for m in muts), f"D179Y not found: {muts}"
        assert any('T243M' in m for m in muts), f"T243M not found: {muts}"


class TestBlaOXA48:
    """blaOXA-48 – class D OXA-type carbapenemase."""

    def test_wildtype_no_mutations(self, gene_seqs):
        check_no_mutations('blaOXA-48', gene_seqs['blaOXA-48'])

    def test_P68A(self, gene_seqs):
        dna = mutate(gene_seqs['blaOXA-48'], 68, 'GCT')  # P(CCT) -> A(GCT)
        check_detected('blaOXA-48', dna, 'P68A')

    def test_Y211S(self, gene_seqs):
        dna = mutate(gene_seqs['blaOXA-48'], 211, 'TCT')  # Y(TAT) -> S(TCT)
        check_detected('blaOXA-48', dna, 'Y211S')

    def test_P68A_Y211S_double(self, gene_seqs):
        """P68A/Y211S double mutant (CAZAVI resistance combination)."""
        dna = mutate(gene_seqs['blaOXA-48'], 68, 'GCT')
        dna = mutate(dna, 211, 'TCT')
        muts = detect_mutations('blaOXA-48', dna, KNOWN_MUTATIONS)
        assert any('P68A' in m for m in muts), f"P68A not found: {muts}"
        assert any('Y211S' in m for m in muts), f"Y211S not found: {muts}"


class TestBlaCMY178:
    """blaCMY-178 – AmpC-type beta-lactamase with N70T CAZAVI mutation."""

    def test_wildtype_no_mutations(self, gene_seqs):
        check_no_mutations('blaCMY', gene_seqs['blaCMY-178'])

    def test_N70T(self, gene_seqs):
        """N70T (Asn70Thr) – blaCMY-178 CAZAVI resistance mutation."""
        dna = mutate(gene_seqs['blaCMY-178'], 70, 'ACT')  # N(AAT) -> T(ACT)
        check_detected('blaCMY', dna, 'N70T')


class TestBlaSHV12:
    """blaSHV-12 – ESBL; G238S/E240K define the SHV-12 variant."""

    def test_wildtype_no_mutations(self, gene_seqs):
        check_no_mutations('blaSHV', gene_seqs['blaSHV-12'])

    def test_G238S(self, gene_seqs):
        dna = mutate(gene_seqs['blaSHV-12'], 238, 'TCT')  # G(GGT) -> S(TCT)
        check_detected('blaSHV', dna, 'G238S')

    def test_E240K(self, gene_seqs):
        dna = mutate(gene_seqs['blaSHV-12'], 240, 'AAA')  # E(GAA) -> K(AAA)
        check_detected('blaSHV', dna, 'E240K')


# ══════════════════════════════════════════════════════════════════════════════
# CAZAVI – Chromosomal genes
# ══════════════════════════════════════════════════════════════════════════════

class TestOmpK36:
    """ompK36 – outer membrane porin K36 (K. pneumoniae)."""

    def test_wildtype_no_mutations(self, gene_seqs):
        check_no_mutations('ompK36', gene_seqs['ompK36'])

    def test_G134D(self, gene_seqs):
        dna = mutate(gene_seqs['ompK36'], 134, 'GAT')  # G -> D
        check_detected('ompK36', dna, 'G134D')

    def test_G134_stop(self, gene_seqs):
        dna = mutate(gene_seqs['ompK36'], 134, 'TAA')  # G -> *
        check_detected('ompK36', dna, 'G134')

    def test_D135_stop(self, gene_seqs):
        dna = mutate(gene_seqs['ompK36'], 135, 'TAA')  # D -> *
        check_detected('ompK36', dna, 'D135')

    def test_G213D(self, gene_seqs):
        dna = mutate(gene_seqs['ompK36'], 213, 'GAT')  # G -> D
        check_detected('ompK36', dna, 'G213D')


class TestOmpK35:
    """ompK35 – outer membrane porin K35 (K. pneumoniae)."""

    def test_wildtype_no_mutations(self, gene_seqs):
        check_no_mutations('ompK35', gene_seqs['ompK35'])

    def test_G134D(self, gene_seqs):
        dna = mutate(gene_seqs['ompK35'], 134, 'GAT')  # G -> D
        check_detected('ompK35', dna, 'G134D')

    def test_D135_stop(self, gene_seqs):
        dna = mutate(gene_seqs['ompK35'], 135, 'TAA')  # D -> *
        check_detected('ompK35', dna, 'D135')

    def test_D181G(self, gene_seqs):
        dna = mutate(gene_seqs['ompK35'], 181, 'GGT')  # D -> G
        check_detected('ompK35', dna, 'D181G')

    def test_D181_stop(self, gene_seqs):
        dna = mutate(gene_seqs['ompK35'], 181, 'TAA')  # D -> *
        check_detected('ompK35', dna, 'D181')


class TestAcrB:
    """acrB – RND efflux pump transporter (K. pneumoniae)."""

    def test_wildtype_no_mutations(self, gene_seqs):
        check_no_mutations('acrB', gene_seqs['acrB'])

    def test_G617D(self, gene_seqs):
        dna = mutate(gene_seqs['acrB'], 617, 'GAT')  # G -> D
        check_detected('acrB', dna, 'G617D')

    def test_G617N(self, gene_seqs):
        dna = mutate(gene_seqs['acrB'], 617, 'AAT')  # G -> N
        check_detected('acrB', dna, 'G617N')

    def test_F626L(self, gene_seqs):
        dna = mutate(gene_seqs['acrB'], 626, 'CTG')  # F -> L
        check_detected('acrB', dna, 'F626L')

    def test_A628T(self, gene_seqs):
        dna = mutate(gene_seqs['acrB'], 628, 'ACT')  # A -> T
        check_detected('acrB', dna, 'A628T')

    def test_A628V(self, gene_seqs):
        dna = mutate(gene_seqs['acrB'], 628, 'GTT')  # A -> V
        check_detected('acrB', dna, 'A628V')


class TestMexR:
    """mexR – efflux pump repressor."""

    def test_wildtype_no_mutations(self, gene_seqs):
        check_no_mutations('mexR', gene_seqs['mexR'])

    def test_W69G(self, gene_seqs):
        dna = mutate(gene_seqs['mexR'], 69, 'GGT')  # W -> G
        check_detected('mexR', dna, 'W69G')

    def test_W69_stop(self, gene_seqs):
        dna = mutate(gene_seqs['mexR'], 69, 'TAA')  # W -> *
        check_detected('mexR', dna, 'W69')

    def test_A75V(self, gene_seqs):
        dna = mutate(gene_seqs['mexR'], 75, 'GTT')  # A -> V
        check_detected('mexR', dna, 'A75V')


class TestNalD:
    """nalD – efflux pump repressor."""

    def test_wildtype_no_mutations(self, gene_seqs):
        check_no_mutations('nalD', gene_seqs['nalD'])

    def test_Q153_stop(self, gene_seqs):
        dna = mutate(gene_seqs['nalD'], 153, 'TAA')  # Q -> *
        check_detected('nalD', dna, 'Q153')

    def test_L174R(self, gene_seqs):
        dna = mutate(gene_seqs['nalD'], 174, 'CGT')  # L -> R
        check_detected('nalD', dna, 'L174R')


class TestFtsI:
    """ftsI / PBP3 – penicillin-binding protein 3."""

    def test_wildtype_no_mutations(self, gene_seqs):
        check_no_mutations('ftsI', gene_seqs['ftsI'])

    def test_A333V(self, gene_seqs):
        dna = mutate(gene_seqs['ftsI'], 333, 'GTT')  # A -> V
        check_detected('ftsI', dna, 'A333V')

    def test_A333T(self, gene_seqs):
        dna = mutate(gene_seqs['ftsI'], 333, 'ACT')  # A -> T
        check_detected('ftsI', dna, 'A333T')

    def test_Y350C(self, gene_seqs):
        dna = mutate(gene_seqs['ftsI'], 350, 'TGT')  # Y -> C
        check_detected('ftsI', dna, 'Y350C')

    def test_Y350S(self, gene_seqs):
        dna = mutate(gene_seqs['ftsI'], 350, 'TCT')  # Y -> S
        check_detected('ftsI', dna, 'Y350S')

    def test_S357N(self, gene_seqs):
        dna = mutate(gene_seqs['ftsI'], 357, 'AAT')  # S -> N
        check_detected('ftsI', dna, 'S357N')


class TestEnvZ:
    """envZ – two-component system sensor histidine kinase."""

    def test_wildtype_no_mutations(self, gene_seqs):
        check_no_mutations('envZ', gene_seqs['envZ'])

    def test_G244S(self, gene_seqs):
        dna = mutate(gene_seqs['envZ'], 244, 'TCT')  # G -> S
        check_detected('envZ', dna, 'G244S')

    def test_G244D(self, gene_seqs):
        dna = mutate(gene_seqs['envZ'], 244, 'GAT')  # G -> D
        check_detected('envZ', dna, 'G244D')

    def test_T324I(self, gene_seqs):
        dna = mutate(gene_seqs['envZ'], 324, 'ATT')  # T -> I
        check_detected('envZ', dna, 'T324I')

    def test_T324A(self, gene_seqs):
        dna = mutate(gene_seqs['envZ'], 324, 'GCT')  # T -> A
        check_detected('envZ', dna, 'T324A')
