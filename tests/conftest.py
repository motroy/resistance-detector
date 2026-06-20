"""
Shared fixtures for resistance-detector test suite.

Synthetic gene sequences are built from scratch using standard codons so that
the reference amino acid is guaranteed to be present at each known mutation
position. This keeps tests completely self-contained – no external tools or
large reference files are required.
"""
import sys
import os
import pytest

# Ensure the project root is on the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# ── codon table ────────────────────────────────────────────────────────────────
AA_TO_CODON = {
    'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT',
    'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'K': 'AAA', 'L': 'CTG',
    'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA', 'R': 'CGT',
    'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT',
    '*': 'TAA',
}


def make_gene(length_aa, ref_positions):
    """
    Build a synthetic CDS (nucleotide string) with 'A' (GCT) at every position
    except where *ref_positions* specifies the desired reference amino acid.
    Position 1 is always M (start codon ATG).

    Parameters
    ----------
    length_aa   : total protein length in amino acids
    ref_positions : {1-indexed_aa_pos: single_letter_aa}
    """
    protein = ['A'] * length_aa
    protein[0] = 'M'
    for pos, aa in ref_positions.items():
        if 1 <= pos <= length_aa:
            protein[pos - 1] = aa
    return ''.join(AA_TO_CODON[aa] for aa in protein)


def mutate(dna, codon_pos, new_codon):
    """Replace the codon at *codon_pos* (1-indexed) in *dna*."""
    start = (codon_pos - 1) * 3
    return dna[:start] + new_codon + dna[start + 3:]


# ── per-gene synthetic reference sequences ─────────────────────────────────────

@pytest.fixture
def gene_seqs():
    """
    Return a dict of {gene_name: reference_dna_string} for every gene listed
    in the gist (CAVI + FOS).  Sequences are synthetic but guarantee the
    correct reference amino acid at each mutation site.
    """
    genes = {
        # ── FOS ──────────────────────────────────────────────────────────────
        # Plasmidic
        'fosA3':  make_gene(192, {90: 'K', 119: 'H'}),
        'fosA4':  make_gene(192, {90: 'K', 119: 'H'}),
        'fosA5':  make_gene(192, {90: 'K', 119: 'H'}),
        'fosA7':  make_gene(192, {90: 'K', 119: 'H'}),
        'fosA11': make_gene(192, {90: 'K', 119: 'H'}),
        # Chromosomal – K. pneumoniae
        'fosAKP': make_gene(150, {91: 'I'}),
        # Chromosomal – E. coli
        'murA':  make_gene(420, {369: 'D', 370: 'L'}),
        'uhpB':  make_gene(500, {350: 'H', 469: 'G'}),
        'uhpC':  make_gene(430, {384: 'F'}),
        'uhpA':  make_gene(200, {54: 'D', 139: 'R'}),
        'uhpT':  make_gene(480, {55: 'G', 198: 'W', 258: 'E', 350: 'W'}),
        'glpT':  make_gene(450, {44: 'E', 88: 'W', 90: 'G', 234: 'W', 362: 'R'}),
        'cyaA':  make_gene(850, {463: 'G'}),
        'ptsI':  make_gene(570, {191: 'H'}),
        'galU':  make_gene(320, {282: 'R'}),
        'lon':   make_gene(600, {558: 'Q'}),
        # ── CAZAVI ───────────────────────────────────────────────────────────
        # Plasmidic – KPC family
        # NB: KNOWN_MUTATIONS uses literal sequential CDS positions (verified
        # against RefSeq WP_004152396.1, KPC-3), which sit one residue before
        # the classic literature/Ambler numbering (D179Y -> seq. pos 178, etc).
        'blaKPC-2':  make_gene(345, {166: 'L', 167: 'E', 168: 'L', 169: 'N', 178: 'D', 239: 'V', 242: 'T'}),
        'blaKPC-3':  make_gene(345, {166: 'L', 167: 'E', 168: 'L', 169: 'N', 178: 'D', 239: 'V', 242: 'T'}),
        'blaKPC-31': make_gene(345, {166: 'L', 167: 'E', 168: 'L', 169: 'N', 178: 'D', 239: 'V', 242: 'T'}),
        'blaKPC-190':make_gene(345, {166: 'L', 167: 'E', 168: 'L', 169: 'N', 178: 'D', 239: 'V', 242: 'T'}),
        # OXA-48
        'blaOXA-48': make_gene(398, {68: 'P', 211: 'Y'}),
        # AmpC CMY
        'blaCMY-178':make_gene(382, {70: 'N'}),
        # SHV
        'blaSHV-12': make_gene(286, {238: 'G', 240: 'E'}),
        # Chromosomal – porins
        'ompK36': make_gene(430, {134: 'G', 135: 'D', 213: 'G'}),
        'ompK35': make_gene(340, {134: 'G', 135: 'D', 181: 'D'}),
        # Chromosomal – efflux
        'acrB':  make_gene(1050, {617: 'G', 626: 'F', 628: 'A'}),
        'mexR':  make_gene(147,  {69: 'W', 75: 'A'}),
        'nalD':  make_gene(215,  {153: 'Q', 174: 'L'}),
        # Chromosomal – PBP / sensor
        'ftsI':  make_gene(750,  {333: 'A', 350: 'Y', 357: 'S'}),
        'envZ':  make_gene(450,  {244: 'G', 324: 'T'}),
    }
    return genes
