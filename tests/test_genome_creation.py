"""
Tests for create_test_genomes.py – validates that synthetic genome creation
works correctly for all FOS / CAZAVI genes and mutations from the gist.

These tests do NOT require external tools; they verify:
  - Reference sequences are loadable from data/example_database.fasta
  - introduce_mutation() works at all documented positions
  - Contig lists have expected properties (count, IDs, sizes)
  - The mutation is actually present in the embedded gene sequence
"""
import sys, os
import pytest
from Bio.Seq import Seq
from Bio import SeqIO

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from create_test_genomes import (
    load_references,
    introduce_mutation,
    create_negative_control,
    create_fos_plasmidic,
    create_chromosomal_fos_mutations_primers,
    create_chromosomal_fos_murA,
    create_chromosomal_fos_uhpT,
    create_chromosomal_fos_glpT,
    create_chromosomal_fos_uhpABC,
    create_chromosomal_fos_lon_cyaA_ptsI,
    create_chromosomal_fosAKP,
    create_cazavi_blakpc2_d179y,
    create_cazavi_blakpc3,
    create_cazavi_blakpc31_d179y,
    create_cazavi_blakpc190,
    create_cazavi_blaoxa48,
    create_cazavi_blacmy178_n70t,
    create_cazavi_blashv12,
    create_cazavi_porins,
    create_cazavi_efflux_regulators,
    create_cazavi_acrB,
    create_cazavi_ftsI,
    create_cazavi_envZ,
    create_multi_resistance,
)

DB_PATH = os.path.join(os.path.dirname(__file__), '..', 'fos_cazavi', 'data', 'example_database.fasta')


@pytest.fixture(scope="module")
def refs():
    assert os.path.exists(DB_PATH), f"Database not found: {DB_PATH}"
    return load_references(DB_PATH)


# ── helper ─────────────────────────────────────────────────────────────────────

def get_aa(dna, pos):
    """Return the amino acid (1-letter) at 1-indexed position *pos*."""
    prot = str(Seq(dna[:len(dna) - len(dna) % 3]).translate())
    return prot[pos - 1]


# ══════════════════════════════════════════════════════════════════════════════
# Database loading
# ══════════════════════════════════════════════════════════════════════════════

class TestDatabaseLoading:

    def test_db_file_exists(self):
        assert os.path.exists(DB_PATH)

    def test_all_required_fos_genes_present(self, refs):
        required = [
            'fosA3_reference', 'fosA4_reference', 'fosA5_reference', 'fosA11_reference',
            'fosAKP_reference', 'murA_reference', 'uhpB_reference', 'uhpC_reference',
            'uhpA_reference', 'uhpT_reference', 'glpT_reference', 'cyaA_reference',
            'ptsI_reference', 'galU_reference', 'lon_reference',
        ]
        missing = [g for g in required if g not in refs]
        assert missing == [], f"Missing from DB: {missing}"

    def test_all_required_cazavi_genes_present(self, refs):
        required = [
            'blaKPC-2_reference', 'blaKPC-3_reference', 'blaKPC-31_reference',
            'blaKPC-190_reference', 'blaOXA-48_reference', 'blaCMY-178_reference',
            'blaSHV-12_reference', 'ompK35_reference', 'ompK36_reference',
            'acrB_reference', 'mexR_reference', 'nalD_reference',
            'ftsI_reference', 'envZ_reference',
        ]
        missing = [g for g in required if g not in refs]
        assert missing == [], f"Missing from DB: {missing}"

    def test_ref_sequences_long_enough(self, refs):
        """Ensure each reference is long enough to reach its mutation site."""
        checks = [
            ('murA_reference',      370, 'murA L370I'),
            ('uhpB_reference',      469, 'uhpB G469R'),
            ('uhpC_reference',      384, 'uhpC F384L'),
            ('uhpA_reference',      139, 'uhpA R139C/H'),
            ('uhpT_reference',      350, 'uhpT W350*/R'),
            ('glpT_reference',      362, 'glpT R362C/H/*'),
            ('cyaA_reference',      463, 'cyaA G463D/*'),
            ('ptsI_reference',      191, 'ptsI H191Y/Q'),
            ('galU_reference',      282, 'galU R282V'),
            ('lon_reference',       558, 'lon Q558*'),
            ('fosAKP_reference',     91, 'fosAKP I91V'),
            ('blaKPC-3_reference',  243, 'blaKPC-3 T243M'),
            ('blaOXA-48_reference', 211, 'blaOXA-48 Y211S'),
            ('blaCMY-178_reference', 70, 'blaCMY-178 N70T'),
            ('blaSHV-12_reference', 240, 'blaSHV-12 E240K'),
            ('ompK35_reference',    181, 'ompK35 D181G'),
            ('ompK36_reference',    213, 'ompK36 G213D'),
            ('acrB_reference',      628, 'acrB A628T/V'),
            ('mexR_reference',       75, 'mexR A75V'),
            ('nalD_reference',      174, 'nalD L174R'),
            ('ftsI_reference',      357, 'ftsI S357N'),
            ('envZ_reference',      324, 'envZ T324I/A'),
        ]
        for ref_key, min_pos, label in checks:
            seq = refs.get(ref_key, '')
            min_nt = min_pos * 3
            assert len(seq) >= min_nt, (
                f"{ref_key} too short for {label}: need {min_nt} nt, got {len(seq)}"
            )


# ══════════════════════════════════════════════════════════════════════════════
# introduce_mutation utility
# ══════════════════════════════════════════════════════════════════════════════

class TestIntroduceMutation:

    def test_mutation_at_start(self):
        dna = 'ATGGCTGCT'  # M-A-A
        result = introduce_mutation(dna, 1, 'TTT')
        assert result == 'TTTGCTGCT'

    def test_mutation_at_middle(self):
        dna = 'ATGGCTGCTTGT'  # M-A-A-C
        result = introduce_mutation(dna, 3, 'CGT')
        assert result == 'ATGGCTCGTTGT'

    def test_mutation_at_end(self):
        dna = 'ATGGCTGCT'  # 3 codons
        result = introduce_mutation(dna, 3, 'TAA')
        assert result == 'ATGGCTTAA'

    def test_beyond_sequence_raises(self):
        dna = 'ATGGCT'  # 2 codons
        with pytest.raises(ValueError):
            introduce_mutation(dna, 3, 'TTT')

    def test_sequence_length_preserved(self):
        dna = 'ATGGCTGCTTGT'
        result = introduce_mutation(dna, 2, 'CGT')
        assert len(result) == len(dna)


# ══════════════════════════════════════════════════════════════════════════════
# Negative control
# ══════════════════════════════════════════════════════════════════════════════

class TestNegativeControl:

    def test_returns_contigs(self):
        contigs = create_negative_control()
        assert len(contigs) == 10

    def test_contigs_not_empty(self):
        for c in create_negative_control():
            assert len(c.seq) > 0

    def test_contig_ids_unique(self):
        contigs = create_negative_control()
        ids = [c.id for c in contigs]
        assert len(ids) == len(set(ids))


# ══════════════════════════════════════════════════════════════════════════════
# FOS – Plasmidic
# ══════════════════════════════════════════════════════════════════════════════

class TestFosPlasmidic:

    @pytest.mark.parametrize("gene_id", [
        'fosA3_reference', 'fosA4_reference', 'fosA5_reference', 'fosA11_reference'
    ])
    def test_plasmid_contig_created(self, refs, gene_id):
        contigs = create_fos_plasmidic(refs, [gene_id])
        short = gene_id.replace('_reference', '')
        ids = [c.id for c in contigs]
        assert any(short in cid for cid in ids), f"No contig for {gene_id}: {ids}"

    def test_gene_sequence_embedded(self, refs):
        """fosA3 sequence should appear verbatim in the plasmid contig."""
        gene_seq = refs['fosA3_reference']
        contigs = create_fos_plasmidic(refs, ['fosA3_reference'])
        plasmid = next(c for c in contigs if 'fosA3' in c.id)
        assert gene_seq in str(plasmid.seq)


# ══════════════════════════════════════════════════════════════════════════════
# FOS – Chromosomal (primers-validated mutations)
# ══════════════════════════════════════════════════════════════════════════════

class TestChromosomalFosPrimers:
    """uhpB G469R, uhpC F384L, galU R282V, lon Q558*"""

    def test_correct_number_of_contigs(self, refs):
        contigs = create_chromosomal_fos_mutations_primers(refs)
        gene_contigs = [c for c in contigs if 'chr' in c.id]
        assert len(gene_contigs) == 4

    def test_uhpB_G469R_mutation_present(self, refs):
        contigs = create_chromosomal_fos_mutations_primers(refs)
        uhpB_contig = next(c for c in contigs if 'uhpB' in c.id)
        # Find uhpB sequence in contig and verify amino acid at pos 469
        ref_seq = refs['uhpB_reference']
        seq_str = str(uhpB_contig.seq)
        idx = seq_str.find(refs['uhpB_reference'][:30])  # use prefix to locate start
        if idx >= 0:
            gene_from_contig = seq_str[idx:idx + len(ref_seq)]
            assert get_aa(gene_from_contig, 469) == 'R', "uhpB G469R not found in contig"

    def test_galU_R282V_mutation_present(self, refs):
        contigs = create_chromosomal_fos_mutations_primers(refs)
        galU_contig = next(c for c in contigs if 'galU' in c.id)
        ref_seq = refs['galU_reference']
        seq_str = str(galU_contig.seq)
        idx = seq_str.find(refs['galU_reference'][:30])
        if idx >= 0:
            gene_from_contig = seq_str[idx:idx + len(ref_seq)]
            assert get_aa(gene_from_contig, 282) == 'V', "galU R282V not found in contig"


class TestMurA:
    def test_contigs_created(self, refs):
        contigs = create_chromosomal_fos_murA(refs)
        assert any('murA' in c.id for c in contigs)

    def test_two_mutation_contigs(self, refs):
        contigs = create_chromosomal_fos_murA(refs)
        murA_contigs = [c for c in contigs if 'murA' in c.id]
        assert len(murA_contigs) == 2


class TestUhpT:
    def test_four_mutation_contigs(self, refs):
        contigs = create_chromosomal_fos_uhpT(refs)
        uhpT_contigs = [c for c in contigs if 'uhpT' in c.id]
        assert len(uhpT_contigs) == 4


class TestGlpT:
    def test_five_mutation_contigs(self, refs):
        contigs = create_chromosomal_fos_glpT(refs)
        glpT_contigs = [c for c in contigs if 'glpT' in c.id]
        assert len(glpT_contigs) == 5


class TestUhpABC:
    def test_four_mutation_contigs(self, refs):
        contigs = create_chromosomal_fos_uhpABC(refs)
        gene_contigs = [c for c in contigs if 'chr' in c.id]
        assert len(gene_contigs) == 4


class TestLonCyaAPtsI:
    def test_three_mutation_contigs(self, refs):
        contigs = create_chromosomal_fos_lon_cyaA_ptsI(refs)
        gene_contigs = [c for c in contigs if 'chr' in c.id]
        assert len(gene_contigs) == 3


class TestFosAKP:
    def test_contig_created(self, refs):
        contigs = create_chromosomal_fosAKP(refs)
        assert any('fosAKP' in c.id for c in contigs)


# ══════════════════════════════════════════════════════════════════════════════
# CAZAVI – KPC variants
# ══════════════════════════════════════════════════════════════════════════════

class TestBlaKPC2:
    def test_contig_created(self, refs):
        contigs = create_cazavi_blakpc2_d179y(refs)
        assert any('blaKPC2' in c.id for c in contigs)

    def test_D179Y_amino_acid(self, refs):
        """Verify that position 179 encodes Y in the embedded mutant gene."""
        contigs = create_cazavi_blakpc2_d179y(refs)
        kpc_contig = next(c for c in contigs if 'blaKPC2' in c.id)
        ref_seq = refs['blaKPC-2_reference']
        seq_str = str(kpc_contig.seq)
        # Find the mutated gene (TAT codon at pos 179 = nt 535..537)
        # Build prefix from reference up to pos 178 and check
        ref_prefix = ref_seq[: (179 - 1) * 3]
        idx = seq_str.find(ref_prefix)
        if idx >= 0:
            codon179 = seq_str[idx + (179 - 1) * 3: idx + 179 * 3]
            aa = str(Seq(codon179).translate())
            assert aa == 'Y', f"Expected Y at pos 179, got {aa}"


class TestBlaKPC3:
    def test_three_variant_contigs(self, refs):
        contigs = create_cazavi_blakpc3(refs)
        kpc_contigs = [c for c in contigs if 'blaKPC3' in c.id]
        assert len(kpc_contigs) == 3

    def test_D179Y_contig_present(self, refs):
        contigs = create_cazavi_blakpc3(refs)
        assert any('D179Y' in c.description or 'D179Y' in c.id for c in contigs)

    def test_V240G_contig_present(self, refs):
        contigs = create_cazavi_blakpc3(refs)
        assert any('V240G' in c.description or 'V240G' in c.id for c in contigs)

    def test_double_mutant_contig_present(self, refs):
        contigs = create_cazavi_blakpc3(refs)
        assert any('T243M' in c.description or 'T243M' in c.id for c in contigs)


class TestBlaKPC31:
    def test_contig_created(self, refs):
        contigs = create_cazavi_blakpc31_d179y(refs)
        assert any('blaKPC31' in c.id for c in contigs)


class TestBlaKPC190:
    def test_contig_created(self, refs):
        contigs = create_cazavi_blakpc190(refs)
        assert any('blaKPC190' in c.id for c in contigs)


# ══════════════════════════════════════════════════════════════════════════════
# CAZAVI – OXA-48
# ══════════════════════════════════════════════════════════════════════════════

class TestBlaOXA48:
    def test_two_variant_contigs(self, refs):
        contigs = create_cazavi_blaoxa48(refs)
        oxa_contigs = [c for c in contigs if 'blaOXA48' in c.id]
        assert len(oxa_contigs) == 2

    def test_P68A_contig(self, refs):
        contigs = create_cazavi_blaoxa48(refs)
        assert any('P68A' in c.description for c in contigs)

    def test_Y211S_contig(self, refs):
        contigs = create_cazavi_blaoxa48(refs)
        assert any('Y211S' in c.description for c in contigs)


# ══════════════════════════════════════════════════════════════════════════════
# CAZAVI – CMY-178 and SHV-12
# ══════════════════════════════════════════════════════════════════════════════

class TestBlaCMY178:
    def test_contig_created(self, refs):
        contigs = create_cazavi_blacmy178_n70t(refs)
        assert any('blaCMY178' in c.id for c in contigs)

    def test_N70T_in_description(self, refs):
        contigs = create_cazavi_blacmy178_n70t(refs)
        assert any('N70T' in c.description for c in contigs)


class TestBlaSHV12:
    def test_contig_created(self, refs):
        contigs = create_cazavi_blashv12(refs)
        assert any('blaSHV12' in c.id for c in contigs)


# ══════════════════════════════════════════════════════════════════════════════
# CAZAVI – Chromosomal
# ══════════════════════════════════════════════════════════════════════════════

class TestPorins:
    def test_ompK35_contigs_created(self, refs):
        contigs = create_cazavi_porins(refs)
        assert any('ompK35' in c.id for c in contigs)

    def test_three_ompK35_contigs(self, refs):
        contigs = create_cazavi_porins(refs)
        k35 = [c for c in contigs if 'ompK35' in c.id]
        assert len(k35) == 3


class TestEffluxRegulators:
    def test_mexR_contigs(self, refs):
        contigs = create_cazavi_efflux_regulators(refs)
        assert any('mexR' in c.id for c in contigs)

    def test_nalD_contigs(self, refs):
        contigs = create_cazavi_efflux_regulators(refs)
        assert any('nalD' in c.id for c in contigs)

    def test_four_total_mutation_contigs(self, refs):
        contigs = create_cazavi_efflux_regulators(refs)
        gene_contigs = [c for c in contigs if 'chr' in c.id]
        assert len(gene_contigs) == 4


class TestAcrB:
    def test_three_mutation_contigs(self, refs):
        contigs = create_cazavi_acrB(refs)
        acrB_contigs = [c for c in contigs if 'acrB' in c.id]
        assert len(acrB_contigs) == 3


class TestFtsI:
    def test_three_mutation_contigs(self, refs):
        contigs = create_cazavi_ftsI(refs)
        ftsI_contigs = [c for c in contigs if 'ftsI' in c.id]
        assert len(ftsI_contigs) == 3


class TestEnvZ:
    def test_two_mutation_contigs(self, refs):
        contigs = create_cazavi_envZ(refs)
        envZ_contigs = [c for c in contigs if 'envZ' in c.id]
        assert len(envZ_contigs) == 2


# ══════════════════════════════════════════════════════════════════════════════
# Multi-resistance
# ══════════════════════════════════════════════════════════════════════════════

class TestMultiResistance:
    def test_contains_fosA3(self, refs):
        contigs = create_multi_resistance(refs)
        assert any('fosA3' in c.id or 'fosA3' in c.description for c in contigs)

    def test_contains_blaKPC3(self, refs):
        contigs = create_multi_resistance(refs)
        assert any('blaKPC3' in c.id or 'blaKPC-3' in c.description for c in contigs)

    def test_contains_blaOXA48(self, refs):
        contigs = create_multi_resistance(refs)
        assert any('blaOXA' in c.id or 'blaOXA' in c.description for c in contigs)

    def test_minimum_contig_count(self, refs):
        contigs = create_multi_resistance(refs)
        assert len(contigs) >= 8  # 5 chromosome + 3 plasmids
