"""Unit tests for the protein-level variant caller."""

from Bio.Seq import Seq

from fos_cazavi.variants import (
    TRUNCATION_TOLERANCE, ambler_to_sequential, call_variants, collapse_runs,
    compare_proteins, extract_gene_span, sequential_to_ambler, translate_cds,
)

CODON = {
    'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT',
    'H': 'CAT', 'I': 'ATT', 'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAT',
    'P': 'CCT', 'Q': 'CAA', 'R': 'CGT', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT',
    'W': 'TGG', 'Y': 'TAT', '*': 'TAA',
}


def cds(protein):
    return ''.join(CODON[aa] for aa in protein) + CODON['*']


class TestAmblerNumbering:
    def test_positions_58_and_253_do_not_exist(self):
        # The standard class A numbering skips these two positions, which is
        # why sequential indices drift from published labels.
        assert 58 not in [sequential_to_ambler(p) for p in range(1, 300)]
        assert 253 not in [sequential_to_ambler(p) for p in range(1, 300)]

    def test_known_anchors(self):
        assert sequential_to_ambler(57) == 57
        assert sequential_to_ambler(58) == 59
        assert sequential_to_ambler(251) == 252
        assert sequential_to_ambler(252) == 254

    def test_round_trip(self):
        for position in range(1, 300):
            ambler = sequential_to_ambler(position)
            assert ambler_to_sequential(ambler) == position

    def test_absent_positions_map_to_none(self):
        assert ambler_to_sequential(58) is None
        assert ambler_to_sequential(253) is None


class TestTranslate:
    def test_normal_cds(self):
        protein, info = translate_cds(cds('MKATYG'))
        assert protein == 'MKATYG'
        assert info['in_frame'] and info['has_stop_codon']

    def test_internal_stop_truncates(self):
        protein, info = translate_cds(cds('MKA') [:-3] + CODON['*'] + cds('TYG'))
        assert protein == 'MKA'
        assert info['has_stop_codon']

    def test_out_of_frame_sequence_is_flagged(self):
        _, info = translate_cds(cds('MKATYG') + 'AT')
        assert not info['in_frame']
        assert info['trailing_bases'] == 2


class TestCompareProteins:
    def test_substitution(self):
        changes = compare_proteins('MKDTYG', 'MKYTYG')
        assert [c['label'] for c in changes] == ['D3Y']

    def test_deletion_does_not_shift_downstream_calls(self):
        # The residue after the deletion is unchanged; a fixed-position lookup
        # would report a substitution for every position after the gap.
        changes = compare_proteins('MKDEFGHIK', 'MKEFGHIK')
        assert [c['label'] for c in changes] == ['D3del']

    def test_insertion_is_reported_as_an_insertion(self):
        changes = compare_proteins('MKDEFG', 'MKDPQEFG')
        assert [c['label'] for c in changes] == ['insPQ@3']

    def test_ambler_numbering_is_applied(self):
        reference = 'M' * 300
        query = 'M' * 178 + 'Y' + 'M' * 121
        changes = compare_proteins(reference, query, numbering='ambler')
        # Sequential residue 179 is Ambler 180; the label must say so.
        assert [c['label'] for c in changes] == ['M180Y']


class TestCollapseRuns:
    def test_adjacent_deletions_merge(self):
        changes = compare_proteins('MKELNSAIP', 'MKNSAIP')
        collapsed = collapse_runs(changes)
        assert [c['label'] for c in collapsed] == ['E3_L4del']

    def test_separate_deletions_stay_separate(self):
        changes = collapse_runs(compare_proteins('MKELNSAIP', 'MKLNSAP'))
        labels = [c['label'] for c in changes]
        assert 'E3del' in labels and 'I8del' in labels


class TestCallVariants:
    def test_wildtype_has_no_changes(self):
        call = call_variants('MKDTYG', cds('MKDTYG'))
        assert call['changes'] == []
        assert not call['loss_of_function']

    def test_nonsense_mutation_is_loss_of_function(self):
        reference = 'M' + 'K' * 199
        query = 'M' + 'K' * 50
        call = call_variants(reference, cds(query))
        assert call['loss_of_function']
        assert call['premature_stop'] == 52   # 51 residues translated, then the stop

    def test_small_in_frame_deletion_is_not_loss_of_function(self):
        # KPC-66 loses two residues and is still a functioning enzyme.
        reference = 'M' + 'K' * 292
        query = 'M' + 'K' * 290
        call = call_variants(reference, cds(query))
        assert not call['loss_of_function']

    def test_frameshift_is_detected(self):
        call = call_variants('MKDTYG', cds('MKDTYG') + 'AT')
        assert call['frameshift']
        assert call['loss_of_function']

    def test_truncation_tolerance_boundary(self):
        reference = 'M' * 100
        just_inside = 'M' * int(100 * (1 - TRUNCATION_TOLERANCE) + 1)
        assert not call_variants(reference, cds(just_inside))['truncated']
        assert call_variants(reference, cds('M' * 50))['truncated']


class TestExtractGeneSpan:
    def test_partial_alignment_is_extended_to_the_full_gene(self):
        gene = cds('MKDTYGHIKLMNP')
        contig = 'AAAA' + gene + 'TTTT'
        # Pretend BLAST aligned only the middle of the gene.
        sequence, complete = extract_gene_span(
            contig, qstart=5 + 6, qend=4 + len(gene) - 6,
            sstart=7, send=len(gene) - 6, reference_length=len(gene))
        assert complete
        assert sequence == gene

    def test_reverse_strand_is_oriented_to_the_reference(self):
        gene = cds('MKDTYGHIKLMNP')
        contig = 'AAAA' + str(Seq(gene).reverse_complement()) + 'TTTT'
        sequence, complete = extract_gene_span(
            contig, qstart=4 + len(gene), qend=5, sstart=1, send=len(gene),
            reference_length=len(gene))
        assert complete
        assert sequence == gene

    def test_gene_running_off_the_contig_is_flagged_incomplete(self):
        gene = cds('MKDTYGHIKLMNP')
        contig = gene[30:]          # first 30 bases missing from the contig
        _, complete = extract_gene_span(
            contig, qstart=1, qend=len(contig), sstart=31, send=len(gene),
            reference_length=len(gene))
        assert not complete
