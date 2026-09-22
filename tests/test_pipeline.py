"""End-to-end tests: build a genome with a known genotype, run the pipeline,
check that the reported genotype and phenotype are the ones that were planted.
"""

import json

import pytest
from Bio.Seq import Seq

from fos_cazavi.acquired import BlastDetector
from fos_cazavi.cli import write_summary
from fos_cazavi.phenotype import predict_phenotypes
from fos_cazavi.variants import ambler_to_sequential

from .conftest import embed, requires_blast, write_genome

pytestmark = requires_blast


def substitute(cds_sequence, position, new_residue, codon_table=None):
    """Replace the codon at 1-based residue ``position``."""
    codons = {'Y': 'TAT', 'P': 'CCT', 'G': 'GGT', 'M': 'ATG', '*': 'TAA',
              'N': 'AAT', 'K': 'AAA', 'A': 'GCT'}
    index = (position - 1) * 3
    return cds_sequence[:index] + codons[new_residue] + cds_sequence[index + 3:]


def delete_residues(cds_sequence, position, count):
    index = (position - 1) * 3
    return cds_sequence[:index] + cds_sequence[index + 3 * count:]


def run(tmp_path, database, contigs, organism=None, min_coverage=80):
    assembly = write_genome(tmp_path / 'genome.fasta', contigs)
    detector = BlastDetector(assembly, database, str(tmp_path / 'out'),
                             min_coverage=min_coverage, organism=organism)
    results = detector.run()
    return assembly, results


def find(results, gene):
    return next((r for r in results if r['gene'] == gene), None)


class TestKpcTyping:
    def test_wildtype_kpc2_is_susceptible(self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database,
                         {'contig1': embed(reference_cds['blaKPC-2'])})
        kpc = find(results, 'blaKPC-2')
        assert kpc['allele'] == 'blaKPC-2'
        assert kpc['changes'] == []
        assert predict_phenotypes(results)['ceftazidime_avibactam']['phenotype'] == 'Susceptible'

    def test_kpc3_background_substitution_is_not_a_resistance_marker(
            self, tmp_path, database, reference_cds):
        # H274Y is what distinguishes KPC-3 from KPC-2. It is not an avibactam
        # escape mutation and must not drive a Resistant call on its own.
        sequence = substitute(reference_cds['blaKPC-2'],
                              ambler_to_sequential(274), 'Y')
        _, results = run(tmp_path, database, {'contig1': embed(sequence)})
        kpc = find(results, 'blaKPC-2')
        assert kpc['allele'] == 'blaKPC-3'
        assert kpc['changes'] == ['H274Y']
        assert kpc['reported_mutations'] == []
        assert predict_phenotypes(results)['ceftazidime_avibactam']['phenotype'] == 'Susceptible'

    def test_d179y_is_called_resistant_and_typed_as_kpc33(
            self, tmp_path, database, reference_cds):
        sequence = substitute(reference_cds['blaKPC-2'],
                              ambler_to_sequential(179), 'Y')
        _, results = run(tmp_path, database, {'contig1': embed(sequence)})
        kpc = find(results, 'blaKPC-2')
        assert kpc['changes'] == ['D179Y']
        assert kpc['allele'] == 'blaKPC-33'
        prediction = predict_phenotypes(results)['ceftazidime_avibactam']
        assert prediction['phenotype'] == 'Resistant'
        assert any('D179Y' in item for item in prediction['evidence'])

    def test_omega_loop_deletion_is_called_resistant(
            self, tmp_path, database, reference_cds):
        # Two residues removed from the Omega loop, as in KPC-66.
        sequence = delete_residues(reference_cds['blaKPC-2'],
                                   ambler_to_sequential(166), 2)
        _, results = run(tmp_path, database, {'contig1': embed(sequence)})
        kpc = find(results, 'blaKPC-2')
        assert kpc['changes'] == ['E166_L167del']
        prediction = predict_phenotypes(results)['ceftazidime_avibactam']
        assert prediction['phenotype'] == 'Resistant'
        assert any('Omega loop' in item for item in prediction['evidence'])

    def test_novel_omega_loop_substitution_is_indeterminate(
            self, tmp_path, database, reference_cds):
        # A change in the hotspot that is not a documented variant must not be
        # asserted either way.
        sequence = substitute(reference_cds['blaKPC-2'],
                              ambler_to_sequential(172), 'P')
        _, results = run(tmp_path, database, {'contig1': embed(sequence)})
        prediction = predict_phenotypes(results)['ceftazidime_avibactam']
        assert prediction['phenotype'] == 'Indeterminate'

    def test_reverse_complemented_gene_gives_the_same_call(
            self, tmp_path, database, reference_cds):
        sequence = substitute(reference_cds['blaKPC-2'],
                              ambler_to_sequential(179), 'Y')
        reverse = str(Seq(embed(sequence)).reverse_complement())
        _, results = run(tmp_path, database, {'contig1': reverse})
        assert find(results, 'blaKPC-2')['changes'] == ['D179Y']


class TestMetalloBetaLactamase:
    def test_ndm_alone_makes_cazavi_resistant(self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database,
                         {'contig1': embed(reference_cds['blaNDM-1'])})
        prediction = predict_phenotypes(results)['ceftazidime_avibactam']
        assert prediction['phenotype'] == 'Resistant'
        assert any('etallo' in item for item in prediction['evidence'])

    def test_ndm_with_wildtype_kpc_still_resistant(self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database, {
            'contig1': embed(reference_cds['blaNDM-1']),
            'contig2': embed(reference_cds['blaKPC-2']),
        })
        assert predict_phenotypes(results)['ceftazidime_avibactam']['phenotype'] == 'Resistant'


class TestFosfomycin:
    def test_acquired_fosa3_is_resistant(self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database, {'contig1': embed(reference_cds['fosA3'])})
        prediction = predict_phenotypes(results)['fosfomycin']
        assert prediction['phenotype'] == 'Resistant'

    def test_intrinsic_fosakp_alone_is_not_resistance(self, tmp_path, database, reference_cds):
        # fosAKP is present in fosfomycin-susceptible K. pneumoniae.
        _, results = run(tmp_path, database, {'contig1': embed(reference_cds['fosAKP'])})
        assert predict_phenotypes(results)['fosfomycin']['phenotype'] == 'Susceptible'

    def test_nonsense_mutation_in_uhpt_is_resistant(self, tmp_path, database, reference_cds):
        # A stop codon a third of the way into uhpT knocks out fosfomycin uptake.
        uhpt = reference_cds['uhpT']
        broken = substitute(uhpt, 150, '*')
        _, results = run(tmp_path, database, {'contig1': embed(broken)},
                         organism='Escherichia')
        gene = find(results, 'uhpT')
        assert gene['loss_of_function']
        prediction = predict_phenotypes(results)['fosfomycin']
        assert prediction['phenotype'] == 'Resistant'
        assert any('uhpT' in item for item in prediction['evidence'])

    def test_intact_uhpt_is_susceptible(self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database,
                         {'contig1': embed(reference_cds['uhpT'])},
                         organism='Escherichia')
        assert not find(results, 'uhpT')['loss_of_function']
        assert predict_phenotypes(results)['fosfomycin']['phenotype'] == 'Susceptible'

    def test_truncation_at_a_contig_boundary_is_indeterminate(
            self, tmp_path, database, reference_cds):
        # The gene is cut by the end of the contig, not by a mutation.
        uhpt = reference_cds['uhpT']
        _, results = run(tmp_path, database, {'contig1': 'ACGT' * 50 + uhpt[:900]},
                         organism='Escherichia', min_coverage=50)
        gene = find(results, 'uhpT')
        assert gene is not None and not gene['complete']
        prediction = predict_phenotypes(results)['fosfomycin']
        # A gene cut off by the assembly is not evidence of a knockout.
        assert prediction['phenotype'] != 'Resistant'
        assert any('contig boundary' in item for item in prediction['evidence'])


class TestOutputs:
    def test_summary_files_are_written_and_parsable(
            self, tmp_path, database, reference_cds):
        assembly, results = run(tmp_path, database,
                                {'contig1': embed(reference_cds['blaKPC-2'])})
        prefix = str(tmp_path / 'report')
        write_summary(prefix, assembly, results, None, None,
                      organism='Klebsiella_pneumoniae')

        with open(f"{prefix}_summary.json") as handle:
            summary = json.load(handle)
        assert summary['organism'] == 'Klebsiella_pneumoniae'
        assert summary['predicted_phenotypes']['ceftazidime_avibactam']['phenotype']

        header = open(f"{prefix}_summary.tsv").readline().rstrip('\n').split('\t')
        assert 'Predicted_Phenotype_Ceftazidime_Avibactam' in header

        text = open(f"{prefix}_summary.txt").read()
        assert 'PREDICTED PHENOTYPES' in text

    def test_copy_number_counts_distinct_loci(self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database, {
            'contig1': embed(reference_cds['blaKPC-2']),
            'contig2': embed(reference_cds['blaKPC-2']),
        })
        kpc = [r for r in results if r['gene'] == 'blaKPC-2']
        assert len(kpc) == 2
        assert all(r['copy_number'] == 2 for r in kpc)


class TestChromosomalMutationsNeedAnOrganism:
    def test_curated_mutations_are_not_reported_without_organism(
            self, tmp_path, database, reference_cds):
        broken = substitute(reference_cds['uhpT'], 150, '*')
        _, results = run(tmp_path, database, {'contig1': embed(broken)})
        gene = find(results, 'uhpT')
        # Loss of function is organism independent and is still detected...
        assert gene['loss_of_function']
        # ...but position-specific curated calls are not made blind.
        assert gene['reported_mutations'] == []
