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
              'N': 'AAT', 'K': 'AAA', 'A': 'GCT', 'C': 'TGT', 'Q': 'CAA',
              'V': 'GTT', 'T': 'ACT', 'I': 'ATT'}
    index = (position - 1) * 3
    return cds_sequence[:index] + codons[new_residue] + cds_sequence[index + 3:]


def ambler_substitute_kpc(cds_sequence, ambler_position, new_residue):
    return substitute(cds_sequence, ambler_to_sequential(ambler_position), new_residue)


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


class TestContributoryCazaviEvidence:
    """Porin, PBP3 and EnvZ changes raise MICs but are not sufficient alone, so
    they must produce Indeterminate rather than Resistant or Susceptible."""

    def test_envz_r397c_gives_indeterminate(self, tmp_path, database, reference_cds):
        # envZ_R397C is curated by AMRFinderPlus for CEFTAZIDIME-AVIBACTAM.
        sequence = substitute(reference_cds['envZ'], 397, 'C')
        _, results = run(tmp_path, database, {'contig1': embed(sequence)},
                         organism='Klebsiella_pneumoniae')
        envz = find(results, 'envZ')
        assert 'envZ_R397C' in envz['reported_mutations']
        prediction = predict_phenotypes(results)['ceftazidime_avibactam']
        assert prediction['phenotype'] == 'Indeterminate'
        assert any('envZ_R397C' in item for item in prediction['evidence'])

    def test_ftsi_l367q_gives_indeterminate(self, tmp_path, database, reference_cds):
        sequence = substitute(reference_cds['ftsI'], 367, 'Q')
        _, results = run(tmp_path, database, {'contig1': embed(sequence)},
                         organism='Klebsiella_pneumoniae')
        assert 'ftsI_L367Q' in find(results, 'ftsI')['reported_mutations']
        assert (predict_phenotypes(results)['ceftazidime_avibactam']['phenotype']
                == 'Indeterminate')

    def test_ompk36_a21v_gives_indeterminate(self, tmp_path, database, reference_cds):
        sequence = substitute(reference_cds['ompK36'], 21, 'V')
        _, results = run(tmp_path, database, {'contig1': embed(sequence)},
                         organism='Klebsiella_pneumoniae')
        assert 'ompK36_A21V' in find(results, 'ompK36')['reported_mutations']
        assert (predict_phenotypes(results)['ceftazidime_avibactam']['phenotype']
                == 'Indeterminate')

    def test_contributory_change_does_not_override_a_resistant_call(
            self, tmp_path, database, reference_cds):
        # With a KPC escape variant present the call stays Resistant, and the
        # porin change appears as supporting context.
        _, results = run(tmp_path, database, {
            'contig1': embed(substitute(reference_cds['envZ'], 397, 'C')),
            'contig2': embed(ambler_substitute_kpc(reference_cds['blaKPC-2'], 179, 'Y')),
        }, organism='Klebsiella_pneumoniae')
        prediction = predict_phenotypes(results)['ceftazidime_avibactam']
        assert prediction['phenotype'] == 'Resistant'
        assert any('D179Y' in item for item in prediction['evidence'])
        assert any('envZ_R397C' in item for item in prediction['evidence'])

    def test_porin_knockout_with_betalactamase_gives_indeterminate(
            self, tmp_path, database, reference_cds):
        broken_porin = substitute(reference_cds['ompK36'], 60, '*')
        _, results = run(tmp_path, database, {
            'contig1': embed(broken_porin),
            'contig2': embed(reference_cds['blaKPC-2']),
        }, organism='Klebsiella_pneumoniae')
        assert find(results, 'ompK36')['loss_of_function']
        prediction = predict_phenotypes(results)['ceftazidime_avibactam']
        assert prediction['phenotype'] == 'Indeterminate'
        assert any('porin' in item for item in prediction['evidence'])

    def test_porin_knockout_without_betalactamase_is_not_scored(
            self, tmp_path, database, reference_cds):
        # Porin loss amplifies a beta-lactamase; on its own it is not evidence
        # of ceftazidime-avibactam resistance.
        broken_porin = substitute(reference_cds['ompK36'], 60, '*')
        _, results = run(tmp_path, database, {'contig1': embed(broken_porin)},
                         organism='Klebsiella_pneumoniae')
        assert find(results, 'ompK36')['loss_of_function']
        assert (predict_phenotypes(results)['ceftazidime_avibactam']['phenotype']
                == 'Susceptible')

    def test_contributory_genes_are_silent_without_organism(
            self, tmp_path, database, reference_cds):
        sequence = substitute(reference_cds['envZ'], 397, 'C')
        _, results = run(tmp_path, database, {'contig1': embed(sequence)})
        assert find(results, 'envZ')['reported_mutations'] == []
        assert (predict_phenotypes(results)['ceftazidime_avibactam']['phenotype']
                == 'Susceptible')


class TestDrugScope:
    """This tool reports on fosfomycin and ceftazidime-avibactam. A mutation
    curated for a different drug must never be presented as a finding for
    either of them."""

    def test_fosmidomycin_mutation_in_cyaa_is_not_fosfomycin_resistance(
            self, tmp_path, database, reference_cds):
        # cyaA_S352T is curated for FOSMIDOMYCIN - a different drug whose name
        # merely looks similar. cyaA is a fosfomycin uptake/regulatory gene, so
        # an unscoped implementation reports this as fosfomycin resistance.
        sequence = substitute(reference_cds['cyaA'], 352, 'T')
        _, results = run(tmp_path, database, {'contig1': embed(sequence)},
                         organism='Escherichia')
        cyaa = find(results, 'cyaA')
        assert 'cyaA_S352T' not in cyaa['reported_mutations']
        assert any('cyaA_S352T' in item for item in cyaa['other_drug_mutations'])
        prediction = predict_phenotypes(results)['fosfomycin']
        assert prediction['phenotype'] == 'Susceptible'

    def test_cephalosporin_mutation_in_galu_is_not_fosfomycin_resistance(
            self, tmp_path, database, reference_cds):
        # galU_R101C is curated for CEPHALOSPORIN; galU is a fosfomycin gene.
        sequence = substitute(reference_cds['galU'], 101, 'C')
        _, results = run(tmp_path, database, {'contig1': embed(sequence)},
                         organism='Pseudomonas_aeruginosa')
        galu = find(results, 'galU')
        assert 'galU_R101C' not in galu['reported_mutations']
        assert predict_phenotypes(results)['fosfomycin']['phenotype'] == 'Susceptible'

    def test_genuine_fosfomycin_mutation_is_still_scored(
            self, tmp_path, database, reference_cds):
        # murA_L370I is curated for FOSFOMYCIN and must still be reported.
        sequence = substitute(reference_cds['murA'], 370, 'I')
        _, results = run(tmp_path, database, {'contig1': embed(sequence)},
                         organism='Escherichia')
        assert 'murA_L370I' in find(results, 'murA')['reported_mutations']
        prediction = predict_phenotypes(results)['fosfomycin']
        assert prediction['phenotype'] == 'Resistant'
        assert any('murA' in item for item in prediction['evidence'])


class TestIntrinsicFosAEnzymes:
    def test_kp_chromosomal_fosa_is_treated_as_intrinsic(
            self, tmp_path, database, reference_cds):
        # The K. pneumoniae chromosomal enzyme (fosAKP, which AMRFinderPlus
        # calls fosA6) is present in fosfomycin-susceptible isolates.
        _, results = run(tmp_path, database, {'contig1': embed(reference_cds['fosAKP'])},
                         organism='Klebsiella_pneumoniae')
        gene = find(results, 'fosAKP')
        assert gene is not None and not gene['acquired']
        assert predict_phenotypes(results)['fosfomycin']['phenotype'] == 'Susceptible'

    def test_acquired_fosa2_is_resistance(self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database, {'contig1': embed(reference_cds['fosA2'])})
        assert predict_phenotypes(results)['fosfomycin']['phenotype'] == 'Resistant'


class TestIntrinsicVersusAcquiredFosA:
    """Every Klebsiella carries a chromosomal fosA. A lone fosA hit with no
    intrinsic copy recognised cannot be told apart from a divergent chromosomal
    enzyme by identity alone, so it must not be asserted as acquired."""

    def test_acquired_fosa_alongside_intrinsic_is_resistant(
            self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database, {
            'contig1': embed(reference_cds['fosA3']),
            'contig2': embed(reference_cds['fosAKP']),
        }, organism='Klebsiella_pneumoniae')
        prediction = predict_phenotypes(results, organism='Klebsiella_pneumoniae')['fosfomycin']
        assert prediction['phenotype'] == 'Resistant'
        assert any('Acquired' in item for item in prediction['evidence'])

    def test_lone_fosa_without_intrinsic_copy_is_indeterminate(
            self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database, {'contig1': embed(reference_cds['fosA3'])},
                         organism='Klebsiella_pneumoniae')
        prediction = predict_phenotypes(results, organism='Klebsiella_pneumoniae')['fosfomycin']
        assert prediction['phenotype'] == 'Indeterminate'
        assert any('divergent chromosomal enzyme' in item
                   for item in prediction['evidence'])

    def test_lone_fosa_in_ecoli_is_resistant(self, tmp_path, database, reference_cds):
        # E. coli has no intrinsic chromosomal fosA, so there is no ambiguity.
        _, results = run(tmp_path, database, {'contig1': embed(reference_cds['fosA3'])},
                         organism='Escherichia')
        prediction = predict_phenotypes(results, organism='Escherichia')['fosfomycin']
        assert prediction['phenotype'] == 'Resistant'


class TestPseudomonasAeruginosa:
    """P. aeruginosa needs species-level handling: it is intrinsically
    fosfomycin-resistant, and its dominant CAZ/AVI mechanism (PDC/AmpC) is not
    assessed by this tool."""

    def test_fosfomycin_is_intrinsically_resistant(
            self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database,
                         {'contig1': embed(reference_cds['fosA_PA1129'])},
                         organism='Pseudomonas_aeruginosa')
        prediction = predict_phenotypes(
            results, organism='Pseudomonas_aeruginosa')['fosfomycin']
        assert prediction['phenotype'] == 'Resistant'
        assert any('intrinsically resistant' in item for item in prediction['evidence'])

    def test_cazavi_without_a_mechanism_is_indeterminate_not_susceptible(
            self, tmp_path, database, reference_cds):
        # PDC derepression drives most CAZ/AVI resistance in this species and is
        # not assessed, so "susceptible" would overstate what was ruled out.
        _, results = run(tmp_path, database,
                         {'contig1': embed(reference_cds['fosA_PA1129'])},
                         organism='Pseudomonas_aeruginosa')
        prediction = predict_phenotypes(
            results, organism='Pseudomonas_aeruginosa')['ceftazidime_avibactam']
        assert prediction['phenotype'] == 'Indeterminate'
        assert any('PDC' in item for item in prediction['evidence'])

    def test_metallo_betalactamase_still_gives_a_definite_resistant_call(
            self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database,
                         {'contig1': embed(reference_cds['blaVIM-2'])},
                         organism='Pseudomonas_aeruginosa')
        prediction = predict_phenotypes(
            results, organism='Pseudomonas_aeruginosa')['ceftazidime_avibactam']
        assert prediction['phenotype'] == 'Resistant'

    def test_enterobacterales_are_unaffected_by_the_pseudomonas_rules(
            self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database, {'contig1': embed(reference_cds['blaKPC-2'])},
                         organism='Klebsiella_pneumoniae')
        phenotypes = predict_phenotypes(results, organism='Klebsiella_pneumoniae')
        assert phenotypes['ceftazidime_avibactam']['phenotype'] == 'Susceptible'
        assert phenotypes['fosfomycin']['phenotype'] == 'Susceptible'


class TestAlleleNamingHonesty:
    def test_exact_match_keeps_the_allele_name(self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database, {'contig1': embed(reference_cds['blaVIM-2'])})
        assert find(results, 'blaVIM-2')['allele'] == 'blaVIM-2'

    def test_inexact_match_reports_the_family_not_a_guessed_allele(
            self, tmp_path, database, reference_cds):
        # blaIMP alleles are up to 99.7% identical to each other, so the closest
        # reference is not evidence of which allele this actually is.
        altered = substitute(reference_cds['blaIMP-1'], 50, 'Y')
        _, results = run(tmp_path, database, {'contig1': embed(altered)})
        hit = next(r for r in results if r['gene'].startswith('blaIMP'))
        assert hit['allele'] == 'blaIMP-like'


class TestPorinAmplification:
    """Porin loss alongside a beta-lactamase avibactam inhibits is a documented
    route to resistance without any carbapenemase (E. coli E2257 in the CREC
    validation set: OmpF truncation + CMY-2, CZA MIC >128)."""

    def test_porin_loss_with_ampc_is_indeterminate_not_susceptible(
            self, tmp_path, database, reference_cds):
        broken_porin = substitute(reference_cds['ompF'], 257, '*')
        _, results = run(tmp_path, database, {
            'contig1': embed(broken_porin),
            'contig2': embed(reference_cds['blaCMY-2']),
        }, organism='Escherichia')
        assert find(results, 'ompF')['loss_of_function']
        prediction = predict_phenotypes(
            results, organism='Escherichia')['ceftazidime_avibactam']
        assert prediction['phenotype'] == 'Indeterminate'
        assert any('ompF' in item and 'blaCMY-2' in item
                   for item in prediction['evidence'])

    def test_porin_loss_without_a_betalactamase_is_not_scored(
            self, tmp_path, database, reference_cds):
        broken_porin = substitute(reference_cds['ompF'], 257, '*')
        _, results = run(tmp_path, database, {'contig1': embed(broken_porin)},
                         organism='Escherichia')
        assert (predict_phenotypes(results, organism='Escherichia')
                ['ceftazidime_avibactam']['phenotype'] == 'Susceptible')

    def test_intact_porin_with_ampc_is_susceptible(
            self, tmp_path, database, reference_cds):
        _, results = run(tmp_path, database, {
            'contig1': embed(reference_cds['ompF']),
            'contig2': embed(reference_cds['blaCMY-2']),
        }, organism='Escherichia')
        assert (predict_phenotypes(results, organism='Escherichia')
                ['ceftazidime_avibactam']['phenotype'] == 'Susceptible')

    def test_acquired_fosa3_in_ecoli_is_unambiguously_resistant(
            self, tmp_path, database, reference_cds):
        # E. coli has no intrinsic chromosomal fosA, so a fosA3 hit is acquired.
        _, results = run(tmp_path, database, {'contig1': embed(reference_cds['fosA3'])},
                         organism='Escherichia')
        prediction = predict_phenotypes(results, organism='Escherichia')['fosfomycin']
        assert prediction['phenotype'] == 'Resistant'
