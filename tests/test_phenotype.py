"""
Unit tests for genotype-to-phenotype prediction (fos_cazavi.phenotype).
"""
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from fos_cazavi.phenotype import predict_fos_phenotype, predict_cazavi_phenotype, predict_phenotypes


def make_blast(gene, mutations='-', identity=100.0, coverage=100.0):
    return {'gene': gene, 'identity': identity, 'coverage': coverage, 'mutations': mutations}


def make_unified(gene, mutation):
    return {'gene': gene, 'mutation': mutation}


class TestFosPhenotype:

    def test_no_genes_susceptible(self):
        result = predict_fos_phenotype([], [])
        assert result['phenotype'] == 'Susceptible'

    def test_acquired_fosa3_resistant(self):
        result = predict_fos_phenotype([make_blast('fosA3')], [])
        assert result['phenotype'] == 'Resistant'
        assert any('fosA3' in e for e in result['evidence'])

    def test_chromosomal_fosakp_i91v_not_resistant(self):
        # fosAKP I91V is intrinsic/near-universal, not a resistance signal
        result = predict_fos_phenotype([make_blast('fosAKP')], [make_unified('fosAKP', 'I91V')])
        assert result['phenotype'] == 'Susceptible'

    def test_transport_gene_stop_codon_resistant(self):
        result = predict_fos_phenotype([], [make_unified('uhpT', 'W198*')])
        assert result['phenotype'] == 'Resistant'

    def test_transport_gene_missense_not_resistant(self):
        result = predict_fos_phenotype([], [make_unified('uhpA', 'D54N')])
        assert result['phenotype'] == 'Susceptible'


class TestCazAviPhenotype:

    def test_no_kpc_susceptible(self):
        result = predict_cazavi_phenotype([], [])
        assert result['phenotype'] == 'Susceptible'

    def test_wildtype_kpc_susceptible(self):
        result = predict_cazavi_phenotype([make_blast('blaKPC-3', mutations='-')], [])
        assert result['phenotype'] == 'Susceptible'

    def test_kpc_xloop_substitution_resistant(self):
        result = predict_cazavi_phenotype(
            [make_blast('blaKPC-3', mutations='L166W (KPC-66 X-loop)')], []
        )
        assert result['phenotype'] == 'Resistant'

    def test_kpc_omega_loop_deletion_resistant(self):
        result = predict_cazavi_phenotype([], [make_unified('blaKPC', 'L166del')])
        assert result['phenotype'] == 'Resistant'

    def test_kpc_unrelated_mutation_not_resistant(self):
        result = predict_cazavi_phenotype([make_blast('blaKPC-3', mutations='A50V')], [])
        assert result['phenotype'] == 'Susceptible'


class TestPredictPhenotypes:

    def test_returns_both_drug_classes_and_disclaimer(self):
        result = predict_phenotypes([], [])
        assert 'fosfomycin' in result
        assert 'ceftazidime_avibactam' in result
        assert 'disclaimer' in result
