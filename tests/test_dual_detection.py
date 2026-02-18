"""
Tests for dual-method (GAMMA + seqkit) resistance mutation detection.

Covers:
  - MutationDetector._normalize_gene_name()
  - MutationDetector._extract_pair_id_from_header()
  - MutationDetector.merge_detection_results() – confidence scoring and method tracking
  - Edge cases: empty inputs, gene name variants, overlapping / non-overlapping results
"""
import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from fos_cazavi.mutations import MutationDetector


# ── Fixture helpers ────────────────────────────────────────────────────────────

def make_detector(gamma_results=None):
    """Create a MutationDetector without running any subprocesses."""
    d = MutationDetector.__new__(MutationDetector)
    d.assembly = "fake.fasta"
    d.output_prefix = "fake_output"
    d.genes_file = None
    d.primers_file = None
    d.gamma_results = gamma_results or []
    d.amplicon_results = []
    d.seqkit_mut_results = []
    d.unified_results = []
    d.primers = {}
    return d


def gamma_hit(protein, mutations, identity=99.0, coverage=100.0):
    """Build a synthetic GAMMA result dict."""
    return {
        'protein': protein,
        'contig': 'contig1',
        'contig_start': 1000,
        'contig_end': 2000,
        'identity': identity,
        'coverage': coverage,
        'mutations': mutations,
        'match_type': 'mutant' if mutations else 'native',
    }


def seqkit_hit(gene, mutations, pair_id='test_pair'):
    """Build a synthetic seqkit mutation result dict."""
    return {
        'gene': gene,
        'pair_id': pair_id,
        'mutations': mutations,
        'method': 'seqkit',
    }


# ══════════════════════════════════════════════════════════════════════════════
# _normalize_gene_name
# ══════════════════════════════════════════════════════════════════════════════

class TestNormalizeGeneName:

    @pytest.mark.parametrize("input_name,expected", [
        # Plain gene names (from primers.tsv)
        ('uhpB',              'uhpb'),
        ('uhpT',              'uhpt'),
        ('galU',              'galu'),
        ('lon',               'lon'),
        ('murA',              'mura'),
        ('acrB',              'acrb'),
        ('ompK36',            'ompk36'),
        ('ompK35',            'ompk35'),
        # bla family – numeric variant suffix stripped
        ('blaKPC',            'blakpc'),
        ('blaKPC-3',          'blakpc'),
        ('blaKPC-31',         'blakpc'),
        ('blaKPC-190',        'blakpc'),
        ('blaOXA-48',         'blaoxa'),
        ('blaCMY-178',        'blacmy'),
        ('blaSHV-12',         'blashv'),
        # GAMMA FASTA ID style (accession after first underscore)
        ('acrB_WP_002892069.1',  'acrb'),
        ('OmpK35_CAA09665.1',    'ompk35'),
        ('OmpK36_ADG56549.1',    'ompk36'),
        ('acrA_WP_002892072.1',  'acra'),
    ])
    def test_normalize(self, input_name, expected):
        assert MutationDetector._normalize_gene_name(input_name) == expected


# ══════════════════════════════════════════════════════════════════════════════
# _extract_pair_id_from_header
# ══════════════════════════════════════════════════════════════════════════════

class TestExtractPairIdFromHeader:

    def _pairs(self):
        return {
            'uhpB_ver': {},
            'Mut3_uhpB_1': {},
            'galU_R282V_1': {},
        }

    def test_pair_id_in_simple_header(self):
        d = make_detector()
        header = 'contig1_uhpB_ver:100-500'
        assert d._extract_pair_id_from_header(header, self._pairs()) == 'uhpB_ver'

    def test_pair_id_in_complex_header(self):
        d = make_detector()
        header = 'NODE_1_Mut3_uhpB_1:2000-2800 len=800'
        assert d._extract_pair_id_from_header(header, self._pairs()) == 'Mut3_uhpB_1'

    def test_no_match_returns_none(self):
        d = make_detector()
        header = 'contig1_unknown_pair:100-500'
        assert d._extract_pair_id_from_header(header, self._pairs()) is None

    def test_empty_header(self):
        d = make_detector()
        assert d._extract_pair_id_from_header('', self._pairs()) is None


# ══════════════════════════════════════════════════════════════════════════════
# merge_detection_results – confidence scoring
# ══════════════════════════════════════════════════════════════════════════════

class TestMergeDetectionResults:

    # ── Basic confidence rules ─────────────────────────────────────────────

    def test_both_methods_gives_100_confidence(self):
        """Mutation detected by GAMMA AND seqkit → 100% confidence."""
        d = make_detector([gamma_hit('blaKPC-3', ['D179Y'])])
        seqkit = [seqkit_hit('blaKPC', ['D179Y'])]
        unified = d.merge_detection_results(seqkit)

        assert len(unified) == 1
        r = unified[0]
        assert r['mutation'] == 'D179Y'
        assert r['confidence'] == 100
        assert 'gamma' in r['methods']
        assert 'seqkit' in r['methods']
        assert r['gamma_detail'] is not None
        assert r['seqkit_detail'] is not None

    def test_gamma_only_gives_50_confidence(self):
        """Mutation detected only by GAMMA → 50% confidence."""
        d = make_detector([gamma_hit('blaKPC-3', ['D179Y'])])
        unified = d.merge_detection_results([])

        assert len(unified) == 1
        r = unified[0]
        assert r['confidence'] == 50
        assert r['methods'] == ['gamma']
        assert r['seqkit_detail'] is None

    def test_seqkit_only_gives_50_confidence(self):
        """Mutation detected only by seqkit → 50% confidence."""
        d = make_detector()
        seqkit = [seqkit_hit('blaKPC', ['D179Y'])]
        unified = d.merge_detection_results(seqkit)

        assert len(unified) == 1
        r = unified[0]
        assert r['confidence'] == 50
        assert r['methods'] == ['seqkit']
        assert r['gamma_detail'] is None

    # ── Empty inputs ───────────────────────────────────────────────────────

    def test_no_mutations_from_either_method(self):
        d = make_detector()
        unified = d.merge_detection_results([])
        assert unified == []

    def test_gamma_wildtype_no_mutations(self):
        """GAMMA hit with empty mutation list contributes nothing."""
        d = make_detector([gamma_hit('blaKPC-3', [])])
        unified = d.merge_detection_results([])
        assert unified == []

    # ── Gene name normalization across methods ─────────────────────────────

    def test_gamma_accession_matches_seqkit_plain_name(self):
        """blaKPC-3 (GAMMA) and blaKPC (seqkit) resolve to same gene."""
        d = make_detector([gamma_hit('blaKPC-3', ['V240G'])])
        seqkit = [seqkit_hit('blaKPC', ['V240G'])]
        unified = d.merge_detection_results(seqkit)

        assert len(unified) == 1
        assert unified[0]['confidence'] == 100

    def test_gamma_full_id_matches_seqkit_plain(self):
        """acrB_WP_002892069.1 (GAMMA) matches acrB (seqkit)."""
        d = make_detector([gamma_hit('acrB_WP_002892069.1', ['G617D'])])
        seqkit = [seqkit_hit('acrB', ['G617D'])]
        unified = d.merge_detection_results(seqkit)

        assert len(unified) == 1
        assert unified[0]['confidence'] == 100

    def test_ompK36_case_insensitive_match(self):
        """OmpK36 (GAMMA FASTA) matches ompK36 (seqkit/primers)."""
        d = make_detector([gamma_hit('OmpK36_ADG56549.1', ['G134D'])])
        seqkit = [seqkit_hit('ompK36', ['G134D'])]
        unified = d.merge_detection_results(seqkit)

        assert len(unified) == 1
        assert unified[0]['confidence'] == 100

    # ── Multiple mutations same gene ───────────────────────────────────────

    def test_multiple_mutations_both_confirmed(self):
        """D179Y + T243M both confirmed by both methods → both 100%."""
        d = make_detector([gamma_hit('blaKPC-3', ['D179Y', 'T243M'])])
        seqkit = [seqkit_hit('blaKPC', ['D179Y', 'T243M'])]
        unified = d.merge_detection_results(seqkit)

        assert len(unified) == 2
        confidences = {r['mutation']: r['confidence'] for r in unified}
        assert confidences['D179Y'] == 100
        assert confidences['T243M'] == 100

    def test_partial_overlap_mixed_confidence(self):
        """GAMMA finds D179Y+T243M; seqkit finds only D179Y → mixed confidence."""
        d = make_detector([gamma_hit('blaKPC-3', ['D179Y', 'T243M'])])
        seqkit = [seqkit_hit('blaKPC', ['D179Y'])]
        unified = d.merge_detection_results(seqkit)

        assert len(unified) == 2
        conf = {r['mutation']: r['confidence'] for r in unified}
        assert conf['D179Y'] == 100   # both methods
        assert conf['T243M'] == 50    # GAMMA only

    def test_non_overlapping_mutations_each_50(self):
        """GAMMA finds D179Y; seqkit finds V240G → both 50%, separate entries."""
        d = make_detector([gamma_hit('blaKPC-3', ['D179Y'])])
        seqkit = [seqkit_hit('blaKPC', ['V240G'])]
        unified = d.merge_detection_results(seqkit)

        assert len(unified) == 2
        conf = {r['mutation']: r['confidence'] for r in unified}
        assert conf['D179Y'] == 50
        assert conf['V240G'] == 50
        methods = {r['mutation']: r['methods'] for r in unified}
        assert methods['D179Y'] == ['gamma']
        assert methods['V240G'] == ['seqkit']

    # ── Multiple genes ─────────────────────────────────────────────────────

    def test_multiple_genes_independent_confidence(self):
        """KPC D179Y confirmed by both (100%); OmpK36 G134D GAMMA only (50%)."""
        d = make_detector([
            gamma_hit('blaKPC-3', ['D179Y']),
            gamma_hit('OmpK36_ADG56549.1', ['G134D']),
        ])
        seqkit = [seqkit_hit('blaKPC', ['D179Y'])]  # OmpK36 not in seqkit
        unified = d.merge_detection_results(seqkit)

        assert len(unified) == 2
        by_mut = {r['mutation']: r for r in unified}
        assert by_mut['D179Y']['confidence'] == 100
        assert by_mut['G134D']['confidence'] == 50

    def test_seqkit_multiple_genes(self):
        """seqkit detects mutations in two genes; GAMMA has no hits → both 50%."""
        d = make_detector()
        seqkit = [
            seqkit_hit('uhpB', ['G469R'], pair_id='Mut3_uhpB_1'),
            seqkit_hit('galU', ['R282V'], pair_id='galU_R282V_1'),
        ]
        unified = d.merge_detection_results(seqkit)

        assert len(unified) == 2
        for r in unified:
            assert r['confidence'] == 50
            assert r['methods'] == ['seqkit']

    # ── Idempotency / result stored on self ───────────────────────────────

    def test_unified_results_stored_on_self(self):
        """merge_detection_results() sets self.unified_results."""
        d = make_detector([gamma_hit('blaKPC-3', ['D179Y'])])
        seqkit = [seqkit_hit('blaKPC', ['D179Y'])]
        result = d.merge_detection_results(seqkit)

        assert result is d.unified_results
        assert len(d.unified_results) == 1

    # ── Sorting ────────────────────────────────────────────────────────────

    def test_results_sorted_by_gene_then_mutation(self):
        """Unified results are sorted lexicographically by (gene, mutation)."""
        d = make_detector([
            gamma_hit('OmpK36_ADG56549.1', ['G213D', 'G134D']),
            gamma_hit('blaKPC-3', ['V240G', 'D179Y']),
        ])
        unified = d.merge_detection_results([])
        keys = [(r['gene'], r['mutation']) for r in unified]
        assert keys == sorted(keys)

    # ── Double-counting prevention ─────────────────────────────────────────

    def test_duplicate_gamma_hits_not_double_counted(self):
        """Two GAMMA hits for the same protein/mutation → only one entry."""
        d = make_detector([
            gamma_hit('blaKPC-3', ['D179Y']),
            gamma_hit('blaKPC-3', ['D179Y']),   # duplicate
        ])
        unified = d.merge_detection_results([])
        # Only one entry per (gene, mutation) key
        muts = [r['mutation'] for r in unified]
        assert muts.count('D179Y') == 1

    def test_duplicate_seqkit_hits_not_double_counted(self):
        """Two seqkit hits for the same gene/mutation → only one entry."""
        d = make_detector()
        seqkit = [
            seqkit_hit('blaKPC', ['D179Y'], pair_id='pair_A'),
            seqkit_hit('blaKPC', ['D179Y'], pair_id='pair_B'),  # duplicate
        ]
        unified = d.merge_detection_results(seqkit)
        muts = [r['mutation'] for r in unified]
        assert muts.count('D179Y') == 1

    # ── OXA-48 example ────────────────────────────────────────────────────

    def test_oxa48_dual_mutation_both_confirmed(self):
        """blaOXA-48 P68A + Y211S, both confirmed by GAMMA and seqkit."""
        d = make_detector([gamma_hit('blaOXA-48', ['P68A', 'Y211S'])])
        seqkit = [seqkit_hit('blaOXA-48', ['P68A', 'Y211S'])]
        unified = d.merge_detection_results(seqkit)

        assert len(unified) == 2
        for r in unified:
            assert r['confidence'] == 100
            assert set(r['methods']) == {'gamma', 'seqkit'}
