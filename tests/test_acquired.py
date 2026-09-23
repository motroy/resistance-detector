"""Unit tests for BLAST hit selection in acquired.py, in particular the
intrinsic-vs-acquired fosA disambiguation.

Bug this guards against: a genome's single, native chromosomal fosA gene can
score a marginally *lower* bitscore against the intrinsic reference (fosAKP,
one strain's sequence) than against a horizontally-acquired allele from the
same family (fosA5, fosA10, ...), purely because that other lineage's native
copy happens to sit a percentage point or two closer to a different catalogued
allele. Plain best-bitscore selection then reports the genome's own gene as an
acquired, resistance-suggestive allele. Confirmed on four real K. pneumoniae
assemblies from the ESKAPE-fosfomycin GOLD validation set (see
bioproject_tests/ESKAPE_fos_GOLD_Kpneumoniae/), each with exactly one fosA
copy, at 94.9-99.1% identity to fosAKP versus 95.0-99.0% to the winning
acquired name - well within the noise of ordinary chromosomal-gene divergence
between K. pneumoniae lineages, not evidence of a second, acquired gene.
"""

from fos_cazavi.acquired import BlastDetector


def make_detector():
    """A BlastDetector with no real assembly/database - only used to call
    _prefer_intrinsic_naming(), which touches no other instance state."""
    return BlastDetector.__new__(BlastDetector)


def hit(gene, query_id='contig1', qstart=1, qend=420, identity=95.0):
    return {'query_id': query_id, 'gene': gene, 'identity': identity,
           'coverage': 100.0, 'qstart': qstart, 'qend': qend,
           'sstart': 1, 'send': qend - qstart + 1, 'bitscore': identity * 7}


class TestPreferIntrinsicNaming:
    def test_close_call_is_renamed_to_the_intrinsic_gene(self):
        # The genome's only fosA copy: fosA5 won the bitscore race by a hair,
        # but fosAKP hit the identical span within tolerance.
        acquired_hit = hit('fosA5', identity=96.19)
        intrinsic_hit = hit('fosAKP', identity=95.71)
        detector = make_detector()

        corrected = detector._prefer_intrinsic_naming(
            all_hits=[acquired_hit, intrinsic_hit], kept_hits=[acquired_hit])

        assert len(corrected) == 1
        assert corrected[0]['gene'] == 'fosAKP'
        assert corrected[0]['identity'] == 95.71

    def test_clear_win_outside_tolerance_is_not_renamed(self):
        # A real acquired copy: near-exact match to fosA5, much worse to
        # fosAKP - a large gap, not ordinary chromosomal divergence.
        acquired_hit = hit('fosA5', identity=99.9)
        intrinsic_hit = hit('fosAKP', identity=90.0)
        detector = make_detector()

        corrected = detector._prefer_intrinsic_naming(
            all_hits=[acquired_hit, intrinsic_hit], kept_hits=[acquired_hit])

        assert corrected[0]['gene'] == 'fosA5'

    def test_no_competing_intrinsic_hit_is_unaffected(self):
        acquired_hit = hit('fosA5', identity=96.19)
        detector = make_detector()

        corrected = detector._prefer_intrinsic_naming(
            all_hits=[acquired_hit], kept_hits=[acquired_hit])

        assert corrected[0]['gene'] == 'fosA5'

    def test_non_overlapping_intrinsic_hit_on_the_same_contig_is_ignored(self):
        # A second, unrelated fosAKP-named locus elsewhere in the assembly
        # must not be treated as competing for this one.
        acquired_hit = hit('fosA5', query_id='contig1', qstart=1, qend=420, identity=96.0)
        elsewhere = hit('fosAKP', query_id='contig1', qstart=50000, qend=50420, identity=95.0)
        detector = make_detector()

        corrected = detector._prefer_intrinsic_naming(
            all_hits=[acquired_hit, elsewhere], kept_hits=[acquired_hit])

        assert corrected[0]['gene'] == 'fosA5'

    def test_intrinsic_hit_on_a_different_contig_is_ignored(self):
        acquired_hit = hit('fosA5', query_id='contig1', identity=96.0)
        other_contig = hit('fosAKP', query_id='contig2', identity=95.0)
        detector = make_detector()

        corrected = detector._prefer_intrinsic_naming(
            all_hits=[acquired_hit, other_contig], kept_hits=[acquired_hit])

        assert corrected[0]['gene'] == 'fosA5'

    def test_hit_already_named_intrinsic_passes_through(self):
        intrinsic_hit = hit('fosAKP', identity=99.0)
        detector = make_detector()

        corrected = detector._prefer_intrinsic_naming(
            all_hits=[intrinsic_hit], kept_hits=[intrinsic_hit])

        assert corrected[0]['gene'] == 'fosAKP'

    def test_unrelated_gene_family_is_unaffected(self):
        # Only families with a curated intrinsic counterpart (fosA) are
        # touched; an unrelated acquired gene must never be renamed.
        blakpc_hit = hit('blaKPC-2', identity=99.0)
        detector = make_detector()

        corrected = detector._prefer_intrinsic_naming(
            all_hits=[blakpc_hit], kept_hits=[blakpc_hit])

        assert corrected[0]['gene'] == 'blaKPC-2'

    def test_picks_the_best_of_several_competing_intrinsic_candidates(self):
        acquired_hit = hit('fosA10', identity=97.0)
        weaker_intrinsic = hit('fosAKP', identity=95.5)
        stronger_intrinsic = hit('fosA_PA1129', identity=96.5)
        detector = make_detector()

        corrected = detector._prefer_intrinsic_naming(
            all_hits=[acquired_hit, weaker_intrinsic, stronger_intrinsic],
            kept_hits=[acquired_hit])

        assert corrected[0]['gene'] == 'fosA_PA1129'

    def test_only_kept_hits_are_returned_in_order(self):
        acquired1 = hit('fosA5', query_id='c1', identity=96.0)
        acquired2 = hit('fosA3', query_id='c2', identity=99.0)  # no competitor
        intrinsic = hit('fosAKP', query_id='c1', identity=95.0)
        detector = make_detector()

        corrected = detector._prefer_intrinsic_naming(
            all_hits=[acquired1, acquired2, intrinsic],
            kept_hits=[acquired1, acquired2])

        assert [c['gene'] for c in corrected] == ['fosAKP', 'fosA3']
