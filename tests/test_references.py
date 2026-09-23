"""Unit tests for the gene-name alias mechanism in references.py.

Some chromosomal genes need more than one nucleotide reference: a single
strain's sequence for the fosfomycin transport/regulatory genes is only
84-89% identical to the K. pneumoniae ortholog (see
bioproject_tests/ESKAPE_fos_GOLD_Kpneumoniae/RESULTS_SUMMARY.md), below the
90% default detection threshold, so those genes were never actually being
located in a K. pneumoniae assembly. The fix stores a second,
K. pneumoniae-specific reference under its own BLAST-unique name
(``<gene>_Kpn``) and maps it back to the canonical gene name everywhere else
- these tests guard that mapping.
"""

from fos_cazavi.references import (
    DATA_DIR, FOS_TRANSPORT_GENES, GENE_ALIASES, canonical_gene_name,
    gene_family,
)


class TestCanonicalGeneName:
    def test_aliased_name_maps_to_canonical(self):
        assert canonical_gene_name('uhpT_Kpn') == 'uhpT'
        assert canonical_gene_name('galU_Kpn') == 'galU'

    def test_unaliased_name_is_unchanged(self):
        assert canonical_gene_name('uhpT') == 'uhpT'
        assert canonical_gene_name('blaKPC-2') == 'blaKPC-2'
        assert canonical_gene_name('fosAKP') == 'fosAKP'

    def test_every_alias_target_is_a_fosfomycin_transport_gene_or_mura(self):
        # canonical_gene_name is only meant to fold organism-specific
        # transport-gene references back to the name FOS_TRANSPORT_GENES (and
        # the phenotype logic's separate murA check) already recognises. An
        # alias pointing anywhere else would silently stop being scored.
        for canonical in GENE_ALIASES.values():
            assert canonical in FOS_TRANSPORT_GENES or canonical == 'murA'

    def test_every_alias_source_ends_in_the_organism_suffix(self):
        for source in GENE_ALIASES:
            assert source.endswith('_Kpn')


class TestGeneFamilyCanonicalises:
    def test_aliased_transport_gene_has_the_canonical_family(self):
        # This is what lets mutations.py group a BLAST hit named ``uhpT``
        # with a GAMMA hit against the same locus named ``uhpT_Kpn`` - GAMMA
        # runs against the raw shared database and never sees the alias.
        assert gene_family('uhpT_Kpn') == gene_family('uhpT') == 'uhpT'
        assert gene_family('galU_Kpn') == gene_family('galU') == 'galU'

    def test_bla_and_fosa_families_are_unaffected(self):
        assert gene_family('blaKPC-33') == 'blaKPC'
        assert gene_family('blaCTX-M-15') == 'blaCTX-M'
        assert gene_family('fosA5') == 'fosA'
        assert gene_family('fosAKP') == 'fosAKP'


class TestAliasedReferencesArePresentInTheBundledDatabase:
    def test_every_alias_source_has_a_nucleotide_reference(self):
        # A dangling alias (a name in GENE_ALIASES with no matching sequence
        # in the shipped database) would mean the K. pneumoniae reference was
        # never actually built into example_database.fasta.
        database = (DATA_DIR / 'example_database.fasta').read_text()
        for source in GENE_ALIASES:
            assert f'>{source}' in database, f"{source} missing from the bundled database"
