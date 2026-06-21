"""
Unit tests for the copy-number aggregation and machine-readable summary
output produced by fos_cazavi.cli.
"""
import json
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from fos_cazavi.cli import _aggregate_genes_by_copy_number, write_machine_summary, write_combined_tsv


def make_blast_result(gene, contig, start, end, identity=99.0, coverage=100.0, mutations='-'):
    return {
        'contig': contig,
        'gene': gene,
        'identity': f"{identity:.2f}",
        'coverage': f"{coverage:.2f}",
        'mutations': mutations,
        'start': start,
        'end': end,
    }


class TestAggregateGenesByCopyNumber:

    def test_single_locus_copy_number_one(self):
        results = [make_blast_result('fosA3', 'contig1', 1, 576)]
        entries = _aggregate_genes_by_copy_number(results)
        assert len(entries) == 1
        assert entries[0]['copy_number'] == 1
        assert len(entries[0]['loci']) == 1

    def test_two_loci_same_gene_copy_number_two(self):
        results = [
            make_blast_result('fosA3', 'contig1', 1, 576),
            make_blast_result('fosA3', 'contig2', 100, 676),
        ]
        entries = _aggregate_genes_by_copy_number(results)
        assert len(entries) == 1
        assert entries[0]['gene'] == 'fosA3'
        assert entries[0]['copy_number'] == 2
        assert len(entries[0]['loci']) == 2

    def test_distinct_genes_kept_separate(self):
        results = [
            make_blast_result('fosA3', 'contig1', 1, 576),
            make_blast_result('blaKPC-3', 'contig1', 2000, 3000),
        ]
        entries = _aggregate_genes_by_copy_number(results)
        assert {e['gene'] for e in entries} == {'fosA3', 'blaKPC-3'}
        assert all(e['copy_number'] == 1 for e in entries)

    def test_empty_results(self):
        assert _aggregate_genes_by_copy_number([]) == []
        assert _aggregate_genes_by_copy_number(None) == []


class TestWriteMachineSummary:

    def test_json_and_tsv_report_copy_numbers(self, tmp_path):
        prefix = str(tmp_path / "sample")
        results = [
            make_blast_result('fosA3', 'contig1', 1, 576, identity=100.0, coverage=100.0),
            make_blast_result('fosA3', 'contig2', 50, 626, identity=95.0, coverage=98.0),
            make_blast_result('blaKPC-3', 'contig1', 2000, 3000, identity=99.0, coverage=100.0),
        ]

        write_machine_summary(prefix, "sample.fasta", results, None, None, None)

        json_path = tmp_path / "sample_summary.json"
        tsv_path = tmp_path / "sample_summary.tsv"
        assert json_path.exists()
        assert tsv_path.exists()

        data = json.loads(json_path.read_text())
        assert data['total_acquired_genes_detected'] == 3
        assert data['unique_acquired_genes'] == 2
        genes_by_name = {g['gene']: g for g in data['acquired_genes']}
        assert genes_by_name['fosA3']['copy_number'] == 2
        assert genes_by_name['blaKPC-3']['copy_number'] == 1

        tsv_lines = tsv_path.read_text().strip().split('\n')
        header = tsv_lines[0].split('\t')
        assert 'Copy_Number' in header
        rows = {line.split('\t')[1]: line.split('\t') for line in tsv_lines[1:]}
        assert rows['fosA3'][header.index('Copy_Number')] == '2'
        assert rows['blaKPC-3'][header.index('Copy_Number')] == '1'

    def test_no_results_writes_empty_summary(self, tmp_path):
        prefix = str(tmp_path / "sample")
        write_machine_summary(prefix, "sample.fasta", [], None, None, None)

        data = json.loads((tmp_path / "sample_summary.json").read_text())
        assert data['total_acquired_genes_detected'] == 0
        assert data['unique_acquired_genes'] == 0
        assert data['acquired_genes'] == []


class TestWriteCombinedTsvCopyNumber:

    def test_combined_tsv_includes_copy_number_column(self, tmp_path):
        prefix = str(tmp_path / "sample")
        results = [
            make_blast_result('fosA3', 'contig1', 1, 576),
            make_blast_result('fosA3', 'contig2', 50, 626),
        ]

        write_combined_tsv(prefix, "sample.fasta", results, None, None, None)

        out_path = tmp_path / "sample_all_results.tsv"
        lines = out_path.read_text().strip().split('\n')
        header = lines[0].split('\t')
        assert 'Copy_Number' in header
        idx = header.index('Copy_Number')
        for line in lines[1:]:
            assert line.split('\t')[idx] == '2'
