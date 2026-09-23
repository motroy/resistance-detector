"""Batch processing, threading, and the combined output tables."""

import csv
import subprocess
import sys
from pathlib import Path

import pytest

from fos_cazavi.batch import (
    ASSEMBLY_SUFFIXES, COMBINED_COLUMNS, GENE_COLUMNS, default_jobs,
    find_assemblies, find_summaries, summarise, write_combined,
)

from .conftest import embed, requires_blast, write_genome

ROOT = Path(__file__).resolve().parent.parent


class TestAssemblyDiscovery:
    def test_finds_assemblies_in_a_directory(self, tmp_path, reference_cds):
        for name in ('a.fasta', 'b.fna', 'c.fa'):
            write_genome(tmp_path / name, {'contig1': reference_cds['blaKPC-2']})
        (tmp_path / 'notes.txt').write_text('ignore me')
        found = find_assemblies([str(tmp_path)])
        assert sorted(p.name for p in found) == ['a.fasta', 'b.fna', 'c.fa']

    def test_accepts_explicit_files(self, tmp_path, reference_cds):
        path = write_genome(tmp_path / 'one.fasta', {'c': reference_cds['blaKPC-2']})
        assert [p.name for p in find_assemblies([path])] == ['one.fasta']

    def test_missing_input_is_skipped_not_fatal(self, tmp_path):
        assert find_assemblies([str(tmp_path / 'nope.fasta')]) == []

    def test_default_jobs_is_at_least_one(self):
        assert default_jobs() >= 1


@requires_blast
class TestBatchRun:
    @pytest.fixture(scope='class')
    @classmethod
    def batch_output(cls, tmp_path_factory, database, reference_cds):
        genomes = tmp_path_factory.mktemp('genomes')
        write_genome(genomes / 'resistant.fasta',
                     {'c1': embed(reference_cds['blaNDM-1']),
                      'c2': embed(reference_cds['fosA3'])})
        write_genome(genomes / 'susceptible.fasta',
                     {'c1': embed(reference_cds['blaKPC-2'])})
        out = tmp_path_factory.mktemp('batch')
        subprocess.run(
            [sys.executable, '-m', 'fos_cazavi.cli', 'batch',
             '-i', str(genomes), '-o', str(out), '-d', database,
             '--organism', 'Escherichia', '-j', '2'],
            check=True, capture_output=True, cwd=ROOT)
        return out

    def test_combined_summary_has_one_row_per_sample(self, batch_output):
        rows = list(csv.DictReader(
            open(batch_output / 'batch_combined_summary.tsv'), delimiter='\t'))
        assert len(rows) == 2
        assert [r['Sample'] for r in rows] == ['resistant.fasta', 'susceptible.fasta']

    def test_combined_summary_columns_are_stable(self, batch_output):
        header = open(batch_output / 'batch_combined_summary.tsv').readline()
        assert header.rstrip('\n').split('\t') == COMBINED_COLUMNS

    def test_phenotypes_are_carried_into_the_combined_table(self, batch_output):
        rows = {r['Sample']: r for r in csv.DictReader(
            open(batch_output / 'batch_combined_summary.tsv'), delimiter='\t')}
        resistant = rows['resistant.fasta']
        assert resistant['Predicted_Phenotype_Fosfomycin'] == 'Resistant'
        assert resistant['Predicted_Phenotype_Ceftazidime_Avibactam'] == 'Resistant'
        assert 'blaNDM-1' in resistant['Carbapenemases']
        assert 'fosA3' in resistant['Fosfomycin_Enzymes']
        assert rows['susceptible.fasta']['Predicted_Phenotype_Ceftazidime_Avibactam'] \
            == 'Susceptible'

    def test_per_gene_table_is_written(self, batch_output):
        rows = list(csv.DictReader(
            open(batch_output / 'batch_combined_genes.tsv'), delimiter='\t'))
        assert rows and set(rows[0]) == set(GENE_COLUMNS)
        assert {r['Sample'] for r in rows} == {'resistant.fasta', 'susceptible.fasta'}

    def test_each_sample_gets_its_own_log(self, batch_output):
        # The logger is a process-wide singleton; without resetting its handlers
        # every sample after the first would write into the first one's file.
        logs = sorted(batch_output.glob('*_analysis.log'))
        assert len(logs) == 2
        for log in logs:
            text = log.read_text()
            assert log.name.split('_analysis')[0] in text

    def test_per_sample_outputs_are_still_written(self, batch_output):
        for sample in ('resistant', 'susceptible'):
            assert (batch_output / f'{sample}_summary.tsv').exists()
            assert (batch_output / f'{sample}_summary.json').exists()
            assert (batch_output / f'{sample}_results.tsv').exists()


@requires_blast
class TestParallelMatchesSerial:
    def test_same_results_either_way(self, tmp_path, database, reference_cds):
        genomes = tmp_path / 'genomes'
        genomes.mkdir()
        for index in range(3):
            write_genome(genomes / f'sample{index}.fasta',
                         {'c1': embed(reference_cds['blaNDM-1'])})

        outputs = {}
        for jobs in (1, 3):
            out = tmp_path / f'out{jobs}'
            subprocess.run(
                [sys.executable, '-m', 'fos_cazavi.cli', 'batch',
                 '-i', str(genomes), '-o', str(out), '-d', database, '-j', str(jobs)],
                check=True, capture_output=True, cwd=ROOT)
            outputs[jobs] = (out / 'batch_combined_summary.tsv').read_text()

        assert outputs[1] == outputs[3]


@requires_blast
class TestCombineExistingResults:
    def test_combine_reads_result_directories(self, tmp_path, database, reference_cds):
        genomes = tmp_path / 'g'
        genomes.mkdir()
        write_genome(genomes / 'x.fasta', {'c1': embed(reference_cds['fosA3'])})
        out = tmp_path / 'run'
        subprocess.run(
            [sys.executable, '-m', 'fos_cazavi.cli', 'batch',
             '-i', str(genomes), '-o', str(out), '-d', database, '-j', '1'],
            check=True, capture_output=True, cwd=ROOT)

        combined = tmp_path / 'again'
        subprocess.run(
            [sys.executable, '-m', 'fos_cazavi.cli', 'combine',
             '-i', str(out), '-o', str(combined)],
            check=True, capture_output=True, cwd=ROOT)

        rows = list(csv.DictReader(
            open(f'{combined}_combined_summary.tsv'), delimiter='\t'))
        assert len(rows) == 1
        assert rows[0]['Predicted_Phenotype_Fosfomycin'] == 'Resistant'

    def test_find_summaries_ignores_other_json(self, tmp_path):
        (tmp_path / 'a_summary.json').write_text('{}')
        (tmp_path / 'other.json').write_text('{}')
        found = find_summaries([str(tmp_path)])
        assert [p.name for p in found] == ['a_summary.json']


@requires_blast
class TestThreads:
    def test_threads_option_is_accepted_and_changes_nothing(
            self, tmp_path, database, reference_cds):
        # Threading is a performance knob: the calls must be identical.
        genome = write_genome(tmp_path / 'g.fasta', {'c1': embed(reference_cds['blaNDM-1'])})
        results = {}
        for threads in (1, 4):
            prefix = tmp_path / f'out{threads}'
            subprocess.run(
                [sys.executable, '-m', 'fos_cazavi.cli', 'fos-cazavi-all',
                 '-a', genome, '-o', str(prefix), '-d', database,
                 '--threads', str(threads)],
                check=True, capture_output=True, cwd=ROOT)
            results[threads] = Path(f'{prefix}_results.tsv').read_text()
        assert results[1] == results[4]
