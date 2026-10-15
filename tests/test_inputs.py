"""Input validation, staging, the BLAST index cache, and how the CLI and batch
mode fail.

A bad input must never produce a clean-looking result: before these checks an
empty file ran "successfully" and reported both drugs Susceptible.
"""

import gzip
import os
import time
from pathlib import Path

import pytest

from fos_cazavi import batch, cli
from fos_cazavi.inputs import (
    InputError, ensure_blast_db, stage_assembly, stage_file,
)

from .conftest import embed, have, requires_blast, write_genome

requires_gamma = pytest.mark.skipif(
    not (have('GAMMA.py') and have('blat')), reason='GAMMA/blat is not installed')


def stage(tmp_path, name, content, mode='w'):
    path = tmp_path / name
    if mode == 'gz':
        with gzip.open(path, 'wt') as handle:
            handle.write(content)
    else:
        with open(path, 'w', newline='') as handle:
            handle.write(content)
    work = tmp_path / 'work'
    work.mkdir(exist_ok=True)
    return stage_assembly(str(path), str(work))


class TestRejectedInputs:
    @pytest.mark.parametrize('name, content, message', [
        ('empty.fasta', '', 'empty'),
        ('blank.fasta', '\n\n', 'not a FASTA'),
        ('text.fasta', 'hello\nworld\n', 'not a FASTA'),
        ('reads.fasta', '@r1\nACGT\n+\nIIII\n', 'FASTQ'),
        ('protein.fasta', '>p\nMSLRRSNTLLLLAAALGAQAAWAHPLLQRSSTMPAEITQ\n', 'protein'),
        ('allN.fasta', '>c\n' + 'N' * 500 + '\n', 'no called bases'),
        ('dup.fasta', '>c1\nACGTACGT\n>c1\nACGTACGT\n', 'duplicate contig names'),
        ('nothing.fasta', '>c1\n\n', 'empty'),
    ])
    def test_bad_input_raises_a_clear_error(self, tmp_path, name, content, message):
        with pytest.raises(InputError, match=message):
            stage(tmp_path, name, content)

    def test_missing_file(self, tmp_path):
        with pytest.raises(InputError, match='not found'):
            stage_assembly(str(tmp_path / 'nope.fasta'), str(tmp_path))

    def test_directory_is_not_an_assembly(self, tmp_path):
        with pytest.raises(InputError, match='not a file'):
            stage_assembly(str(tmp_path), str(tmp_path))


class TestAcceptedInputs:
    def test_gzip_crlf_and_lowercase_are_normalised(self, tmp_path):
        staged, contigs, length = stage(
            tmp_path, 'g.fasta.gz', '>c1 some description\r\nacgtacgtNN\r\n>c2\r\nGGCC\r\n',
            mode='gz')
        assert (contigs, length) == (2, 14)
        text = Path(staged).read_text()
        assert text == '>c1 some description\nACGTACGTNN\n>c2\nGGCC\n'

    def test_ambiguity_codes_are_allowed(self, tmp_path):
        _, contigs, _ = stage(tmp_path, 'iupac.fasta', '>c\n' + 'ACGT' * 50 + 'RYKM\n')
        assert contigs == 1

    def test_stage_file_gives_a_safe_name(self, tmp_path):
        source = tmp_path / 'a name; with $(odd) chars.fasta'
        source.write_text('>c\nACGT\n')
        staged = stage_file(str(source), str(tmp_path), 'assembly.fasta')
        assert Path(staged).name == 'assembly.fasta'
        assert Path(staged).read_text() == '>c\nACGT\n'


@requires_blast
class TestBlastIndexCache:
    def test_fresh_index_beside_the_fasta_is_used(self, database):
        assert ensure_blast_db(database) == database

    def test_missing_index_is_built_in_the_cache(self, tmp_path, monkeypatch):
        monkeypatch.setenv('FOS_CAZAVI_CACHE', str(tmp_path / 'cache'))
        fasta = tmp_path / 'db.fasta'
        fasta.write_text('>g1\n' + 'ACGTTGCA' * 30 + '\n')
        built = ensure_blast_db(str(fasta))
        assert built.startswith(str(tmp_path / 'cache'))
        assert Path(f"{built}.nsq").exists()
        # Nothing is written next to the original (it may be read-only).
        assert not list(tmp_path.glob('db.fasta.n*'))
        # A second call reuses the same index.
        assert ensure_blast_db(str(fasta)) == built

    def test_edited_database_never_uses_a_stale_index(self, tmp_path, monkeypatch):
        monkeypatch.setenv('FOS_CAZAVI_CACHE', str(tmp_path / 'cache'))
        fasta = tmp_path / 'db.fasta'
        fasta.write_text('>g1\n' + 'ACGTTGCA' * 30 + '\n')
        first = ensure_blast_db(str(fasta))
        time.sleep(0.01)
        fasta.write_text('>g2\n' + 'TTGGCCAA' * 30 + '\n')
        second = ensure_blast_db(str(fasta))
        assert first != second

    def test_index_older_than_the_fasta_is_ignored(self, tmp_path, monkeypatch):
        monkeypatch.setenv('FOS_CAZAVI_CACHE', str(tmp_path / 'cache'))
        fasta = tmp_path / 'db.fasta'
        fasta.write_text('>g1\n' + 'ACGTTGCA' * 30 + '\n')
        for extension in ('nhr', 'nin', 'nsq'):
            Path(f"{fasta}.{extension}").write_text('stale')
        old = time.time() - 3600
        for extension in ('nhr', 'nin', 'nsq'):
            os.utime(f"{fasta}.{extension}", (old, old))
        assert ensure_blast_db(str(fasta)) != str(fasta)


class TestCliErrors:
    @pytest.mark.parametrize('command', ['fos-cazavi-all', 'fos-cazavi-acquired',
                                         'fos-cazavi-mutations'])
    def test_empty_assembly_fails_with_exit_code_2(self, tmp_path, capsys, command):
        empty = tmp_path / 'empty.fasta'
        empty.write_text('')
        code = cli.main([command, '-a', str(empty), '-o', str(tmp_path / 'out')])
        assert code == 2
        assert 'ERROR: assembly is empty' in capsys.readouterr().err
        assert not (tmp_path / 'out_summary.txt').exists()

    def test_missing_assembly_fails_cleanly(self, tmp_path, capsys):
        code = cli.main(['fos-cazavi-all', '-a', str(tmp_path / 'nope.fasta'),
                         '-o', str(tmp_path / 'out')])
        assert code == 2
        assert 'Traceback' not in capsys.readouterr().err

    def test_version(self, capsys):
        with pytest.raises(SystemExit) as exit_info:
            cli.main(['--version'])
        assert exit_info.value.code == 0
        assert 'AMRFinderPlus' in capsys.readouterr().out


@requires_blast
class TestAwkwardPaths:
    def test_output_directory_is_created(self, tmp_path, reference_cds):
        assembly = write_genome(tmp_path / 'g.fasta', {'c1': embed(reference_cds['fosA3'])})
        prefix = tmp_path / 'not' / 'yet' / 'there' / 'sample'
        assert cli.main(['fos-cazavi-acquired', '-a', assembly, '-o', str(prefix),
                         '--organism', 'Escherichia']) == 0
        assert Path(f"{prefix}_summary.txt").exists()

    def test_gzipped_assembly_gives_the_same_call(self, tmp_path, reference_cds):
        plain = write_genome(tmp_path / 'g.fasta', {'c1': embed(reference_cds['blaNDM-1'])})
        gz = tmp_path / 'g.fasta.gz'
        with open(plain, 'rb') as source, gzip.open(gz, 'wb') as target:
            target.write(source.read())
        assert cli.main(['fos-cazavi-acquired', '-a', str(gz), '-o', str(tmp_path / 'z')]) == 0
        assert 'Ceftazidime-Avibactam (CAZ/AVI): Resistant' in \
            (tmp_path / 'z_summary.txt').read_text()

    @requires_gamma
    def test_shell_metacharacters_in_paths_are_not_executed(
            self, tmp_path, monkeypatch, reference_cds):
        # GAMMA builds a shell command from its arguments; a file name must
        # never be run as one, and must not stop GAMMA from running either.
        directory = tmp_path / 'with space'
        directory.mkdir()
        assembly = write_genome(directory / 'x;touch INJECTED;.fasta',
                                {'c1': embed(reference_cds['blaKPC-2'])})
        monkeypatch.chdir(directory)
        assert cli.main(['fos-cazavi-all', '-a', assembly,
                         '-o', str(directory / 'out')]) == 0
        assert not (directory / 'INJECTED').exists()
        assert not Path('INJECTED').exists()
        assert (directory / 'out_gamma.gamma').exists()


@requires_blast
class TestBatchIsolation:
    def test_one_bad_assembly_does_not_stop_the_batch(self, tmp_path, reference_cds):
        genomes = tmp_path / 'genomes'
        genomes.mkdir()
        write_genome(genomes / 'good.fasta', {'c1': embed(reference_cds['blaNDM-1'])})
        (genomes / 'bad.fasta').write_text('')
        out = tmp_path / 'out'
        code = cli.main(['batch', '-i', str(genomes), '-o', str(out), '-j', '1'])
        assert code == 1
        summary = (out / 'batch_combined_summary.tsv').read_text()
        assert 'good' in summary and 'bad' not in summary

    def test_gzipped_assemblies_are_found(self, tmp_path):
        (tmp_path / 'a.fasta.gz').write_bytes(b'')
        (tmp_path / 'b.fna').write_text('')
        (tmp_path / 'c.txt.gz').write_bytes(b'')
        found = batch.find_assemblies([str(tmp_path)])
        assert sorted(p.name for p in found) == ['a.fasta.gz', 'b.fna']
        assert batch._sample_name('x/a.fasta.gz') == 'a'

    def test_clashing_sample_names_are_refused(self, tmp_path, capsys, reference_cds):
        write_genome(tmp_path / 's1.fasta', {'c1': embed(reference_cds['blaNDM-1'])})
        write_genome(tmp_path / 's1.fna', {'c1': embed(reference_cds['blaNDM-1'])})
        code = cli.main(['batch', '-i', str(tmp_path), '-o', str(tmp_path / 'out')])
        assert code == 2
        assert "sample 's1'" in capsys.readouterr().err
