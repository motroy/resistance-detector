"""Shared fixtures.

The tests build small synthetic genomes and run the real pipeline against the
real bundled reference data, so they exercise the same code path a user does.
Tests needing BLAST/GAMMA/seqkit are skipped when those tools are absent.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fos_cazavi.references import DEFAULT_DATABASE  # noqa: E402


def have(tool):
    return shutil.which(tool) is not None


requires_blast = pytest.mark.skipif(
    not (have('blastn') and have('makeblastdb')),
    reason='BLAST+ is not installed')


@pytest.fixture(scope='session')
def reference_cds():
    """{gene: nucleotide CDS} from the bundled reference database."""
    return {record.id: str(record.seq)
            for record in SeqIO.parse(DEFAULT_DATABASE, 'fasta')}


@pytest.fixture(scope='session')
def database(tmp_path_factory):
    """A BLAST database built from the bundled reference sequences."""
    directory = tmp_path_factory.mktemp('db')
    path = directory / 'reference.fasta'
    shutil.copy(DEFAULT_DATABASE, path)
    if have('makeblastdb'):
        subprocess.run(['makeblastdb', '-in', str(path), '-dbtype', 'nucl'],
                       check=True, capture_output=True)
    return str(path)


def write_genome(path, contigs):
    """Write {contig_name: sequence} to a FASTA file and return the path."""
    with open(path, 'w') as handle:
        for name, sequence in contigs.items():
            handle.write(f">{name}\n")
            for index in range(0, len(sequence), 70):
                handle.write(sequence[index:index + 70] + '\n')
    return str(path)


FILLER = 'ATGCCTAGCTAGGCTTACGATCCGATCGGATCAGCTTAACGGATCCATGCATGCTAGCTAAC' * 20


def embed(sequence, flank=FILLER):
    """Put a gene inside neutral flanking sequence, as it sits in a contig."""
    return flank + sequence + flank
