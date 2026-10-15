"""Input validation, staging, and the BLAST index.

Everything the external tools read goes through here first, for three reasons:

* A bad input must fail loudly. An empty file, a FASTQ, a protein FASTA or an
  all-N assembly would otherwise run "successfully" and report both drugs as
  Susceptible - a clean-looking result for a sample that was never analysed.
* GAMMA builds its ``blat`` command as a shell string, so a path containing a
  space breaks it and a path containing ``;`` would be executed. Staging every
  input under a fixed, safe name in a private directory removes both problems.
* The BLAST index must be buildable when the package is installed somewhere
  read-only (a system-wide or admin-owned environment), and must not be
  corrupted when several runs start at once.
"""

import contextlib
import gzip
import hashlib
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO

# Fraction of a sequence's residues that must be nucleotide codes (ACGTN, any
# case) for the file to be treated as DNA rather than protein.
_MIN_NUCLEOTIDE_FRACTION = 0.9
_BLAST_INDEX_EXTENSIONS = ('nhr', 'nin', 'nsq')


class InputError(Exception):
    """An input file that cannot be analysed. The message is user-facing."""


def _open_text(path):
    with open(path, 'rb') as handle:
        magic = handle.read(2)
    if magic == b'\x1f\x8b':
        return gzip.open(path, 'rt')
    return open(path)


def stage_assembly(assembly, directory):
    """Validate an assembly and write a clean, uncompressed copy into
    ``directory``.

    Returns ``(staged_path, n_contigs, total_length)``. Raises
    :class:`InputError` with a message naming the problem when the file is
    missing, empty, not FASTA, not nucleotide, has duplicate contig names, or
    contains no called bases.
    """
    path = Path(assembly)
    if not path.exists():
        raise InputError(f"assembly not found: {assembly}")
    if not path.is_file():
        raise InputError(f"assembly is not a file: {assembly}")
    if path.stat().st_size == 0:
        raise InputError(f"assembly is empty: {assembly}")

    try:
        with _open_text(path) as handle:
            first = ''
            for line in handle:
                if line.strip():
                    first = line.lstrip()
                    break
            if first.startswith('@'):
                raise InputError(f"{assembly} looks like FASTQ (reads); an assembly "
                                 f"in FASTA format is required")
            if not first.startswith('>'):
                raise InputError(f"{assembly} is not a FASTA file (first line does "
                                 f"not start with '>')")
            handle.seek(0)
            records = list(SeqIO.parse(handle, 'fasta'))
    except (OSError, EOFError, UnicodeDecodeError) as error:
        raise InputError(f"could not read {assembly}: {error}") from error
    except ValueError as error:
        raise InputError(f"{assembly} is not valid FASTA: {error}") from error

    if not records:
        raise InputError(f"no sequences in {assembly}")

    seen = set()
    duplicates = []
    for record in records:
        if record.id in seen:
            duplicates.append(record.id)
        seen.add(record.id)
    if duplicates:
        shown = ', '.join(sorted(set(duplicates))[:5])
        raise InputError(f"{assembly} has duplicate contig names ({shown}); "
                         f"contig names must be unique")

    total = 0
    nucleotide = 0
    called = 0
    for record in records:
        sequence = str(record.seq).upper()
        total += len(sequence)
        bases = sum(sequence.count(base) for base in 'ACGT')
        called += bases
        nucleotide += bases + sequence.count('N') + sequence.count('U')

    if total == 0:
        raise InputError(f"all sequences in {assembly} are empty")
    if nucleotide < _MIN_NUCLEOTIDE_FRACTION * total:
        raise InputError(f"{assembly} does not look like a nucleotide assembly "
                         f"(only {100 * nucleotide / total:.0f}% of residues are "
                         f"A/C/G/T/N); is it a protein FASTA?")
    if called == 0:
        raise InputError(f"{assembly} contains no called bases (only N)")

    staged = Path(directory) / 'assembly.fasta'
    with open(staged, 'w') as handle:
        for record in records:
            # Upper-case so every downstream tool sees the same sequence; the
            # header is kept verbatim so contig names in the output match the
            # user's file.
            handle.write(f">{record.description}\n")
            sequence = str(record.seq).upper()
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + '\n')
    return str(staged), len(records), total


def stage_file(source, directory, name):
    """Make ``source`` available under a safe name in ``directory``."""
    target = Path(directory) / name
    try:
        os.symlink(Path(source).resolve(), target)
    except OSError:
        shutil.copyfile(source, target)
    return str(target)


# --------------------------------------------------------------- BLAST index

def _index_is_fresh(fasta):
    fasta = Path(fasta)
    index = [Path(f"{fasta}.{extension}") for extension in _BLAST_INDEX_EXTENSIONS]
    if not all(path.exists() for path in index):
        return False
    source_time = fasta.stat().st_mtime
    return all(path.stat().st_mtime >= source_time for path in index)


def cache_dir():
    """Where indexes for read-only or shared reference files are kept."""
    configured = os.environ.get('FOS_CAZAVI_CACHE')
    if configured:
        return Path(configured)
    base = os.environ.get('XDG_CACHE_HOME') or os.path.join(os.path.expanduser('~'), '.cache')
    return Path(base) / 'fos-cazavi'


def _file_digest(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()[:20]


def ensure_blast_db(fasta):
    """Return a FASTA path that has an up-to-date BLAST index beside it.

    An index next to ``fasta`` is used when it exists and is newer than the
    FASTA. Otherwise one is built in a cache directory keyed by the FASTA's
    content, so a read-only install works, an edited database is never matched
    against a stale index, and concurrent runs cannot corrupt each other: each
    builds in a private directory and the first to finish is moved into place
    atomically.
    """
    if _index_is_fresh(fasta):
        return str(fasta)

    candidates = [cache_dir(), Path(tempfile.gettempdir()) / f"fos-cazavi-{os.getuid()}"]
    digest = _file_digest(fasta)
    last_error = None
    for root in candidates:
        final = root / 'blastdb' / digest
        target = final / 'reference.fasta'
        if _index_is_fresh(target):
            return str(target)
        try:
            final.parent.mkdir(parents=True, exist_ok=True)
            building = Path(tempfile.mkdtemp(prefix=f"{digest}.", dir=final.parent))
        except OSError as error:
            last_error = error
            continue
        try:
            shutil.copyfile(fasta, building / 'reference.fasta')
            completed = subprocess.run(
                ['makeblastdb', '-in', 'reference.fasta', '-dbtype', 'nucl'],
                cwd=building, capture_output=True, text=True)
            if completed.returncode != 0:
                raise InputError(f"makeblastdb failed on {fasta}: "
                                 f"{completed.stderr.strip() or completed.stdout.strip()}")
            try:
                os.rename(building, final)
            except OSError:
                # Another run finished first; its index is equivalent.
                if not _index_is_fresh(target):
                    raise
        finally:
            shutil.rmtree(building, ignore_errors=True)
        return str(target)

    raise InputError(f"could not create a BLAST index for {fasta}: no writable "
                     f"cache directory ({last_error}). Set FOS_CAZAVI_CACHE to a "
                     f"writable directory.")


class ToolError(Exception):
    """An external tool is missing or failed. The message is user-facing."""


@contextlib.contextmanager
def staged_assembly(assembly):
    """Validate ``assembly`` and yield the path of a clean staged copy, which
    is removed again afterwards."""
    with tempfile.TemporaryDirectory(prefix='fos-cazavi-') as work:
        staged, contigs, length = stage_assembly(assembly, work)
        print(f"Assembly {assembly}: {contigs} contig(s), {length:,} bp")
        yield staged
