"""Reference database creation.

The database is built by :mod:`fos_cazavi.build_data`, which derives
every sequence, every mutation position and the whole blaKPC allele table from
one AMRFinderPlus release and validates them against each other.  This module
is a thin wrapper so ``fos-cazavi create-db`` and the packaged data are always
produced by exactly the same code - there is no second, divergent definition of
what the reference data should contain.
"""

import subprocess
import sys
from pathlib import Path

def create_db(email, output_dir, output_prefix='resistance_db'):
    """Rebuild the reference data into ``output_dir``.

    ``email`` is accepted for backwards compatibility with the previous NCBI
    Entrez-based builder; the current builder downloads the AMRFinderPlus
    release files directly and does not need it.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    command = [
        sys.executable, '-m', 'fos_cazavi.build_data',
        '--amr-dir', str(output_dir / 'amrfinder_data'),
        '--out-dir', str(output_dir),
    ]
    print(f"Building reference data into {output_dir} ...")
    completed = subprocess.run(command)
    if completed.returncode != 0:
        print('ERROR: reference data build failed', file=sys.stderr)
        return completed.returncode

    database = output_dir / 'example_database.fasta'
    if output_prefix != 'example_database' and database.exists():
        target = output_dir / f"{output_prefix}.fasta"
        target.write_bytes(database.read_bytes())
        print(f"Reference database available as {target}")
    print('Done.')
    return 0
