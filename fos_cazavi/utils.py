import csv
import datetime
import logging
import subprocess
import sys
from pathlib import Path

def setup_logger(output_prefix, args, console=True):
    """Set up logging to the sample's log file, and optionally to the console.

    Batch workers pass ``console=False``: their log file already records
    everything, and fifty samples logging to one terminal at once is unreadable.
    """
    log_file = f"{output_prefix}_analysis.log"

    # Create logger.  Handlers from a previous sample are closed and dropped
    # first: this logger is a process-wide singleton, so leaving them attached
    # would send every later sample's log into the first sample's file.
    logger = logging.getLogger('fos_cazavi')
    logger.setLevel(logging.INFO)
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()

    # File handler
    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(logging.INFO)

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    if console:
        logger.addHandler(ch)

    # Log run details
    logger.info("=" * 60)
    logger.info("FOS-CAZAVI Resistance Detector Analysis Log")
    logger.info("=" * 60)
    logger.info(f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Command: {' '.join(sys.argv)}")
    logger.info("-" * 60)
    logger.info("Parameters:")
    for arg, value in vars(args).items():
        if arg != 'func':
            logger.info(f"  {arg}: {value}")
    logger.info("-" * 60)

    return logger

def get_tool_version(tool):
    """Get version of external tool"""
    try:
        # Try version (seqkit)
        result = subprocess.run([tool, 'version'], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.split('\n')[0]

        # Try -version first (BLAST style)
        result = subprocess.run([tool, '-version'], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.split('\n')[0]

        # Try --version (others)
        result = subprocess.run([tool, '--version'], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.split('\n')[0]

        return "Unknown"
    except:
        return "Not found"

def log_tool_versions(logger):
    """Log versions of dependencies"""
    logger.info("External Tools:")
    for tool in ['blastn', 'GAMMA.py', 'seqkit']:
        version = get_tool_version(tool)
        logger.info(f"  {tool}: {version}")
    logger.info("-" * 60)

def check_dependencies(tools):
    """Check if required tools are installed"""
    missing = []

    for tool in tools:
        try:
            subprocess.run([tool, '-version'],
                         capture_output=True,
                         check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            # Try --version instead
            try:
                subprocess.run([tool, '--version'],
                             capture_output=True,
                             check=True)
            except (subprocess.CalledProcessError, FileNotFoundError):
                missing.append(tool)

    if missing:
        print(f"ERROR: Missing required tools: {', '.join(missing)}",
              file=sys.stderr)
        return False
    return True

def load_primers(primers_file):
    """Load primers from TSV file"""
    primers = {}
    if not primers_file or not Path(primers_file).exists():
        return primers

    print(f"Loading primers from {primers_file}...")

    try:
        with open(primers_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                name = row.get('Primer', row.get('name'))
                seq = row.get('Nucleotide_sequence', row.get('seq', row.get('sequence')))
                purpose = row.get('Purpose', row.get('purpose', ''))
                mutation = row.get('Mutation', row.get('mutation', ''))
                gene = row.get('Gene', row.get('gene', ''))
                pair_id = row.get('Pair_ID', row.get('pair_id', ''))

                if name and seq:
                    primers[name] = {
                        'seq': seq.strip(),
                        'purpose': purpose.strip(),
                        'mutation': mutation.strip() if mutation and mutation.strip() != '-' else None,
                        'gene': gene.strip() if gene and gene.strip() != '-' else None,
                        'pair_id': pair_id.strip() if pair_id and pair_id.strip() != '-' else None
                    }

        print(f"Loaded {len(primers)} primers")
        return primers

    except Exception as e:
        print(f"ERROR loading primers file: {e}", file=sys.stderr)
        return {}
