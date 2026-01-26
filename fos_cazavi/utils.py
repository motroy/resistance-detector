import subprocess
import sys
from Bio.Seq import Seq
from pathlib import Path
from collections import defaultdict
import csv

# Default known resistance mutations (fallback)
KNOWN_MUTATIONS = {
    'blaKPC': {
        179: {'ref': 'D', 'variants': ['Y', 'N'], 'name': 'D179Y/N'},
        240: {'ref': 'V', 'variants': ['G'], 'name': 'V240G'},
        243: {'ref': 'T', 'variants': ['M'], 'name': 'T243M'}
    },
    'blaOXA-48': {
        68: {'ref': 'P', 'variants': ['A'], 'name': 'P68A'},
        211: {'ref': 'Y', 'variants': ['S'], 'name': 'Y211S'}
    },
    'murA': {
        369: {'ref': 'D', 'variants': ['N'], 'name': 'D369N'},
        370: {'ref': 'L', 'variants': ['I'], 'name': 'L370I'}
    },
    'uhpT': {
        55: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G55D/*'},
        198: {'ref': 'W', 'variants': ['*', 'R'], 'name': 'W198*/R'},
        258: {'ref': 'E', 'variants': ['*', 'K'], 'name': 'E258*/K'},
        350: {'ref': 'W', 'variants': ['*', 'R'], 'name': 'W350*/R'}
    },
    'glpT': {
        44: {'ref': 'E', 'variants': ['*', 'K'], 'name': 'E44*/K'},
        88: {'ref': 'W', 'variants': ['*', 'R', 'G'], 'name': 'W88*/R/G'},
        90: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G90D/*'},
        234: {'ref': 'W', 'variants': ['*', 'R'], 'name': 'W234*/R'},
        362: {'ref': 'R', 'variants': ['C', 'H', '*'], 'name': 'R362C/H/*'}
    },
    'uhpA': {
        54: {'ref': 'D', 'variants': ['N', 'A'], 'name': 'D54N/A'},
        139: {'ref': 'R', 'variants': ['C', 'H'], 'name': 'R139C/H'}
    },
    'uhpB': {
        469: {'ref': 'G', 'variants': ['R'], 'name': 'G469R'},
        350: {'ref': 'H', 'variants': ['Y', 'Q'], 'name': 'H350Y/Q'}
    },
    'uhpC': {
        384: {'ref': 'F', 'variants': ['L'], 'name': 'F384L'}
    },
    'cyaA': {
        463: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G463D/*'}
    },
    'ptsI': {
        191: {'ref': 'H', 'variants': ['Y', 'Q'], 'name': 'H191Y/Q'}
    },
    'galU': {
        282: {'ref': 'R', 'variants': ['V'], 'name': 'R282V'}
    },
    'lon': {
        558: {'ref': 'Q', 'variants': ['*'], 'name': 'Q558*'}
    },
    'fosAKP': {
        91: {'ref': 'I', 'variants': ['V'], 'name': 'I91V'}
    },
    'fosA': {
        90: {'ref': 'K', 'variants': ['E', 'Q'], 'name': 'K90E/Q'},
        119: {'ref': 'H', 'variants': ['Q', 'R'], 'name': 'H119Q/R'}
    },
    'acrB': {
        617: {'ref': 'G', 'variants': ['D', 'N'], 'name': 'G617D/N'},
        626: {'ref': 'F', 'variants': ['L'], 'name': 'F626L'},
        628: {'ref': 'A', 'variants': ['T', 'V'], 'name': 'A628T/V'}
    },
    'ompK36': {
        134: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G134D/*'},
        135: {'ref': 'D', 'variants': ['*', 'N'], 'name': 'D135*/N'},
        213: {'ref': 'G', 'variants': ['D', '*'], 'name': 'G213D/*'}
    },
    'ftsI': {
        333: {'ref': 'A', 'variants': ['V', 'T'], 'name': 'A333V/T'},
        350: {'ref': 'Y', 'variants': ['C', 'S'], 'name': 'Y350C/S'},
        357: {'ref': 'S', 'variants': ['N'], 'name': 'S357N'}
    },
    'envZ': {
        244: {'ref': 'G', 'variants': ['S', 'D'], 'name': 'G244S/D'},
        324: {'ref': 'T', 'variants': ['I', 'A'], 'name': 'T324I/A'}
    }
}

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

def detect_mutations(gene_name, sequence, mutation_db=None):
    """Detect known mutations in a gene sequence"""
    mutations_found = []

    # Use provided DB or fallback to default
    db = mutation_db if mutation_db else KNOWN_MUTATIONS

    # Get base gene name (remove variant numbers)
    # This logic matches ResistanceDetector.detect_mutations
    base_gene = gene_name.split('-')[0]
    if base_gene.startswith('bla'):
        base_gene = 'bla' + base_gene[3:].rstrip('0123456789')

    # Try exact match first
    if gene_name in db:
        mutation_dict = db[gene_name]
    elif base_gene in db:
        mutation_dict = db[base_gene]
    else:
        return mutations_found

    # Translate to protein
    try:
        # Ensure sequence length is multiple of 3
        if len(sequence) % 3 != 0:
            sequence = sequence[:-(len(sequence) % 3)]

        protein = str(Seq(sequence).translate())

        # Check each known mutation position
        for pos, mut_info in mutation_dict.items():
            if pos <= len(protein):
                observed_aa = protein[pos-1]
                ref_aa = mut_info['ref']

                if observed_aa in mut_info['variants']:
                    mutations_found.append(mut_info['name'])
                elif observed_aa != ref_aa:
                    # Novel mutation at known position
                    mutations_found.append(f"{ref_aa}{pos}{observed_aa}")

    except Exception as e:
        print(f"Warning: Could not translate {gene_name}: {e}",
              file=sys.stderr)

    return mutations_found

def load_mutation_db(mutation_file):
    """Load mutations from TSV file"""
    if not mutation_file or not Path(mutation_file).exists():
        return KNOWN_MUTATIONS

    print(f"Loading mutations from {mutation_file}...")

    new_db = defaultdict(dict)

    try:
        with open(mutation_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene = row['Gene']
                pos = int(row['Position'])
                ref = row['Ref']
                var = row['Variant']

                vars_list = var.split('/')

                if gene not in new_db:
                    new_db[gene] = {}

                if pos not in new_db[gene]:
                    new_db[gene][pos] = {
                        'ref': ref,
                        'variants': set()
                    }

                new_db[gene][pos]['variants'].update(vars_list)

        # Post-process to format names
        for gene in new_db:
            for pos in new_db[gene]:
                entry = new_db[gene][pos]
                variants_str = '/'.join(sorted(entry['variants']))
                entry['name'] = f"{entry['ref']}{pos}{variants_str}"
                entry['variants'] = list(entry['variants'])

        print(f"Loaded mutation definitions for {len(new_db)} genes")
        return new_db

    except Exception as e:
        print(f"ERROR loading mutation file: {e}", file=sys.stderr)
        return KNOWN_MUTATIONS

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
