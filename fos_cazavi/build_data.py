#!/usr/bin/env python3
"""Rebuild the bundled reference data from the AMRFinderPlus database.

This module is the single source of truth for everything in
``fos_cazavi/data/``.  It is deliberately reproducible and offline-checkable:
it downloads the AMRFinderPlus release files, then derives

  * ``example_database.fasta``    - nucleotide CDS used to locate genes (BLAST/GAMMA),
                                    built from AMRFinderPlus alleles plus the
                                    committed chromosomal reference CDS
  * ``reference_proteins.faa``    - per gene+organism protein references that
                                    define residue numbering for point mutations
  * ``point_mutations.tsv``       - curated point mutations, positions validated
                                    against ``reference_proteins.faa``
  * ``blaKPC_alleles.tsv``        - every known blaKPC allele expressed as its
                                    amino-acid differences from KPC-2 in
                                    standardised Ambler numbering

Every emitted mutation position is checked against the reference protein it is
numbered in; a row whose reference residue does not match is dropped and
reported, so a silently wrong position cannot reach the bundled data.

Usage:
    python3 -m fos_cazavi.build_data --amr-dir amrfinder_data
"""

import argparse
import csv
import re
import sys
import urllib.request
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .variants import collapse_runs, compare_proteins

AMRFINDER_BASE = (
    'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/'
    'AMRFinderPlus/database/latest/'
)
AMRFINDER_FILES = [
    'ReferenceGeneCatalog.txt',
    'AMR_CDS.fa',
    'AMRProt.fa',
    'AMRProt-mutation.tsv',
    'version.txt',
]

# ---------------------------------------------------------------------------
# What we track
# ---------------------------------------------------------------------------

# Acquired genes pulled straight out of AMR_CDS.fa by allele name.  These are
# horizontally acquired and species independent, so a single reference allele
# per family is enough to find them and to type them by alignment.
ACQUIRED_ALLELES = [
    # Class A carbapenemase - the ceftazidime-avibactam-relevant family
    'blaKPC-2',                     # canonical KPC reference for allele typing
    # Class D carbapenemases - inhibited by avibactam
    'blaOXA-48', 'blaOXA-181', 'blaOXA-232',
    # ESBL / AmpC context.  blaPDC is the P. aeruginosa chromosomal AmpC: every
    # isolate has it, so its presence is context, not a finding - what matters
    # clinically is expression and PDC variant, neither of which is assessed.
    'blaCTX-M-15', 'blaSHV-12', 'blaCMY-2', 'blaPDC-1',
    # Acquired fosfomycin-modifying enzymes.  fosA6 and fosA_PA1129 are
    # deliberately absent here and listed as intrinsic below.
    'fosA', 'fosA2', 'fosA3', 'fosA4', 'fosA5', 'fosA7', 'fosA8', 'fosA9',
    'fosA10', 'fosA11', 'fosA12', 'fosA13',
    'fosC2', 'fosB', 'fosL1', 'fosL2',
]

# Fosfomycin-modifying enzymes that are a normal part of a species' chromosome.
# They are detected so the user sees them, but they are present in
# fosfomycin-susceptible isolates and are never scored as acquired resistance:
# fosA_PA1129 is the P. aeruginosa chromosomal enzyme.  The K. pneumoniae one,
# which AMRFinderPlus calls fosA6, is not listed here because it is already
# carried as `fosAKP` in the chromosomal reference set and encodes an identical
# protein; having both would split hits between two names for one gene.
INTRINSIC_ALLELES = ['fosA_PA1129']

# Families where EVERY known allele is included rather than a representative.
# These are highly diverse - 62 of the 108 blaIMP alleles are under 90% identity
# to blaIMP-1 - so a couple of references would miss most of the family at the
# default identity threshold.  Their presence is what drives a resistant
# ceftazidime-avibactam call, so family-level sensitivity matters more than
# keeping the database small.
ACQUIRED_FAMILIES = [
    # Metallo-beta-lactamases - NOT inhibited by avibactam
    'blaNDM', 'blaVIM', 'blaIMP', 'blaSPM', 'blaGIM', 'blaSIM',
    # Class A carbapenemases/ESBLs that avibactam does inhibit
    'blaGES',
]


def alleles_in_families(amr_cds, families):
    """Every `<family>-<number>` allele present in AMR_CDS.fa."""
    pattern = re.compile(r'^(' + '|'.join(families) + r')-\d+$')
    return sorted(name for name in amr_cds if pattern.match(name))


# Chromosomal genes kept in the nucleotide database purely so the locus can be
# found in an assembly (BLAST/GAMMA/amplicons).  They are carried over from the
# previous build rather than re-fetched, and they never define residue
# numbering - that comes from reference_proteins.faa.
CHROMOSOMAL_GENES = [
    'murA', 'uhpT', 'glpT', 'uhpA', 'uhpB', 'uhpC', 'cyaA', 'ptsI', 'galU',
    'lon', 'acrB', 'ompK36', 'ompK35', 'ftsI', 'envZ', 'mexR', 'nalD', 'fosAKP',
    # E. coli porins: loss of either, alongside an AmpC or ESBL, is a documented
    # route to carbapenem and ceftazidime-avibactam resistance without any
    # carbapenemase being present.
    'ompC', 'ompF',
]

# Point-mutation reference proteins, taken from AMRFinderPlus.  The positions in
# AMRProt-mutation.tsv are defined against exactly these proteins, so using them
# as the numbering reference makes every position correct by construction.
# gene -> {organism: protein accession}
POINT_REFERENCES = {
    'murA':   {'Escherichia': 'WP_000357259.1'},
    'uhpT':   {'Escherichia': 'WP_000879194.1'},
    'uhpA':   {'Escherichia': 'WP_000633668.1'},
    'cyaA':   {'Escherichia': 'WP_000281668.1'},
    'ptsI':   {'Escherichia': 'WP_000623140.1'},
    'lon':    {'Escherichia': 'WP_001295325.1'},
    'acrB':   {'Escherichia': 'WP_086259262.1'},
    'ftsI':   {'Escherichia': 'WP_000642196.1',
               'Klebsiella_pneumoniae': 'WP_016532398.1',
               'Pseudomonas_aeruginosa': 'WP_003094139.1'},
    'ompK36': {'Klebsiella_pneumoniae': 'WP_002913005.1'},
    'ompK35': {'Klebsiella_pneumoniae': 'WP_004141771.1'},
    'envZ':   {'Klebsiella_pneumoniae': 'WP_317867775.1'},
    'glpT':   {'Pseudomonas_aeruginosa': 'WP_003096330.1'},
    'galU':   {'Pseudomonas_aeruginosa': 'WP_003088639.1'},
    'mexR':   {'Pseudomonas_aeruginosa': 'WP_003114897.1'},
    'nalD':   {'Pseudomonas_aeruginosa': 'WP_003092152.1'},
}


def download(amr_dir):
    amr_dir.mkdir(parents=True, exist_ok=True)
    for name in AMRFINDER_FILES:
        destination = amr_dir / name
        if destination.exists():
            continue
        print(f"Downloading {name} ...")
        urllib.request.urlretrieve(AMRFINDER_BASE + name, destination)
    version = (amr_dir / 'version.txt').read_text().strip()
    print(f"AMRFinderPlus database version: {version}")
    return version


def load_amr_cds(amr_dir):
    """{allele_name: SeqRecord} for every nucleotide CDS in AMR_CDS.fa."""
    by_allele = {}
    for record in SeqIO.parse(amr_dir / 'AMR_CDS.fa', 'fasta'):
        fields = record.description.split('|')
        if len(fields) > 5:
            by_allele.setdefault(fields[4], record)
    return by_allele


def load_amr_proteins(amr_dir):
    proteins = {}
    for record in SeqIO.parse(amr_dir / 'AMRProt.fa', 'fasta'):
        accession = record.description.split('|')[0]
        proteins.setdefault(accession, str(record.seq).rstrip('*'))
    return proteins


def load_chromosomal(chromosomal_fasta):
    """Load the curated chromosomal CDS shipped with the package.

    These were fetched once from curated RefSeq/GenBank loci and are committed
    as ``data/chromosomal_reference_cds.fasta`` so a rebuild is reproducible and
    does not depend on NCBI being reachable.  They are re-validated below.
    """
    records = []
    if not Path(chromosomal_fasta).exists():
        print(f"WARNING: {chromosomal_fasta} not found; no chromosomal CDS included",
              file=sys.stderr)
        return records
    keep = set(CHROMOSOMAL_GENES)
    for record in SeqIO.parse(chromosomal_fasta, 'fasta'):
        if record.id in keep:
            records.append(record)
    found = {record.id for record in records}
    for gene in CHROMOSOMAL_GENES:
        if gene not in found:
            print(f"WARNING: chromosomal gene {gene} missing from {chromosomal_fasta}",
                  file=sys.stderr)
    return records


def validate_cds(record):
    """Return a list of problems with a nucleotide CDS record."""
    problems = []
    sequence = str(record.seq).upper()
    if len(sequence) % 3:
        problems.append(f"length {len(sequence)} is not a multiple of 3")
    if sequence[:3] not in ('ATG', 'GTG', 'TTG'):
        problems.append(f"does not start with a start codon ({sequence[:3]})")
    protein = str(record.seq.translate())
    if not protein.endswith('*'):
        problems.append('does not end with a stop codon')
    if '*' in protein[:-1]:
        problems.append(f"internal stop codon at residue {protein.index('*') + 1}")
    return problems


# This tool reports on two drugs.  Every curated mutation is tagged with which
# of them it belongs to, so a mutation curated for tigecycline or carbapenems
# can never be presented as a fosfomycin or ceftazidime-avibactam finding.
def drug_scope(class_field, subclass_field):
    """Classify a curated mutation against this tool's scope."""
    combined = f"{class_field}/{subclass_field}".upper()
    if 'FOSFOMYCIN' in combined:
        return 'fosfomycin'
    if 'CEFTAZIDIME-AVIBACTAM' in combined:
        return 'ceftazidime-avibactam'
    if 'AVIBACTAM' in combined:
        # Ceftibuten- or aztreonam-avibactam: the inhibitor is shared but the
        # partner drug is not, so these are context, not in-scope findings.
        return 'avibactam-combination'
    return 'other'


def build_point_mutations(amr_dir, proteins):
    """Curated point mutations, each validated against its reference protein."""
    rows, dropped = [], []
    definitions = list(csv.DictReader(open(amr_dir / 'AMRProt-mutation.tsv'),
                                      delimiter='\t'))

    wanted = {}
    for gene, per_organism in POINT_REFERENCES.items():
        for organism, accession in per_organism.items():
            wanted[(organism, accession)] = gene

    for row in definitions:
        key = (row['#taxgroup'], row['accession_version'])
        if key not in wanted:
            continue
        gene = wanted[key]
        organism, accession = key
        symbol = row['standard_mutation_symbol']
        change = symbol.split('_', 1)[1] if '_' in symbol else symbol

        match = re.match(r'^([A-Z])(\d+)([A-Z]+|Ter|del)$', change)
        if not match:
            dropped.append((gene, symbol, 'unparsable symbol'))
            continue
        reference_aa, position, variant = match.group(1), int(match.group(2)), match.group(3)

        protein = proteins.get(accession)
        if protein is None:
            dropped.append((gene, symbol, f'reference protein {accession} not found'))
            continue
        if position > len(protein):
            dropped.append((gene, symbol, f'position {position} beyond protein length'))
            continue
        if protein[position - 1] != reference_aa:
            dropped.append((gene, symbol,
                            f'reference residue mismatch: table says {reference_aa}, '
                            f'{accession} has {protein[position - 1]}'))
            continue

        # Normalise the variant: Ter is a premature stop; a multi-residue string
        # at one position is an insertion/duplication, kept verbatim.
        if variant == 'Ter':
            normalised, kind = '*', 'nonsense'
        elif variant == 'del':
            normalised, kind = 'del', 'deletion'
        elif len(variant) == 1:
            normalised, kind = variant, 'substitution'
        else:
            normalised, kind = variant, 'insertion'

        rows.append({
            'Gene': gene,
            'Organism': organism,
            'Reference_Protein': accession,
            'Position': position,
            'Ref': reference_aa,
            'Variant': normalised,
            'Kind': kind,
            'Label': symbol,
            'Class': row['class'],
            'Subclass': row['subclass'],
            'Drug_Scope': drug_scope(row['class'], row['subclass']),
            'Source': 'AMRFinderPlus',
        })

    rows.sort(key=lambda row: (row['Gene'], row['Organism'], row['Position'], row['Variant']))
    return rows, dropped


def describe_changes(reference, query):
    """Amino-acid differences of ``query`` from ``reference``, in Ambler numbering.

    Uses the same comparison code as the detector itself, so an allele's stored
    definition and a called variant are guaranteed to be written identically.
    """
    changes = collapse_runs(compare_proteins(reference, query, numbering='ambler'))
    return [change['label'] for change in changes]


def build_kpc_allele_table(amr_cds, catalog):
    """Express every blaKPC allele as its differences from KPC-2."""
    kpc = {}
    for allele, record in amr_cds.items():
        if allele.startswith('blaKPC-'):
            kpc[allele] = str(record.seq.translate()).rstrip('*')

    reference = kpc.get('blaKPC-2')
    if reference is None:
        raise SystemExit('blaKPC-2 not found in AMR_CDS.fa; cannot build allele table')

    def allele_sort_key(name):
        suffix = name.split('-', 1)[1]
        return (0, int(suffix)) if suffix.isdigit() else (1, 0)

    rows = []
    for allele in sorted(kpc, key=allele_sort_key):
        changes = describe_changes(reference, kpc[allele])
        entry = catalog.get(allele, {})
        product = entry.get('product_name', '')
        rows.append({
            'Allele': allele,
            'Changes_Ambler': ';'.join(changes),
            'NCBI_subclass': entry.get('subclass', ''),
            'NCBI_product_name': product,
            'Inhibitor_resistant': 'yes' if 'inhibitor-resistant' in product else 'no',
        })
    return rows


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--amr-dir', default='amrfinder_data',
                        help='Directory holding (or receiving) the AMRFinderPlus files')
    parser.add_argument('--out-dir', default=str(Path(__file__).parent / 'data'),
                        help='Directory to write the bundled reference data into')
    args = parser.parse_args()

    amr_dir = Path(args.amr_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    version = download(amr_dir)
    amr_cds = load_amr_cds(amr_dir)
    proteins = load_amr_proteins(amr_dir)
    catalog = {row['allele']: row for row in
               csv.DictReader(open(amr_dir / 'ReferenceGeneCatalog.txt'), delimiter='\t')}

    # ---- nucleotide database ------------------------------------------------
    records, missing = [], []
    family_alleles = alleles_in_families(amr_cds, ACQUIRED_FAMILIES)
    print(f"Including {len(family_alleles)} alleles from the diverse families: "
          f"{', '.join(ACQUIRED_FAMILIES)}")
    for allele in ACQUIRED_ALLELES + INTRINSIC_ALLELES + family_alleles:
        record = amr_cds.get(allele)
        if record is None:
            missing.append(allele)
            continue
        record.id = allele
        record.description = f"{allele} acquired resistance gene (AMRFinderPlus {version})"
        records.append(record)
    if missing:
        print(f"WARNING: not found in AMR_CDS.fa: {', '.join(missing)}", file=sys.stderr)

    records.extend(load_chromosomal(out_dir / 'chromosomal_reference_cds.fasta'))

    # AMRFinderPlus ships its CDS lowercase while the chromosomal references are
    # uppercase.  Nothing downstream depends on case, but normalising it keeps
    # the file consistent and stops case from masking sequence comparisons.
    for record in records:
        record.seq = record.seq.upper()

    problems = 0
    for record in records:
        for problem in validate_cds(record):
            print(f"WARNING: {record.id}: {problem}", file=sys.stderr)
            problems += 1

    database = out_dir / 'example_database.fasta'
    SeqIO.write(records, database, 'fasta')
    print(f"Wrote {len(records)} reference CDS to {database} ({problems} warnings)")

    # ---- point-mutation reference proteins ---------------------------------
    protein_records = []
    for gene, per_organism in sorted(POINT_REFERENCES.items()):
        for organism, accession in sorted(per_organism.items()):
            sequence = proteins.get(accession)
            if sequence is None:
                print(f"WARNING: protein {accession} for {gene} not in AMRProt.fa",
                      file=sys.stderr)
                continue
            protein_records.append(SeqRecord(
                Seq(sequence),
                id=f"{gene}|{organism}|{accession}",
                description=f"{gene} point-mutation reference protein ({organism})",
            ))
    reference_proteins = out_dir / 'reference_proteins.faa'
    SeqIO.write(protein_records, reference_proteins, 'fasta')
    print(f"Wrote {len(protein_records)} reference proteins to {reference_proteins}")

    # ---- point mutations ----------------------------------------------------
    rows, dropped = build_point_mutations(amr_dir, proteins)
    mutations_file = out_dir / 'point_mutations.tsv'
    fieldnames = ['Gene', 'Organism', 'Reference_Protein', 'Position', 'Ref',
                  'Variant', 'Kind', 'Label', 'Class', 'Subclass', 'Drug_Scope',
                  'Source']
    with open(mutations_file, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    in_scope = sum(1 for row in rows
                   if row['Drug_Scope'] in ('fosfomycin', 'ceftazidime-avibactam'))
    print(f"Wrote {len(rows)} validated point mutations to {mutations_file} "
          f"({in_scope} curated for fosfomycin or ceftazidime-avibactam, "
          f"{len(rows) - in_scope} for other drugs, kept as context)")
    for gene, symbol, reason in dropped:
        print(f"  dropped {gene} {symbol}: {reason}", file=sys.stderr)

    # ---- blaKPC allele table ------------------------------------------------
    kpc_rows = build_kpc_allele_table(amr_cds, catalog)
    alleles_file = out_dir / 'blaKPC_alleles.tsv'
    with open(alleles_file, 'w', newline='') as handle:
        writer = csv.DictWriter(
            handle, delimiter='\t',
            fieldnames=['Allele', 'Changes_Ambler', 'NCBI_subclass',
                        'NCBI_product_name', 'Inhibitor_resistant'])
        writer.writeheader()
        writer.writerows(kpc_rows)
    print(f"Wrote {len(kpc_rows)} blaKPC alleles to {alleles_file}")

    version_file = out_dir / 'DATA_VERSION.txt'
    version_file.write_text(
        f"AMRFinderPlus database version: {version}\n"
        f"Generated by fos_cazavi/build_data.py\n"
    )
    print(f"Wrote {version_file}")


if __name__ == '__main__':
    main()
