#!/usr/bin/env python3
"""Build synthetic test genomes with a known, declared genotype.

Each scenario is a small assembly carrying reference genes, some of them
deliberately mutated, together with the result the pipeline is expected to
produce.  The expectations are written next to the genomes as
``expected_results.tsv``, so a run can be checked against ground truth rather
than against a previous run of the same code.

Positions are given in the numbering the tool reports: standardised Ambler
numbering for class A beta-lactamases, sequential residue numbering elsewhere.

Usage:
    python3 create_test_genomes.py [output_directory]
"""

import csv
import random
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fos_cazavi.references import DEFAULT_DATABASE          # noqa: E402
from fos_cazavi.variants import ambler_to_sequential        # noqa: E402

random.seed(42)

AA_TO_CODON = {
    'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT',
    'H': 'CAT', 'I': 'ATT', 'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAT',
    'P': 'CCT', 'Q': 'CAA', 'R': 'CGT', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT',
    'W': 'TGG', 'Y': 'TAT', '*': 'TAA',
}


def load_references():
    return {record.id: str(record.seq) for record in SeqIO.parse(DEFAULT_DATABASE, 'fasta')}


def random_sequence(length):
    return ''.join(random.choices('ACGT', k=length))


def embed(sequence, flank=800):
    """Place a gene inside random flanking sequence, as it sits in a contig."""
    return random_sequence(flank) + sequence + random_sequence(flank)


def substitute(sequence, residue_position, new_residue):
    """Replace the codon at a 1-based sequential residue position."""
    index = (residue_position - 1) * 3
    return sequence[:index] + AA_TO_CODON[new_residue] + sequence[index + 3:]


def delete(sequence, residue_position, count=1):
    index = (residue_position - 1) * 3
    return sequence[:index] + sequence[index + 3 * count:]


def insert(sequence, residue_position, residues):
    index = residue_position * 3
    codons = ''.join(AA_TO_CODON[residue] for residue in residues)
    return sequence[:index] + codons + sequence[index:]


def ambler_substitute(sequence, ambler_position, new_residue):
    return substitute(sequence, ambler_to_sequential(ambler_position), new_residue)


def ambler_delete(sequence, ambler_position, count=1):
    return delete(sequence, ambler_to_sequential(ambler_position), count)


# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------
# Each scenario is (name, builder, expected fosfomycin call, expected CAZ/AVI
# call, organism, note).  The builder returns {contig_name: sequence}.

def scenarios(references):
    kpc2 = references['blaKPC-2']

    def genome(**genes):
        return {f"contig_{name}": embed(sequence) for name, sequence in genes.items()}

    return [
        (
            'susceptible_kpc2',
            lambda: genome(blaKPC=kpc2, fosAKP=references['fosAKP'],
                           uhpT=references['uhpT']),
            'Susceptible', 'Susceptible', 'Klebsiella_pneumoniae',
            'Wild-type KPC-2 is inhibited by avibactam; intrinsic fosAKP is not '
            'acquired resistance',
        ),
        (
            'kpc3_background_only',
            lambda: genome(blaKPC=ambler_substitute(kpc2, 274, 'Y')),
            'Susceptible', 'Susceptible', 'Klebsiella_pneumoniae',
            'H274Y defines KPC-3 and is not an avibactam escape mutation',
        ),
        (
            'kpc33_d179y',
            lambda: genome(blaKPC=ambler_substitute(kpc2, 179, 'Y')),
            'Susceptible', 'Resistant', 'Klebsiella_pneumoniae',
            'D179Y is the most frequently reported ceftazidime-avibactam '
            'resistance substitution',
        ),
        (
            'kpc31_d179y_on_kpc3',
            lambda: genome(blaKPC=ambler_substitute(
                ambler_substitute(kpc2, 274, 'Y'), 179, 'Y')),
            'Susceptible', 'Resistant', 'Klebsiella_pneumoniae',
            'KPC-31 = KPC-3 + D179Y',
        ),
        (
            'kpc66_omega_loop_deletion',
            lambda: genome(blaKPC=ambler_delete(
                ambler_substitute(kpc2, 274, 'Y'), 166, 2)),
            'Susceptible', 'Resistant', 'Klebsiella_pneumoniae',
            'In-frame Omega-loop deletion (KPC-66)',
        ),
        (
            'kpc_v240g',
            lambda: genome(blaKPC=ambler_substitute(kpc2, 240, 'G')),
            'Susceptible', 'Resistant', 'Klebsiella_pneumoniae',
            'V240G in the active-site region',
        ),
        (
            'kpc_novel_omega_loop_change',
            lambda: genome(blaKPC=ambler_substitute(kpc2, 172, 'P')),
            'Susceptible', 'Indeterminate', 'Klebsiella_pneumoniae',
            'Undocumented change in a hotspot must not be asserted either way',
        ),
        (
            'metallo_betalactamase_ndm1',
            lambda: genome(blaNDM=references['blaNDM-1'], blaKPC=kpc2),
            'Susceptible', 'Resistant', 'Klebsiella_pneumoniae',
            'Avibactam does not inhibit metallo-beta-lactamases',
        ),
        (
            'oxa48_only',
            lambda: genome(blaOXA=references['blaOXA-48']),
            'Susceptible', 'Susceptible', 'Klebsiella_pneumoniae',
            'OXA-48 is inhibited by avibactam',
        ),
        (
            'acquired_fosa3',
            lambda: genome(fosA3=references['fosA3'], blaKPC=kpc2),
            'Resistant', 'Susceptible', 'Escherichia',
            'Acquired fosA3 inactivates fosfomycin',
        ),
        (
            'uhpt_nonsense',
            lambda: genome(uhpT=substitute(references['uhpT'], 150, '*'),
                           glpT=references['glpT']),
            'Resistant', 'Susceptible', 'Escherichia',
            'Premature stop in uhpT abolishes fosfomycin uptake',
        ),
        (
            'glpt_frameshift',
            lambda: genome(glpT=references['glpT'][:600] + 'A' + references['glpT'][600:],
                           uhpT=references['uhpT']),
            'Resistant', 'Susceptible', 'Escherichia',
            'Single-base insertion breaks the glpT reading frame',
        ),
        (
            'wildtype_transporters',
            lambda: genome(uhpT=references['uhpT'], glpT=references['glpT'],
                           murA=references['murA']),
            'Susceptible', 'Susceptible', 'Escherichia',
            'Intact uptake genes and no acquired enzyme',
        ),
        (
            'multi_mechanism',
            lambda: genome(blaKPC=ambler_substitute(kpc2, 179, 'Y'),
                           fosA3=references['fosA3'],
                           fosAKP=references['fosAKP'],
                           blaCTX_M=references['blaCTX-M-15']),
            'Resistant', 'Resistant', 'Klebsiella_pneumoniae',
            'Both drugs compromised in one isolate: acquired fosA3 alongside '
            'the intrinsic chromosomal fosA, plus a KPC escape variant',
        ),
        (
            'lone_fosa_no_intrinsic_copy',
            lambda: genome(fosA3=references['fosA3']),
            'Indeterminate', 'Susceptible', 'Klebsiella_pneumoniae',
            'A single fosA hit with no intrinsic chromosomal copy recognised '
            'cannot be told apart from a divergent chromosomal enzyme',
        ),
    ]


def main():
    output_dir = Path(sys.argv[1] if len(sys.argv) > 1 else 'test_data')
    output_dir.mkdir(parents=True, exist_ok=True)

    references = load_references()
    rows = []

    for name, builder, fos_call, cazavi_call, organism, note in scenarios(references):
        contigs = builder()
        records = [SeqRecord(Seq(sequence), id=contig, description=f"scenario={name}")
                   for contig, sequence in contigs.items()]
        path = output_dir / f"{name}.fasta"
        SeqIO.write(records, path, 'fasta')
        rows.append({
            'Genome': path.name,
            'Organism': organism,
            'Expected_Fosfomycin': fos_call,
            'Expected_Ceftazidime_Avibactam': cazavi_call,
            'Rationale': note,
        })
        print(f"Wrote {path} ({len(records)} contigs)")

    expected = output_dir / 'expected_results.tsv'
    with open(expected, 'w', newline='') as handle:
        writer = csv.DictWriter(
            handle, delimiter='\t',
            fieldnames=['Genome', 'Organism', 'Expected_Fosfomycin',
                        'Expected_Ceftazidime_Avibactam', 'Rationale'])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {expected} with {len(rows)} expected results")


if __name__ == '__main__':
    main()
