"""Protein-level variant calling against a reference coding sequence.

The functions here answer one question: given a stretch of assembly that
aligns to a reference gene, what does the encoded protein actually look like
compared with that reference?

They are deliberately explicit about the things that make naive mutation
callers wrong:

* A reference codon is compared with the query codon the *alignment* pairs it
  with, so an indel upstream does not shift every downstream call.
* Insertions, deletions, premature stop codons and frameshifts are reported as
  such instead of being silently turned into substitutions.
* Class A beta-lactamase positions are reported in standardised Ambler
  numbering, which omits positions 58 and 253, so the labels this tool prints
  match the KPC literature and NCBI allele definitions.
"""

from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices

# Class A beta-lactamase (ABL) standard numbering omits residues 58 and 253.
# Sequence index -> Ambler position.
_ABL_GAPS = (58, 253)

# How much shorter than the reference a product may be before it is called
# truncated.  5% leaves room for in-frame indels and for reference alleles that
# differ slightly in length, while still catching real nonsense mutations.
TRUNCATION_TOLERANCE = 0.05


def sequential_to_ambler(position):
    """Map a 1-based sequential residue index in a class A beta-lactamase to
    its standardised Ambler position."""
    if position < 58:
        return position
    if position < 252:
        return position + 1
    return position + 2


def ambler_to_sequential(position):
    """Inverse of :func:`sequential_to_ambler`.

    Returns ``None`` for the two Ambler positions that do not exist in the
    sequence (58 and 253)."""
    if position in _ABL_GAPS:
        return None
    if position < 58:
        return position
    if position < 253:
        return position - 1
    return position - 2


def _protein_aligner():
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = 'global'
    return aligner


def translate_cds(nucleotide_seq):
    """Translate a CDS and describe how it terminates.

    Returns ``(protein, info)`` where ``protein`` is the translation up to (not
    including) the first stop codon, and ``info`` records whether the sequence
    was in frame, whether translation hit a stop codon, and where.
    """
    sequence = str(nucleotide_seq).upper().replace('-', '')
    remainder = len(sequence) % 3
    trimmed = sequence[:len(sequence) - remainder] if remainder else sequence

    full = str(Seq(trimmed).translate()) if trimmed else ''
    stop_index = full.find('*')
    protein = full[:stop_index] if stop_index >= 0 else full

    return protein, {
        'in_frame': remainder == 0,
        'trailing_bases': remainder,
        'has_stop_codon': stop_index >= 0,
        'stop_position': stop_index + 1 if stop_index >= 0 else None,
        'translated_length': len(protein),
    }


def compare_proteins(reference_protein, query_protein, numbering='sequential'):
    """Amino-acid differences of ``query_protein`` relative to ``reference_protein``.

    Returns a list of dicts, each with ``kind`` (substitution / deletion /
    insertion), ``position`` (in the requested numbering), and ``label`` (e.g.
    ``D179Y``, ``E166del``, ``insPNK@269``).
    """
    if not reference_protein or not query_protein:
        return []

    renumber = sequential_to_ambler if numbering == 'ambler' else (lambda p: p)

    alignment = _protein_aligner().align(reference_protein, query_protein)[0]
    reference_row, query_row = str(alignment[0]), str(alignment[1])

    changes = []
    reference_position = 0
    insertion = []
    insertion_after = 0

    def flush_insertion():
        if insertion:
            residues = ''.join(insertion)
            changes.append({
                'kind': 'insertion',
                'position': renumber(insertion_after) if insertion_after else 0,
                'ref': '',
                'alt': residues,
                'label': f"ins{residues}@{renumber(insertion_after) if insertion_after else 0}",
            })
            insertion.clear()

    for reference_aa, query_aa in zip(reference_row, query_row):
        if reference_aa == '-':
            insertion.append(query_aa)
            continue
        flush_insertion()
        reference_position += 1
        insertion_after = reference_position
        if query_aa == reference_aa:
            continue
        position = renumber(reference_position)
        if query_aa == '-':
            changes.append({
                'kind': 'deletion',
                'position': position,
                'ref': reference_aa,
                'alt': '',
                'label': f"{reference_aa}{position}del",
            })
        else:
            changes.append({
                'kind': 'substitution',
                'position': position,
                'ref': reference_aa,
                'alt': query_aa,
                'label': f"{reference_aa}{position}{query_aa}",
            })
    flush_insertion()
    return changes


def collapse_runs(changes):
    """Merge consecutive single-residue deletions into one range label.

    ``E166del`` + ``L167del`` becomes ``E166_L167del``, which is how such
    variants are named in the literature (e.g. KPC-66).
    """
    collapsed = []
    run = []

    def flush():
        if not run:
            return
        if len(run) == 1:
            collapsed.append(run[0])
        else:
            first, last = run[0], run[-1]
            collapsed.append({
                'kind': 'deletion',
                'position': first['position'],
                'ref': ''.join(c['ref'] for c in run),
                'alt': '',
                'label': f"{first['ref']}{first['position']}_{last['ref']}{last['position']}del",
            })
        run.clear()

    for change in changes:
        if change['kind'] == 'deletion':
            if run and change['position'] == run[-1]['position'] + 1:
                run.append(change)
                continue
            flush()
            run.append(change)
        else:
            flush()
            collapsed.append(change)
    flush()
    return collapsed


def call_variants(reference_protein, query_nucleotides, numbering='sequential'):
    """Full variant call for one gene copy.

    ``query_nucleotides`` is the assembly sequence spanning the gene, oriented
    in the reference's sense and starting at the reference's first codon.

    Returns a dict with the protein changes plus the loss-of-function signals
    (premature stop, frameshift, truncation) that matter for resistance.
    """
    protein, translation = translate_cds(query_nucleotides)
    changes = collapse_runs(compare_proteins(reference_protein, protein, numbering))

    reference_length = len(reference_protein)
    # A product shortened by less than this is not treated as truncated: an
    # in-frame deletion of a few residues (KPC-66 loses two) legitimately
    # yields a slightly shorter protein and must not be called a knockout.
    minimum_length = reference_length * (1 - TRUNCATION_TOLERANCE)

    truncated = bool(protein) and len(protein) < minimum_length
    # A length difference that is not a whole number of codons means the
    # reading frame is broken somewhere in this gene copy.
    frameshift = not translation['in_frame']

    premature_stop = None
    if translation['has_stop_codon'] and len(protein) < minimum_length:
        premature_stop = len(protein) + 1

    return {
        'protein': protein,
        'protein_length': len(protein),
        'reference_length': reference_length,
        'changes': changes,
        'frameshift': frameshift,
        'premature_stop': premature_stop,
        'truncated': truncated,
        'loss_of_function': bool(frameshift or premature_stop or truncated),
        'translation': translation,
    }


def loss_of_function_label(call, numbering='sequential'):
    """A short human-readable label for a loss-of-function call, or None."""
    if call['premature_stop']:
        renumber = sequential_to_ambler if numbering == 'ambler' else (lambda p: p)
        return (f"premature stop at residue {renumber(call['premature_stop'])} "
                f"({call['protein_length']}/{call['reference_length']} aa)")
    if call['frameshift']:
        return (f"frameshift ({call['translation']['trailing_bases']} bp out of frame "
                f"relative to the reference coding sequence)")
    if call['truncated']:
        return (f"truncated product ({call['protein_length']}/"
                f"{call['reference_length']} aa)")
    return None


def extract_gene_span(contig_sequence, qstart, qend, sstart, send, reference_length):
    """Recover the full reference-length span of a gene from a contig.

    BLAST reports the aligned region only.  When the alignment starts at
    reference base ``sstart`` > 1, the corresponding assembly bases upstream are
    added so translation begins at the gene's first codon; the same is done at
    the 3' end.  The returned sequence is oriented in the reference's sense.

    Returns ``(sequence, complete)`` where ``complete`` is False when the contig
    runs out before the gene does (a gene broken by a contig boundary).
    """
    query_forward = qstart <= qend
    subject_forward = sstart <= send
    same_strand = query_forward == subject_forward

    query_low, query_high = min(qstart, qend), max(qstart, qend)
    subject_low, subject_high = min(sstart, send), max(sstart, send)

    missing_left = subject_low - 1
    missing_right = reference_length - subject_high

    if same_strand:
        start = query_low - 1 - missing_left
        end = query_high + missing_right
    else:
        start = query_low - 1 - missing_right
        end = query_high + missing_left

    complete = start >= 0 and end <= len(contig_sequence)
    start = max(start, 0)
    end = min(end, len(contig_sequence))

    sequence = contig_sequence[start:end]
    if not same_strand:
        sequence = str(Seq(str(sequence)).reverse_complement())
    return str(sequence), complete
