"""Gene detection and variant calling from an assembly, via BLAST.

Each BLAST hit is turned into a full protein-level call: the gene span is
recovered from the contig (including the parts of the gene the alignment did
not cover), translated, and compared with the reference protein.  What comes
out is what the assembly actually encodes - substitutions, indels, premature
stops, frameshifts and truncations - rather than a lookup at fixed positions.
"""

import subprocess
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import betalactamase
from .references import (
    ACQUIRED_PREFIXES, INTRINSIC_GENES, IN_SCOPE_DRUGS, canonical_gene_name,
    gene_family, is_acquired_gene, load_point_mutations, load_reference_proteins,
    mutation_lookup_key, mutation_scope,
)
from .utils import check_dependencies
from .variants import call_variants, extract_gene_span, loss_of_function_label

# BLAST output columns, in the order requested below.
_BLAST_FIELDS = ('qseqid sseqid pident length mismatch gapopen qstart qend '
                 'sstart send evalue bitscore slen qlen')


class BlastDetector:
    """Find reference genes in an assembly and call their variants."""

    def __init__(self, assembly, database, output_prefix, min_identity=90,
                 min_coverage=80, mutation_db_file=None, organism=None, threads=1):
        self.assembly = assembly
        self.database = database
        self.output_prefix = output_prefix
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.organism = organism
        self.threads = max(1, int(threads or 1))
        self.results = []
        self.detected_genes = []

        self.reference_cds = {
            record.id: str(record.seq)
            for record in SeqIO.parse(self.database, 'fasta')
        }
        self.contigs = {
            record.id: str(record.seq)
            for record in SeqIO.parse(self.assembly, 'fasta')
        }
        self.point_mutations = load_point_mutations(mutation_db_file)
        self.reference_proteins = load_reference_proteins()

    # ------------------------------------------------------------------ BLAST

    def prepare_database(self):
        db_files = [f"{self.database}.{extension}" for extension in ('nhr', 'nin', 'nsq')]
        if all(Path(path).exists() for path in db_files):
            return
        print(f"Creating BLAST database from {self.database}...")
        try:
            subprocess.run(['makeblastdb', '-in', self.database, '-dbtype', 'nucl'],
                           check=True, capture_output=True)
        except subprocess.CalledProcessError as error:
            print(f"ERROR creating database: {error.stderr.decode()}", file=sys.stderr)
            sys.exit(1)

    def run_blast(self):
        if not check_dependencies(['blastn', 'makeblastdb']):
            sys.exit(1)
        self.prepare_database()

        print(f"Running BLAST search (min_id={self.min_identity}%, "
              f"min_cov={self.min_coverage}%)...")

        command = [
            'blastn',
            '-query', self.assembly,
            '-db', self.database,
            '-outfmt', f'6 {_BLAST_FIELDS}',
            '-evalue', '1e-20',
            # One reference gene can have many near-identical relatives in the
            # database (every fosA, every MBL).  Keeping all of them and
            # choosing by bitscore afterwards is what makes allele assignment
            # correct; capping the list here would hide the true best match.
            '-max_target_seqs', '5000',
            '-perc_identity', str(max(self.min_identity - 5, 70)),
            '-num_threads', str(self.threads),
        ]
        try:
            completed = subprocess.run(command, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as error:
            print(f"ERROR running BLAST: {error.stderr}", file=sys.stderr)
            sys.exit(1)

        with open(f"{self.output_prefix}_blast.txt", 'w') as handle:
            handle.write(completed.stdout)

        hits = self.parse_blast_output(completed.stdout)
        kept = self.filter_redundant_hits(hits)
        return self._prefer_intrinsic_naming(hits, kept)

    def parse_blast_output(self, blast_output):
        hits = []
        for line in blast_output.strip().split('\n'):
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 14:
                continue
            (query_id, subject_id, identity, length, _mismatch, _gapopen,
             qstart, qend, sstart, send, _evalue, bitscore, slen, _qlen) = fields[:14]

            subject_length = int(slen)
            coverage = (abs(int(send) - int(sstart)) + 1) / subject_length * 100

            if float(identity) < self.min_identity or coverage < self.min_coverage:
                continue

            hits.append({
                'query_id': query_id,
                'gene': subject_id,
                'identity': float(identity),
                'coverage': coverage,
                'length': int(length),
                'qstart': int(qstart),
                'qend': int(qend),
                'sstart': int(sstart),
                'send': int(send),
                'slen': subject_length,
                'bitscore': float(bitscore),
            })

        print(f"Found {len(hits)} gene hits passing thresholds")
        return hits

    def filter_redundant_hits(self, hits):
        """Keep the single best reference match per genomic locus.

        Alleles of one family are near-identical, so a locus hits many of them.
        Ranking by bitscore (then identity, then coverage) and dropping
        overlapping lower-scoring hits leaves one call per locus against the
        closest reference - which is what allele assignment needs.
        """
        if not hits:
            return []

        hits.sort(key=lambda hit: (hit['bitscore'], hit['identity'], hit['coverage']),
                  reverse=True)

        kept = []
        for hit in hits:
            start, end = min(hit['qstart'], hit['qend']), max(hit['qstart'], hit['qend'])
            redundant = False
            for previous in kept:
                if hit['query_id'] != previous['query_id']:
                    continue
                previous_start = min(previous['qstart'], previous['qend'])
                previous_end = max(previous['qstart'], previous['qend'])
                overlap = max(0, min(end, previous_end) - max(start, previous_start) + 1)
                if overlap / (end - start + 1) * 100 > 50:
                    redundant = True
                    break
            if not redundant:
                kept.append(hit)

        print(f"Filtered to {len(kept)} hits after redundancy check")
        return kept

    # Reference identity gap, in percentage points, within which a locus's
    # winning acquired-family hit and a competing intrinsic-gene hit at the
    # same span are treated as indistinguishable (see _prefer_intrinsic_naming).
    _INTRINSIC_NAMING_TOLERANCE = 3.0

    def _prefer_intrinsic_naming(self, all_hits, kept_hits):
        """Correct a specific, verified mis-naming: a locus's single copy of a
        species' intrinsic gene (fosAKP/fosA6/fosA_PA1129) can score a
        marginally lower bitscore than a horizontally-acquired allele from the
        same family (fosA5, fosA10, ...), because the intrinsic reference is
        one strain's sequence and other lineages' native copy is sometimes, by
        chance, a percentage point or two closer to a different catalogued
        allele than to that one reference strain.

        Plain best-bitscore selection then reports the genome's own
        chromosomal gene as an ambiguous or acquired allele - which matters a
        great deal here, because that naming is exactly what the phenotype
        logic uses to tell an acquired enzyme from the intrinsic one. Confirmed
        on real assemblies: four genomes in the ESKAPE-fosfomycin GOLD
        validation set had their single fosA copy (Copy_Number 1, i.e. there is
        nothing else it could be) named fosA5-like or fosA10-like this way,
        each within 2 percentage points of the intrinsic fosAKP reference at
        the identical query span.

        For every kept hit that lost this way, if a competing intrinsic-gene
        hit exists at (near enough) the same span and within
        ``_INTRINSIC_NAMING_TOLERANCE`` percentage points of identity, that
        intrinsic hit is substituted in place of the acquired one - restoring
        both the correct name and the correct reference for translation and
        variant calling. This never *invents* resistance calls; it only ever
        prevents an ordinary chromosomal-allele match from masquerading as one.
        """
        intrinsic_hits = [hit for hit in all_hits if hit['gene'] in INTRINSIC_GENES]
        if not intrinsic_hits:
            return kept_hits

        intrinsic_families = {gene_family(g) for g in INTRINSIC_GENES}

        corrected = []
        for hit in kept_hits:
            # Only hits already named as *acquired* members of a family that
            # has a curated intrinsic counterpart are in scope (currently:
            # fosA). A hit already named as the intrinsic gene, or belonging
            # to an unrelated family, passes through untouched.
            if hit['gene'] in INTRINSIC_GENES or gene_family(hit['gene']) not in intrinsic_families:
                corrected.append(hit)
                continue

            start, end = min(hit['qstart'], hit['qend']), max(hit['qstart'], hit['qend'])
            best_intrinsic = None
            for candidate in intrinsic_hits:
                if candidate['query_id'] != hit['query_id']:
                    continue
                candidate_start = min(candidate['qstart'], candidate['qend'])
                candidate_end = max(candidate['qstart'], candidate['qend'])
                overlap = max(0, min(end, candidate_end) - max(start, candidate_start) + 1)
                if overlap / (end - start + 1) * 100 <= 50:
                    continue
                if hit['identity'] - candidate['identity'] > self._INTRINSIC_NAMING_TOLERANCE:
                    continue
                if best_intrinsic is None or candidate['identity'] > best_intrinsic['identity']:
                    best_intrinsic = candidate

            if best_intrinsic is not None:
                print(f"  {hit['query_id']}: {hit['gene']} ({hit['identity']:.2f}%) "
                      f"renamed to intrinsic {best_intrinsic['gene']} "
                      f"({best_intrinsic['identity']:.2f}%) - same locus, within "
                      f"{self._INTRINSIC_NAMING_TOLERANCE} points, and this is the "
                      f"genome's only copy at this locus")
                corrected.append(best_intrinsic)
            else:
                corrected.append(hit)

        return corrected

    # ------------------------------------------------------------- variant call

    def analyze_hits(self, hits):
        print("Analyzing hits and calling variants...")

        for hit in hits:
            # The BLAST subject ID (raw_gene) is what the reference CDS is
            # actually stored under - for a gene with more than one
            # organism-specific reference (see GENE_ALIASES) that is not the
            # same as the canonical name every other check uses.
            raw_gene = hit['gene']
            gene = canonical_gene_name(raw_gene)
            family = gene_family(gene)
            contig = self.contigs.get(hit['query_id'], '')
            reference_cds = self.reference_cds.get(raw_gene, '')
            if not contig or not reference_cds:
                continue

            reference_protein = str(Seq(reference_cds).translate()).rstrip('*')
            sequence, complete = extract_gene_span(
                contig, hit['qstart'], hit['qend'], hit['sstart'], hit['send'],
                len(reference_cds))

            numbering = 'ambler' if family in ('blaKPC', 'blaSHV', 'blaOXA', 'blaCTX-M') else 'sequential'
            call = call_variants(reference_protein, sequence, numbering=numbering)
            change_labels = [change['label'] for change in call['changes']]

            allele_row = None
            allele_name = gene
            if family == 'blaKPC':
                allele_row = betalactamase.identify_kpc_allele(change_labels)
                allele_name = betalactamase.describe_kpc_result(change_labels, allele_row)
            elif is_acquired_gene(gene) and change_labels:
                # Allele-level typing is only done for blaKPC.  Elsewhere the
                # matched reference is merely the closest of many near-identical
                # alleles - blaIMP-18 and blaIMP-99 are 99.7% identical - so
                # naming a specific allele would assert more than the data
                # supports.  The family is reported instead; the Gene column
                # still records which reference was closest.
                allele_name = f"{family}-like"

            # Curated chromosomal mutations are numbered against the organism's
            # own reference protein, so they are matched against a second call
            # made against that protein rather than against the nucleotide
            # database entry, whose numbering can differ.
            curated_call = call
            curated_reference = self.reference_proteins.get((gene, self.organism))
            if curated_reference and not is_acquired_gene(gene):
                curated_call = call_variants(curated_reference[1], sequence)

            reported, unreported, curated_rows = self._classify_changes(
                gene, family, call, curated_call)

            result = {
                'contig': hit['query_id'],
                'gene': gene,
                'family': family,
                'allele': allele_name,
                'allele_row': allele_row,
                'identity': f"{hit['identity']:.2f}",
                'coverage': f"{hit['coverage']:.2f}",
                'complete': complete,
                'acquired': is_acquired_gene(gene),
                'changes': change_labels,
                'reported_mutations': reported,
                'curated_mutations': curated_rows,
                'other_drug_mutations': [
                    f"{row['Label']} ({row['Subclass'] or row['Class']})"
                    for row in curated_rows
                    if mutation_scope(row) not in IN_SCOPE_DRUGS
                ],
                'other_changes': unreported,
                'loss_of_function': call['loss_of_function'],
                'lof_description': loss_of_function_label(call, numbering),
                'mutations': ','.join(reported) if reported else '-',
                'sequence': sequence,
                'start': hit['qstart'],
                'end': hit['qend'],
                'call': call,
            }
            self.results.append(result)
            self.detected_genes.append(SeqRecord(
                Seq(sequence),
                id=f"{hit['query_id']}_{gene}",
                description=(f"allele={allele_name} identity={hit['identity']:.2f}% "
                             f"coverage={hit['coverage']:.2f}% "
                             f"changes={';'.join(change_labels) if change_labels else 'none'}"),
            ))

        copy_counts = Counter(result['gene'] for result in self.results)
        for result in self.results:
            result['copy_number'] = copy_counts[result['gene']]

    def _classify_changes(self, gene, family, call, curated_call=None):
        """Split called changes into ones this tool is willing to report as
        resistance mutations, and everything else.

        For an acquired beta-lactamase every change is meaningful, because the
        reference is the canonical allele of that same gene.  For a chromosomal
        gene a difference from the reference is usually just natural sequence
        variation, so only positions curated for the sample's organism are
        reported - and only when the organism was actually declared.

        Returns ``(reported_labels, other_labels, curated_rows)``.  The curated
        rows carry the drug class each mutation was curated for, which is what
        lets the phenotype logic tell a ceftazidime-avibactam mutation from,
        say, a tigecycline one.
        """
        labels = [change['label'] for change in call['changes']]

        if family == 'blaKPC':
            # Only the changes the ceftazidime-avibactam assessor recognises are
            # reported as resistance mutations; the rest (e.g. H274Y, which
            # simply distinguishes KPC-3 from KPC-2) stay in the full change list.
            assessment = betalactamase.assess_kpc_changes(labels)
            flagged = [label for label in labels
                       if any(label in item for item in assessment['evidence'])]
            return flagged, [label for label in labels if label not in flagged], []

        if is_acquired_gene(gene):
            return labels, [], []

        curated = self.point_mutations.get((gene, self.organism), {}) if self.organism else {}
        source = curated_call or call
        reported, other, rows = [], [], []
        for change in source['changes']:
            entry = curated.get(mutation_lookup_key(change))
            if entry is None:
                other.append(change['label'])
                continue
            rows.append(entry)
            if mutation_scope(entry) in IN_SCOPE_DRUGS:
                reported.append(entry['Label'])
            else:
                # Curated, but for a drug this tool does not report on (e.g.
                # tigecycline, carbapenems).  Kept visible as context rather
                # than presented as a fosfomycin/ceftazidime-avibactam finding.
                other.append(change['label'])
        return reported, other, rows

    # ----------------------------------------------------------------- reports

    def write_report(self):
        report_file = f"{self.output_prefix}_results.tsv"
        print(f"Writing results to {report_file}...")

        with open(report_file, 'w') as handle:
            handle.write('\t'.join([
                'Contig', 'Gene', 'Allele', 'Identity%', 'Coverage%', 'Complete',
                'Reported_Mutations', 'Other_Drug_Mutations',
                'All_Protein_Changes', 'Loss_Of_Function',
                'Method', 'Copy_Number',
            ]) + '\n')

            for result in self.results:
                handle.write('\t'.join([
                    result['contig'],
                    result['gene'],
                    result['allele'],
                    result['identity'],
                    result['coverage'],
                    'yes' if result['complete'] else 'no (contig boundary)',
                    result['mutations'],
                    ';'.join(result['other_drug_mutations'])
                    if result['other_drug_mutations'] else '-',
                    ';'.join(result['changes']) if result['changes'] else '-',
                    result['lof_description'] or '-',
                    'BLAST',
                    str(result['copy_number']),
                ]) + '\n')

        print(f"Detected {len(self.results)} genes")

    def write_sequences(self):
        fasta_file = f"{self.output_prefix}_genes.fasta"
        if self.detected_genes:
            print(f"Writing {len(self.detected_genes)} gene sequences to {fasta_file}...")
            SeqIO.write(self.detected_genes, fasta_file, 'fasta')
        else:
            print("No genes detected to write")

    def run(self):
        hits = self.run_blast()
        if hits:
            self.analyze_hits(hits)
        self.write_report()
        self.write_sequences()
        return self.results


def run_acquired_detection(assembly, database, output, min_id, min_cov,
                           mutation_db=None, organism=None, threads=1):
    detector = BlastDetector(assembly, database, output, min_id, min_cov,
                             mutation_db, organism, threads)
    return detector.run()
