#!/usr/bin/env python3
"""
Create synthetic test genomes for resistance detection validation.
Simulates E. coli / K. pneumoniae assemblies with resistance genes and mutations.

Genes and mutations tested (from gist motroy/c10bfb64aa1f8a5a86bb2d3986c8c724):

FOS (Fosfomycin) resistance:
  Plasmidic:
    - fosA3, fosA4, fosA5, fosA11
  Chromosomal:
    - fosAKP (I91V)
    - murA (D369N, L370I)
    - uhpB (G469R), uhpC (F384L), uhpA, uhpT, glpT
    - lon (Q558*), cyaA, ptsI
    - galU (R282V)

CAZAVI (Ceftazidime-Avibactam) resistance:
  Plasmidic:
    - blaKPC-2 (D179Y)
    - blaKPC-3 (D179Y, V240G, D179Y/T243M)
    - blaKPC-31 (D179Y)
    - blaKPC-190
    - blaOXA-48 (P68A, P68A/Y211S)
    - blaCMY-178 (N70T)
    - blaSHV-12
  Chromosomal:
    - ompK36, ompK35 (porin mutations)
    - acrB, mexR, nalD (efflux pump regulators)
    - ftsI/PBP3
    - envZ
"""
import random
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Seed for reproducibility
random.seed(42)

# Standard codon per amino acid
AA_TO_CODON = {
    'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT',
    'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'K': 'AAA', 'L': 'CTG',
    'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA', 'R': 'CGT',
    'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT',
    '*': 'TAA'
}


def random_seq(length):
    """Generate random DNA sequence"""
    return ''.join(random.choices('ACGT', k=length))


def create_contig(seq, name, description=""):
    """Create a SeqRecord for a contig"""
    return SeqRecord(Seq(seq), id=name, description=description)


def introduce_mutation(seq, codon_pos, new_codon):
    """
    Introduce a codon mutation at specified position (1-indexed codon).
    codon_pos: 1-indexed codon position
    new_codon: the mutant codon to insert
    """
    start_nt = (codon_pos - 1) * 3
    if start_nt + 3 > len(seq):
        raise ValueError(
            f"Codon position {codon_pos} (nt {start_nt}-{start_nt+3}) "
            f"is beyond sequence length {len(seq)}"
        )
    return seq[:start_nt] + new_codon + seq[start_nt + 3:]


def load_references(db_path):
    """Load reference sequences from FASTA database"""
    refs = {}
    for record in SeqIO.parse(db_path, "fasta"):
        refs[record.id] = str(record.seq)
        # Short name: first part before underscore, e.g. 'fosA3'
        gene_name = record.id.split('_')[0]
        refs[gene_name] = str(record.seq)
        # Also store by full id minus '' suffix
        clean = record.id
        refs[clean] = str(record.seq)
    return refs


def embed_in_contig(gene_seq, padding=50000):
    """Embed gene sequence in random flanking sequence"""
    return random_seq(padding) + gene_seq + random_seq(padding)


def write_genome(contigs, filename):
    """Write contigs to FASTA file"""
    with open(filename, "w") as f:
        SeqIO.write(contigs, f, "fasta")
    total_bp = sum(len(c.seq) for c in contigs)
    print(f"  Created {filename}: {len(contigs)} contigs, {total_bp:,} bp")


# ─── Negative control ────────────────────────────────────────────────────────

def create_negative_control():
    """Create E. coli K-12 simulation without resistance genes"""
    contigs = []
    for i in range(1, 11):
        size = random.randint(400000, 600000)
        seq = random_seq(size)
        contigs.append(create_contig(seq, f"contig_{i}", f"length={size}"))
    return contigs


# ─── FOS plasmidic genes ─────────────────────────────────────────────────────

def create_fos_plasmidic(refs, gene_ids):
    """
    Create genome with one or more plasmid-mediated fosA-family genes.
    gene_ids: list of database keys, e.g. ['fosA3', 'fosA4']
    """
    contigs = []
    for i in range(1, 6):
        size = random.randint(400000, 600000)
        contigs.append(create_contig(random_seq(size), f"contig_{i}", f"length={size}"))
    for gene_id in gene_ids:
        seq = refs.get(gene_id)
        if seq is None:
            print(f"  WARNING: {gene_id} not found in database, skipping")
            continue
        plasmid_seq = embed_in_contig(seq)
        short = gene_id
        contigs.append(create_contig(plasmid_seq, f"contig_plasmid_{short}",
                                     f"contains={short}"))
    return contigs


# ─── FOS chromosomal mutations ───────────────────────────────────────────────

def _make_fos_chromo_contigs(refs, mutations):
    """
    Build a set of contigs each containing one mutated chromosomal gene.
    mutations: list of (ref_key, codon_pos, new_codon, label)
    """
    contigs = []
    for i in range(1, 5):
        size = random.randint(400000, 600000)
        contigs.append(create_contig(random_seq(size), f"contig_{i}", f"length={size}"))

    for ref_key, codon_pos, new_codon, label in mutations:
        seq = refs.get(ref_key)
        if seq is None:
            print(f"  WARNING: {ref_key} not found in database, skipping")
            continue
        try:
            mut_seq = introduce_mutation(seq, codon_pos, new_codon)
            contig_seq = embed_in_contig(mut_seq)
            safe = label.replace('/', '_').replace('*', 'stop')
            gene_name = ref_key
            contigs.append(create_contig(
                contig_seq,
                f"contig_chr_{gene_name}_{safe}",
                f"contains={label}"
            ))
        except ValueError as e:
            print(f"  WARNING: {label} - {e}")
    return contigs


def create_chromosomal_fos_mutations_primers(refs):
    """
    FOS chromosomal mutations validated by primers.tsv:
    uhpB G469R, uhpC F384L, galU R282V, lon Q558*
    """
    mutations = [
        ('uhpB',  469, 'CGT', 'uhpB_G469R'),
        ('uhpC',  384, 'CTG', 'uhpC_F384L'),
        ('galU',  282, 'GTG', 'galU_R282V'),
        ('lon',   558, 'TAG', 'lon_Q558stop'),
    ]
    return _make_fos_chromo_contigs(refs, mutations)


def create_chromosomal_fos_murA(refs):
    """murA D369N and L370I (fosfomycin target)"""
    mutations = [
        ('murA', 369, 'AAT', 'murA_D369N'),
        ('murA', 370, 'ATT', 'murA_L370I'),
    ]
    return _make_fos_chromo_contigs(refs, mutations)


def create_chromosomal_fos_uhpT(refs):
    """uhpT transport mutations (fosfomycin uptake)"""
    mutations = [
        ('uhpT',  55, 'GAT', 'uhpT_G55D'),
        ('uhpT', 198, 'TAA', 'uhpT_W198stop'),
        ('uhpT', 258, 'TAA', 'uhpT_E258stop'),
        ('uhpT', 350, 'TAA', 'uhpT_W350stop'),
    ]
    return _make_fos_chromo_contigs(refs, mutations)


def create_chromosomal_fos_glpT(refs):
    """glpT transport mutations (fosfomycin uptake)"""
    mutations = [
        ('glpT',  44, 'TAA', 'glpT_E44stop'),
        ('glpT',  88, 'TAA', 'glpT_W88stop'),
        ('glpT',  90, 'GAT', 'glpT_G90D'),
        ('glpT', 234, 'TAA', 'glpT_W234stop'),
        ('glpT', 362, 'TGT', 'glpT_R362C'),
    ]
    return _make_fos_chromo_contigs(refs, mutations)


def create_chromosomal_fos_uhpABC(refs):
    """uhpA, uhpB, uhpC two-component system mutations"""
    mutations = [
        ('uhpA',  54, 'AAT', 'uhpA_D54N'),
        ('uhpA', 139, 'TGT', 'uhpA_R139C'),
        ('uhpB', 350, 'TAT', 'uhpB_H350Y'),
        ('uhpC', 384, 'CTG', 'uhpC_F384L'),
    ]
    return _make_fos_chromo_contigs(refs, mutations)


def create_chromosomal_fos_lon_cyaA_ptsI(refs):
    """lon Q558*, cyaA G463D, ptsI H191Y"""
    mutations = [
        ('lon',  558, 'TAG', 'lon_Q558stop'),
        ('cyaA', 463, 'GAT', 'cyaA_G463D'),
        ('ptsI', 191, 'TAT', 'ptsI_H191Y'),
    ]
    return _make_fos_chromo_contigs(refs, mutations)


def create_chromosomal_fosAKP(refs):
    """fosAKP I91V (chromosomal Klebsiella fosfomycin resistance)"""
    mutations = [
        ('fosAKP', 91, 'GTT', 'fosAKP_I91V'),
    ]
    return _make_fos_chromo_contigs(refs, mutations)


# ─── CAZAVI plasmidic genes ──────────────────────────────────────────────────

def create_cazavi_blakpc2_d179y(refs):
    """blaKPC-2 with D179Y"""
    contigs = []
    for i in range(1, 6):
        contigs.append(create_contig(random_seq(random.randint(400000, 600000)),
                                     f"contig_{i}", ""))
    seq = refs.get('blaKPC-2')
    if seq:
        mut = introduce_mutation(seq, 179, 'TAT')   # D179Y: GAT->TAT
        contigs.append(create_contig(embed_in_contig(mut),
                                     'contig_plasmid_blaKPC2',
                                     'contains=blaKPC-2_D179Y'))
    return contigs


def create_cazavi_blakpc3(refs):
    """blaKPC-3 with D179Y, V240G, and D179Y/T243M"""
    contigs = []
    for i in range(1, 4):
        contigs.append(create_contig(random_seq(random.randint(400000, 600000)),
                                     f"contig_{i}", ""))
    seq = refs.get('blaKPC-3')
    if seq:
        # D179Y alone
        mut_d179y = introduce_mutation(seq, 179, 'TAT')
        contigs.append(create_contig(embed_in_contig(mut_d179y),
                                     'contig_plasmid_blaKPC3_D179Y',
                                     'contains=blaKPC-3_D179Y'))
        # V240G alone
        mut_v240g = introduce_mutation(seq, 240, 'GGT')
        contigs.append(create_contig(embed_in_contig(mut_v240g),
                                     'contig_plasmid_blaKPC3_V240G',
                                     'contains=blaKPC-3_V240G'))
        # D179Y + T243M double mutant
        mut_double = introduce_mutation(seq, 179, 'TAT')
        mut_double = introduce_mutation(mut_double, 243, 'ATG')
        contigs.append(create_contig(embed_in_contig(mut_double),
                                     'contig_plasmid_blaKPC3_D179Y_T243M',
                                     'contains=blaKPC-3_D179Y_T243M'))
    return contigs


def create_cazavi_blakpc31_d179y(refs):
    """blaKPC-31 with D179Y"""
    contigs = []
    for i in range(1, 4):
        contigs.append(create_contig(random_seq(random.randint(400000, 600000)),
                                     f"contig_{i}", ""))
    seq = refs.get('blaKPC-31')
    if seq:
        mut = introduce_mutation(seq, 179, 'TAT')
        contigs.append(create_contig(embed_in_contig(mut),
                                     'contig_plasmid_blaKPC31',
                                     'contains=blaKPC-31_D179Y'))
    return contigs


def create_cazavi_blakpc190(refs):
    """blaKPC-190 (wildtype for detection, no specific resistance mutation)"""
    contigs = []
    for i in range(1, 4):
        contigs.append(create_contig(random_seq(random.randint(400000, 600000)),
                                     f"contig_{i}", ""))
    seq = refs.get('blaKPC-190')
    if seq:
        contigs.append(create_contig(embed_in_contig(seq),
                                     'contig_plasmid_blaKPC190',
                                     'contains=blaKPC-190'))
    return contigs


def create_cazavi_blaoxa48(refs):
    """blaOXA-48 with P68A and P68A/Y211S"""
    contigs = []
    for i in range(1, 4):
        contigs.append(create_contig(random_seq(random.randint(400000, 600000)),
                                     f"contig_{i}", ""))
    seq = refs.get('blaOXA-48')
    if seq:
        # P68A alone
        mut_p68a = introduce_mutation(seq, 68, 'GCT')
        contigs.append(create_contig(embed_in_contig(mut_p68a),
                                     'contig_plasmid_blaOXA48_P68A',
                                     'contains=blaOXA-48_P68A'))
        # P68A + Y211S double
        mut_double = introduce_mutation(seq, 68, 'GCT')
        mut_double = introduce_mutation(mut_double, 211, 'TCT')
        contigs.append(create_contig(embed_in_contig(mut_double),
                                     'contig_plasmid_blaOXA48_P68A_Y211S',
                                     'contains=blaOXA-48_P68A_Y211S'))
    return contigs


def create_cazavi_blacmy178_n70t(refs):
    """blaCMY-178 with N70T (Asn70Thr)"""
    contigs = []
    for i in range(1, 4):
        contigs.append(create_contig(random_seq(random.randint(400000, 600000)),
                                     f"contig_{i}", ""))
    seq = refs.get('blaCMY-178')
    if seq:
        mut = introduce_mutation(seq, 70, 'ACT')   # N70T: AAT->ACT
        contigs.append(create_contig(embed_in_contig(mut),
                                     'contig_plasmid_blaCMY178',
                                     'contains=blaCMY-178_N70T'))
    return contigs


def create_cazavi_blashv12(refs):
    """blaSHV-12 (wildtype reference for gene detection)"""
    contigs = []
    for i in range(1, 4):
        contigs.append(create_contig(random_seq(random.randint(400000, 600000)),
                                     f"contig_{i}", ""))
    seq = refs.get('blaSHV-12')
    if seq:
        contigs.append(create_contig(embed_in_contig(seq),
                                     'contig_plasmid_blaSHV12',
                                     'contains=blaSHV-12'))
    return contigs


# ─── CAZAVI chromosomal mutations ────────────────────────────────────────────

def create_cazavi_porins(refs):
    """ompK35 and ompK36 porin mutations"""
    mutations = [
        ('ompK35', 134, 'GAT', 'ompK35_G134D'),
        ('ompK35', 135, 'TAA', 'ompK35_D135stop'),
        ('ompK35', 181, 'GGT', 'ompK35_D181G'),
        ('ompK36', 134, 'GAT', 'ompK36_G134D') if 'ompK36' in refs else None,
    ]
    mutations = [m for m in mutations if m is not None]
    return _make_fos_chromo_contigs(refs, mutations)


def create_cazavi_efflux_regulators(refs):
    """mexR and nalD efflux pump regulator mutations"""
    mutations = [
        ('mexR', 69, 'GGT', 'mexR_W69G'),
        ('mexR', 75, 'GTT', 'mexR_A75V'),
        ('nalD', 153, 'TAA', 'nalD_Q153stop'),
        ('nalD', 174, 'CGT', 'nalD_L174R'),
    ]
    return _make_fos_chromo_contigs(refs, mutations)


def create_cazavi_acrB(refs):
    """acrB efflux pump mutations"""
    mutations = [
        ('acrB', 617, 'GAT', 'acrB_G617D'),
        ('acrB', 626, 'CTG', 'acrB_F626L'),
        ('acrB', 628, 'ACT', 'acrB_A628T'),
    ]
    return _make_fos_chromo_contigs(refs, mutations)


def create_cazavi_ftsI(refs):
    """ftsI/PBP3 penicillin-binding protein mutations"""
    mutations = [
        ('ftsI', 333, 'GTT', 'ftsI_A333V'),
        ('ftsI', 350, 'TGT', 'ftsI_Y350C'),
        ('ftsI', 357, 'AAT', 'ftsI_S357N'),
    ]
    return _make_fos_chromo_contigs(refs, mutations)


def create_cazavi_envZ(refs):
    """envZ two-component sensor mutations"""
    mutations = [
        ('envZ', 244, 'TCT', 'envZ_G244S'),
        ('envZ', 324, 'ATT', 'envZ_T324I'),
    ]
    return _make_fos_chromo_contigs(refs, mutations)


# ─── Combined multi-resistance ────────────────────────────────────────────────

def create_multi_resistance(refs):
    """Create genome with multiple resistance genes and mutations"""
    contigs = []

    fosa3 = refs.get('fosA3', refs.get('fosA3'))
    blakpc3 = refs.get('blaKPC-3', refs.get('blaKPC-3'))
    blaoxa48 = refs.get('blaOXA-48', refs.get('blaOXA-48'))

    for i in range(1, 6):
        size = random.randint(400000, 600000)
        contigs.append(create_contig(random_seq(size), f"contig_{i}", f"length={size}"))

    if fosa3:
        p1_seq = embed_in_contig(fosa3, padding=20000)
        contigs.append(create_contig(p1_seq, "contig_plasmid1_fosA3", "contains=fosA3"))

    if blakpc3:
        blakpc3_mut = introduce_mutation(blakpc3, 179, "TAT")   # D179Y
        blakpc3_mut = introduce_mutation(blakpc3_mut, 243, "ATG")  # T243M
        p2_seq = embed_in_contig(blakpc3_mut, padding=25000)
        contigs.append(create_contig(p2_seq, "contig_plasmid2_blaKPC3",
                                     "contains=blaKPC-3_D179Y_T243M"))

    if blaoxa48:
        p3_seq = embed_in_contig(blaoxa48, padding=30000)
        contigs.append(create_contig(p3_seq, "contig_plasmid3_blaOXA48",
                                     "contains=blaOXA-48"))

    return contigs


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    os.makedirs("test_genomes", exist_ok=True)

    db_path = "fos_cazavi/data/example_database.fasta"
    print(f"Loading references from {db_path}...")
    refs = load_references(db_path)
    print(f"  Loaded {len(refs)} reference sequences\n")

    print("Creating synthetic test genomes...")
    print("=" * 60)

    scenarios = [
        # (label, function, output_filename, expected_detection)
        ("1.  Negative control",
         create_negative_control,
         "ecoli_negative_control.fasta",
         "No resistance genes"),

        # FOS plasmidic
        ("2.  fosA3 (plasmidic)",
         lambda: create_fos_plasmidic(refs, ['fosA3']),
         "ecoli_fosA3.fasta",
         "fosA3"),

        ("3.  fosA4 (plasmidic)",
         lambda: create_fos_plasmidic(refs, ['fosA4']),
         "ecoli_fosA4.fasta",
         "fosA4"),

        ("4.  fosA5 (plasmidic)",
         lambda: create_fos_plasmidic(refs, ['fosA5']),
         "ecoli_fosA5.fasta",
         "fosA5"),

        ("5.  fosA11 (plasmidic)",
         lambda: create_fos_plasmidic(refs, ['fosA11']),
         "ecoli_fosA11.fasta",
         "fosA11"),

        # FOS chromosomal
        ("6.  uhpB G469R + uhpC F384L + galU R282V + lon Q558*",
         lambda: create_chromosomal_fos_mutations_primers(refs),
         "ecoli_chromosomal_fos_primers.fasta",
         "uhpB G469R, uhpC F384L, galU R282V, lon Q558*"),

        ("7.  murA D369N + L370I",
         lambda: create_chromosomal_fos_murA(refs),
         "ecoli_murA.fasta",
         "murA D369N, murA L370I"),

        ("8.  uhpT transport mutations",
         lambda: create_chromosomal_fos_uhpT(refs),
         "ecoli_uhpT.fasta",
         "uhpT G55D, W198*, E258*, W350*"),

        ("9.  glpT transport mutations",
         lambda: create_chromosomal_fos_glpT(refs),
         "ecoli_glpT.fasta",
         "glpT E44*, W88*, G90D, W234*, R362C"),

        ("10. uhpA + uhpB + uhpC mutations",
         lambda: create_chromosomal_fos_uhpABC(refs),
         "ecoli_uhpABC.fasta",
         "uhpA D54N, uhpA R139C, uhpB H350Y, uhpC F384L"),

        ("11. lon Q558* + cyaA G463D + ptsI H191Y",
         lambda: create_chromosomal_fos_lon_cyaA_ptsI(refs),
         "ecoli_lon_cyaA_ptsI.fasta",
         "lon Q558*, cyaA G463D, ptsI H191Y"),

        ("12. fosAKP I91V (chromosomal Klebsiella)",
         lambda: create_chromosomal_fosAKP(refs),
         "kpneumo_fosAKP.fasta",
         "fosAKP I91V"),

        # CAZAVI plasmidic
        ("13. blaKPC-2 D179Y",
         lambda: create_cazavi_blakpc2_d179y(refs),
         "kpneumo_blaKPC2_D179Y.fasta",
         "blaKPC-2 D179Y"),

        ("14. blaKPC-3 D179Y / V240G / D179Y+T243M",
         lambda: create_cazavi_blakpc3(refs),
         "kpneumo_blaKPC3_variants.fasta",
         "blaKPC-3 D179Y, V240G, D179Y/T243M"),

        ("15. blaKPC-31 D179Y",
         lambda: create_cazavi_blakpc31_d179y(refs),
         "kpneumo_blaKPC31_D179Y.fasta",
         "blaKPC-31 D179Y"),

        ("16. blaKPC-190",
         lambda: create_cazavi_blakpc190(refs),
         "kpneumo_blaKPC190.fasta",
         "blaKPC-190"),

        ("17. blaOXA-48 P68A + P68A/Y211S",
         lambda: create_cazavi_blaoxa48(refs),
         "kpneumo_blaOXA48_variants.fasta",
         "blaOXA-48 P68A, P68A/Y211S"),

        ("18. blaCMY-178 N70T",
         lambda: create_cazavi_blacmy178_n70t(refs),
         "kpneumo_blaCMY178_N70T.fasta",
         "blaCMY-178 N70T"),

        ("19. blaSHV-12",
         lambda: create_cazavi_blashv12(refs),
         "kpneumo_blaSHV12.fasta",
         "blaSHV-12"),

        # CAZAVI chromosomal
        ("20. ompK35 + ompK36 porin mutations",
         lambda: create_cazavi_porins(refs),
         "kpneumo_porins.fasta",
         "ompK35 G134D, D135*, D181G"),

        ("21. acrB efflux mutations",
         lambda: create_cazavi_acrB(refs),
         "kpneumo_acrB.fasta",
         "acrB G617D, F626L, A628T"),

        ("22. mexR + nalD efflux regulator mutations",
         lambda: create_cazavi_efflux_regulators(refs),
         "kpneumo_efflux_regulators.fasta",
         "mexR W69G, A75V; nalD Q153*, L174R"),

        ("23. ftsI/PBP3 mutations",
         lambda: create_cazavi_ftsI(refs),
         "kpneumo_ftsI.fasta",
         "ftsI A333V, Y350C, S357N"),

        ("24. envZ two-component sensor mutations",
         lambda: create_cazavi_envZ(refs),
         "kpneumo_envZ.fasta",
         "envZ G244S, T324I"),

        # Multi-resistance
        ("25. Multi-resistance (fosA3 + blaKPC-3 D179Y/T243M + blaOXA-48)",
         lambda: create_multi_resistance(refs),
         "ecoli_multi_resistance.fasta",
         "fosA3, blaKPC-3 D179Y/T243M, blaOXA-48"),
    ]

    print("\nExpected detection results:")
    for label, func, fname, expected in scenarios:
        print(f"\n{label}")
        try:
            contigs = func()
            write_genome(contigs, f"test_genomes/{fname}")
            print(f"  Expected: {expected}")
        except Exception as e:
            print(f"  ERROR: {e}")

    print("\n" + "=" * 60)
    print(f"Test genomes created in test_genomes/")


if __name__ == "__main__":
    main()
