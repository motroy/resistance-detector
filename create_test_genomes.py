#!/usr/bin/env python3
"""
Create synthetic test genomes for resistance detection validation.
Simulates E. coli assemblies with resistance genes and chromosomal mutations.

Mutations tested (from primers.tsv):
- uhpB G469R (fosfomycin resistance)
- uhpC F384L (fosfomycin resistance)
- galU R282V (fosfomycin resistance)
- lon Q558* (fosfomycin resistance - stop codon)
- blaKPC-3 D179Y, T243M (ceftazidime-avibactam resistance)
"""
import random
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Seed for reproducibility
random.seed(42)

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
    return seq[:start_nt] + new_codon + seq[start_nt + 3:]

# Load reference sequences from the database
def load_references(db_path):
    """Load reference sequences from FASTA database"""
    refs = {}
    for record in SeqIO.parse(db_path, "fasta"):
        # Extract gene name (first part before underscore)
        gene_name = record.id.split('_')[0]
        refs[record.id] = str(record.seq)
        refs[gene_name] = str(record.seq)  # Also index by short name
    return refs

def create_negative_control():
    """Create E. coli K-12 simulation without resistance genes"""
    contigs = []
    for i in range(1, 11):
        size = random.randint(400000, 600000)
        seq = random_seq(size)
        contigs.append(create_contig(seq, f"contig_{i}", f"length={size}"))
    return contigs

def create_positive_fosa3(refs):
    """Create genome with fosA3 plasmid-mediated gene"""
    contigs = []
    fosa3 = refs.get('fosA3_reference', refs.get('fosA3'))

    for i in range(1, 8):
        size = random.randint(400000, 600000)
        contigs.append(create_contig(random_seq(size), f"contig_{i}", f"length={size}"))

    # Plasmid with fosA3
    plasmid_size = 50000
    insert_pos = 25000
    plasmid_seq = random_seq(insert_pos) + fosa3 + random_seq(plasmid_size - insert_pos - len(fosa3))
    contigs.append(create_contig(plasmid_seq, "contig_plasmid_fosA3", "contains=fosA3"))
    return contigs

def create_positive_blakpc3_d179y(refs):
    """Create genome with blaKPC-3 containing D179Y mutation"""
    contigs = []
    blakpc3 = refs.get('blaKPC-3_reference', refs.get('blaKPC'))

    for i in range(1, 8):
        size = random.randint(400000, 600000)
        contigs.append(create_contig(random_seq(size), f"contig_{i}", f"length={size}"))

    # D179Y mutation: codon 179 GAT (D) -> TAT (Y)
    blakpc3_mut = introduce_mutation(blakpc3, 179, "TAT")

    plasmid_size = 60000
    insert_pos = 30000
    plasmid_seq = random_seq(insert_pos) + blakpc3_mut + random_seq(plasmid_size - insert_pos - len(blakpc3_mut))
    contigs.append(create_contig(plasmid_seq, "contig_plasmid_blaKPC3", "contains=blaKPC-3_D179Y"))
    return contigs

def create_chromosomal_fos_mutations(refs):
    """
    Create genome with chromosomal fosfomycin resistance mutations.
    Mutations from primers.tsv:
    - uhpB G469R
    - uhpC F384L
    - galU R282V
    - lon Q558*
    """
    contigs = []

    # Get reference sequences
    uhpB = refs.get('uhpB_reference', refs.get('uhpB'))
    uhpC = refs.get('uhpC_reference', refs.get('uhpC'))
    galU = refs.get('galU_reference', refs.get('galU'))
    lon = refs.get('lon_reference', refs.get('lon'))

    # Introduce mutations
    # uhpB G469R: GGT (G) -> CGT (R) at position 469
    uhpB_mut = introduce_mutation(uhpB, 469, "CGT") if uhpB and len(uhpB) >= 469*3 else uhpB

    # uhpC F384L: TTC (F) -> CTG (L) at position 384
    uhpC_mut = introduce_mutation(uhpC, 384, "CTG") if uhpC and len(uhpC) >= 384*3 else uhpC

    # galU R282V: CGT (R) -> GTG (V) at position 282
    galU_mut = introduce_mutation(galU, 282, "GTG") if galU and len(galU) >= 282*3 else galU

    # lon Q558*: CAG (Q) -> TAG (stop) at position 558
    lon_mut = introduce_mutation(lon, 558, "TAG") if lon and len(lon) >= 558*3 else lon

    # Create chromosome-like contigs with mutant genes embedded
    for i in range(1, 5):
        size = random.randint(400000, 600000)
        contigs.append(create_contig(random_seq(size), f"contig_{i}", f"length={size}"))

    # Contig with uhpB G469R
    if uhpB_mut:
        seq = random_seq(50000) + uhpB_mut + random_seq(50000)
        contigs.append(create_contig(seq, "contig_chr_uhpB", "contains=uhpB_G469R"))

    # Contig with uhpC F384L
    if uhpC_mut:
        seq = random_seq(50000) + uhpC_mut + random_seq(50000)
        contigs.append(create_contig(seq, "contig_chr_uhpC", "contains=uhpC_F384L"))

    # Contig with galU R282V
    if galU_mut:
        seq = random_seq(50000) + galU_mut + random_seq(50000)
        contigs.append(create_contig(seq, "contig_chr_galU", "contains=galU_R282V"))

    # Contig with lon Q558*
    if lon_mut:
        seq = random_seq(50000) + lon_mut + random_seq(50000)
        contigs.append(create_contig(seq, "contig_chr_lon", "contains=lon_Q558*"))

    return contigs

def create_multi_resistance(refs):
    """Create genome with multiple resistance genes and mutations"""
    contigs = []

    fosa3 = refs.get('fosA3_reference', refs.get('fosA3'))
    blakpc3 = refs.get('blaKPC-3_reference', refs.get('blaKPC'))
    blaoxa48 = refs.get('blaOXA-48_reference', refs.get('blaOXA'))

    for i in range(1, 6):
        size = random.randint(400000, 600000)
        contigs.append(create_contig(random_seq(size), f"contig_{i}", f"length={size}"))

    # Plasmid 1: fosA3
    p1_seq = random_seq(20000) + fosa3 + random_seq(25000)
    contigs.append(create_contig(p1_seq, "contig_plasmid1_fosA3", "contains=fosA3"))

    # Plasmid 2: blaKPC-3 with D179Y and T243M
    blakpc3_mut = introduce_mutation(blakpc3, 179, "TAT")  # D179Y
    blakpc3_mut = introduce_mutation(blakpc3_mut, 243, "ATG")  # T243M
    p2_seq = random_seq(25000) + blakpc3_mut + random_seq(30000)
    contigs.append(create_contig(p2_seq, "contig_plasmid2_blaKPC3", "contains=blaKPC-3_D179Y_T243M"))

    # Plasmid 3: blaOXA-48
    p3_seq = random_seq(30000) + blaoxa48 + random_seq(35000)
    contigs.append(create_contig(p3_seq, "contig_plasmid3_blaOXA48", "contains=blaOXA-48"))

    return contigs

def write_genome(contigs, filename):
    """Write contigs to FASTA file"""
    with open(filename, "w") as f:
        SeqIO.write(contigs, f, "fasta")
    total_bp = sum(len(c.seq) for c in contigs)
    print(f"  Created {filename}: {len(contigs)} contigs, {total_bp:,} bp")

def main():
    os.makedirs("test_genomes", exist_ok=True)

    # Load reference sequences from database
    db_path = "data/example_database.fasta"
    print(f"Loading references from {db_path}...")
    refs = load_references(db_path)
    print(f"  Loaded {len(refs)//2} reference sequences\n")

    print("Creating synthetic test genomes...")
    print("=" * 60)

    # 1. Negative control
    print("\n1. Negative control (no resistance genes)")
    negative = create_negative_control()
    write_genome(negative, "test_genomes/ecoli_negative_control.fasta")

    # 2. fosA3 positive
    print("\n2. Plasmid-mediated fosA3")
    fosa3 = create_positive_fosa3(refs)
    write_genome(fosa3, "test_genomes/ecoli_fosA3_positive.fasta")

    # 3. blaKPC-3 D179Y
    print("\n3. blaKPC-3 with D179Y mutation")
    kpc = create_positive_blakpc3_d179y(refs)
    write_genome(kpc, "test_genomes/ecoli_blaKPC3_D179Y.fasta")

    # 4. Chromosomal FOS mutations (from primers.tsv)
    print("\n4. Chromosomal fosfomycin resistance mutations")
    print("   (uhpB G469R, uhpC F384L, galU R282V, lon Q558*)")
    chromosomal = create_chromosomal_fos_mutations(refs)
    write_genome(chromosomal, "test_genomes/ecoli_chromosomal_fos.fasta")

    # 5. Multi-resistance
    print("\n5. Multi-resistance (fosA3 + blaKPC-3 D179Y/T243M + blaOXA-48)")
    multi = create_multi_resistance(refs)
    write_genome(multi, "test_genomes/ecoli_multi_resistance.fasta")

    print("\n" + "=" * 60)
    print("Test genomes created successfully!")
    print("\nExpected detection results:")
    print("  1. ecoli_negative_control.fasta  -> No resistance genes")
    print("  2. ecoli_fosA3_positive.fasta    -> fosA3")
    print("  3. ecoli_blaKPC3_D179Y.fasta     -> blaKPC-3 + D179Y mutation")
    print("  4. ecoli_chromosomal_fos.fasta   -> uhpB, uhpC, galU, lon with mutations")
    print("  5. ecoli_multi_resistance.fasta  -> fosA3, blaKPC-3, blaOXA-48")

if __name__ == "__main__":
    main()
