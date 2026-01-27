#!/usr/bin/env python3
"""
Script to download real ESBL genomes and simulate acquired genes/mutations
to test fos-cazavi detection capabilities.
"""

import os
import sys
import subprocess
import random
import time
from pathlib import Path
from Bio import Entrez, SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

# Configure Entrez
Entrez.email = "tool_test@example.com"

def setup_logging():
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler("simulation.log")
        ]
    )
    return logging.getLogger("simulation")

logger = setup_logging()

def run_command(cmd, desc):
    logger.info(f"Running: {desc}")
    try:
        subprocess.run(cmd, check=True, shell=isinstance(cmd, str))
        logger.info(f"  Success: {desc}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"  Failed: {desc} (Error: {e})")
        return False

def generate_database():
    """Run fos-cazavi create-db to get full reference database"""
    output_prefix = "resistance_db"
    if Path(f"{output_prefix}.fasta").exists():
        logger.info(f"Database {output_prefix}.fasta already exists, skipping creation.")
        return output_prefix

    logger.info("Generating resistance database (this may take a while)...")
    cmd = [
        sys.executable, "-m", "fos_cazavi.cli", "create-db",
        "--email", Entrez.email,
        "--output", output_prefix
    ]
    if run_command(cmd, "fos-cazavi create-db"):
        return output_prefix
    else:
        logger.error("Failed to create database. Cannot proceed.")
        sys.exit(1)

def search_and_download_genomes():
    """Download 2 E. coli and 2 K. pneumoniae genomes"""
    Path("downloaded_genomes").mkdir(exist_ok=True)
    
    targets = [
        {
            "organism": "Escherichia coli",
            "query": 'Escherichia coli[Organism] AND "complete genome"[Title] AND ("blaCTX-M" OR "blaTEM")',
            "count": 2,
            "prefix": "ecoli"
        },
        {
            "organism": "Klebsiella pneumoniae",
            "query": 'Klebsiella pneumoniae[Organism] AND "complete genome"[Title] AND ("blaKPC" OR "blaNDM")',
            "count": 2,
            "prefix": "kpneumo"
        }
    ]
    
    downloaded_files = []
    
    for target in targets:
        logger.info(f"Searching for {target['organism']} genomes...")
        try:
            handle = Entrez.esearch(db="nucleotide", term=target["query"], retmax=target["count"], sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            
            ids = record["IdList"]
            logger.info(f"  Found IDs: {ids}")
            
            for i, uid in enumerate(ids):
                filename = f"downloaded_genomes/{target['prefix']}_{i+1}.fasta"
                if Path(filename).exists():
                    logger.info(f"  File {filename} exists, skipping download.")
                    downloaded_files.append(filename)
                    continue
                
                logger.info(f"  Downloading {uid} to {filename}...")
                try:
                    handle = Entrez.efetch(db="nucleotide", id=uid, rettype="fasta", retmode="text")
                    with open(filename, "w") as f:
                        f.write(handle.read())
                    handle.close()
                    downloaded_files.append(filename)
                    # Be nice to NCBI
                    time.sleep(1)
                except Exception as e:
                    logger.error(f"  Error downloading {uid}: {e}")
                    
        except Exception as e:
            logger.error(f"Search failed for {target['organism']}: {e}")
            
    return downloaded_files

def load_sequences(fasta_file):
    """Load sequences from FASTA into a dict"""
    seqs = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Key by simple name (e.g., blaKPC-3) if possible
        name = record.id
        if "reference sequence" in record.description:
             # Try to extract the gene name we want
             parts = record.description.split()
             if parts:
                 name = parts[0]
        
        # Also store by ID just in case
        seqs[record.id] = str(record.seq)
        seqs[name] = str(record.seq)
    return seqs

def introduce_mutation(seq, pos, new_aa, ref_aa):
    """
    Introduce a mutation at amino acid position `pos` (1-based).
    This assumes `seq` is a coding sequence (nucleotides).
    """
    seq_obj = Seq.Seq(seq)
    protein = seq_obj.translate()
    
    if pos > len(protein):
        logger.warning(f"  Position {pos} out of range for protein length {len(protein)}")
        return seq
    
    original_aa = protein[pos-1]
    if original_aa != ref_aa:
        # Just log a warning but proceed (sequence might differ slightly from ref)
        pass # logger.warning(f"  Ref AA at {pos} is {original_aa}, expected {ref_aa}")

    # Determine codon index
    codon_start = (pos - 1) * 3
    current_codon = seq[codon_start:codon_start+3]
    
    # Simple strategy: find a codon for new_aa that is closest to current?
    # Or just pick a common one.
    from Bio.Data import CodonTable
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    
    # Find a codon for new_aa
    target_codons = [c for c, aa in standard_table.forward_table.items() if aa == new_aa]
    if new_aa == '*':
        target_codons = standard_table.stop_codons
        
    if not target_codons:
        logger.warning(f"  No codon found for {new_aa}")
        return seq
        
    new_codon = target_codons[0] # Pick first one
    
    mutated_seq = seq[:codon_start] + new_codon + seq[codon_start+3:]
    return mutated_seq

def create_random_seq(length):
    return ''.join(random.choices('ACGT', k=length))

def simulate_genomes(genome_files, db_prefix):
    """Inject genes and mutations into downloaded genomes"""
    logger.info("Simulating genes and mutations...")
    
    db_fasta = f"{db_prefix}.fasta"
    refs = load_sequences(db_fasta)
    
    simulated_files = []
    
    # Define scenarios for the 4 genomes
    scenarios = [
        {
            # Genome 1 (E. coli): fosA3 + uhpB G469R
            "file_idx": 0,
            "genes": ["fosA3"],
            "mutations": [("uhpB", 469, "R", "G")]
        },
        {
            # Genome 2 (E. coli): blaKPC-3 (D179Y) + galU R282V
            "file_idx": 1,
            "genes": [], # blaKPC is added with mutation
            "mutations": [("blaKPC-3", 179, "Y", "D"), ("galU", 282, "V", "R")]
        },
        {
            # Genome 3 (K. pneumo): blaOXA-48 + lon Q558*
            "file_idx": 2,
            "genes": ["blaOXA-48"],
            "mutations": [("lon", 558, "*", "Q")]
        },
        {
            # Genome 4 (K. pneumo): blaKPC-3 (T243M) + uhpC F384L
            "file_idx": 3,
            "genes": [],
            "mutations": [("blaKPC-3", 243, "M", "T"), ("uhpC", 384, "L", "F")]
        }
    ]
    
    Path("simulated_genomes").mkdir(exist_ok=True)

    for i, scenario in enumerate(scenarios):
        if scenario["file_idx"] >= len(genome_files):
            logger.warning(f"Not enough downloaded genomes for scenario {i+1}")
            continue
            
        original_file = genome_files[scenario["file_idx"]]
        base_name = Path(original_file).stem
        output_file = f"simulated_genomes/{base_name}_simulated.fasta"
        
        logger.info(f"Processing {original_file} -> {output_file}")
        
        # Read original contigs
        contigs = list(SeqIO.parse(original_file, "fasta"))
        
        # Add genes
        for gene in scenario["genes"]:
            if gene in refs:
                seq = refs[gene]
                # Embed in a small contig
                contig_seq = create_random_seq(2000) + seq + create_random_seq(2000)
                contigs.append(SeqRecord(
                    Seq.Seq(contig_seq),
                    id=f"simulated_{gene}",
                    description=f"contains {gene}"
                ))
                logger.info(f"  Added {gene}")
            else:
                logger.warning(f"  Gene {gene} not found in database")
                
        # Add mutations (as new contigs with the mutated gene)
        for gene, pos, new_aa, ref_aa in scenario["mutations"]:
            # Check if we have the reference
            ref_gene = gene
            if gene not in refs:
                # Try simple name (e.g. blaKPC-3 -> blaKPC)
                if gene.startswith("bla") and gene not in refs:
                     # This logic is tricky, assume full name provided matches DB
                     pass

            if ref_gene in refs:
                original_seq = refs[ref_gene]
                mutated_seq = introduce_mutation(original_seq, pos, new_aa, ref_aa)
                
                contig_seq = create_random_seq(2000) + mutated_seq + create_random_seq(2000)
                contigs.append(SeqRecord(
                    Seq.Seq(contig_seq),
                    id=f"simulated_{gene}_mut_{ref_aa}{pos}{new_aa}",
                    description=f"contains {gene} {ref_aa}{pos}{new_aa}"
                ))
                logger.info(f"  Added {gene} with mutation {ref_aa}{pos}{new_aa}")
            else:
                 logger.warning(f"  Gene {ref_gene} not found in database for mutation")

        # Write output
        SeqIO.write(contigs, output_file, "fasta")
        simulated_files.append(output_file)
        
    return simulated_files

def run_detection(simulated_files, db_prefix):
    """Run fos-cazavi-all on simulated genomes"""
    logger.info("Running detection on simulated genomes...")
    Path("results").mkdir(exist_ok=True)
    
    # We also need proteins file for miniprot if we want to test that
    # It usually comes with the repo in data/cazavi_proteins.fasta
    proteins_file = "data/cazavi_proteins.fasta"
    primers_file = "data/primers.tsv"
    
    for fasta in simulated_files:
        base_name = Path(fasta).stem
        output_base = f"results/{base_name}"
        
        cmd = [
            sys.executable, "-m", "fos_cazavi.cli", "fos-cazavi-all",
            "--assembly", fasta,
            "--database", f"{db_prefix}.fasta",
            "--output", output_base,
            "--mutations", f"{db_prefix}_mutations.tsv"
        ]
        
        if Path(proteins_file).exists():
            cmd.extend(["--proteins", proteins_file])
        if Path(primers_file).exists():
            cmd.extend(["--primers", primers_file])
            
        run_command(cmd, f"Analysis of {base_name}")

def main():
    # 1. Generate Database
    db_prefix = generate_database()
    
    # 2. Download Genomes
    genomes = search_and_download_genomes()
    if not genomes:
        logger.error("No genomes downloaded. Exiting.")
        sys.exit(1)
        
    # 3. Simulate
    simulated = simulate_genomes(genomes, db_prefix)
    
    # 4. Test Detection
    run_detection(simulated, db_prefix)
    
    logger.info("Done! Check 'results/' directory.")

if __name__ == "__main__":
    main()
