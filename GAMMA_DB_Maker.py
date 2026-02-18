import sys
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def NonStandard(input_fasta, output_nonstandard_fasta, output_NonStop, output_correct_fasta):
    """Makes a new fasta with all the truncated or nonstop genes"""
    genes = list(SeqIO.parse(input_fasta, 'fasta'))
    Out = open(output_nonstandard_fasta, 'w')
    Out_NonStop = open(output_NonStop, 'w')
    Out_Correct = open(output_correct_fasta, 'w')
    for gene in genes:
        pro = str(gene.seq.translate())
##        if len(pro) == 0:
##            print(gene.id)
##        if pro[-1] != '*':
##            SeqIO.write(gene, Out, 'fasta')
        if ('*' in pro[0:-1]) == True:
            SeqIO.write(gene, Out, 'fasta')
        elif pro[-1] != '*':
            SeqIO.write(gene, Out_NonStop, 'fasta')
        else:
            SeqIO.write(gene, Out_Correct, 'fasta')
    Out.close()
    Out_NonStop.close()
    Out_Correct.close()

def NonStandard_Rename(input_fasta, output_correct_fasta, output_incorrect):
    """Reads in a nonstandard fasta and outputs a corrected version"""
    genes = list(SeqIO.parse(input_fasta, 'fasta'))
    Out = open(output_correct_fasta, 'w')
    Out_Incorrect = open(output_incorrect, 'w')
    for gene in genes:
        pro = str(gene.seq.translate())
        pro1 = str(gene.seq[1:].translate())
        pro2 = str(gene.seq[2:].translate())
        gene_rev = gene.seq.reverse_complement()
        gene1_rev = gene.seq[0:-1].reverse_complement()
        gene2_rev = gene.seq[0:-2].reverse_complement()
        pro_rev = str(gene_rev.translate())
        pro1_rev = str(gene1_rev.translate())
        pro2_rev = str(gene2_rev.translate())
        if ('*' in pro[0:-1]) == False and pro[-1] == '*':
            print(gene.id + '\t' + 'Correct')
        if ('*' in pro1[0:-1]) == False and pro1[-1] == '*':
            print(gene.id + '\tframeshift 1')
            Gene = gene[1:]
            SeqIO.write(Gene, Out, 'fasta')
        elif ('*' in pro2[0:-1]) == False and pro2[-1] == '*':
            print(gene.id + '\tframeshift 2')
            Gene = gene[2:]
            SeqIO.write(Gene, Out, 'fasta')
        elif ('*' in pro_rev[0:-1]) == False and pro_rev[-1] == '*':
            print(gene.id + '\tReverse Complement')
            Gene = gene.seq.reverse_complement()
            Gene = SeqRecord(Gene, id=gene.id, description=gene.description)
            SeqIO.write(Gene, Out, 'fasta')
        elif ('*' in pro1_rev[0:-1]) == False and pro1_rev[-1] == '*':
            print(gene.id + '\tReverse Complement, frameshift 1')
            Gene = gene.seq[0:-1].reverse_complement()
            Gene = SeqRecord(Gene, id=gene.id, description=gene.description)
            SeqIO.write(Gene, Out, 'fasta')
        elif ('*' in pro2_rev[0:-1]) == False and pro2_rev[-1] == '*':
            print(gene.id + '\tReverse Complement, frameshift 2')
            Gene = gene.seq[0:-2].reverse_complement()
            Gene = SeqRecord(Gene, id=gene.id, description=gene.description)
            SeqIO.write(Gene, Out, 'fasta')
        else:
            print(gene.id +'\tNo correction')
            SeqIO.write(gene, Out_Incorrect, 'fasta')
    Out.close()
    Out_Incorrect.close()

def NonStandard_Rename_non_nonstop(input_fasta, output_correct_fasta, output_nonstop, output_incorrect):
    """Reads in a nonstandard fasta and outputs a corrected version"""
    genes = list(SeqIO.parse(input_fasta, 'fasta'))
    Out = open(output_correct_fasta, 'w')
    Out_Incorrect = open(output_incorrect, 'w')
    Out_NonStop = open(output_nonstop, 'w')
    for gene in genes:
        pro = str(gene.seq.translate())
        pro1 = str(gene.seq[1:].translate())
        pro2 = str(gene.seq[2:].translate())
        gene_rev = gene.seq.reverse_complement()
        gene1_rev = gene.seq[0:-1].reverse_complement()
        gene2_rev = gene.seq[0:-2].reverse_complement()
        pro_rev = str(gene_rev.translate())
        pro1_rev = str(gene1_rev.translate())
        pro2_rev = str(gene2_rev.translate())
        if ('*' in pro[0:-1]) == False:
            #print(gene.id + '\t' + 'Correct')
            if  pro[-1] == '*':
                SeqIO.write(gene, Out, 'fasta')
            else:
                SeqIO.write(gene, Out_NonStop, 'fasta')
        if ('*' in pro1[0:-1]) == False:
            #print(gene.id + '\tframeshift 1')
            Gene = gene[1:]
            if pro1[-1] == '*':
                SeqIO.write(Gene, Out, 'fasta')
            else:
                SeqIO.write(Gene, Out_NonStop, 'fasta')
        elif ('*' in pro2[0:-1]) == False:
            #print(gene.id + '\tframeshift 2')
            Gene = gene[2:]
            if pro2[-1] == '*':
                SeqIO.write(Gene, Out, 'fasta')
            else:
                SeqIO.write(Gene, Out_NonStop, 'fasta')
        elif ('*' in pro_rev[0:-1]) == False:
            #print(gene.id + '\tReverse Complement')
            Gene = gene.seq.reverse_complement()
            Gene = SeqRecord(Gene, id=gene.id, description=gene.description)
            if pro_rev[-1] == '*':
                SeqIO.write(Gene, Out, 'fasta')
            else:
                SeqIO.write(Gene, Out_NonStop, 'fasta')
        elif ('*' in pro1_rev[0:-1]) == False:
            #print(gene.id + '\tReverse Complement, frameshift 1')
            Gene = gene.seq[0:-1].reverse_complement()
            Gene = SeqRecord(Gene, id=gene.id, description=gene.description)
            if pro1_rev[-1] == '*':
                SeqIO.write(Gene, Out, 'fasta')
            else:
                SeqIO.write(Gene, Out_NonStop, 'fasta')
        elif ('*' in pro2_rev[0:-1]) == False:
            #print(gene.id + '\tReverse Complement, frameshift 2')
            Gene = gene.seq[0:-2].reverse_complement()
            Gene = SeqRecord(Gene, id=gene.id, description=gene.description)
            if pro2_rev[-1] == '*':
                SeqIO.write(Gene, Out, 'fasta')
            else:
                SeqIO.write(Gene, Out_NonStop, 'fasta')
        else:
            #print(gene.id +'\tNo correction')
            SeqIO.write(gene, Out_Incorrect, 'fasta')
    Out.close()
    Out_NonStop.close()
    Out_Incorrect.close()

def NonStandard_Rename_nonstop(input_fasta, output_correct_fasta, output_nonstop, output_incorrect):
    """Reads in a nonstandard fasta and outputs a corrected version"""
    genes = list(SeqIO.parse(input_fasta, 'fasta'))
    Out = open(output_correct_fasta, 'w')
    Out_Incorrect = open(output_incorrect, 'w')
    Out_NonStop = open(output_nonstop, 'w')
    for gene in genes:
        pro = str(gene.seq.translate())
        pro1 = str(gene.seq[1:].translate())
        pro2 = str(gene.seq[2:].translate())
        gene_rev = gene.seq.reverse_complement()
        gene1_rev = gene.seq[0:-1].reverse_complement()
        gene2_rev = gene.seq[0:-2].reverse_complement()
        pro_rev = str(gene_rev.translate())
        pro1_rev = str(gene1_rev.translate())
        pro2_rev = str(gene2_rev.translate())
        if ('*' in pro[0:-1]) == False:
            #print(gene.id + '\t' + 'Correct')
##            if  pro[-1] == '*':
            SeqIO.write(gene, Out, 'fasta')
            if pro[-1] != '*':
                SeqIO.write(gene, Out_NonStop, 'fasta')
        elif ('*' in pro1[0:-1]) == False:
            #print(gene.id + '\tframeshift 1')
            Gene = gene[1:]
##            if pro1[-1] == '*':
            SeqIO.write(Gene, Out, 'fasta')
            if pro1[-1] != '*':
                SeqIO.write(Gene, Out_NonStop, 'fasta')
        elif ('*' in pro2[0:-1]) == False:
            #print(gene.id + '\tframeshift 2')
            Gene = gene[2:]
##            if pro2[-1] == '*':
            SeqIO.write(Gene, Out, 'fasta')
            if pro2[-1] != '*':
                SeqIO.write(Gene, Out_NonStop, 'fasta')
        elif ('*' in pro_rev[0:-1]) == False:
            #print(gene.id + '\tReverse Complement')
            Gene = gene.seq.reverse_complement()
            Gene = SeqRecord(Gene, id=gene.id, description=gene.description)
##            if pro_rev[-1] == '*':
            SeqIO.write(Gene, Out, 'fasta')
            if pro_rev[-1] != '*':
                SeqIO.write(Gene, Out_NonStop, 'fasta')
        elif ('*' in pro1_rev[0:-1]) == False:
            #print(gene.id + '\tReverse Complement, frameshift 1')
            Gene = gene.seq[0:-1].reverse_complement()
            Gene = SeqRecord(Gene, id=gene.id, description=gene.description)
##            if pro1_rev[-1] == '*':
            SeqIO.write(Gene, Out, 'fasta')
            if pro_rev[-1] != '*':
                SeqIO.write(Gene, Out_NonStop, 'fasta')
        elif ('*' in pro2_rev[0:-1]) == False:
            #print(gene.id + '\tReverse Complement, frameshift 2')
            Gene = gene.seq[0:-2].reverse_complement()
            Gene = SeqRecord(Gene, id=gene.id, description=gene.description)
##            if pro2_rev[-1] == '*':
            SeqIO.write(Gene, Out, 'fasta')
            if pro2_rev[-1] != '*':
                SeqIO.write(Gene, Out_NonStop, 'fasta')
        else:
            print(gene.id +'\tNo correction')
            SeqIO.write(gene, Out_Incorrect, 'fasta')
    Out.close()
    Out_NonStop.close()
    Out_Incorrect.close()

def NonStop_Pro(input_fasta):
    """Identifies nonstop genes"""
    genes = list(SeqIO.parse(input_fasta, 'fasta'))
    for gene in genes:
        pro = str(gene.seq.translate())
        if pro[-1] != '*':
            Add = 0
            for gene2 in genes:
                pro2 = str(gene2.seq.translate())
                if gene.id == gene2.id:
                    continue
                elif (pro in pro2) == True:
                    print(gene.id + '\t' + str(len(pro)) + '\t' + gene2.id + '\t' + str(len(pro2))) 
                    Add = 1
            if Add == 0:
                print(gene.id +'\tNo match')
        
def NonStop(input_fasta):
    """Identifies nonstop genes"""
    genes = list(SeqIO.parse(input_fasta, 'fasta'))
    for gene in genes:
        pro = str(gene.seq.translate())
        if pro[-1] != '*':
            Add = 0
            for gene2 in genes:
                pro2 = str(gene2.seq.translate())
                if gene.id == gene2.id:
                    continue
                elif (str(gene.seq) in str(gene2.seq)) == True:
                    print(gene.id + '\t' + gene2.id)
                    Add = 1
##            if Add == 0:
##                print(gene.id +'\tNo match')      

def NonStop_Internal(correct_fasta, nonstop_fasta):
    """Determines if nonstop fasta are internal to other genes"""
    genes = list(SeqIO.parse(correct_fasta, 'fasta'))
    nonstops = list(SeqIO.parse(nonstop_fasta, 'fasta'))
    for ns_gene in nonstops:
        for gene in genes:
            if (str(ns_gene.seq) in str(gene.seq)) == True:
                print(ns_gene.id)

def NonStop_Internal_Write(correct_fasta, nonstop_fasta, output_fasta):
    """Determines if nonstop fasta are internal to other genes"""
    genes = list(SeqIO.parse(correct_fasta, 'fasta'))
    nonstops = list(SeqIO.parse(nonstop_fasta, 'fasta'))
    Out = open(output_fasta, 'w')
    for gene1 in genes:
        SeqIO.write(gene1, Out, 'fasta')
    for ns_gene in nonstops:
        Add = 1
        for gene in genes:
            if (str(ns_gene.seq) in str(gene.seq)) == True:
##                print(ns_gene.id)
                Add = 0
        if Add == 1:
            SeqIO.write(ns_gene, Out, 'fasta')
    Out.close()

def Internal_Remover(input_fasta, output_fasta):
    genes = list(SeqIO.parse(input_fasta, 'fasta'))
    Out = open(output_fasta, 'w')
    Sequences = []
    for gene in genes:
        Add = 1
        for entries in Sequences:
            if (str(gene.seq.upper()) in entries) == True:
                Add = 0
        if Add == 1:
            SeqIO.write(gene, Out, 'fasta')
            Sequences.append(str(gene.seq).upper())
    Out.close()
            
def Nonbase_Remover(input_fasta, output_fasta, output_incorrect):
    genes = list(SeqIO.parse(input_fasta, 'fasta'))
    Out = open(output_fasta, 'w')
    Out_Incorrect = open(output_incorrect, 'w')
    for gene in genes:
        Bad = 0
        for bases in str(gene.seq).upper():
            if bases != 'C' and bases != 'A' and bases != 'G' and bases != 'T':
                Bad = 1
        if Bad == 0:
            SeqIO.write(gene, Out, 'fasta')
        elif Bad == 1:
            SeqIO.write(gene, Out_Incorrect, 'fasta')
    Out.close()
    Out_Incorrect.close()

Nonbase_Remover(sys.argv[1], sys.argv[1][0:-6] + '_rb.fasta', sys.argv[1][0:-6] + '_nsb.fasta')
NonStandard_Rename_nonstop(sys.argv[1][0:-6] + '_rb.fasta', sys.argv[1][0:-6] + '_correct.fasta', sys.argv[1][0:-6] + '_nonstop.fasta', sys.argv[1][0:-6] + '_nonstandard.fasta')
Internal_Remover(sys.argv[1][0:-6] + '_correct.fasta', sys.argv[1][0:-6] + '_deduplicated.fasta')
