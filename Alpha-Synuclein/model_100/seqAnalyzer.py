#! /usr/bin/env python3
# Prints the residue sequence of the protein

from Bio import SeqIO 
import Bio.SeqUtils as su

MWATER = 18.01528
def get_res_wt(residue):
    return su.molecular_weight(residue, 'protein') - 2*MWATER

def get_res_charge(res):
    # converts letter rep. of AA residue to charge
    if res == 'R' or res == 'K':
        return 1
    elif res == 'D' or res == 'E':
        return -1
    else:
        return 0

def extract_protein_seq(fasta_file):
    fasta_seq = SeqIO.parse(open(fasta_file), 'fasta')
    for fasta in fasta_seq:
        sequence = fasta.seq
    return sequence
