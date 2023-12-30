#!/usr/bin/env python3
# Creates a file with a polymer of size = AA size

import seqAnalyzer as san
import numpy as np
import Bio.SeqUtils as su

fasta_file = "./1BA9.fasta"

seq = san.extract_protein_seq(fasta_file)

natoms = seq.__len__()
mol_wt = su.molecular_weight(seq, 'protein')
atom_wt = mol_wt/natoms

A = (su.molecular_weight('A', 'protein'))
C = (su.molecular_weight('AAA', 'protein'))
MH2O = (3*A-C)/2
print(MH2O)
D = (su.molecular_weight('L', 'protein'))
E = (su.molecular_weight('LLA', 'protein'))

# print(natoms)
