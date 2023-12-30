#! /usr/bin/env python3
import seqAnalyzer as san
import matplotlib.pyplot as plt
import numpy as np

fasta_file = "./4BCY.fasta"

seq = san.extract_protein_seq(fasta_file)

charge_list = np.array([san.res_to_charge(X) for X in seq])
plt.imshow(charge_list[np.newaxis,:],
           aspect=10, cmap='RdBu',
           vmax=1.5, vmin=-1.5)
plt.colorbar(orientation='horizontal')
plt.show()
