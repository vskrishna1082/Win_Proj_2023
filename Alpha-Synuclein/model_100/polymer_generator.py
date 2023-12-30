#! /usr/bin/env python3
# Writes a LAMMPS data file for a single polymer

import mylmp as ml
import seqAnalyzer as san
import numpy as np

TEMP = 300
SIGMA = 0.19
EPS = 0.25
BOND_K = 1.2111

fasta_file = "../1XQ8.fasta"
OUTFILE = "polymer.data"
seq = san.extract_protein_seq(fasta_file)

bond_data = np.loadtxt("../bond_mean_stdev.dat")
angle_data = np.loadtxt("../angle_mean_stdev.dat")
dihedral_data = np.loadtxt("../dihedral_mean_stdev.dat")
bond_data = bond_data*np.array([10,100]) # convert to Angstorm
angle_data = angle_data*np.array([1,10]) # convert to Angstorm

sim_box = ml.SimulationBox(750.0,750.0,750.0)

INIT_POS = np.array(sim_box.array/4)
AVG_BL = np.average(bond_data[:,0])
KB = 0.001985875 # boltzmann constant in kCal/mol.T
KBT = KB*TEMP
def stdev_to_energy(stdev):
    return KBT/(2*stdev*stdev)

print(stdev_to_energy(bond_data[0][1]))

atoms = [] # list for storing atoms - 1 res = 1 atom
bond_types = []
bonds = []
angle_types = []
angles = []
dihedral_types = []
dihedrals = []

for residue in seq:
    cg_atom = ml.Atom(len(atoms)+1)
    cg_atom.mass = san.get_res_wt(residue)
    cg_atom.charge = san.get_res_charge(residue)
    cg_atom.atom_type = cg_atom.atom_id
    # add the avg bond length of atom for new pos
    cg_atom.pos = INIT_POS + (len(atoms))*AVG_BL/(np.sqrt(3))
    cg_atom.vel = (np.random.rand(3)-0.5)*(np.sqrt(12*KBT/cg_atom.mass))
    atoms.append(cg_atom)

while len(bonds) < len(atoms)-1:
    id = len(bond_types) + 1
    myBondType = ml.BondType(id,
                             BOND_K,
                             2*SIGMA)
    myBond = ml.Bond(id, id, id+1)
    bond_types.append(myBondType)
    bonds.append(myBond)

for angle in angle_data[:]:
    id = len(angle_types) + 1
    myAngleType = ml.AngleType(id,
                             stdev_to_energy(angle[1]),
                             angle[0])
    myAngle = ml.Angle(id, id, id+1, id+2)
    angle_types.append(myAngleType)
    angles.append(myAngle)

for dihedral in dihedral_data[:]:
    id = len(dihedral_types) + 1
    myDihedralType = ml.AngleType(id,
                             stdev_to_energy(dihedral[1]),
                             dihedral[0])
    myDihedral = ml.Dihedral(id, id, id+1, id+2, id+3)
    dihedral_types.append(myDihedralType)
    dihedrals.append(myDihedral)

with open(OUTFILE, 'w') as fdata:
    fdata.write("LAMMPS Data File for polymers\n\n")
    # Header
    # detail numbers
    fdata.write(f"{len(atoms)} atoms\n")
    fdata.write(f"{len(bonds)} bonds\n")
    fdata.write(f"{len(angles)} angles\n")
    fdata.write(f"{len(dihedrals)} dihedrals\n")
    fdata.write(f"{len(atoms)} atom types\n")
    fdata.write(f"{len(bond_types)} bond types\n\n")
    fdata.write(f"{len(angle_types)} angle types\n\n")
    fdata.write(f"{len(dihedral_types)} dihedral types\n\n")

    # box dimensions
    fdata.write(sim_box.lammps_string())
    fdata.write('\n')

    # masses
    fdata.write('Masses\n\n')
    for idx, atom in enumerate(atoms):
        fdata.write(f"{idx+1} {atom.mass}\n")
    fdata.write('\n')

    # bond coeffs
    fdata.write('Bond Coeffs\n\n')
    for bond_type in bond_types:
        fdata.write(bond_type.lammps_string())
    fdata.write('\n')

    # angle coeffs
    fdata.write('Angle Coeffs\n\n')
    for angle_type in angle_types:
        fdata.write(angle_type.lammps_string())
    fdata.write('\n')

    # dihedral coeffs
    fdata.write('Dihedral Coeffs\n\n')
    for dihedral_type in dihedral_types:
        fdata.write(dihedral_type.lammps_string())
    fdata.write('\n')

    # atoms section
    fdata.write('Atoms\n\n')
    for atom in atoms:
        fdata.write(atom.lammps_string())
    fdata.write('\n')

    # velocity section
    fdata.write('Velocities\n\n')
    for atom in atoms:
        fdata.write('{} {} {} {}\n'.format(atom.atom_id, *atom.vel))
    fdata.write('\n')

    # bonds section
    fdata.write('Bonds\n\n')
    for bond in bonds:
        fdata.write(bond.lammps_string())

    # angles section
    fdata.write('Angles\n\n')
    for angle in angles:
        fdata.write(angle.lammps_string())

    # dihedrals section
    fdata.write('Dihedrals\n\n')
    for dihedral in dihedrals:
        fdata.write(dihedral.lammps_string())
