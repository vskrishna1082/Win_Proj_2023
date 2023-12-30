#! /usr/bin/env python3
""" Module containing objects for LAMMPS data file """
import numpy as np

class SimulationBox:
    """ A class for the simulationBox parameters """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.array = np.array([x,y,z])
    def lammps_string(self):
        """ writes out a lammps string for the sim box """
        return (f"0.0 {self.x} xlo xhi\n"
                f"0.0 {self.y} ylo yhi\n"
                f"0.0 {self.z} zlo zhi\n")

class Atom:
    """ A class for atom properties
    @params pos,vel,mass,mol_id,id,charge,atom_type
    @func set_pos(),set_vel(),write_stype()
    """
    def __init__(self,
                 atom_id,
                 pos = np.zeros(3),
                 vel = np.zeros(3),
                 mass = 1,
                 mol_id = 1,
                 charge = 0,
                 atom_type = 1):
        self.atom_id = atom_id
        self.pos = pos
        self.vel = vel
        self.mass = mass
        self.mol_id = mol_id
        self.charge = charge
        self.atom_type = atom_type
    def set_pos(self, x,y,z):
        """ sets position of atom """
        self.pos = [x,y,z]
    def set_vel(self, vx,vy,vz):
        """ sets velocity of atom """
        self.vel= [vx,vy,vz]
    def lammps_string(self,style='full'):
        """ writes out an Atoms section line for LAMMPS"""
        match style:
            case 'full':
                return "{} {} {} {} {} {} {}\n".format(self.atom_id,
                                                       self.mol_id,
                                                       self.atom_type,
                                                       self.charge,
                                                       *self.pos)

class Bond:
    """ A class for bonds - bond has id, type
    and atoms it connects.
    (but not properties - that is in BondType)
    """
    def __init__(self,
                 bond_id,
                 atom1_id,
                 atom2_id,
                 bond_type=1):
        self.bond_id = bond_id
        self.atom1_id = atom1_id
        self.atom2_id = atom2_id
        self.bond_type = bond_type
    def connecting_atoms(self):
        """ get the atoms making up the bond """
        return [self.atom1_id, self.atom2_id]
    def lammps_string(self):
        """ writes out a Bonds section line for LAMMPS"""
        return "{} {} {} {}\n".format(self.bond_id,
                                      self.bond_type,
                                      self.atom1_id,
                                      self.atom2_id)

class Angle:
    """ A class for angles - an angle has id, type
    and atoms it connects - Angle L123
    (but not properties - that is in AngleType)
    """
    def __init__(self,
                 angle_id,
                 atom1_id,
                 atom2_id,
                 atom3_id,
                 angle_type=1):
        self.angle_id = angle_id
        self.atom1_id = atom1_id
        self.atom2_id = atom2_id
        self.atom3_id = atom3_id
        self.angle_type = angle_type
    def connecting_atoms(self):
        """ get the atoms making up the angle """
        return [self.atom1_id, self.atom2_id, self.atom3_id]
    def lammps_string(self):
        """ writes out a Angles section line for LAMMPS"""
        return "{} {} {} {} {}\n".format(self.angle_id,
                                      self.angle_type,
                                      self.atom1_id,
                                      self.atom2_id,
                                      self.atom3_id)

class Dihedral:
    """ A class for dihesra;s - a dihedral has id, type
    and atoms it connects - 2-3 is the central bond
    """
    def __init__(self,
                 angle_id,
                 atom1_id,
                 atom2_id,
                 atom3_id,
                 atom4_id,
                 dihedral_type=1):
        self.angle_id = angle_id
        self.atom1_id = atom1_id
        self.atom2_id = atom2_id
        self.atom3_id = atom3_id
        self.atom4_id = atom4_id
        self.dihedral_type = dihedral_type
    def connecting_atoms(self):
        """ get the atoms making up the angle """
        return [self.atom1_id, self.atom2_id, self.atom3_id, self.atom4_id]
    def lammps_string(self):
        """ writes out a Dihedrals section line for LAMMPS"""
        return "{} {} {} {} {} {}\n".format(self.angle_id,
                                      self.dihedral_type,
                                      self.atom1_id,
                                      self.atom2_id,
                                      self.atom3_id,
                                      self.atom4_id)

class BondType:
    """ Holds a bond type's parameters
    Currently, harmonic bonds supported
    >>> myBondType = BondType(1,10,100)
    >>> myBondType.lammps_string()
    '1 10 100\\n'
    """
    def __init__(self, bondtype_id, k, l_0):
        self.bondtype_id = bondtype_id
        self.k = k
        self.l_0 = l_0
    def lammps_string(self):
        """ writes out a Bond Types section line for LAMMPS"""
        return f"{self.bondtype_id} {self.k} {self.l_0}\n"

class AngleType:
    """ Holds a angle type's parameters
    >>> myAngleType = AngleType(1,10,100)
    >>> myAngleType.lammps_string()
    '1 10 100\\n'
    """
    def __init__(self, angletype_id, k, theta_0):
        self.angletype_id = angletype_id
        self.k = k
        self.theta_0 = theta_0 # theta in degrees
    def lammps_string(self):
        """ writes out a Angle Types section line for LAMMPS"""
        return f"{self.angletype_id} {self.k} {self.theta_0}\n"

class DihedralType:
    """ Holds a dihedral angle type's parameters
    """
    def __init__(self, dihedraltype_id, k, theta_0):
        self.dihedraltype_id = dihedraltype_id
        self.k = k
        self.theta_0 = theta_0 # theta in degrees
    def lammps_string(self):
        """ writes out a Angle Types section line for LAMMPS"""
        return f"{self.dihedraltype_id} {self.k} {self.theta_0}\n"
