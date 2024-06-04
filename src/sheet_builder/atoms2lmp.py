# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
June 4th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.masses as masses
mass_map = masses.mass_map 


##########################
# Function to guess mass #
##########################
def guess_mass(atomtype, mass_map, masses):
    # Setup reduced_mass_map to only get elements of interest
    elements = ['C', 'H', 'N', 'B']
    reduced_mass_map = {i:mass_map[i] for i in mass_map if i in elements}
    
    # Set defaults
    mass = 0; element = '';
    
    # Try getting from "masses" first
    if atomtype in masses:
        mass, element = masses[atomtype]
    else: # else attempt setting mass by checking first characters of atomtype
        for element in reduced_mass_map:
            capital = element[0].capitalize()
            lower = element[0].lower()
            if atomtype.startswith(capital) or atomtype.startswith(lower):
                mass = mass_map[element][0]    
    return mass, element

# Class to create object like read_lmp
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz .comment
class Bond: pass  # .type .atomids = [atom1id, atom2id]
class Angle: pass  # .type .atomids = [atom1id, atom2id, atom3id]
class Dihedral: pass  # .type .atomids = [atom1id, atom2id, atom3id, atom4id]
class Improper: pass  # .type .atomids = [atom1,atom2,atom3,atom4]
class Coeff_class: pass  # .type .coeffs = []
class Molecule_File:
    def __init__(self, basename, atoms, bonds, box, masses, header):
        self.filename = basename
        
        # Structure info
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.angles = {}  # {angle number : angle object}
        self.dihedrals = {}  # {dihedral number : dihedral object}
        self.impropers = {}  # {improper number : improper object}
        self.velocities = {}  # {atom number : tuple of velocities}

        # Parameters
        self.masses = {}  # {atom type : list of coeffs}
        self.pair_coeffs = {}  # {atom type : list of coeffs}
        self.bond_coeffs = {}   # {bond type : list of coeffs}  {1: [340,1.5], 2: [450,1.2], ...}
        self.angle_coeffs = {}  # {angle type : list of coeffs}
        self.dihedral_coeffs = {}  # {dihedral type : list of coeffs}
        self.improper_coeffs = {}  # {improper type : list of coeffs}
        self.bondbond_coeffs = {} # {cross-term number : list of coeffs}, angles
        self.bondangle_coeffs = {} # {cross-term number : list of coeffs}, angles
        self.angleangle_coeffs = {} # {cross-term number : list of coeffs}, impropers
        self.angleangletorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.endbondtorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.middlebondtorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.bondbond13_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.angletorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        
        # number of quantities
        self.natoms = 0
        self.natomtypes = 0
        self.nbonds = 0
        self.nbondtypes = 0
        self.nangles = 0
        self.nangletypes = 0
        self.ndihedrals = 0
        self.ndihedraltypes = 0
        self.nimpropers = 0
        self.nimpropertypes = 0
        self.nbondbond = 0
        self.nbondangle = 0
        self.nangleangle = 0
        self.nangleangletorsion = 0
        self.nendbondtorsion = 0
        self.nmiddlebondtorsion = 0
        self.nbondbond13 = 0
        self.nangletorsion = 0
        
        # header, box, and misc
        self.triclinicbox_line = ''
        self.xlo = box['xlo']
        self.xhi = box['xhi']
        self.ylo = box['ylo']
        self.yhi = box['yhi']
        self.zlo = box['zlo']
        self.zhi = box['zhi']
        self.lx = box['xhi'] - box['xlo']
        self.ly = box['yhi'] - box['ylo']
        self.lz = box['zhi'] - box['zlo']
        self.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(self.xlo, self.xhi, 'xlo', 'xhi')
        self.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(self.ylo, self.yhi, 'ylo', 'yhi')
        self.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(self.zlo, self.zhi, 'zlo', 'zhi')
        self.xy = 0
        self.xz = 0
        self.yz = 0
        self.header = header
        
        # Generate masses
        unique_types = sorted({atoms[ID].atomtype for ID in atoms})
        self.natomtypes = len(unique_types)
        types_forward = {} # { atomtype : atomID }
        types_reverse = {} # { atomID : atomtype  }
        for ID, atomtype in enumerate(unique_types, 1):
            types_forward[atomtype] = ID
            types_reverse[ID] = atomtype
            mass, element = guess_mass(atomtype, mass_map, masses)
            c = Coeff_class()
            c.coeffs = [mass]
            c.type = atomtype
            c.element = element
            self.masses[ID] = c
            
        # Generate bonds and bond coeffs
        if bonds:
            self.bond_coeffs_style_hint = 'N/A'
            self.nbonds = len(bonds)
            self.nbondtypes = 1
            c = Coeff_class()
            c.coeffs = [0, 0]
            c.type = 'N/A'
            self.bond_coeffs[1] = c
            for ID, bond in enumerate(bonds, 1):
                B = Bond()
                B.type = 1
                B.atomids = bond
                self.bonds[ID] = B
                
        # file in atoms
        self.natoms = len(atoms)
        for ID in atoms:
            atom =  atoms[ID]
            a = Atom()
            a.molid = atom.molid
            a.type = types_forward[atom.atomtype]
            a.charge = 0.0
            a.x = atom.x
            a.y = atom.y
            a.z = atom.z
            a.ix = atom.ix
            a.iy = atom.iy
            a.iz = atom.iz
            a.comment = atom.atomtype
            self.atoms[ID] = a