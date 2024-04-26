# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 29th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    *********************************************************
    * Requirements:                                         *
    *   python 3.7+                                         *
    *                                                       *   
    * Dependencies:                                         *
    *   python rdkit module:                                *
    *    - pip3 install rdkit (if pip manager is installed) *
    *********************************************************
"""


###############################################################
### Class to manipulate rdkit opperations to interface with ### 
### atom_typing code data structs/classes. This is still    ###
### under construction, but supports basic interfacing with ###
### rdkit and can be added onto as needed                   ###
###############################################################
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

# Function to generate gasteiger from rdkit
def compute_gasteiger_charges(mol):
    n = mol.GetNumAtoms(); charge_dict = {} # {atomid: partial charge}
    maxIter = 100*n # Set maximum interations based on system size (Over Shoot by a large amount)
    AllChem.ComputeGasteigerCharges(mol, nIter=maxIter, throwOnParamFailure=True)
    for i in range(n):
        a = mol.GetAtomWithIdx(i)
        charge = float(a.GetProp('_GasteigerCharge'))
        #hcharge = float(a.GetProp('_GasteigerHCharge'))
        charge_dict[i+1] = charge
    return charge_dict

# Class to help interface rdkit with rest of other codes
class rdkit2lmp:
    def __init__(self, name, addHs=True, gen_mclass=False):
        self.charges = {} # {atomid: partial charge} Initialize as empty but update later on
        self.m = None # Initialize as None but update later on if desired
        
        ################################################
        # Generate rdkit mol class from name extension #
        ################################################
        # Use .mol file to generate mol class from rdkit
        if name.endswith('mol'):
            mol = Chem.MolFromMolFile(name)
        
        # Use .mol2 file to generate mol class from rdkit
        elif name.endswith('mol2'):
            mol = Chem.MolFromMol2File(name)
            
        # Use .smiles file to generate mol class from rdkit
        elif name.endswith('smiles'):
            smiles = name[:name.rfind('.')] # strip .smiles ending
            mol = Chem.MolFromSmiles(smiles)

        
        # Else Warn user and exit
        else:
            print(f'ERROR Input file extension: {name} not currently supported'); sys.exit();
            
        #############################################
        # If addHs add hydrogens, Default is to add #
        #############################################
        if addHs:
            mol = Chem.AddHs(mol)
            
        # Make molecule 3D after Hs have been added
        # https://stackoverflow.com/questions/70335013/why-i-cant-get-3d-mol-structure-in-rdkit
        AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
            
        ######################################################################
        # Find m class like other read methods from rdkit, Default is not to #
        #####################################################################
        if gen_mclass:
            # Get information from mol_block
            mol_block = Chem.MolToMolBlock(mol)
            string = ''; lines = [];
            for i in mol_block:
                string += i # sring together lines
                # If next line character is encountered append
                # string to lines and reset string as new string
                if '\n' in i:
                    lines.append(string); string = '';
                    
            # Generate m class
            self.m = mol_block2lmp(lines)
                    
        ####################################################################
        # Compute Gasteiger charges and update self.charge_dict dictionary #
        ####################################################################
        self.charges = compute_gasteiger_charges(mol)

            
#########################################################    
### Class for converting lines in MolBlock to m class ###
#########################################################
class Atom_mol:
    pass # .element .x .y .z

class Bond_mol:
    pass  # .type .atomids = [atom1id, atom2id]
            
# for reading .mol file
class mol_block2lmp:    
    def __init__(self, lines):        
        # Set system information
        self.filename = 'rdkit.txt'
        self.natoms = 0 # total atoms in .mol file
        self.nbonds = 0 # total bonds in .mol file
        self.natomtypes = 0
        self.nbondtypes = 0
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.xbox_line = ''
        self.ybox_line = ''
        self.zbox_line = ''
        self.xy = 0;
        self.xz = 0;
        self.yz = 0;
        self.header = '{}'.format('HEADER, mol_block2lmp from rdkit')
        self.velocities = {} # { atomid : tuple(velx, vely, velz)}
        
        ####################################
        # Testing if integer function for  #
        # finding bonds line in .mol file. #
        ####################################
        def is_integer(n):
            try:
                float(n)
            except ValueError:
                return False
            else:
                return float(n).is_integer()
        
        ############################################
        ### Looping through lines from mol_block ###
        ############################################
        # Initializing flags
        info_flag = False; atoms_flag = False;
        bonds_flag = False; atom_id = 0; bond_id = 0;
        
            
        # Looping through each line of the lines in file
        for n, line in enumerate(lines):
            
            # Split every 3 chars to get natoms
            char_lst = [line[i:i+3] for i in range(0, len(line), 3)]
            
            # Setting flags
            if char_lst == ['\n']:
                info_flag = False; atoms_flag = False;
                bonds_flag = False; end_flag = False; 
            elif 'V2000' in line:
                info_flag = True
                natoms = int(char_lst[0].strip())
            elif len(char_lst) > 11 and n >= 4:
                info_flag = False
                atoms_flag = True
            elif line != '' and is_integer(char_lst[0].strip()) and n > 3+natoms:
                atoms_flag = False
                bonds_flag = True
            elif 'END' in line or 'M' in line:
                atoms_flag = False
                bonds_flag = False
                end_flag = True

            
            # Finding atom and bond info
            if info_flag:
                # Find char_lst be splitting every 3-characters and then split whitespace from needed info
                char_lst = [line[i:i+3] for i in range(0, len(line), 3)] # Split every 3 chars to get natoms and nbonds
                self.natoms = int(char_lst[0].strip())
                self.nbonds = int(char_lst[1].strip())
            elif atoms_flag:
                atom_id += 1
                a = Atom_mol()
                char_coords = [line[i:i+10] for i in range(0, len(line), 10)] # Split every 10 chars to get coords
                char_element = line[30:33] # Split line at 30-33 chars to get element
                a.x = float(char_coords[0].strip())
                a.y = float(char_coords[1].strip())
                a.z = float(char_coords[2].strip())
                a.element = str(char_element.strip())
                a.charge = 0
                a.molid = 1
                self.atoms[atom_id] = a
            elif bonds_flag:
                char_bonds = [line[i:i+3] for i in range(0, len(line), 3)] # Split every 3 chars to get bonds
                bond_id += 1; b = Bond_mol();
                b.atomids = [int(char_bonds[0].strip()), int(char_bonds[1].strip())]
                b.type = int(char_bonds[2].strip())
                self.bonds[bond_id] = b
            elif end_flag:
                pass

            
        ###############################################
        # Find new atoms x, y, z position after shift #  
        # to build lammps box size to produce image   #
        # flags in x, y, and z direction as zeros     #
        ###############################################
        x = []; y = []; z = [];
        for i in self.atoms:
            x.append(self.atoms[i].x)
            y.append(self.atoms[i].y)
            z.append(self.atoms[i].z)
            
        
        ###################################################
        # Find x, y, z box dims (search for min/max and   #
        # then oversize slightly in each direction also   #
        # if certain dimensions are zero set default dim) #
        ###################################################
        oversize = 0.25 # default oversize of 0.25 angtroms in each value (total over size = 0.25*2 = 0.5 angstroms)
        zero_dim_default = 0.5 # default +- value of box dimension is zero
        xlo = min(x)-oversize; xhi = max(x)+oversize;
        ylo = min(y)-oversize; yhi = max(y)+oversize;
        zlo = min(z)-oversize; zhi = max(z)+oversize;
        
        # if xlo and xhi == 0 reset to +- zero_dim_default value
        if xlo == 0 and xhi == 0 or xlo == -oversize and xhi == oversize:
            xlo = -zero_dim_default
            xhi = zero_dim_default
            
        # if ylo and yhi == 0 reset to +- zero_dim_default value
        if ylo == 0 and yhi == 0 or ylo == -oversize and yhi == oversize:
            ylo = -zero_dim_default
            yhi = zero_dim_default
            
        # if zlo and zhi == 0 reset to +- zero_dim_default value
        if zlo == 0 and zhi == 0 or zlo == -oversize and zhi == oversize:
            zlo = -zero_dim_default
            zhi = zero_dim_default
        
        # Set box dimensions string
        self.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(xlo, xhi, 'xlo', 'xhi')
        self.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(ylo, yhi, 'ylo', 'yhi')
        self.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(zlo, zhi, 'zlo', 'zhi')