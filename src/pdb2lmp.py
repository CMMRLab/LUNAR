# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
Decemeber 12th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

This script allows for a .pdb file to be read in
and an Molecule class 'm' to be created on the info
in the .pdb file, such that the class is exactly the
same that read_lmp creates, when reading in a LAMMPS
data file. This allows for a .pdb file to have the exact
same class type objects that exists when reading a LAMMPS
data file and allows the set of codes to read a .pdb file
as if it was a LAMMPS data file.

The atom cooridinates will be centered around the geometrical
center of the read in molecule and then the box will be centered
around the newly re-centered moleucle. The box will be 2.0 angstroms
larger in the x, y, and z directions and all image flags will be set
to zero. The atom charge will be set to zero and updated later on and 
the atom molid will be set as 1 and will stay as 1 throughout the
rest of the codes that use the information found by this code.
      
This script could serve as a starting point to create other
_ _ _2lmp scripts to have the entire set of scripts be able
to read in other molecule file types and be able to convert
a large number of file types to a lammps datafile.  
"""


###################################    
### Class for reading pdb files ###
###################################  
class Atom_mol: pass # .element .x .y .z .charge
class Bond_mol: pass  # .type .atomids = [atom1id, atom2id]
class read_pdb:    
    def __init__(self, inmolfile):
        self.natoms = 0 # total atoms in .pdb file
        self.nbonds = 0 # total bonds in .pdb file
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.graph = {}  # {atomID : set(of bondedIDs)}
        
        
        ############################################
        ### Opening and reading the vmd.mol file ###
        ############################################
        supported_record_types = ['ATOM', 'HETATM', 'CONECT', 'END', 'REMARK', '#',
                                  'MASTER', 'COMPND', 'AUTHOR', 'HEADER', 'TITLE']
        print(f'\nReading .pdb file, currently supported record types: {", ".join(supported_record_types)}')
        with open(inmolfile, 'r') as f:
            # Looping through each line of the reaxc file
            for line_number, line in enumerate(f, 1):

                # Try finding line record type in 1st position
                record_type = ''; nchars = 0;
                for record in supported_record_types:
                    nchars = len(record)
                    char_lst = [line[i:i+nchars] for i in range(0, len(line), nchars)]
                    if char_lst:
                        if char_lst[0] in supported_record_types:
                            record_type = char_lst[0]; break;
                if record_type == '':
                    print(f'record on line#: {line_number} not supported. line: {line.strip()}')

                
                # Start logging info based record_type (add onto this section
                # as more record types are added to supported_record_type lst)
                if record_type in ['ATOM', 'HETATM']:
                    try:
                        atom_serial_number = int(line[6:11].strip())
                        try: atom_name = int(line[12:16].strip())
                        except: atom_name = line[12:16].strip()
                        location_indicator = line[16].strip()
                        residue_name = line[17:20].strip()
                        chain_identifier = line[21].strip()
                        residue_sequence = line[22:26].strip()
                        residue_insertion_code = line[26].strip()
                        x = float(line[30:38].strip()) # limit of 8-digits
                        y = float(line[38:46].strip()) # limit of 8-digits
                        z = float(line[46:54].strip()) # limit of 8-digits
                        try: occupancy = float(line[54:60].strip()) # limit of 6-digits
                        except: occupancy = float(0)
                        try: temp_factor = float(line[60:66].strip()) # limit of 6-digits
                        except: temp_factor = float(0)
                        segment_identifier = line[72:76].strip() # limit of 4-chars
                        element = line[76:78].strip() # limit of 2-chars
                        try: charge = float(line[78:80].strip()) # limit of 2-digits
                        except: charge = float(0)
                        
                        # Add info to self.atoms
                        atom_id = atom_serial_number
                        a = Atom_mol()
                        a.atom_name = atom_name
                        a.type = atom_name
                        a.location_indicator = location_indicator
                        a.residue_name = residue_name
                        a.chain_identifier = chain_identifier
                        a.residue_sequence = residue_sequence
                        a.residue_insertion_code = residue_insertion_code
                        a.x = x
                        a.y = y
                        a.z = z
                        a.occupancy = occupancy
                        a.temp_factor = temp_factor
                        a.segment_identifier = segment_identifier
                        a.element = element
                        a.charge = charge
                        self.atoms[atom_id] = a
                        self.natoms += 1
                        self.graph[atom_serial_number] = set()
                    except: pass
                if record_type in ['CONECT']:
                    try:
                        atom_serial_number = int(line[6:11].strip())
                        maxIDs = 10 # maximum number of IDs to read (max 1st neighs)
                        start_index = 11; char_span = 5; # set start index and character spans
                        IDs = [0 for i in range(maxIDs)]
                        for ID_local in range(maxIDs):
                            lo_index = int(ID_local*char_span + start_index)
                            hi_index = int(ID_local*char_span + start_index + char_span)
                            try: ID = int(line[lo_index:hi_index].strip())
                            except: ID = int(0)
                            IDs[ID_local] = ID
                            
                        # Add atoms to graph
                        for ID in IDs:
                            if ID != 0:
                                if atom_serial_number in self.graph:
                                    self.graph[atom_serial_number].add(ID)
                                else:
                                    self.graph[atom_serial_number] = set(ID)
                                    print(f'WARNING atom serial number {atom_serial_number} was not inserted into graph during reading atoms')
                    except: pass
                            
        # Generate bonds from graph
        bonds = set() # { tuple(sorted(ID1, ID2)) }
        for ID1 in self.graph:
            for ID2 in self.graph[ID1]:
                bond = tuple(sorted([ID1, ID2]))
                bonds.add(bond)
        bonds = sorted(bonds)
                
        # Find insert bonds into a common data structure
        for bond_id, (ID1, ID2) in enumerate(bonds, 1):
            b = Bond_mol()
            b.atomids = [ID1, ID2]
            #b.type = line_split[3]
            self.bonds[bond_id] = b
        
        # Update natoms and nbonds
        self.natoms = len(self.atoms)
        self.nbonds = len(self.bonds)



###################################################
### Class for converting all info in .pdb file  ### 
### into a class that is exactly the same as if ###
### comes from read_lmpo to be able to use      ###
### use .pdb files in remaining of the codes    ###
###################################################
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz
class Bond: pass  # .type .atomids = [atom1id, atom2id]
import os
class Molecule_File:
    def __init__(self, inmolfile):

        # Read SYBYL pdb file
        molecule = read_pdb(inmolfile) 

        # Set filename
        self.filename = inmolfile
        
        # Set system information
        self.natoms = molecule.natoms # total atoms in .pdb file
        self.nbonds = molecule.nbonds # total bonds in .pdb file
        self.natomtypes = 0
        self.nbondtypes = 0
        self.atoms = {}  # {atom number : atom object}
        self.bonds = molecule.bonds  # {bond number : bond object}
        self.xbox_line = ''
        self.ybox_line = ''
        self.zbox_line = ''
        self.xy = 0;
        self.xz = 0;
        self.yz = 0;
        self.header = '{} {} {}'.format('HEADER, ', os.path.basename(inmolfile), ' read w/ pdb2lmp')
        self.velocities = {} # { atomid : tuple(velx, vely, velz)}
        self.natomtypes = 0
        self.bond_coeffs = {}
        
        #########################################################
        # Find center of atoms and shift needed to center atoms #
        #########################################################
        x = []; y = []; z = [];
        for i in molecule.atoms:
            x.append(molecule.atoms[i].x)
            y.append(molecule.atoms[i].y)
            z.append(molecule.atoms[i].z)
        
        # Finding geometric center
        x_center = sum(x)/len(x); y_center = sum(y)/len(y); z_center = sum(z)/len(z)
        
        # Finding needed shift to center around geometric center
        x_shift = (0 - x_center); y_shift = (0 - y_center); z_shift = (0 - z_center);

        
        #####################################
        # Shifting atoms and building atoms #
        # object exactly as read_lmp does   #
        #####################################
        atomtypes = set()
        for i in molecule.atoms:
            atom = molecule.atoms[i]
            # adding in new atom types and shifting atoms to be centered around geometric center
            a = Atom()
            a.type = atom.type
            a.molid = 1  # Set as 1 (option in main code to update later on)
            try: a.charge = molecule.atoms[i].charge # try getting charge
            except: a.charge = 0 # Set as 0 (option in main code to update later on)
            atomtypes.add(atom.type)
            a.atom_name = atom.atom_name
            a.location_indicator = atom.location_indicator
            a.residue_name = atom.residue_name
            a.chain_identifier = atom.chain_identifier
            a.residue_sequence = atom.residue_sequence
            a.residue_insertion_code = atom.residue_insertion_code
            a.x = atom.x + x_shift
            a.y = atom.y + y_shift
            a.z = atom.z + z_shift
            a.occupancy = atom.occupancy
            a.temp_factor = atom.temp_factor
            a.segment_identifier = atom.segment_identifier
            a.element = atom.element
            a.ix = 0 # box size will be found such that image flags from .pdb files will be zero
            a.iy = 0 # box size will be found such that image flags from .pdb files will be zero
            a.iz = 0 # box size will be found such that image flags from .pdb files will be zero
            self.atoms[i] = a
            self.velocities[i] = (0, 0, 0)
        self.natomtypes = len(atomtypes)
            
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
        oversize = 0.5 # default oversize of 0.5 angtroms in each value (total over size = 0.5*2 = 1.0 angstroms)
        zero_dim_default = 1.0 # default +- value of box dimension is zero
        xlo = min(x)-oversize; xhi = max(x)+oversize;
        ylo = min(y)-oversize; yhi = max(y)+oversize;
        zlo = min(z)-oversize; zhi = max(z)+oversize;
        
        # if xlo and xhi == 0 reset to +- zero_dim_default value
        if xlo == 0 and xhi == 0 or xlo == -oversize and xhi == oversize:
            xlo = -zero_dim_default; xhi = zero_dim_default;
            
        # if ylo and yhi == 0 reset to +- zero_dim_default value
        if ylo == 0 and yhi == 0 or ylo == -oversize and yhi == oversize:
            ylo = -zero_dim_default; yhi = zero_dim_default;
            
        # if zlo and zhi == 0 reset to +- zero_dim_default value
        if zlo == 0 and zhi == 0 or zlo == -oversize and zhi == oversize:
            zlo = -zero_dim_default; zhi = zero_dim_default;
        
        # Set box dimensions string
        self.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(xlo, xhi, 'xlo', 'xhi')
        self.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(ylo, yhi, 'ylo', 'yhi')
        self.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(zlo, zhi, 'zlo', 'zhi')
        

######################
# Test Read Function #
######################
if __name__ =="__main__":
    pdb_file = '../EXAMPLES/pdb_file_support/detda.pdb' 
    m = Molecule_File(pdb_file)
    
    # Check atoms
    for i in m.atoms:
        atom = m.atoms[i]
        print(i, atom.x, atom.y, atom.z, atom.element, atom.charge)