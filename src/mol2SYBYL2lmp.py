# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
October 16th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

This script allows for a VMD MDL .mol2 file to be read in
and an Molecule class 'm' to be created on the info
in the .mol2 file, such that the class is exactly the
same that read_lmp creates, when reading in a LAMMPS
data file. This allows for a .mol2 file to have the exact
same class type objects that exists when reading a LAMMPS
data file and allows the set of codes to read a .mol2 file
as if it was a LAMMPS data file.

The atom cooridinates will be centered around the geometrical
center of the read in molecule and then the box will be centered
around the newly re-centered moleucle. The box will be 1.0 angstroms
larger in the x, y, and z directions and all image flags will be set
to zero. The atom charge will be set to zero and updated later on and 
the atom molid will be set as 1 and will stay as 1 throughout the
rest of the codes that use the information found by this code.

This function has been testing on .mol file created by chemdraw
and .mol2 files that have been download from the NIST Chemistry
webook site. If you have a .mol file that causes an error with
this script, two option exist:
    - address the issue and re-work the read_mol class
    - open .mol file in chemdraw and save from chemdraw to create
      new file that this script should be able to read in.
      
This script could serve as a starting point to create other
_ _ _2lmp scripts to have the entire set of scripts be able
to read in other molecule file types and be able to convert
a large number of file types to a lammps datafile.  
"""


########################################    
### Class for reading chemdraw files ###
########################################  
class Atom_mol: pass # .element .x .y .z
class Bond_mol: pass  # .type .atomids = [atom1id, atom2id]
# for reading .mol file
class read_mol2:    
    def __init__(self, inmolfile):
        self.natoms = 0 # total atoms in .mol file
        self.nbonds = 0 # total bonds in .mol file
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        
        
        ############################################
        ### Opening and reading the vmd.mol file ###
        ############################################
        with open(inmolfile, 'r') as f:
            
            # Initializing flags
            atoms_flag = False; bonds_flag = False;
            
            # Looping through each line of the reaxc file
            for n, line in enumerate(f):
                
                # strip and split line
                str_line = line
                line = line.strip();  line_split = line.split();
                
                # Setting flags
                if line == '': # skip empty lines
                    atoms_flag = False; bonds_flag = False;
                elif 'ATOM' in line and line[0] == '@':
                    atoms_flag = True
                    continue
                elif 'BOND' in line and line[0] == '@':
                    atoms_flag = False
                    bonds_flag = True
                    continue
                elif line[0] == '@': # split sections when encountered
                    atoms_flag = False; bonds_flag = False;
                    
                # Find atoms
                # https://chemicbook.com/2021/02/20/mol2-file-format-explained-for-beginners-part-2.html
                if atoms_flag and len(line_split) >= 5:
                    try:
                        atom_id = int(line_split[0])
                        try: charge = float(line_split[-1])
                        except: charge = 0.0
                        
                        # Get element symbol
                        element = ''
                        for char in line_split[1]:
                            if char.isalpha(): element += char
                            else: break
                            
                        # Save atom info
                        a = Atom_mol()
                        a.element = element
                        a.x = float(line_split[2])
                        a.y = float(line_split[3])
                        a.z = float(line_split[4])
                        a.charge = charge
                        self.atoms[atom_id] = a
                    except: print(f'WARNING can not read line number {n} -> ', str_line)
                    
                # Find bonds
                elif bonds_flag and len(line_split) >= 4:
                    try:
                        bond_id = int(line_split[0])
                        b = Bond_mol()
                        b.atomids = [int(line_split[1]), int(line_split[2])]
                        #b.type = line_split[3]
                        self.bonds[bond_id] = b
                    except: print(f'WARNING can not read line number {n} -> ', str_line)
        
        # Update natoms and nbonds
        self.natoms = len(self.atoms)
        self.nbonds = len(self.bonds)


#################################################
# function to find image flag and atom position #
#################################################
def wrap_atom_compute_iflag(pos, domain, li):
    # pos = position
    # domain = ilo or ihi; were i = x, y, or z
    # li = li; were i = x, y, or z
    
    # Intialize image as zero and pos as pos update based on domain
    image = 0; pos = pos; max_images = 3 # set max images to break while loop
    
    # If domain is positve
    if domain > 0:
        while pos > domain:
            pos = pos - li; image += 1;
            if abs(image) == max_images: break
    # If domain is negative
    elif domain < 0:
        while pos < domain:
            pos = pos + li; image -= 1;  
            if abs(image) == max_images: break
    # if domain is zero or 2D
    else:
        pos = pos; image = 0;
    return pos, image


###################################################
### Class for converting all info in VMD MDL    ### 
### .mol2 into a class that is exactly the same ###
### as it comes from read_lmp to be able to     ###
### use .mol2 files in remaining of the codes   ###
###################################################
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz
class Bond: pass  # .type .atomids = [atom1id, atom2id]
import os
class Molecule_File:
    def __init__(self, inmolfile):

        # Read SYBYL mol2 file
        molecule = read_mol2(inmolfile) 

        # Set filename
        self.filename = inmolfile
        
        # Set system information
        self.natoms = molecule.natoms # total atoms in .mol file
        self.nbonds = molecule.nbonds # total bonds in .mol file
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
        self.header = '{} {} {}'.format('HEADER, ', os.path.basename(inmolfile), ' read w/ mol2SYBYL2lmp')
        self.velocities = {} # { atomid : tuple(velx, vely, velz)}
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
        for i in molecule.atoms:
            # adding in new atom types and shifting atoms to be centered around geometric center
            a = Atom()
            a.type = 1 # Set as 1 this information is not used by all2lmp
            a.molid = 1  # Set as 1 (option in main code to update later on)
            try: a.charge = molecule.atoms[i].charge # try getting charge
            except: a.charge = 0 # Set as 0 (option in main code to update later on)
            a.element = molecule.atoms[i].element
            a.x = molecule.atoms[i].x + x_shift
            a.y = molecule.atoms[i].y + y_shift
            a.z = molecule.atoms[i].z + z_shift
            a.ix = 0 # box size will be found such that image flags from .mol2 files will be zero
            a.iy = 0 # box size will be found such that image flags from .mol2 files will be zero
            a.iz = 0 # box size will be found such that image flags from .mol2 files will be zero
            self.atoms[i] = a
            self.velocities[i] = (0, 0, 0)
            
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
        zero_dim_default = 0.5 # default +- value of box dimension is zero
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
        
        
        #########################################################
        # Check for filename of:  ${filename}_VMD_graphite.mol2 #
        # to perform proper method for setting simulaiton cell  #
        # and re-wraping of atoms between certain sheets.       #
        #########################################################
        if inmolfile.endswith('VMD_graphite.mol2'):
            print('\n\nRead-in .mol2 file ends with "VMD_graphite.mol2". Setting proper simulation cell size as graphene sheets are staggered')
            # Generate graph
            graph = {i:[] for i in self.atoms}
            for i in self.bonds:
                id1, id2 = self.bonds[i].atomids
                graph[id1].append(id2)
                graph[id2].append(id1)
            
            # Find clusters
            checked = {ID:False for ID in self.atoms}; clusters = set([]);
            for ID in graph:
                if checked[ID]: continue
                visited=set([ID]); queue=[ID];
                while queue:
                    s = queue.pop(0) 
                    for neighbor in graph[s]:
                        if checked[neighbor]: continue
                        visited.add(neighbor)
                        queue.append(neighbor)
                        checked[neighbor]=True
                clusters.add( tuple(sorted(visited)) )
            clusters = sorted(clusters, key=lambda x: x[0]) # Sort all clusters based on 1st atomID in cluster
            clusters = sorted(clusters, key=len, reverse=True) # Sort all clusters by number of atoms
            print('  {}  graphene sheets found'.format(len(clusters)))
            
            # Only procede if graphene clusters are greater then 1
            if len(clusters) > 1:
                print('  Reseting box size based on first graphene sheet.')
                x = []; y = []; d002 = 3.354
                for i in clusters[0]:
                    x.append(self.atoms[i].x)
                    y.append(self.atoms[i].y)
                xlo = min(x)-oversize; xhi = max(x)+oversize;
                ylo = min(y)-oversize; yhi = max(y)+oversize;
                zlo = min(z)-d002/2; zhi = max(z)+d002/2;
                print('  Originally constructed simulation cell:')
                print('    {}'.format(self.xbox_line))
                print('    {}'.format(self.ybox_line))
                print('    {}'.format(self.zbox_line))
                
                # Set box dimensions string
                self.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(xlo, xhi, 'xlo', 'xhi')
                self.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(ylo, yhi, 'ylo', 'yhi')
                self.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(zlo, zhi, 'zlo', 'zhi')
                
                print('  Redefined simulation cell:')
                print('    {}'.format(self.xbox_line))
                print('    {}'.format(self.ybox_line))
                print('    {}'.format(self.zbox_line))
                
                # Shart shifting any atom that is outside of the box
                print('  Re-wrapping atoms based on new simulation cell.')
                lx = xhi - xlo; ly = yhi - ylo; lz = zhi - zlo; nwrapped = 0;
                for i in self.atoms:
                    atom = self.atoms[i]
                    x = atom.x; y = atom.y; z = atom.z;
                    # Move x atoms and update ix flag
                    if x > xhi:
                        x, ix = wrap_atom_compute_iflag(x, xhi, lx)
                    elif x < xlo:
                        x, ix = wrap_atom_compute_iflag(x, xlo, lx)
                    else:
                        x = x # else use non-periodic position
                        ix = 0 # else set flag a zero
                        
                    # Move y atoms and update iy flag
                    if y > yhi:
                        y, iy = wrap_atom_compute_iflag(y, yhi, ly)
                    elif y < ylo:
                        y, iy = wrap_atom_compute_iflag(y, ylo, ly)
                    else:
                        y = y # else use non-periodic position
                        iy = 0 # else set flag a zero
                        
                    # Move z atoms and update iz flag
                    if z > zhi:
                        z, iz = wrap_atom_compute_iflag(z, zhi, lz)
                    elif z < zlo:
                        z, iz = wrap_atom_compute_iflag(z, zlo, lz)
                    else:
                        z = z # else use non-periodic position
                        iz = 0 # else set flag a zero
                    
                    # Update atom info
                    atom.x = x; atom.ix = ix;
                    atom.y = y; atom.iy = iy;
                    atom.z = z; atom.iz = iz;
                    if ix != 0 or iy != 0 or iz != 0:
                        nwrapped += 1
                print(f'    Re-wrapped {nwrapped} atoms.')
            else: print('  Skipping "VMD_graphite.mol2" option as system is a single graphene sheet')
            print()        