# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
January 14th, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

This script allows for a material studio .car and .mdf file to be read in
and an Molecule class 'm' to be created on the info from the .car and .mdf
file, such that the class is exactly the same that read_lmp creates, when
reading in a LAMMPS data file. This allows for a .car and .mdf file to have
the exact same class type objects that exists when reading a LAMMPS data
file and allows the set of codes to read a .car and .mdf file as if it was
a LAMMPS data file.

The atom cooridinates will be centered around the geometrical
center of the read in molecule and then the box will be centered
around the newly re-centered moleucle.

IF PBC=NO in the .car file
   The box will be 0.5 angstroms larger in the x, y, and z directions and
   all image flags will be set to zero. The atom charge will be set to 
   zero and updated later on and the atom molid will be set as 1 and will
   stay as 1 throughout the rest of the codes that use the information
   found by this code.
   
IF PBC=YES and Space Group is P1 in the .car file
   The box will be set by the values in the .car file. This code is also
   equiped to handle triclinic cells and atom wrapping/image flag updates.
      
This script could serve as a starting point to create other
_ _ _2lmp scripts to have the entire set of scripts be able
to read in other molecule file types and be able to convert
a large number of file types to a lammps datafile.  
"""

##############################
# Import Necessary Libraries #
##############################
import math
import sys
import os


############################################
# Class to read and organize msi .car file #
############################################
class unit_cell: pass #  [PBC=ON .a .b .c .alpha .beta .gamma .space_group] [PBC=2D .l .k .gamma .plane_group]
class Atom_Car: pass # .type .molid .charge .x .y .z .atom_type .element .ix .iy .iz (iflags will be zeroes and updated later)
class read_car:
    def __init__(self, incarfile):
        self.filename = incarfile # Filename
        self.header = ''
        self.pbc_condition = 'OFF' # update later 
        self.title = '' # set as empty string and update later 
        self.date = '' # set as empty string and update later 
        self.atoms = {}  # { atomid : atom object }
        self.atoms_map = {} # { msi atomid : lmp numeric id}
        self.nta = {} # { lmp number atomid: atom type}
        
        # Initialize unit cell
        uc = unit_cell()
        uc.a = 1 # Default box size in 3D
        uc.b = 1 # Default box size in 3D
        uc.c = 1 # Default box size in 3D
        uc.l = 1 # Default box size in 2D
        uc.k = 1 # Default box size in 2D
        uc.alpha = 90 # Default tilt factor
        uc.beta =  90 # Default tilt factor
        uc.gamma = 90 # Default tilt factor
        uc.space_group = ''
        uc.plane_group = ''
        self.cell = uc
        
        # Open and read file
        with open(incarfile, 'r') as f:
            
            # Itializatizing flags
            header_flag = False
            pbc_flag = False
            title_flag = False
            date_flag = False
            cell_flag = False
            atoms_flag = False
            atoms_lst = [] # list to add atom lines to for futher spliting and parsing
            molecule_count = 1 # To count molecules (start at 1 since updated at the end of each molecule via the 'end' keyword)
            
            # Looping through each line of the file
            for n, line in enumerate(f):
                
                # split and strip line
                line_strip = line.strip()
                line_split = line_strip.split()
                
                # Go through and flag lines
                if line_strip == '':
                    header_flag = False
                    pbc_flag = False
                    title_flag = False
                    date_flag = False
                    cell_flag = False
                    atoms_flag = False
                elif line_strip == 'end':
                    header_flag = False
                    pbc_flag = False
                    title_flag = False
                    date_flag = False
                    cell_flag = False
                    atoms_flag = False
                    molecule_count += 1 # if end is found update molecule_count
                elif n == 0:
                    header_flag = True
                elif 'PBC=' in line_strip and n < 3:
                    header_flag = False
                    pbc_flag = True
                elif n == 2:
                    pbc_flag = False
                    title_flag = True
                elif line_split[0] == '!DATE':
                    title_flag = False
                    date_flag = True
                elif 'PBC' in line_strip and 'PBC=' not in line_strip:
                    date_flag = False
                    cell_flag = True
                elif len(line_split) >= 9:
                    cell_flag = False
                    atoms_flag = True
                    
                # Find header info
                if header_flag:
                    self.header = ' '.join(line_split)
                
                # Find pbc flag and set findings
                if pbc_flag:
                    # Find pbc condition and set it
                    if 'PBC=' in line_strip:
                        self.pbc_condition = line_split[0].split('=')[-1]
                
                # Find title information
                if title_flag: self.title = line_strip
                 
                # Find date information
                if date_flag: self.date = line_strip
                    
                # Find cell information
                if cell_flag:
                    # Find unit cell if PBC=ON
                    if self.pbc_condition == 'ON':
                        # 2D specific and update with true parms if found
                        k = float(0); l = float(0); gamma = float(0); plane_group = ''
                        
                        # ON specific
                        a = float(line_split[1])
                        b = float(line_split[2])
                        c = float(line_split[3])
                        alpha = float(line_split[4])
                        beta = float(line_split[5])
                        gamma = float(line_split[6])
                        space_group = line_split[7].replace('(', '').replace(')', '')
                       
                    # Find unit cell if PBC=2D (THIS WONT BE USED YET....)
                    elif self.pbc_condition == '2D':
                        # ON specific and update with true parms if found
                        a = float(0); b = float(0); c = float(0);
                        alpha = float(0); beta = float(0); space_group = '';
                        
                        # 2D specific
                        k = float(line_split[1])
                        l = float(line_split[2])
                        gamma = float(line_split[3])
                        plane_group = line_split[4].replace('(', '').replace(')', '')
                        
                    # else zeros and empyt strings for all
                    else:
                        # ON specific
                        a = float(0); b = float(0); c = float(0);
                        alpha = float(0); beta = float(0); gamma = float(0);
                        
                        # 2D specific
                        k = float(0); l = float(0); gamma = float(0); plane_group = ''
                        
                    # Save info
                    self.cell.a = a
                    self.cell.b = b
                    self.cell.c = c
                    self.cell.l = l
                    self.cell.k = k
                    self.cell.alpha = alpha
                    self.cell.beta = beta
                    self.cell.gamma = gamma
                    self.cell.space_group = space_group
                    self.cell.plane_group = plane_group


                # Find atoms information and append to list to deal with outside of this for loop
                # Adding molecule count to line_split to keep track of molecules outside of this loop
                if atoms_flag:
                    # recreate line_split and add molecule count to end of list
                    atom = [i for i in line_split] + [molecule_count]
                    atoms_lst.append(atom)

                    
            # Find atoms information (setting numerical id based in index in atoms_lst)
            for atomid, lst in enumerate(atoms_lst):
                lmpatomid = atomid + 1 # add 1 to atomid so it starts at 1
                molid = lst[-1] # molecule id was inserted at the end of the list
                msiatomid = lst[0]
                x = float(lst[1])
                y = float(lst[2])
                z = float(lst[3])
                residue_type = lst[4]
                residue_name = lst[5]
                residue = '{}_{}'.format(residue_type, residue_name)
                atom_type = lst[6]
                element = lst[7]
                charge = float(lst[8])
                
                # Add atoms to atoms instance
                a = Atom_Car()
                a.type = 1
                a.molid = molid
                a.charge = charge
                a.x = x
                a.y = y
                a.z = z
                a.atom_type = atom_type
                a.element = element
                a.ix = 0
                a.iy = 0
                a.iz = 0
                self.atoms[lmpatomid] = a
                
                # Find atom unique name by stringing together msiatomid+residue+molid with '-'
                # character, this will set the unqiue string to be used to map from msi atomid's
                # to LAMMPS numeric atomid's. Will be done this way in read_mdf as well.
                unique_name = '{}-{}-{}'.format(str(msiatomid), str(residue), str(molid))
                self.atoms_map[unique_name] = lmpatomid
                
                # create nta dictionary
                self.nta[lmpatomid] = atom_type
                
            
############################################
# Class to read and organize msi .mdf file #
############################################
class read_mdf:
    def __init__(self, inmdffile, msi2lmp_atomid_map):
        self.filename = inmdffile # Filename        
        self.title =  '' # set as empty string and update later on
        
        # Initialize all column indexes as -1 and update later (use -1
        # and check for positve number later to determine file setup)
        self.element = -1          # column info for element
        self.atom_type = -1        # column info for atom_type
        self.charge_group = -1     # column info for charge_group
        self.isotope = -1          # column info for isotope
        self.formal_charge = -1    # column info for formal_charge
        self.charge = -1           # column info for charge
        self.charge_group = -1     # column info for charge_group
        self.switching_atom = -1   # column info for switching_atom
        self.oop_flag = -1         # column info for oop_flag
        self.chirality_flag = -1   # column info for chirality_flag
        self.occupancy = -1        # column info for occupancy
        self.xray_temp_factor = -1 # column info for xray_temp_factor

        self.molecule = [] # set as empty list to append to
        self.periodicity = '' # set as empty string and update later on
        self.group = '' # set as empty string and update later on
        self.bonds = set([]) #  set of sorted of tuple bonding atomids

        
        # Open and read file
        with open(inmdffile, 'r') as f:
            
            # Itializatizing flags
            header_flag = False
            topology_flag = False
            molecule_flag = False
            bonding_flag = False
            periodicity_flag = False
            group_flag = False
            end_flag = False
            molecule_count = 0 # To count molecules (start at 0 since updated before reading molecule)
            
            # intialize indexes with defaults and update later on (use -1 and check for positve number later)
            column_index = {'element': -1, 'atom_type': -1, 'charge_group': -1, 'charge': -1, 'isotope': -1,
                            'formal_charge': -1, 'switching_atom': -1, 'oop_flag': -1, 'chirality_flag': -1,
                            'occupancy': -1, 'xray_temp_factor': -1, 'connections': -1}

            # Looping through each line of the file
            for n, line in enumerate(f):
                
                # Strip comment's
                line = line.split('!')[0]
                line = line.rstrip()

                # split and strip line
                line_strip = line.strip()
                line_split = line_strip.split()
                
                # Go through and flag lines
                if line_strip == '' or 'end' in line_split:
                    header_flag = False
                    molecule_flag = False
                    bonding_flag = False
                    group_flag = False
                elif n == 0:
                    header_flag = True
                elif '#topology' in line_strip or '@column' in line_strip:
                    header_flag = False
                    topology_flag = True
                elif '@molecule' in line_strip:
                    topology_flag = False
                    molecule_flag = True
                    molecule_count += 1
                elif len(line_split) > 5 and not end_flag:
                    molecule_flag = False
                    bonding_flag = True                    
                elif '@periodicity' in line_strip:
                    bonding_flag = False
                    periodicity_flag = True
                elif '@group' in line_strip:
                    periodicity_flag = False
                    group_flag = True
                elif '!' in line_strip or '#' in line_strip or '@' in line_strip or 'end' in line_strip:
                    end_flag = True
                    
                    
                # Find title info from header
                if header_flag:
                    self.title = line_strip
                    
                # Find topology indexes
                elif topology_flag:
                    # skip over empty spaces
                    if line_strip != '':
                        # Find column type index and update column_index dict
                        if '@column' in line_strip:
                            column_type = line_split[-1]
                            ind = line_split[-2]
                            if column_type in column_index:
                                column_index[column_type] = int(ind)
                
                # Find molecule name
                elif molecule_flag:
                    self.molecule.append(line_split[-1])
                
                # Find bonding information
                elif bonding_flag:

                    # Find atom unique name by stringing together msiatomid+residue+molid with '-' character, this
                    # will set the unqiue string to be used to map from msi atomid's to LAMMPS numeric atomid's.
                    split_zero = line_split[0].split(':')
                    if len(split_zero) == 2:
                        msiatomid = split_zero[-1]
                        residue = split_zero[0]
                        unique_name = '{}-{}-{}'.format(str(msiatomid), str(residue), str(molecule_count))
    
                        
                        # Just in case resname's are inconsitant between .mdf and .car file
                        try:
                            lmpatomid = msi2lmp_atomid_map[unique_name]
                        except:
                            print(f'resname is inconsistant bewteen .car and .mdf file. resname in .mdf file trying to be mapped onto .car file {residue} for atom: {msiatomid}')
                            sys.exit()
                        
                        # Find element if updated index is positive else set as empty string
                        if column_index['element'] >= 0:
                            self.element = line_split[column_index['element']]
                            
                        # Find atom type if updated index is positive else set as empty string
                        if column_index['atom_type'] >= 0:
                            self.atom_type = line_split[column_index['atom_type']]
      
                        # Find charge group if updated index is positive else set as empty string
                        if column_index['charge_group'] >= 0:
                            self.charge_group = line_split[column_index['charge_group']]
                            
                        # Find isotope if updated index is positive else set as empty string
                        if column_index['isotope'] >= 0:
                            self.isotope = line_split[column_index['isotope']]
                            
                        # Find formal charge if updated index is positive else set as empty string
                        if column_index['formal_charge'] >= 0:
                            self.formal_charge = line_split[column_index['formal_charge']]
    
                        # Find charge if updated index is positive else set as empty string
                        if column_index['charge'] >= 0:
                            self.charge = line_split[column_index['charge']]
                            
                        # Find switching atom if updated index is positive else set as empty string
                        if column_index['switching_atom'] >= 0:
                            self.switching_atom = line_split[column_index['switching_atom']]
                            
                        # Find oop flag if updated index is positive else set as empty string
                        if column_index['oop_flag'] >= 0:
                            self.oop_flag = line_split[column_index['oop_flag']]
                            
                        # Find chirality flag if updated index is positive else set as empty string
                        if column_index['chirality_flag'] >= 0:
                            self.chirality_flag = line_split[column_index['chirality_flag']]
                            
                        # Find occupancy if updated index is positive else set as empty string
                        if column_index['occupancy'] >= 0:
                            self.occupancy = line_split[column_index['occupancy']]
                            
                        # Find xray_temp_factor if updated index is positive else set as empty string
                        if column_index['xray_temp_factor'] >= 0:
                            self.xray_temp_factor = line_split[column_index['xray_temp_factor']]
                            
                        # Find connections if updated index is positive else set as empty string
                        if column_index['connections'] >= 0:
                            ind = column_index['connections']
                            connections = [line_split[i] for i in range(ind, len(line_split), 1)]
                        else:
                            connections = []
                        
                        # Split connections that have the msi record keeping tags, except the residue ':' character to split later on
                        split_connections = []; Connectivity_Record_Items = ['%', '#', '/', ','] # pg 39 Product Name - Current File Formats .pdf
                        for bonded in connections:
                            # Check that symmetry opertion is 1, if not raise Excetion
                            try: symop = int(bonded.split('#')[1]) # try getting
                            except: symop = 1 # If symop is non existent assume to be 1
                            if symop != 1:
                                print(f'ERROR this tool is not equipped to handle symmetry operations that are not 1. symop in file: {symop}')
                                sys.exit()
                            
                            # find msiatomid and string together
                            msiatomid_bonded = ''
                            for k in bonded:
                                if k not in Connectivity_Record_Items:
                                  msiatomid_bonded += k
                                else: break
                        
                            # Find residue_bonded if ':' in msiatomid_bonded
                            if ':' in msiatomid_bonded:
                                res_split = msiatomid_bonded.split(':')
                                residue_bonded = res_split[0]
                                msiatomid_bonded = res_split[1]
                             # else set residue_bonded from residue above
                            else: residue_bonded = residue
                            
                            # Find atom unique name by stringing together msiatomid+residue+molid with '-' character, this
                            # will set the unqiue string to be used to map from msi atomid's to LAMMPS numeric atomid's.
                            unique_name = '{}-{}-{}'.format(str(msiatomid_bonded), str(residue_bonded), str(molecule_count))
                            split_connections.append(unique_name)
    
                        # Find lmpatomids from msi connections to build bonds list
                        lmpbondingids = [msi2lmp_atomid_map[i] for i in split_connections]
                        
                        # create bonds from lmpbondingids and lmpatomid
                        for lmpbondingatomid in lmpbondingids:
                            bond = tuple(sorted([lmpbondingatomid, lmpatomid]))
                            self.bonds.add(bond)
                    
                # Find periodicity information
                elif periodicity_flag:
                    self.periodicity = line_split[1:]

                # group information
                elif group_flag:
                    self.group = line_split[-1].replace('(', '').replace(')', '')
                    
        # sorting bonds to give nice ordering of bonds
        self.bonds = sorted(self.bonds)

        
#######################################################
### Class for converting all info material studio   ### 
### .car and .mdf file into a class that is exactly ###
### the same as it comes from read_lmp to be able   ###
### to use .car and .mdf files in remaining of the  ###
### codes.                                          ###
#######################################################
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz
class Bond: pass  # .type .atomids = [atom1id, atom2id]
class Molecule_File:
    def __init__(self, incarfile, inmdffile, wrap_atoms, no_center):
        ###############################
        # Reading .car and .mdf files #
        ###############################
        # Read car file and find space group
        car = read_car(incarfile) 
        
        # Read mdf file with atoms map from car file for setting consitant numeric atomids
        mdf = read_mdf(inmdffile, car.atoms_map)
        
        ################################
        # Error and exiting conditions #
        ################################
        # Create error and exit is space group is not P1
        if car.pbc_condition == 'ON':
            print('PBC conditions were turned ON. The simulation cell will be set based the info in')
            print(f'{incarfile} file.')
            if car.cell.space_group != 'P1':
                print('ERROR this tool only works for space group P1 (no other symmetry operations are currently supported)')
                sys.exit()
        else:
            print('PBC conditions were turned OFF. The simulation cell will be orthogonal and set')
            print('based on the extent of the atoms +- 0.5 angstroms in the 3 principal directions')
            
        # Create warning if mdf has no connections column (not exiting, just in case for reaxFF conversions)
        if len(mdf.bonds) == 0:
            print('WARNING mdf file does not contain connections column. Only valid use for this .mdf file is to convert to reaxff')
        
        ##############################
        # Create/intialize instances #
        ##############################        
        # Set combined file names
        self.filename = car.filename + '_+_' + mdf.filename # Filename
        
        # set nta dictionary
        self.nta = car.nta
        self.name = {i:self.nta[i] for i in self.nta}
        
        # Set system information
        self.natoms = len(car.atoms) # total atoms in .car file
        self.nbonds = len(mdf.bonds) # total bonds in .mdf file
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
        self.header = '{} {} {} {}'.format('HEADER, ', os.path.basename(inmdffile), os.path.basename(incarfile),' read w/ msi2lmp')
        self.velocities = {} # { atomid : tuple(velx, vely, velz)}


        ####################################################################
        # Find center of atoms and shift needed to center atoms (will only #
        # be applied on orthogonal boxes and not triclinic for now ...     #
        ####################################################################
        x = []; y = []; z = [];
        for i in car.atoms:
            atom = car.atoms[i]
            x.append(atom.x)
            y.append(atom.y)
            z.append(atom.z)
    
        # Finding geometric center
        x_center = sum(x)/len(x); y_center = sum(y)/len(y); z_center = sum(z)/len(z)
        
        # Finding needed shift to center around geometric center
        x_center = (0 - x_center); y_center = (0 - y_center); z_center = (0 - z_center);
        
        # If no_center is True set i_centers to zero; i = x, y, z
        if no_center:
            x_center = 0; y_center = 0; z_center = 0;

        ###########################################################################################
        # find if read in .car coordinates are periodic, if they are perform periodic operations. #
        # else search min/max values and determine simulation bix dimensions in that fashion      #
        ###########################################################################################
        if car.pbc_condition == 'ON':
            a = car.cell.a
            b = car.cell.b
            c = car.cell.c
            alpha = car.cell.alpha
            beta = car.cell.beta
            gamma = car.cell.gamma
            triclinic_flag = False
            
            ######################################################################################################
            # If alpha, beta, or gamma does not equal 90 the box is triclinic and perform triclinic oppertations #
            ######################################################################################################
            if alpha != 90 or beta != 90 or gamma != 90:
                triclinic_flag = True # Update triclinic flag
                
                # set center to zeros for triclinic systems and use min/max values (If you want to center
                # the box comment out the values below which will enforce the true center values from above)
                x_center = 0; y_center = 0; z_center = 0;

                # Find  xy, xz, yz tilt factors to handle triclinic cells
                # https://docs.lammps.org/Howto_triclinic.html
                xy = self.xy; xz = self.xz; yz = self.yz; # intialize all as zeros from above and update later
                lx = a
                xy = b*math.cos(math.radians(gamma))
                ly = math.sqrt(b**2 - xy**2)
                xz = c*math.cos(math.radians(beta))
                if abs(math.sin(math.radians(gamma))) > 0.0001:
                    yz = ((b*c)*math.cos(math.radians(alpha)) - (xy*xz))/ly
                else:
                    yz = 0.0
                lz = math.sqrt(c**2 - xz**2 - yz**2)
                                
                # Update self instances
                self.xy = xy; self.xz = xz; self.yz = yz;
                
                # Generate h and h_inv vector like LAMMPS does
                # Taken from LAMMPS domain.cpp file (lammps-23Jun2022)
                h = 6*[0]; h_inv = 6*[0] # set defaults as zeros and then update later on
                boxlo = 3*[-0.5] # intialize lo dims like LAMMPS
                
                # Set h box
                h[0] = lx; h[1] = ly; h[2] = lz;
                h[3] = yz; h[4] = xz; h[5] = xy;

                # Set h-inverse box
                h_inv[0] = 1/h[0]; h_inv[1] = 1/h[1]; h_inv[2] = 1/h[2]
                h_inv[3] = -h[3] / (h[1]*h[2]);
                h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2]);
                h_inv[5] = -h[5] / (h[0]*h[1]);
                
                # set ilo/ihi based on min li's and i_centers; were i = x, y, z 
                xlo = min(x) + x_center - 0.5; xhi = lx + xlo
                ylo = min(y) + y_center - 0.5; yhi = ly + ylo 
                zlo = min(z) + z_center - 0.5; zhi = lz + zlo
                boxlo = [xlo, ylo, zlo]


                # ----------------------------------------------------------------------
                #   convert triclinic 0-1 lamda coords to box coords for one atom
                #   x = H lamda + x0;
                #   lamda and x can point to same 3-vector
                #-----------------------------------------------------------------------
                # converted from c to python from msi2lmp orginal code
                def lamda2pos(lamda, h, boxlo):
                    posx = h[0]*lamda[0] + h[5]*lamda[1] + h[4]*lamda[2] + boxlo[0]
                    posy = h[1]*lamda[1] + h[3]*lamda[2] + boxlo[1]
                    posz = h[2]*lamda[2] + boxlo[2]
                    return [posx, posy, posz]
                
                
                # ----------------------------------------------------------------------
                #   convert box coords to triclinic 0-1 lamda coords for one atom
                #   lamda = H^-1 (x - x0)
                #   x and lamda can point to same 3-vector
                #-----------------------------------------------------------------------
                # converted from c to python from msi2lmp orginal code
                def pos2lamda(pos, h_inv, boxlo):
                    delta = [0 for _ in range(3)] # [0, 0, 0]
                    delta[0] = pos[0] - boxlo[0]
                    delta[1] = pos[1] - boxlo[1]
                    delta[2] = pos[2] - boxlo[2]

                    lamdax = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2]
                    lamday = h_inv[1]*delta[1] + h_inv[3]*delta[2]
                    lamdaz = h_inv[2]*delta[2]
                    return [lamdax, lamday, lamdaz]

                # Find atoms info
                for i in car.atoms:
                    atom = car.atoms[i]
                    
                    # Find atom coords and image flags
                    x = atom.x + x_center;
                    y = atom.y + y_center;
                    z = atom.z + z_center;
                    ix = atom.ix; iy = atom.iy; iz = atom.iz;
                    
                    # If wrap_atoms perfrom wrapping task, otherwise use defaults from above
                    if wrap_atoms:
                        # build pos-vector
                        pos = [atom.x, atom.y, atom.z]
    
                        # perform triclinic wrapping like msi2lmp                    
                        lamda = pos2lamda(pos, h_inv, boxlo)
                        images = 3*[0] # intialize as zero's
                        for k in range(0, 3):
                            tmp = math.floor(lamda[k])
                            images[k] = tmp
                            lamda[k] -= tmp
                        pos = lamda2pos(lamda, h, boxlo)
    
                        # Extract out new positions and image flags
                        x, y, z = pos
                        ix, iy, iz = images
                    
                    # Update atoms class with new info
                    a = Atom()
                    a.type = 1 # Set as 1 this information is not used by all2lmp 
                    a.molid = atom.molid
                    a.charge = atom.charge
                    a.element = atom.element
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    self.atoms[i] = a
                    self.velocities[i] = (0, 0, 0)


            ###############################################################
            # if box is orthoganl set lx, ly, and lz from a, b, c vectors #
            ###############################################################
            if not triclinic_flag:                
                # Assume molecule will be centered about 0, 0, 0 with i_center; i = x, y, z
                lx = a; ly = b; lz = c;
                xlo = -lx/2; xhi = lx/2;
                ylo = -ly/2; yhi = ly/2;
                zlo = -lz/2; zhi = lz/2;
                
                # function to find image flag
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
                        
                # Find atoms info
                for i in car.atoms:
                    atom = car.atoms[i]
                    
                    # Find atom coords and image flags
                    x = atom.x + x_center;
                    y = atom.y + y_center;
                    z = atom.z + z_center;
                    ix = atom.ix; iy = atom.iy; iz = atom.iz;
                    
                    # If wrap_atoms perfrom wrapping task, otherwise use defaults from above
                    if wrap_atoms:
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
                        
                    # Update atoms class with new info
                    a = Atom()
                    a.type = 1 # Set as 1 this information is not used by all2lmp 
                    a.molid = atom.molid
                    a.charge = atom.charge
                    a.element = atom.element
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    self.atoms[i] = a
                    self.velocities[i] = (0, 0, 0)

                   
        ###############################################     
        # If pbc flag is not on use min max searching #
        ###############################################
        else:
            ###################################################################
            # Find new atoms x, y, z position after shift to build lammps box #
            # size to produce image flags in x, y, and z direction as zeros   #
            ###################################################################
            x = []; y = []; z = [];
            for i in car.atoms:
                atom = car.atoms[i]                        
                # Update atoms class with new info
                a = Atom()
                a.type = 1 # Set as 1 this information is not used by all2lmp 
                a.molid = atom.molid
                a.charge = atom.charge
                a.element = atom.element
                a.x = atom.x + x_center
                a.y = atom.y + y_center
                a.z = atom.z + z_center
                x.append(atom.x + x_center)
                y.append(atom.y + y_center)
                z.append(atom.z + z_center)
                a.ix = atom.ix
                a.iy = atom.iy
                a.iz = atom.iz
                self.atoms[i] = a
                self.velocities[i] = (0, 0, 0)

                
            ###########################################################################
            # Find x, y, z box dims (search for min/max and then oversize slightly in #
            # each direction also  if certain dimensions are zero set default dim)    #
            ###########################################################################
            oversize = 0.5 # default oversize of 0.5 angtroms in each value (total over size = 0.5*2 = 1 angstroms)
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
        
        # Find bond information like read_lmp
        for n, atomids in enumerate(mdf.bonds):
            b = Bond()
            b.atomids = list(sorted(atomids))
            b.type = 1 # set as 1 (wont be using for anything any ways)
            self.bonds[n+1] = b