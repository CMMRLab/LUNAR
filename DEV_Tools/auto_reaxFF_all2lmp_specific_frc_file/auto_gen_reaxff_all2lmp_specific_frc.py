#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
December 19th, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

This code is meant to auto generate reaxff specific .frc file
to work with all2lmp. This file will be used to identify certain
masses from different force fields so all2lmp can back out element
types based on masses. Thus this code will find all masses for all
elements that exists in all .frc files that are in the current 
working directory
"""


##############################
# Import Necessary Libraries #
##############################
import read_frc as read_frc
from datetime import date
import os


##########
# Inputs #
##########
# Set frc file name
reaxff_frc = 'reaxff_all2lmp.frc'


# Set version of frc file
version = '1.0'


# set periodic_table_masses based on element to use as formal mass. IF element is not in periodic_table_masses
# the most frequent mass from the .frc file will be used. https://ptable.com/?lang=en#Properties
reference = {'ref#': '1', 'ref location': 'https://ptable.com/?lang=en#Properties'} # set reference info
periodic_table_masses = {'H': 1.008, 'Li': 6.94, 'Na': 22.990, 'K': 39.098, 'Rb': 85.468, 'Cs': 132.91, 'Fr': 223,  # Group 1
                         'Be': 9.0122, 'Mg':24.305, 'Ca': 40.078, 'Sr': 87.62, 'Ba': 137.33, 'Ra': 226,             # Group 2
                         'B': 10.81, 'Al': 26.982, 'Ga': 69.723, 'In': 114.82, 'TI': 204.38, 'Nh': 286,             # Group 13
                         'C': 12.00, 'Si': 28.085, 'Ge': 72.630, 'Sn': 118.71, 'Pb': 207.2, 'FI': 289,             # Group 14
                         'N': 14.000, 'P': 30.974, 'As': 74.922, 'Sb': 121.76, 'Bi': 208.98, 'Mc': 290,             # Group 15
                         'O': 15.999, 'S': 32.06, 'Se': 78.971, 'Te': 127.60, 'Po': 209, 'Lv': 293,                 # Group 16
                         'F': 18.998, 'Cl': 35.45, 'Br': 79.904, 'I': 126.90, 'At': 210, 'Ts': 294,                 # Group 17
                         'He': 4.0026, 'Ne': 20.180, 'Ar': 39.948, 'Kr': 83.798, 'Xe': 131.29, 'Rn': 222, 'Og': 294 # Group 18
                         }



# Find todays date
today = date.today()


# Find present working directory
pwd = os.getcwd()

# Find all .frc files that are in pwd      
frcfiles = [file for file in os.listdir(pwd) if file.endswith('.frc') if file != reaxff_frc]
print('\n\n\nRead in .frc files: ')
for i in frcfiles:
    print(i)
    

###############################################
# Read in all .frc files and save in frc_dict #
###############################################
frc_dict = {} # {.frc file name : frc object }
for i in frcfiles:
    frc = read_frc.forcefield_file(i)
    frc_dict[i] = frc
    
#####################
# Find all elements #
#####################
elements = set([])
for i in frc_dict:
    frc = frc_dict[i]
    
    # Loop through atom types
    for j in frc.atom_types:
        atom_type = frc.atom_types[j]
        element = atom_type.element
        elements.add(element)
# sort found elements
elements = list(elements) + list(periodic_table_masses.keys()) # add elements from periodic_table_masses
elements = sorted(elements)
print('n\n\nElements found in all files')
print(elements)

#########################################
# Create elements_dict to add masses to #
#########################################
elements_dict = {i:set([]) for i in elements}
for i in frc_dict:
    frc = frc_dict[i]
    
    # Loop through atom types
    for j in frc.atom_types:
        atom_type = frc.atom_types[j]
        element = atom_type.element
        mass = atom_type.mass
        elements_dict[element].add(mass)

# Add info from periodic_table_masses
for element in periodic_table_masses:
    elements_dict[element].add(periodic_table_masses[element])
#print(elements_dict)


###################
# Write .frc file #
###################
# Function to find most frequent occurance in list
def most_frequent(List):
    counter = 0; num = List[0];
    for i in List:
        curr_frequency = List.count(i)
        if(curr_frequency> counter):
            counter = curr_frequency
            num = i
    return num

# Writing new file with new atoms
with open(reaxff_frc,'w') as f: 
    # Write header
    f.write('!all2lmp reaxFF specific forcefield for mass definition and mass mapping from all .frc files used to generate this file.\n') 
    f.write('!The "!" character still acts as a comment character and all2lmp will disregard and information after the "!" character\n\n') 
    
    # Write date
    f.write(f'#version {reaxff_frc}     {version}     {today.strftime("%b-%d-%Y")}\n')
    f.write('#define all2lmp reaxFF .frc file\n')
    
    # Write atom types
    f.write('#atom_types     elements for reaxFF\n\n')
    f.write('> Atom type definitions for most variants of reaxFF potential (Ref column is for where Formal Mass comes from)\n')
    f.write(f'> Masses from {"   ".join(frcfiles)} frc files\n\n')
    f.write('> Formal Mass will be used to set mass of each ReaxFF element that gets set as an atom type (UPDATE AS NEEDED)\n')
    f.write('> Masses will be used to find element types for specific read in files to the code like LAMMPS .data files that do not have element symbols (UPDATE AS NEEDED)\n\n')

    f.write('!Ver   Ref    Element      Formal Mass      Masses\n')
    f.write('!---   -----  -------      -----------      ----------------------------------------\n')
    for i in elements_dict:
        # Find masses list
        masses = sorted(list(elements_dict[i]))
        
        # Set formal_mass from periodic_table_masses if element is in dict
        if i in periodic_table_masses:
            formal_mass = periodic_table_masses[i]
            ref = reference['ref#']
        # else set as most frequent for .frc files
        else:
            formal_mass = most_frequent(masses)
            ref = '0'
            
        # write to .frc file
        f.write('{:<6} {:<6} {:<12} {:<16} {:<18}\n'.format(version, ref, i, formal_mass, '   '.join([str(i) for i in masses])))
        
    
    # Write simulation cell options
    f.write('\n\n\n#simulation_cell USER DEFINED (UPDATE AS NEEDED)\n')
    f.write('Nx: 3 !(sets Ntimes as large of box in xlo and xhi direction - 1 sets box at min/max, 2 sets double lx, ...)\n')
    f.write('Ny: 3 !(sets Ntimes as large of box in ylo and yhi direction - 1 sets box at min/max, 2 sets double ly, ...)\n')
    f.write('Nz: 3 !(sets Ntimes as large of box in zlo and zhi direction - 1 sets box at min/max, 2 sets double lz, ...)\n')
    f.write('Ax: 0 !(adds Alength to box in xlo and xhi directions in angstroms - setting Nx: 1 as 1 will allow for control from min/max values)\n')
    f.write('Ay: 0 !(adds Alength to box in ylo and yhi directions in angstroms - setting Ny: 1 as 1 will allow for control from min/max values)\n')
    f.write('Az: 0 !(adds Alength to box in zlo and zhi directions in angstroms - setting Nz: 1 as 1 will allow for control from min/max values)\n')
    f.write('zero_dim_buffer: 3 !(adds dimension to zero dimension boxs in lo and hi directions - reaxFF will detect bonds across PBC during intialization that is not desired)\n')
    
    
    # Write extend elements options
    f.write('\n\n#extend_elements  USER DEFINED (UPDATE AS NEEDED)\n')
    f.write('!will add element to list of atom types in all2lmp internal data struture for atom type numbering - makes atom types consistent between different files for LAMMPS reaxff\n') 
    f.write('!pair_coeff command the "Add:" keyword must be below #extend_elements header! Elements will be seperated by white space. If no elements are listed after Add: keyword no\n')
    f.write('!elements will be extended. EXAMPLES:\n')
    f.write('!Add:     (only keyword specified to allow the code to read this file - no elements will be used to extend atom typing that code performs)')
    f.write('!Add: S   (will add S to all2lmp interanl data structure to keep atom types consistant between a file that has S and one that does not, but has every other element similar)\n')
    f.write('!Add: S N (will add S and N to all2lmp interanl data structure to keep atom types consistant between a file that one has S and the other has N, but has every other element similar)\n')
    f.write('Add: S\n')
    
    # Write references
    f.write('\n\n\n#reference 0\n')
    f.write('@Author tester\n')
    f.write(f'@Date {today.strftime("%b-%d-%Y")}\n')
    f.write(f'{"   ".join(frcfiles)} frc files\n')
    
    f.write('\n\n#reference 1\n')
    f.write('@Author tester\n')
    f.write(f'@Date {today.strftime("%b-%d-%Y")}\n')
    f.write(f'{reference["ref location"]}\n')