#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.4
June 5th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
# Function to check if variable is a float
def check_float(variable):
    try:
        float(variable)
        return_boolean = True
    except:
        return_boolean = False
    return return_boolean

# Open and read .frc file and store info in dictionaries and classes
class Atom_types: pass # .type .mass, .element, .masses
class forcefield_file:
    def __init__(self, frc_file):
        # set frc file type
        self.type = 'interatomic'
        
        # atom types info/bondics
        self.atom_types = {}                # reaxff {atom type : atom object}
        self.extend_elements = []           # List of elements to extend for atom typing purposes
        self.bond_increments = {}           # class1/class2 but will be empty for reaxFF (using reset_charges will then zero charges)
        self.equivalences = {}              # class1/class2 but will be empty for reaxFF
        self.auto_equivalences = {}         # class1/class2 but will be empty for reaxFF
        self.pair_coeffs_9_6 = {}           # class1/class2 but will be empty for reaxF
        self.pair_coeffs_12_6 = {}          # class1/class2 but will be empty for reaxF
        self.filename = frc_file
        
        # Box controls (intialize as defaults and update if found)
        self.nx = 1; self.ny = 1; self.nz = 1;
        self.ax = 0.5; self.ay = 0.5; self.az = 0.5;
        self.zero_dim_buffer = 0.5;

        # Opening and reading the frc file
        with open(frc_file, 'r') as f:
            
            # Initializing flags    
            atomtypes_flag = False
            cell_flag = False
            extend_flag = False
            
            
            # Looping through each line of the frc file
            for line in f:
                
                # Strip comment's
                line = line.split('!')[0]
                line = line.rstrip()

                # Split line and strip line
                line = line.strip()
                line = line.split()
                
                # Setting flags if '#' character is in string of 1st element in list
                # as False to break each section flag from previous iteration
                if len(line) > 0 and line[0][0] == '#':
                    # atom types info/bondics
                    atomtypes_flag = False
                    cell_flag = False
                    extend_flag = False
                    
                # atom types info
                if '#atom_types' in line:
                    atomtypes_flag = True
                elif '#simulation_cell' in line:
                    cell_flag = True
                elif '#extend_elements' in line:
                    extend_flag = True

                    
                # Finding atom types information    
                if atomtypes_flag:
                    try:
                        # Check if len(line) >= 5 and if line[0] is a float
                        if len(line) >= 5 and check_float(line[0]):
                            #print(line)
                            a = Atom_types()
                            ver = float(line[0])
                            ref = float(line[1])
                            Type = line[2]
                            a.type = Type
                            a.mass = float(line[3])
                            a.element = Type
                            a.masses = [float(i) for i in line[4:]]
                            a.ver = ver
                            a.ref = ref
                            self.atom_types[Type] = a
                    except: # update self.type and perform error handling in all2lmp file
                        self.type = 'Inconsistent force field file'
                        pass

                        
                        
                # Find cell options and update
                elif cell_flag and len(line) >= 2:
                    try:
                        if line[0] == 'Nx:':
                            self.nx = float(line[1])
                        elif line[0] == 'Ny:':
                            self.ny = float(line[1])
                        elif line[0] == 'Nz:':
                            self.nz = float(line[1])
                        elif line[0] == 'Ax:':
                            self.ax = float(line[1])
                        elif line[0] == 'Ay:':
                            self.ay = float(line[1])
                        elif line[0] == 'Az:':
                            self.az = float(line[1])
                        elif line[0] == 'zero_dim_buffer:':
                            self.zero_dim_buffer = float(line[1])
                    except: # update self.type and perform error handling in all2lmp file
                        self.type = 'Inconsistent force field file'
                        pass
                    
                # Find extend elements flag
                elif extend_flag:
                    try:
                        if line[0] == 'Add:' and len(line) > 1:
                            self.extend_elements = line[1:]
                    except:
                        pass
                
                        
# Function to generate nta file for reaxFF conversions                     
def generate_nta(m, frc, reset_charges, log):
    nta = {} # { atomid : new atom type}
    name = {} # {atomid : nta:NAME}
    edge = {} # return empty since interatomic force fields does not use this option
    charges = {} # return empty since reaxff does not use this option
    neutralize = {'all': False, 'bond-inc': False, 'user-defined': False, 'zero': False} # Set as False since ReaxFF wont use this option
    remove = {'angle-nta':[], 'dihedral-nta':[], 'improper-nta':[], 'angle-ID':[], 'dihedral-ID':[], 'improper-ID':[], 'zero':{'angle':False, 'dihedral':False, 'improper':False}}
    
    # Function to find elements
    def find_element(atomtype, log):
        element = ''; # set element as empty and update later on
        mass = m.masses[atomtype].coeffs[0]
        for i in frc.atom_types:
            if mass in frc.atom_types[i].masses:
                element = frc.atom_types[i].element
                break
        # warn and exit if element is still empty
        if element == '': log.error(f'\nERROR mass {mass} is not in {frc.filename}. Add {mass} into Masses column')
        return element
    
    # Find nta atom types for interatomic force fields and find x, y, z pos for box positioning
    x = []; y = []; z = []; system_periodicity = [];
    for i in m.atoms:
        atom = m.atoms[i]
        
        # reset charges if user wants
        if reset_charges:
            atom.charge = 0 # zero charge for reaxFF simulations
        
        # Try getting element
        try: element = atom.element
            
        # except try to find element from mass from m and masses from frc
        except: element = find_element(atom.type, log)

        # Add element as nta
        nta[i] = element
        name[i] = element
        
        # Save x, y, and z data for box scaling
        x.append(atom.x); y.append(atom.y); z.append(atom.z);
        system_periodicity.append([atom.ix, atom.iy, atom.iz]) 
        
    # add extend_elements if they are found in read in .frc file
    if frc.extend_elements:
        # print warning
        log.warn(f'\nWARNING interatomic force field conversion is using {frc.filename} with #extend_elements Add: {" ".join(frc.extend_elements)}. This option is meant to make')
        log.out(f'atom types consistant between files. If this is not desired delete {" ".join(frc.extend_elements)} after Add: the keyword.')
        
        # Add elements to nta dict with 'extend'n key, where n is index from element in extend_elements. All atom types is a integer thus keys will never
        # be mistaken for atom typing each atom when to code performs that analysis. The find_BADI.atom_types(nta) will code will only loop through nta dict
        for n, i in enumerate(frc.extend_elements):
            extend_key = 'extend' + str(n)
            nta[extend_key] = i

    ####################################################
    # reset box dimensions based on info from frc file #
    # - if and only if the system has zero image flags # 
    # - else set box dims based on read in box dims    #
    ####################################################
    def check_images(system_periodicity):
        periodicity = False
        for iflags in system_periodicity:
            if any(iflags) != 0: periodicity = True; break
        return periodicity  
    # check periodicity and if not periodic reset box dims
    if not check_images(system_periodicity):
        # Find x, y, z box dims based on info from frc file
        xlo = frc.nx*min(x); xhi = frc.nx*max(x);
        ylo = frc.ny*min(y); yhi = frc.ny*max(y);
        zlo = frc.nz*min(z); zhi = frc.nz*max(z);
        
        # Add ai's to each direction
        if xlo <  0: xlo -= frc.ax
        if xhi >= 0: xhi += frc.ax 
        if ylo <  0: ylo -= frc.ax
        if yhi >= 0: yhi += frc.ax
        if zlo <  0: zlo -= frc.ax
        if zhi >= 0: zhi += frc.ax
     
        # Check if any values are zero if so reset to +- buffer
        # update buffer if zero and override with 0.5 angstroms, becuase zero volume is not desired
        if frc.zero_dim_buffer == 0: 
            buffer = 0.5; log.out(f'Interatomic_force_fields zero_dim_buffer was zero, resulting in a {buffer} angstrom over ride')
        else: buffer = frc.zero_dim_buffer

        if xlo == 0:
            xlo = xlo - buffer;
        if xhi == 0:
            xhi = xhi + buffer;
        if ylo == 0:
            ylo = ylo - buffer;
        if yhi == 0:
            yhi = yhi + buffer;
        if zlo == 0:
            zlo = zlo - buffer;
        if zhi == 0:
            zhi = zhi + buffer;
        
        # Update m.box lines
        m.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(xlo, xhi, 'xlo', 'xhi')
        m.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(ylo, yhi, 'ylo', 'yhi')
        m.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(zlo, zhi, 'zlo', 'zhi')
    return nta, name, edge, charges, neutralize, remove