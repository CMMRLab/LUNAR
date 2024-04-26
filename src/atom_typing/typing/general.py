# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
February 6th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
from collections import OrderedDict
import math
import os


######################################################################
# Function to perform atom-typing for PCFF-IFF forcefield atom types #
######################################################################
def nta(mm, basename, ff_name):
    tally = {'found':0, 'assumed':0, 'failed':0} # To tally findings
    generate_small_lmp = True # Option to generate "small" LAMMPS data files to show atom type topology
    if len(mm.atoms) <= 100: generate_small_lmp = False # comment/uncomment, if the file has few atoms there really is no need to generate small lmp files
    
    ######################################
    # Set general informations and flags #
    ######################################
    # Log supported types
    supported_types = {'general':      ['has', 'no', 'specific', 'atom', 'types']}
    
    ######################################
    # Function to string neigh list info #
    ######################################
    def string_neigh_info(lst):
        complete = [tuple(i) for i in lst]
        unique = set([tuple(i) for i in lst])
        count = {i:complete.count(i) for i in unique}
        count = dict(OrderedDict( sorted(count.items(), key=lambda x:x[0][2]) )) # [0=keys;1=values][index position in key tuple lst]
        count = dict(OrderedDict( sorted(count.items(), key=lambda x:x[0][1]) )) # [0=keys;1=values][index position in key tuple lst]
        count = dict(OrderedDict( sorted(count.items(), key=lambda x:x[0][0]) )) # [0=keys;1=values][index position in key tuple lst]
        string = ''
        for i in count:
            string += '{};'.format(count[i])
            string += ''.join([str(j) for j in i])
            string += ','
        return string[:-1]
    

    ########################################################
    # Opening files to dump assumed/failed atom types into #
    ########################################################
    a = open(basename+'_assumed.txt', 'w')
    f = open(basename+'_failed.txt', 'w')


    ########################################################    
    # Loop through mm.atoms and start assigning atom-types #
    ########################################################
    for i in mm.atoms:
        # Set atom instance and pull out desired info
        atom = mm.atoms[i]#; formula = atom.molecule.formula;
        nb = int(atom.nb); ring = int(atom.ring); element = atom.element;
        info = '{}{}{}'.format(element, ring, nb)
        neighbor_info = atom.neighbor_info
        neigh1_info = string_neigh_info(neighbor_info[1])
        neigh2_info = string_neigh_info(neighbor_info[2])
        neigh3_info = string_neigh_info(neighbor_info[3])
        neigh4_info = string_neigh_info(neighbor_info[4])
        
        
        # Set intial .nta and .nta_comments and update later on if found
        atom.nta_type = '{}-type-yourself'.format(element)
        atom.nta_info = 'FAILED TO BE TYPED:  element: {}, ring: {}, nb: {}'.format(element, ring, nb)

        # Set new type and comment if ff_name = general:0
        if ff_name == 'general:0':
            atom.nta_type = info
            tally['found'] += 1;
            atom.nta_info = '{} type; where type = elementRINGnb'.format(ff_name)
            
        # Set new type and comment if ff_name = general:1
        if ff_name == 'general:1':
            atom.nta_type = '{}-({})'.format(info, neigh1_info)
            tally['found'] += 1;
            atom.nta_info = '{} type; where type = elementRINGnb-(1st-neighs) and (1st-neighs) = (count;elementRINGnb,...)'.format(ff_name)
            
        # Set new type and comment if ff_name = general:2
        if ff_name == 'general:2':
            atom.nta_type = '{}-({})-({})'.format(info, neigh1_info, neigh2_info)
            tally['found'] += 1;
            atom.nta_info = '{} type; where type = elementRINGnb-(1st-neighs)-(2nd-neighs) and (ith-neighs) = (count;elementRINGnb,...)'.format(ff_name)
            
        # Set new type and comment if ff_name = general:3
        if ff_name == 'general:3':
            atom.nta_type = '{}-({})-({})-({})'.format(info, neigh1_info, neigh2_info, neigh3_info)
            tally['found'] += 1;
            atom.nta_info = '{} type; where type = elementRINGnb-(1st-neighs)-(2nd-neighs)-(3rd-neighs) and (ith-neighs) = (count;elementRINGnb,...)'.format(ff_name)
            
        # Set new type and comment if ff_name = general:4
        if ff_name == 'general:4':
            atom.nta_type = '{}-({})-({})-({})-({})'.format(info, neigh1_info, neigh2_info, neigh3_info, neigh4_info)
            tally['found'] += 1;
            atom.nta_info = '{} type; where type = elementRINGnb-(1st-neighs)-(2nd-neighs)-(3rd-neighs)-(4th-neighs) and (ith-neighs) = (count;elementRINGnb,...)'.format(ff_name)
                
        
    ####################################################################################
    # Sort final mm.atoms[ID] to ensure written files have atomids in asscending order #
    ####################################################################################
    mm.atoms = dict(OrderedDict(sorted(mm.atoms.items())))
    
    
    ##################################################################
    # Add to the following instances to mm for data logging purposes #
    ##################################################################
    mm.supported_types = supported_types; mm.tally = tally; mm.assumed = {};
    
    
    ##########################################################################
    # Close the open files and check if they are zero in size; if so delete #
    #########################################################################
    a.close(); f.close(); asize = os.path.getsize(a.name); fsize = os.path.getsize(f.name);
    try:
        if asize == 0: os.remove(a.name);
        if fsize == 0: os.remove(f.name);
    except: a.close(); f.close();
    
    #######################################################################################
    # Generate "small" LAMMPS data files to show the topologies of each general atom type #
    #######################################################################################
    # Function to find bonds in "small" fragment
    def reduce_atoms_bonds(m, atoms2keep):
        bonds = [];
        for i in m.bonds:
            id1, id2 = m.bonds[i].atomids
            if id1 in atoms2keep and id2 in atoms2keep:
                bonds.append([id1, id2])
        return bonds
    
    # Function to compute distances
    def compute_distance(x1, y1, z1, x2, y2, z2):
        dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
        return math.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Function to get sign of a number
    def get_sign(number):
       if number >= 0:
          return 1
       elif number < 0:
          return -1
       else:
          return 0
    
    # Function to unwrap atoms via anchor method
    def unwrap_atoms_anchor(mm, atoms2keep, lx, ly, lz, cx, cy, cz):        
        # set max_x, max_y, max_z w/ minimum image convention
        max_x = lx/2; max_y = ly/2; max_z = lz/2;
        
        # find distances to center of box
        distances = {} # { atomID : dist2center }
        for i in atoms2keep:
            atom = mm.atoms[i]
            distances[i] = compute_distance(atom.x, atom.y, atom.z, cx, cy, cz)
            
        # find closets atom to center
        closests2center_ID = min(distances, key=distances.get)
        closests2center_atom = mm.atoms[closests2center_ID]
        basex = closests2center_atom.x
        basey = closests2center_atom.y
        basez = closests2center_atom.z
        
        # unwrap atoms
        atoms = {} # { atomID : [unwrapped-x, unwrapped-y, unwrapped-z]}
        for i in atoms2keep:
            atom = mm.atoms[i]
            x = atom.x; y = atom.y; z = atom.z;
            
            # Find vectors in each direction
            diffx = basex - x
            diffy = basey - y
            diffz = basez - z
            
            # Apply minimum image convention
            if abs(diffx) > max_x:
                x += get_sign(diffx)*lx
            if abs(diffy) > max_y:
                y += get_sign(diffy)*ly
            if abs(diffz) > max_z:
                z += get_sign(diffz)*lz
            atoms[i] = [x, y, z]
        return atoms
    
    # Function to unwrap atoms using image flags
    def unwrap_atoms_iflags(mm, atoms2keep, lx, ly, lz):
        atoms = {} # { atomID : [unwrapped-x, unwrapped-y, unwrapped-z]}
        for i in atoms2keep:
            atom = mm.atoms[i]
            x = atom.x + atom.ix*lx
            y = atom.y + atom.iy*ly
            z = atom.z + atom.iz*lz
            atoms[i] = [x, y, z]
        return atoms
    
    # Function to write "small" lmp data file
    def write_small_lmp(mm, atoms2keep, bonds2keep, smallfile, center=True):
        #--------------#
        # unwrap atoms #
        #--------------#
        # Find box dimensions
        xline = mm.xbox_line.split();
        yline = mm.ybox_line.split();
        zline = mm.zbox_line.split();
        
        # Lx, Ly, and Ly box dimensions
        lx = float(xline[1])-float(xline[0])
        ly = float(yline[1])-float(yline[0])
        lz = float(zline[1])-float(zline[0])
        
        # Find cx, cy, cz box centers
        cx = (float(xline[1])+float(xline[0]))/2
        cy = (float(yline[1])+float(yline[0]))/2
        cz = (float(zline[1])+float(zline[0]))/2
        
        # Get unwrapped atoms (anchoring is the best method since it does not require consistent iflags)
        #atoms = unwrap_atoms_iflags(mm, atoms2keep, lx, ly, lz)
        atoms = unwrap_atoms_anchor(mm, atoms2keep, lx, ly, lz, cx, cy, cz)

        #----------------------------#
        # Find shift to center atoms #
        #----------------------------#
        if center:
            x = []; y = []; z = [];
            for i in atoms:
                x.append(atoms[i][0])
                y.append(atoms[i][1])
                z.append(atoms[i][2])
            
            # Finding geometric center
            x_center = sum(x)/len(x); y_center = sum(y)/len(y); z_center = sum(z)/len(z)
            
            # Finding needed shift to center around geometric center
            x_shift = (0 - x_center); y_shift = (0 - y_center); z_shift = (0 - z_center);
                        
            #-------------------------------------------------------------------------#
            # Find x, y, z box dims (search for min/max and then oversize slightly in #
            # each direction also if certain dimensions are zero set default dim)     #
            #-------------------------------------------------------------------------#
            # shift positions as needed
            x = [i + x_shift for i in x]
            y = [i + y_shift for i in y]
            z = [i + z_shift for i in z]
            
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
            xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(xlo, xhi, 'xlo', 'xhi')
            ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(ylo, yhi, 'ylo', 'yhi')
            zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(zlo, zhi, 'zlo', 'zhi')
        else:
            xbox_line = mm.xbox_line
            ybox_line = mm.ybox_line
            zbox_line = mm.zbox_line
        
        #-----------------------------#
        # write small lammps datafile #
        #-----------------------------#
        with open(smallfile,'w') as f: 
            header = '{} > atom_typing {}'.format(mm.header, ff_name)
            f.write(f'{header[-220:len(header)]}\n\n') # Make max header length of 220 characters
            f.write(f'{len(atoms2keep)} atoms\n')
            f.write(f'{len(bonds2keep)} bonds\n\n')
            f.write(f'{mm.natomtypes} atom types\n')
            f.write('1 bond types\n\n')
            f.write(f'{xbox_line}\n{ybox_line}\n{zbox_line}\n')
            if mm.xy != 0 or mm.xz != 0 or mm.yz != 0:
                f.write('{:>12.9f} {:^9.9f} {:^9.9f} {} {} {}\n'.format(mm.xy, mm.xz, mm.yz, 'xy', 'xz', 'yz'))
                
            # Write masses
            f.write('\nMasses\n\n')
            for i in mm.masses: 
                comment = '{:^2} {:5}'.format('#', mm.masses[i].type)
                mass = mm.masses[i].coeffs[0]
                f.write('{:^3} {:^10.5f} {:^2}\n'.format(i, mass, comment))
                
            # Write atoms    
            f.write('\nAtoms # full\n\n')
            for i in atoms2keep:
                atom = mm.atoms[i]
                # Try to get comment data
                try:  comment = '{:^5} {:>3} {:>10}     avg-angle: {:>10.6f}      type: {}'.format('#', atom.comment, atom.hybridization, atom.avg_angle, atom.nta_type)
                except: comment = ''
                
                # Find atom postions
                x = atoms[i][0]; y = atoms[i][1]; z = atoms[i][2];
                if center:
                    x += x_shift
                    y += y_shift
                    z += z_shift
                
                # trying to write image flags if they exist (Skipping Iflags since VMD struggles with ringed systems and image flags)
                try: f.write('{:^3} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2} {:^5}\n'.format(i, atom.molid, atom.type, atom.charge, x, y, z, atom.ix, atom.iy, atom.iz, comment)) 
                except: f.write('{:^3} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^5}\n'.format(i, atom.molid, atom.type, atom.charge, x, y, z, comment)) 
                    
            # Write bonds
            f.write('\nBonds\n\n')
            for i, (id1, id2) in enumerate(bonds2keep, 1):
                f.write('{:^2} {:^2} {:^2} {:^2}\n'.format(i, 1, id1, id2)) 
        return
        
    # Find atoms 2 write and write them
    if generate_small_lmp:
        # Find every atomID and map it to the equiv string
        equivs = sorted({mm.atoms[i].nta_type for i in mm.atoms})
        equiv2atomIDs = {i:[] for i in equivs}
        for i in mm.atoms:
            atom = mm.atoms[i]
            equiv2atomIDs[atom.nta_type].append(i)
        
        # Find the atoms2keep for each atom type string
        for string in equiv2atomIDs:
            atomIDs = equiv2atomIDs[string]
            if atomIDs:
                atomID = min(atomIDs)
                neighs = mm.atoms[atomID].neighbor_ids
                atoms2keep = {ID for depth in neighs for ID in neighs[depth]}
                atoms2keep.add(atomID)
                bonds2keep = reduce_atoms_bonds(mm, atoms2keep)
                smallfile = '{}-atomtype={}-atomID={}.data'.format(basename, string, atomID)
                write_small_lmp(mm, atoms2keep, bonds2keep, smallfile, center=True)    
    return mm