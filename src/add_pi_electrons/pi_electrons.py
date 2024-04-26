# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
Oct 3rd, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import numpy as np


# Function to update dictionary of coeffs
def update_coeff_dict(dict2update, typeID, coeffs, typestring):
    coeff = Coeff_class()
    coeff.coeffs = coeffs
    coeff.type = typestring
    dict2update.update({typeID:coeff})
    return dict2update

    
# Function to shift atom via box dim (only shift atom2)
def shift_pbc_atom(x1, x2, lx, max_x):
    diffx = x1 - x2; x22 = x2;
    if diffx > max_x:
        x22 = x2 + lx 
    elif diffx < -max_x:
        x22 = x2 - lx
    return x1, x22


# Function to create atom object
def create_atoms(x, y, z, ix, iy, iz, molid, comment, charge, Type):
    a = Atom()
    a.x = x
    a.y = y
    a.z = z
    a.ix = ix
    a.iy = iy
    a.iz = iz
    a.molid = molid
    a.comment = comment
    a.charge = charge
    a.type = Type
    return a


# Function for stringing together parameter types
def string_parameter_type(parameter_type):
    string = ''; str_buffer = 6 # Set string buffer size
    for n, i in enumerate(parameter_type):
        if n == 0: string += '{:<{str_buffer}} '.format(i, str_buffer=str_buffer)
        elif n < len(parameter_type)-1: string += ' {:^{str_buffer}} '.format(i, str_buffer=str_buffer+2)
        else: string += ' {:>{str_buffer}} {:^2}'.format(i, '', str_buffer=str_buffer)
    return string


# Function to fit plane to points
# https://stackoverflow.com/questions/35070178/fit-plane-to-a-set-of-points-in-3d-scipy-optimize-minimize-vs-scipy-linalg-lsts
def fitplane(coords):
    # barycenter of the points (centered coordinates)
    coords = np.array(coords)
    c = coords.sum(axis=0) / coords.shape[0]
    
    # run SVD
    u, s, vh = np.linalg.svd(coords - c)
    
    # unitary normal vector
    normal = vh[2, :]
    return c, normal


# Function to create pi-electrons and add info to m class    
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz .comment
class Bond: pass  # .type .atomids = (atom1id, atom2id)
class Angle: pass  # .type .atomids = (atom1id, atom2id, atom3id)
class Coeff_class: pass  # .type .coeffs = []   
def add(m, types2convert, charges, graph, log, neighbor_charge_constraint, reset_box_dims=True):
    #######################################################
    # Exit if neighbor_charge_constraint is not supported #
    #######################################################
    supported = ['check-neighbors', 'accumulate-neighbor', 'accumulate-carbon', 'accumulate-pi-electron', 'none']
    if neighbor_charge_constraint not in supported:
        log.error(f'ERROR neighbor_charge_constraint {neighbor_charge_constraint} not supported.')
    
    ##############################################################
    # Finding orginal number of atoms, bonds, angles, types, etc #
    ##############################################################
    atom_types_count = len(m.masses); atoms_count = len(m.atoms);
    bond_types_count = len(m.bond_coeffs); bonds_count = len(m.bonds);
    angle_types_count = len(m.angle_coeffs); angles_count = len(m.angles);

    
    #################################################
    # Adding new masses/pair coeffs for cg1 and cge #
    #################################################
    cg1 = atom_types_count + 1 # new TypeID for cg1 atom
    cge = atom_types_count + 2 # new TypeID for cge atom
    m.cg1_type = cg1; m.cge_type = cge # Set ***_type attributes to access for charge nuetrality
    
    # Update m.masses with cg1 and cge info
    m.masses = update_coeff_dict(m.masses, cg1, [10.011150], 'cg1') # cg1 mass
    m.masses = update_coeff_dict(m.masses, cge, [1.0000000], 'cge') # cge mass
    
    # Update m.pair_coeffs with cg1 and cge info
    m.pair_coeffs = update_coeff_dict(m.pair_coeffs, cg1, [0.06200, 3.93200],  'cg1') # cg1 pair
    m.pair_coeffs = update_coeff_dict(m.pair_coeffs, cge, [0.00001, 0.00001],  'cge') # cge pair
    
    # update m.natomtypes
    m.natomtypes = len(m.masses)
    
    
    #########################################
    # Adding new bond coeff for cg1 and cge #
    #########################################
    cg1_cge = bond_types_count + 1 # new TypeID for cg1-cge bond
    m.bond_coeffs = update_coeff_dict(m.bond_coeffs, cg1_cge, [0.65, 250.0, 0.0, 0.0], string_parameter_type(('cg1', 'cge'))) # cg1-cge bond
    
    # update m.nbondtypes
    m.nbondtypes = len(m.bond_coeffs)
    
    
    ##########################################
    # Adding new ange coeffs for cg1 and cge #
    ##########################################
    cge_cg1_cg1 = angle_types_count + 1 # new TypeID for cge-cg1-cg1 bond
    cge_cg1_cge = angle_types_count + 2 # new TypeID for cge-cg1-cge bond
    
    # Update m.angle_coeffs with cge-cg1-cg1 and cge-cg1-cge info
    m.angle_coeffs = update_coeff_dict(m.angle_coeffs, cge_cg1_cg1, [90.0, 50.0, 0.0, 0.0],  string_parameter_type(('cge', 'cg1', 'cg1'))) # cge-cge-cg1 angle
    m.angle_coeffs = update_coeff_dict(m.angle_coeffs, cge_cg1_cge, [180.0, 50.0, 0.0, 0.0], string_parameter_type(('cge', 'cg1', 'cge'))) # cge-cg1-cge angle
    
    # Update m.bondbond_coeffs with cge-cg1-cg1 and cge-cg1-cge info
    m.bondbond_coeffs = update_coeff_dict(m.bondbond_coeffs, cge_cg1_cg1, [0.0, 0.0, 0.0], string_parameter_type(('cge', 'cg1', 'cg1'))) # cge-cge-cg1 bondbond
    m.bondbond_coeffs = update_coeff_dict(m.bondbond_coeffs, cge_cg1_cge, [0.0, 0.0, 0.0], string_parameter_type(('cge', 'cg1', 'cge'))) # cge-cg1-cge bondbond
    
    # Update m.bondangle_coeffs with cge-cg1-cg1 and cge-cg1-cge info
    m.bondangle_coeffs = update_coeff_dict(m.bondangle_coeffs, cge_cg1_cg1, [0.0, 0.0, 0.0, 0.0], string_parameter_type(('cge', 'cg1', 'cg1'))) # cge-cge-cg1 bondangle
    m.bondangle_coeffs = update_coeff_dict(m.bondangle_coeffs, cge_cg1_cge, [0.0, 0.0, 0.0, 0.0], string_parameter_type(('cge', 'cg1', 'cge'))) # cge-cg1-cge bondangle

    # update m.nangletypes
    m.nangletypes = len(m.angle_coeffs)
        
        
    ##############################################################
    # Find box dimensions to remove periodic boundary conditions #
    ##############################################################
    xline = m.xbox_line.split()
    yline = m.ybox_line.split()
    zline = m.zbox_line.split()
    xlo = float(xline[0]); xhi = float(xline[1]);
    ylo = float(yline[0]); yhi = float(yline[1]);
    zlo = float(zline[0]); zhi = float(zline[1]);
    lx = xhi - xlo; ly = yhi - ylo; lz = zhi - zlo;
    
    # set max_x, max_y, max_z w/ minimum image convention
    max_x = lx/2; max_y = ly/2; max_z = lz/2;
    
    
    #########################################################################################################
    # Finding atoms to add pi electrons to and then creating new atoms/converting orginal atom types to cg1 #
    #########################################################################################################
    # Creating list to hold sub list of flagged angles that have already been found
    flagged_angles = set([])
    
    # Looping through atoms in read in data file
    atoms2add = {} # { atomID : atom object }
    nb_criteria = [2,3] # only add pi-electrons to atoms that have 2 or 3 bonded neighbors
    accetable_neightypes = types2convert + [cg1, cge]
    for id1 in m.atoms:
        atom1 = m.atoms[id1]
        
        # Get atom1 postion
        x1 = atom1.x; y1 = atom1.y; z1 = atom1.z;
        
        # Get neighbor types and check if they are all in accetable_neightypes
        if neighbor_charge_constraint == 'check-neighbors':
            neightypes = [m.atoms[i].type for i in graph[id1]]
            neightype_check = all(i in accetable_neightypes for i in neightypes)
        else: neightype_check = True
        
        # If atom.type is specified in atom_type_list, add pi electrons to
        if atom1.type in types2convert and len(graph[id1]) in nb_criteria and neightype_check:
            
            # Set charges and update for compound_neighbor_charge_constraints
            cg1_q = charges['cg1']; cge_q = charges['cge'];
            if neighbor_charge_constraint == 'accumulate-carbon':
                cg1_q = charges['cg1'] + atom1.charge
            if neighbor_charge_constraint == 'accumulate-pi-electron':
                cge_q = charges['cge'] + atom1.charge/2
            if neighbor_charge_constraint == 'accumulate-neighbor' and len(graph[id1]) == 3:
                neighs = [i for i in graph[id1] if m.atoms[i].type not in accetable_neightypes]
                for j in neighs:
                    m.atoms[j].charge = m.atoms[j].charge + m.atoms[id1].charge
                        
            # Convert atom to 'cg1' by updating charge and setting new comment
            atom1.charge = cg1_q        # update charge
            atom1.comment = 'cg1/C'     # update comment w/bond_react_merge comptable comment
            atom1.type = cg1            # update type
            
            # Loop through 1st-neighs to get xyz data (remove PBCs as well)
            xyz = []
            for id2 in graph[id1]:
                atom2 = m.atoms[id2]
                x2 = atom2.x; y2 = atom2.y; z2 = atom2.z;
                
                # Shift x, y and z-direction if periodic
                x11, x22 = shift_pbc_atom(x1, x2, lx, max_x)
                y11, y22 = shift_pbc_atom(y1, y2, ly, max_y)
                z11, z22 = shift_pbc_atom(z1, z2, lz, max_z)
                
                # Log I22-postions
                xyz.append([x22, y22, z22])
                
            # If len(xyz) <= 2 add I1-postion to xyz lst (to add more atoms to fit plane to)
            if len(xyz) <= 2: xyz.append([x1, y1, z1])
            
            # Fitting plane to 1st neighbors non-periodic positions
            c, normal = fitplane(xyz)
            
            #--------------------------------------------------#
            # Create virtual atom in positive normal direction #
            #--------------------------------------------------#
            atoms_count += 1
            pi1 = atoms_count # For finding cge-cg1-cge angle
            x1 = atom1.x + 0.65*normal[0] # cg1-cge bond length = 0.6500
            y1 = atom1.y + 0.65*normal[1] # cg1-cge bond length = 0.6500
            z1 = atom1.z + 0.65*normal[2] # cg1-cge bond length = 0.6500
            atoms2add[pi1] = create_atoms(x1, y1, z1, atom1.ix, atom1.iy, atom1.iz, atom1.molid, 'cge/C', cge_q, cge)
            
            # Creating virtual atom bonds
            bonds_count += 1
            b = Bond()
            b.type = cg1_cge
            b.atomids = (id1, pi1)
            m.bonds[bonds_count] = b
            
            # Finding new angles for cge-cg1-cg1 atom pi1
            for i in graph[id1]:  
                #             [cge, cg1, cg1]
                angle = tuple([pi1, id1, i])
                if angle not in flagged_angles:
                    angles_count += 1
                    ang = Angle()
                    ang.type = cge_cg1_cg1                        
                    ang.atomids = angle
                    m.angles[angles_count] = ang
                    flagged_angles.add(angle)
            

            #--------------------------------------------------#
            # Create virtual atom in negative normal direction #
            #--------------------------------------------------#
            atoms_count += 1
            pi2 = atoms_count # For finding cge-cg1-cge angle
            x1 = atom1.x - 0.65*normal[0] # cg1-cge bond length = 0.6500
            y1 = atom1.y - 0.65*normal[1] # cg1-cge bond length = 0.6500
            z1 = atom1.z - 0.65*normal[2] # cg1-cge bond length = 0.6500
            atoms2add[pi2] = create_atoms(x1, y1, z1, atom1.ix, atom1.iy, atom1.iz, atom1.molid, 'cge/C', cge_q, cge)
            
            # Creating virtual atom bonds
            bonds_count += 1
            b = Bond()
            b.type = cg1_cge
            b.atomids = (id1, pi2)
            m.bonds[bonds_count] = b
            
            # Finding new angles for cge-cg1-cg1 atom pi2
            for i in graph[id1]:  
                #             [cge, cg1, cg1]
                angle = tuple([pi2, id1, i])
                if angle not in flagged_angles:
                    angles_count += 1
                    ang = Angle()
                    ang.type = cge_cg1_cg1                        
                    ang.atomids = angle
                    m.angles[angles_count] = ang
                    flagged_angles.add(angle)
                  
            #------------------------------------------------------------#
            # Finding new angle for cge-cg1-cge atoms, where cge is from #
            # pos and neg direction and cg1 is from converted atom type  #
            #------------------------------------------------------------#
            angles_count += 1
            ang = Angle()
            ang.type = cge_cg1_cge
            #                   [cge, cg1, cge]
            ang.atomids = tuple([pi1, id1, pi2])
            m.angles[angles_count] = ang
            
    ####################################################
    # Print number of created atoms, bonds, and angles #
    ####################################################
    log.out('    {:>6} atoms'.format(atoms_count-m.natoms))
    log.out('    {:>6} bonds'.format(bonds_count-m.nbonds))
    log.out('    {:>6} angles'.format(angles_count-m.nangles))

            
            
    #####################################################################
    # Add atoms to m.atoms and update m.natoms, m.nbonds, and m.nangles #
    #####################################################################
    for i in atoms2add:
        m.atoms[i] = atoms2add[i]
    m.natoms = len(m.atoms)
    m.nbonds = len(m.bonds)
    m.nangles = len(m.angles)
    system_charge = [m.atoms[i].charge for i in m.atoms]
    net_charge = round(sum(system_charge)/len(system_charge), 6)
    log.out('Used neighbor_charge_constraint: {}'.format(neighbor_charge_constraint))
    log.out('Current system charge: {}'.format(net_charge))
    
    
    ##########################
    # Reset box dims if flag #
    ##########################
    # Function to adjust box (only grows box, will never let box shrink)
    def adjust_cell(box_corner, dim, add2box):
        if dim < 0: 
            newdim = dim - add2box
            if newdim > box_corner:
                newdim = dim
        else:
            newdim = dim + add2box
            if newdim < box_corner:
                newdim = dim
        return newdim
    if reset_box_dims and m.xy == 0 or m.xz == 0 or m.yz == 0:
        log.out('Resetting simulation cell dimensions after adding pi-electrons.')
        x = []; y = []; z = [];
        for i in m.atoms:
            atom = m.atoms[i]
            x.append(atom.x)
            y.append(atom.y)
            z.append(atom.z)
            
        # Find lo/hi and ln values
        xlo_pos = min(x); xhi_pos = max(x);
        ylo_pos = min(y); yhi_pos = max(y);
        zlo_pos = min(z); zhi_pos = max(z);
        
        # adjust cell directions by plus_minus
        add2box = 0.5 # add extra to each direction
        xlo = adjust_cell(xlo, xlo_pos, add2box)
        xhi = adjust_cell(xhi, xhi_pos, add2box)
        ylo = adjust_cell(ylo, ylo_pos, add2box)
        yhi = adjust_cell(yhi, yhi_pos, add2box)
        zlo = adjust_cell(zlo, zlo_pos, add2box)
        zhi = adjust_cell(zhi, zhi_pos, add2box)
    
        # Create box lines
        m.xbox_line = '{:>12.9f} {:^9.9f} {} {}'.format(xlo, xhi, 'xlo', 'xhi')
        m.ybox_line = '{:>12.9f} {:^9.9f} {} {}'.format(ylo, yhi, 'ylo', 'yhi')
        m.zbox_line = '{:>12.9f} {:^9.9f} {} {}'.format(zlo, zhi, 'zlo', 'zhi')
    return m