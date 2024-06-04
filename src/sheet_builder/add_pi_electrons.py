# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
May 29th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import numpy as np


######################################
# Function to shift atom via box dim #
######################################
def shift_pbc_atom(x1, x2, lx, max_x):
    diffx = x1 - x2; x11 = x1; x22 = x2;
    if diffx > max_x:
        x11 = x1 + lx
        x22 = x2 + lx 
    elif diffx < -max_x:
        x11 = x1 - lx
        x22 = x2 - lx
    return x11, x22


###################################
# Function to fit plane to points #
###################################
def fitplane(coords):
    # barycenter of the points (centered coordinates)
    coords = np.array(coords)
    c = coords.sum(axis=0) / coords.shape[0]
    
    # run SVD
    u, s, vh = np.linalg.svd(coords - c)
    
    # unitary normal vector
    normal = vh[2, :]
    return c, normal


###################################
# Function to create atoms object #
###################################
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz .comment .atomtype
def create_atoms(molid, typeint, atomtype, x, y, z, ix, iy, iz):
    a = Atom()
    a.molid = molid
    a.type = typeint
    a.atomtype = atomtype
    a.comment = atomtype
    a.charge = 0.0
    a.x = x
    a.y = y
    a.z = z
    a.ix = ix
    a.iy = iy
    a.iz = iz
    return a


################################
# Function to add pi-electrons #
################################
def add(atoms, bonds, box, pi_electrons, log):
    # Find orginal number of atoms and bonds
    natoms = len(atoms); nbonds = len(bonds);
    
    # Find simulation cell size
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']
    
    # set max_x, max_y, max_z w/ minimum image convention
    max_x = lx/2; max_y = ly/2; max_z = lz/2;
    
    # Generate graph
    graph = {i:[] for i in atoms}
    for id1, id2 in bonds:     
        graph[id1].append(id2)
        graph[id2].append(id1)
        
    # Start adding pi-electrons
    atoms_count = len(atoms); bonds_count = len(bonds); 
    bond_length = 0.65 # cg1-cge bond length = 0.6500
    for id1 in graph:
        atom1 = atoms[id1]
        electron = pi_electrons[atom1.type]
        
        # Get atom1 postion
        x1 = atom1.x; y1 = atom1.y; z1 = atom1.z;
        
        # Only add pi-electron if it is not an empty string and the neighor count of the atoms is 2 or 3
        if electron != '' and len(graph[id1]) in [2, 3]:
            
            # Loop through 1st-neighs to get xyz data (remove PBCs as well)
            xyz = []
            for id2 in graph[id1]:
                atom2 = atoms[id2]
            
                # Get atom2 postion
                x2 = atom2.x; y2 = atom2.y; z2 = atom2.z;
                
                # Shift x, y and z-direction if periodic
                x11, x22 = shift_pbc_atom(x1, x2, lx, max_x)
                y11, y22 = shift_pbc_atom(y1, y2, ly, max_y)
                z11, z22 = shift_pbc_atom(z1, z2, lz, max_z)
                
                # Log I22-postions
                xyz.append([x22, y22, z22])
                
            # If len(xyz) <= 2 add I11-postion to xyz lst (to add more atoms to fit plane to)
            if len(xyz) <= 2: xyz.append([x11, y11, z11])
            
            # Fitting plane to 1st neighbors non-periodic positions
            c, normal = fitplane(xyz)
            
            
            # Create virtual atom in positive normal direction
            atoms_count += 1
            pi_x1 = atom1.x + bond_length*normal[0] 
            pi_y1 = atom1.y + bond_length*normal[1]
            pi_z1 = atom1.z + bond_length*normal[2]
            molid = atom1.molid
            typeint = atom1.type + 4
            atoms[atoms_count] = create_atoms(molid, typeint, electron, pi_x1, pi_y1, pi_z1, atom1.ix, atom1.iy, atom1.iz)
            bonds_count += 1
            bonds.append(tuple(sorted([id1, atoms_count])))
            
            # Create virtual atom in negative normal direction
            atoms_count += 1
            pi_x1 = atom1.x - bond_length*normal[0] 
            pi_y1 = atom1.y - bond_length*normal[1]
            pi_z1 = atom1.z - bond_length*normal[2]
            molid = atom1.molid
            typeint = atom1.type + 4
            atoms[atoms_count] = create_atoms(molid, typeint, electron, pi_x1, pi_y1, pi_z1, atom1.ix, atom1.iy, atom1.iz)
            bonds_count += 1
            bonds.append(tuple(sorted([id1, atoms_count])))
            
    # Print number of created atoms and bonds
    log.out('  {:>6} atoms'.format(atoms_count - natoms))
    log.out('  {:>6} bonds'.format(bonds_count - nbonds))
    return atoms, bonds