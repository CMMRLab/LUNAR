# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
April 10th, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


Inconsistent image flags:
    The image flags for a pair on bonded atoms appear to be inconsistent. Inconsistent means
    that when the coordinates of the two atoms are unwrapped using the image flags, the two
    atoms are far apart. Specifically they are further apart than half a periodic box length.
    Or they are more than a box length apart in a non-periodic dimension. This is usually due
    to the initial data file not having correct image flags for the 2 atoms in a bond that straddles
    a periodic boundary. They should be different by 1 in that case. This is a warning because
    inconsistent image flags will not cause problems for dynamics or most LAMMPS simulations.
    However they can cause problems when such atoms are used with the fix rigid or replicate commands.
    Note that if you have an infinite periodic crystal with bonds then it is impossible to have
    fully consistent image flags, since some bonds will cross periodic boundaries and connect two
    atoms with the same image flag.
"""

##############################
# Import Necessary Libraries #
##############################
import math




###############################################
# Function to find molecules and resid molids #
###############################################  
def update_molids(m):
    # Generate graph
    graph = {i:[] for i in m.atoms}
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        graph[id1].append(id2)
        graph[id2].append(id1)
    
    # Initialize clusters and molids
    checked = {ID:False for ID in m.atoms}; clusters = set([]);
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

    # Sort clusters in a very unique fashion
    clusters = list(clusters)
    clusters = sorted(clusters, key=lambda x: x[0]) # Sort all clusters based on 1st atomID in cluster
    clusters = sorted(clusters, key=len, reverse=True) # Sort all clusters by number of atoms
    
    # Update atom molids
    for molID, cluster in enumerate(clusters, 1):
        for ID in cluster:
            m.atoms[ID].molid = molID
    return m, clusters


###################################################################################
# Function to compute distance when math.dist is not available. math.dist will be #
# quicker but for some reason not all python has math.dist .... compute_distance  #
# will be time optimized but still will be slower then math.dist.                 #
###################################################################################
def compute_distance(x1, y1, z1, x2, y2, z2):
    dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
    return math.sqrt(dx*dx + dy*dy + dz*dz)


#################################################################
# Convert triclinic 0-1 lamda coords to box coords for one atom #
#   x = H lamda + x0                                            #
#################################################################
def lamda2pos(lamda, h, boxlo):
    posx = h[0]*lamda[0] + h[5]*lamda[1] + h[4]*lamda[2] + boxlo[0]
    posy = h[1]*lamda[1] + h[3]*lamda[2] + boxlo[1]
    posz = h[2]*lamda[2] + boxlo[2]
    return [posx, posy, posz]


#################################################################
# Convert box coords to triclinic 0-1 lamda coords for one atom #
#   lamda = H^-1 (x - x0)                                       #
#################################################################
def pos2lamda(pos, h_inv, boxlo):
    dx = pos[0] - boxlo[0]
    dy = pos[1] - boxlo[1]
    dz = pos[2] - boxlo[2]
    delta = [dx, dy, dz]

    lamdax = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2]
    lamday = h_inv[1]*delta[1] + h_inv[3]*delta[2]
    lamdaz = h_inv[2]*delta[2]
    return [lamdax, lamday, lamdaz]


###################################################
# Get box parameters like LAMMPS and msi2lmp does #
###################################################
def get_box_parameters(m):
    # Find set simulataion cell parameters
    xline = m.xbox_line.split(); yline = m.ybox_line.split(); zline = m.zbox_line.split()
    xlo = float(xline[0]); xhi = float(xline[1])
    ylo = float(yline[0]); yhi = float(yline[1])
    zlo = float(zline[0]); zhi = float(zline[1])
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo
    yz = m.yz
    xz = m.xz
    xy = m.xy
    
    # Generate h and h_inv vector like LAMMPS does
    # Taken from LAMMPS domain.cpp file (lammps-23Jun2022)
    h = [lx, ly, lz, yz, xz, xy]
    boxlo = [xlo, ylo, zlo]
    boxhi = [xhi, yhi, zhi]

    # Set h-inverse box
    h_inv = 6*[0]
    h_inv[0] = 1/h[0]
    h_inv[1] = 1/h[1]
    h_inv[2] = 1/h[2]
    h_inv[3] = -h[3] / (h[1]*h[2])
    h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2])
    h_inv[5] = -h[5] / (h[0]*h[1])
    
    return h, h_inv, boxlo, boxhi


###############################################################################
# Function to update image flags in an orthogonal or triclinc simulation cell #
###############################################################################
def update_image_flags(m):
    
    # Generate h and h_inv vector like LAMMPS does
    h, h_inv, boxlo, boxhi = get_box_parameters(m)
    cx = boxhi[0] - boxlo[0]
    cy = boxhi[1] - boxlo[1]
    cz = boxhi[2] - boxlo[2]

    # Find clusters and update molIDs
    m, clusters = update_molids(m)    
    
    # Find mapping of atomIDs to bondIDs
    atoms2bondIDs = {i:set() for i in m.atoms} # {atomIDs : set(bondIDs)}
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        atoms2bondIDs[id1].add(i)
        atoms2bondIDs[id2].add(i)
        
    # Map bonds to cluster index
    clusters_bonds = []
    for cluster in clusters:
        bondIDs = set()
        for i in cluster:
            bondIDs.update(atoms2bondIDs[i])
        clusters_bonds.append(bondIDs)
    
    # Go through and update image flags in each cluster
    half_box_lambda_space = 0.499
    for cluster, bonds in zip(clusters, clusters_bonds):
        # Find lowest coordinate in each direction
        xs = []; ys = []; zs = []
        for i in cluster:
            atom = m.atoms[i]
            xs.append(atom.x)
            ys.append(atom.y)
            zs.append(atom.z)
        cx = min(xs)
        cy = min(ys)
        cz = min(zs)
        
        # Find closest atom to lowest cooridates
        dist2anchor = {} # {atomID:distance-to-anchor position}
        for i in cluster:
            atom = m.atoms[i]
            dist2anchor[i] = compute_distance(cx, cy, cz, atom.x, atom.y, atom.z)
        anchorID = min(dist2anchor, key=dist2anchor.get)
        anchor_atom = m.atoms[anchorID]
        anchor_pos = [anchor_atom.x, anchor_atom.y, anchor_atom.z]
        anchor_lamda = pos2lamda(anchor_pos, h_inv, boxlo)
        
        # Update anchoring atom image flags to be zeros
        anchor_atom.ix = 0
        anchor_atom.iy = 0
        anchor_atom.iz = 0
                
        # Find directions that are currently periodically bonded
        periodic_flags = [False, False, False]
        for i in bonds:
            id1, id2 = m.bonds[i].atomids
            atom1 = m.atoms[id1]
            atom2 = m.atoms[id2]
            pos1 = [atom1.x, atom1.y, atom1.z]
            pos2 = [atom2.x, atom2.y, atom2.z]
            lamda1 = pos2lamda(pos1, h_inv, boxlo)
            lamda2 = pos2lamda(pos2, h_inv, boxlo)
            dx = lamda1[0] - lamda2[0]
            dy = lamda1[1] - lamda2[1]
            dz = lamda1[2] - lamda2[2]
            if abs(dx) >= half_box_lambda_space: periodic_flags[0] = True
            if abs(dy) >= half_box_lambda_space: periodic_flags[1] = True
            if abs(dz) >= half_box_lambda_space: periodic_flags[2] = True 

        
        # Reset image flags based on anchoring atom position
        for i in cluster:
            if i == anchorID: continue
            atom = m.atoms[i]
            pos = [atom.x, atom.y, atom.z]
            lamda = pos2lamda(pos, h_inv, boxlo)
            
            # Check if differents lamba space from anchoringID to current ID is greater
            # then half of the simulation cell dimensions in lamda space and check to be
            # sure the cluster is periodically bonded in that direction (to avoid issues
            # with small simulation cells of single molecules).                 
            images = [0, 0, 0]
            for k in range(0, 3):
                if not periodic_flags[k]: continue
                diff_lamda = lamda[k] - anchor_lamda[k]
                image = 0
                if diff_lamda >= half_box_lambda_space:
                    image = -1
                if diff_lamda <= -half_box_lambda_space:
                    image = 1
                images[k] = image
                lamda[k] -= image

            # Update atoms image flags
            atom.ix = images[0]
            atom.iy = images[1]
            atom.iz = images[2]
    return m


##########################################################################################
# Function to unwrap atoms with image flags in an orthogonal or triclinc simulation cell #
##########################################################################################
def unwrap_atoms_with_iflags(m):
    # Update image flags
    update_image_flags(m)
    
    # Generate h and h_inv vector like LAMMPS does
    h, h_inv, boxlo, boxhi = get_box_parameters(m)
    
    for i in m.atoms:
        atom = m.atoms[i]
        pos = [atom.x, atom.y, atom.z]
        images = [atom.ix, atom.iy, atom.iz]
        lamda = pos2lamda(pos, h_inv, boxlo)
        
        # scale lamda with image flags
        lamda[0] += images[0]
        lamda[1] += images[1]
        lamda[2] += images[2]
        
        # Unwrap atoms if desired
        pos = lamda2pos(lamda, h, boxlo)
        atom.x = pos[0] 
        atom.y = pos[1] 
        atom.z = pos[2] 
        # atom.ix -= images[0]
        # atom.iy -= images[1]
        # atom.iz -= images[2]
    return m