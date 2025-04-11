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
    
    https://docs.lammps.org/Howto_triclinic.html
"""

##############################
# Import Necessary Libraries #
##############################
import numpy as np



###################################################
# Get box parameters like LAMMPS and msi2lmp does #
###################################################
def get_box_parameters(m):
    # Find set simulataion cell parameters from LAMMPS datafile
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
    
    # Generate transformation matrix to convert
    # to and from scaled or Cartesian coordinates
    a = np.array([lx, 0.0, 0.0])
    b = np.array([xy, ly, 0.0])
    c = np.array([xz, yz, lz])
    h = np.vstack([a, b, c])
    h_inv = np.linalg.inv(h)
    box = [xlo, xhi,
           ylo, yhi,
           zlo, zhi,
           xy, xz, yz]
    
    return h, h_inv, box


##########################################################################################
# Function to unwrap atoms with image flags in an orthogonal or triclinc simulation cell #
##########################################################################################
def unwrap_atoms(m):    
    # Generate transformation matrix to convert
    # to and from scaled or Cartesian coordinates
    h, h_inv, box = get_box_parameters(m)
    
    # Generate graph
    graph = {i:[] for i in m.atoms}
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        graph[id1].append(id2)
        graph[id2].append(id1)
    
    # Iterative unwrapping
    checked = {ID:False for ID in m.atoms}
    for ID in graph:
        if checked[ID]:
            continue
        queue = [ID]
        while queue:
            id1 = queue.pop(0) 
            for id2 in graph[id1]:
                if checked[id2]: continue
                atom1 = m.atoms[id1]
                atom2 = m.atoms[id2]
                pos1 = np.array([atom1.x, atom1.y, atom1.z])
                pos2 = np.array([atom2.x, atom2.y, atom2.z])
                delta = pos2 - pos1
                
                # Apply minimum image convention
                frac = delta @ h_inv
                shift = np.round(frac).astype(int)
                delta_mic = (frac - shift) @ h
                unwrapped = pos1 + delta_mic

                # Update atom parameters
                atom2.x = unwrapped[0]
                atom2.y = unwrapped[1]
                atom2.z = unwrapped[2]
                atom2.ix = -shift[0]
                atom2.iy = -shift[1]
                atom2.iz = -shift[2]
                
                queue.append(id2)
                checked[id2] = True
    #m = wrap_atoms(m)
    return m


def wrap_atoms(m): 
    
    # Generate transformation matrix to convert
    # to and from scaled or Cartesian coordinates
    h, h_inv, box = get_box_parameters(m)
    
    # Iterative wrapping
    for i in m.atoms:
        atom = m.atoms[i]
        pos = np.array([atom.x, atom.y, atom.z])
        image = np.array([atom.ix, atom.iy, atom.iz])
        wrapped = pos
        
        print()
        print(image)
        

        #-----------------------------#
        # Option1: Seems to work okay #
        #-----------------------------#
        # # Reset image flags
        # pos = pos - pos1
        # frac = pos @ h_inv
        # image = np.floor(frac)
        # print(image)
                
        # Move atoms using image flags and box vectors
        # displacement = a*image[0] + b*image[1] + c*image[2]
        displacement =  image @ h
        wrapped = pos - displacement
        
        print(image)

        
        #------------------------------#
        # Option 2: Doesnt really work #
        #------------------------------#
        # # Move atoms using MIC
        # frac = pos @ h_inv
        # image_flags = np.floor(frac).astype(int)
        # frac_wrapped = frac - image_flags
        # wrapped = frac_wrapped @ h
        

        
        # Update atom parameters
        atom.x = wrapped[0]
        atom.y = wrapped[1]
        atom.z = wrapped[2]
        atom.ix = image[0]
        atom.iy = image[1]
        atom.iz = image[2]
    
    return m

    
    
