# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
April 25, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import numpy as np



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


################################################################################
# Function to find normal vector from atoms (assumes molecule is non-periodic) #
################################################################################
def find_normal_vector(m):
    normal = [0, 0, 1]
    coords = []
    for i in m.atoms:
        atom = m.atoms[i]
        coords.append([atom.x, atom.y, atom.z])
    
    c, normal = fitplane(coords)
    return normal 
    


#################################################################################################################
# Define a function that can find the rotational matrix that rotates the unit vector a onto the unit vector b.  #
# This function requires that unit vector a points toward unit vector b. The math was taken from:               #
# https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d #
#################################################################################################################
def rotational_matrix(a, b):
    v = np.cross(a, b) # axis of rotation
    s = np.linalg.norm(v) # sine of angle
    c = np.dot(a, b) # cosine of angle
    
    # Generate the skew-symmetric cross-product matric (and its square)
    if c != -1:
        vx = np.array([[0,    -v[2],  v[1]],
                       [v[2],     0, -v[0]],
                       [-v[1], v[0],    0]])
        vx2 = np.linalg.matrix_power(vx, 2)
        
        # Generate the final rotational matrix
        r = np.eye(3) + vx + vx2*((1 - c)/s**2)
    else:
        r = np.eye(3)
        print('WARNING can not generate rotational matrix as vector a and b point in opposite directions. Setting rotational matrix as the identity matrix')
    return r


#############################################################
# Function to rotate a molecule at centered about (0, 0, 0) #
#############################################################
def rotate_molecule_from_a_to_b(m, a, b):
    # Update simulation cell
    xline = m.xbox_line.split(); yline = m.ybox_line.split(); zline = m.zbox_line.split();
    xlo = float(xline[0]); xhi = float(xline[1]);
    ylo = float(yline[0]); yhi = float(yline[1]);
    zlo = float(zline[0]); zhi = float(zline[1]);
    
    # Find rotational matrix
    r = rotational_matrix(a, b)
    old_pos = {'x':set(), 'y':set(), 'z':set()}
    new_pos = {'x':set(), 'y':set(), 'z':set()}
    for i in m.atoms:
        atom = m.atoms[i]
        pos = np.array([atom.x, atom.y, atom.z])
        old_pos['x'].add(atom.x)
        old_pos['y'].add(atom.y)
        old_pos['z'].add(atom.z)
        
        # Rotate around atoms
        x, y, z = r @ pos

        # Save data into m.atoms
        atom.x = x
        atom.y = y
        atom.z = z
        new_pos['x'].add(atom.x)
        new_pos['y'].add(atom.y)
        new_pos['z'].add(atom.z)
        

    # Find dilo and diho buffers from orginal system
    if xlo < 0: dxlo = xlo - min(old_pos['x'])
    else: dxlo = min(old_pos['x']) - xlo
    if xhi > 0: dxhi = xhi - max(old_pos['x'])
    else: dxhi = xhi - max(old_pos['x'])
    
    if ylo < 0: dylo = ylo - min(old_pos['y'])
    else: dylo = min(old_pos['y']) - ylo
    if yhi > 0: dyhi = yhi - max(old_pos['y'])
    else: dyhi = yhi - max(old_pos['y'])
    
    if zlo < 0: dzlo = zlo - min(old_pos['z'])
    else: dzlo = min(old_pos['z']) - zlo
    if zhi > 0: dzhi = zhi - max(old_pos['z'])
    else: dzhi = zhi - max(old_pos['z'])
    
    # reset box dims  
    #print(dxlo, dxhi)
    #print(dylo, dyhi)
    #print(dzlo, dzhi)
    xlo = min(new_pos['x']) + dxlo; xhi = max(new_pos['x']) + dxhi;
    ylo = min(new_pos['y']) + dylo; yhi = max(new_pos['y']) + dyhi;
    zlo = min(new_pos['z']) + dzlo; zhi = max(new_pos['z']) + dzhi;
    m.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(xlo, xhi, 'xlo', 'xhi')
    m.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(ylo, yhi, 'ylo', 'yhi')
    m.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(zlo, zhi, 'zlo', 'zhi')
    return m

