# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
April 13th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

%matplotlib auto
%matplotlib qt
"""



##############################
# Import Necessary Libraries #
##############################
import matplotlib.pyplot as plt
import numpy as np
import math

# Set script preferences (this may only work for Spyder IDE, not sure yet?)
# %matplotlib auto
# %matplotlib qt


##########################################################################################################################
# cp-cp-cp-cp  Dihedral Coeff:     8.36670000    0.00000000      1.19320000      0.00000000   0.0000  0.00               #
# cp-cp-cp-cp  endbondtorsion:    -0.11850000    6.32040000      0.00000000      -0.11850000  6.3204  0.00  1.417  1.417 #
# cp-cp Class2 Bond Coeff:         1.41700000    470.83610000   -627.61790000    1327.63450000                           #
# cp-cp Morse  Bond Coeff:         29.63670000   2.90000000      1.41700000                                              #
##########################################################################################################################
# bond radius to search          start,  end,  inc    decimal place round (Angstrom)
#radius = list(np.around(np.arange(0.8,   4,   0.05), 2))
radius = list(np.around(np.arange(0.8,   8,   0.1), 2))
print('Total radius data points: ', len(radius))

# theta angle to search        start,  end,   inc   decimal place round (Degrees)
phi = list(np.around(np.arange(-180,   180,   10), 1))
print('Total angle data points: ', len(phi))


# Generate class2_bond_coeff data from radius
r0 = 1.41700000 
k2 = 470.83610000 
k3 = -627.61790000
k4 = 1327.63450000
def compute_class2_bond(r, rot):
    class2 = []
    for i in r:
        class2.append(k2*(i-r0)**2 + k3*(i-r0)**3 + k4*(i-r0)**4)
    e = np.array(class2)*np.ones(len(rot))
    e[e > 4*d] = 4*d # Will set anything >600 to 600
    return  e

    
# Generate class2_bond_coeff data from radius
r0 = 1.41700000 
a = 2.90000000 
d = 29.63670000
def compute_morse_bond(r, rot):
    morse = []
    for i in r:
        morse.append(d*( 1 - math.exp(-a*(i-r0) ) )**2)
    e = np.array(morse)*np.ones(len(rot))
    e[e > 4*d] = 4*d # Will set anything >600 to 600
    return  e


    
# Generate endbondtorsion surface data from radius and phi
r1 = 1.41700000
r3 = 1.41700000
b1 = -0.1185
b2 = 6.3204
b3 = 0.0
c1 = -0.1185
c2 = 6.3204
c3 = 0.0
def compute_endbondtorsion(r, rot):
    p = np.array([math.radians(i) for i in rot])
    e = (r-r1)*( b1*np.cos(p) + b2*np.cos(2*p) + b3*np.cos(3*p) ) + (r-r3)*( c1*np.cos(p) + c2*np.cos(2*p) + c3*np.cos(3*p) )
    #e[e > 4*d] = 4*d # Will set anything >600 to 600
    #e[e < 0] = 0
    return e


#################
# Create figure #
#################
fig = plt.figure()
ax = plt.axes(projection='3d')
X, Y = np.meshgrid(radius, phi)

# compute class2bond, morsebond, and bondangle
class2bond = np.array(compute_class2_bond(np.ravel(X), np.ravel(Y))).reshape(X.shape)
morsebond = np.array(compute_morse_bond(np.ravel(X), np.ravel(Y))).reshape(X.shape)
endbondtorsion = np.array(compute_endbondtorsion(np.ravel(X), np.ravel(Y))).reshape(X.shape)


# Plot surfaces
t = 0.5
ax.plot_surface(X, Y, class2bond, color='tab:orange', alpha=t, label='class2 bond')
ax.plot_surface(X, Y, morsebond, color='tab:green', alpha=t, label='morse bond')
ax.plot_surface(X, Y, endbondtorsion, color='tab:blue', alpha=t, label='endbondtorsion')


# Adjust axis and set label axis
ax.set_zlim3d(-1*d, 2*d)
ax.set_xlabel('Bond Radius (A)')
ax.set_ylabel('Torsion (Degrees)')
ax.set_zlabel('Energy (Kcal/mol)')



        