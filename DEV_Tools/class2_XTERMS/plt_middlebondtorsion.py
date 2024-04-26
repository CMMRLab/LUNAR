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


##############################################################################################################
# cp-cp-cp-cp  Dihedral Coeff:     8.36670000    0.00000000      1.19320000      0.00000000   0.0000  0.00   #
# cp-cp-cp-cp  middlebondtorsion:  27.59890000   -2.31200000     0.00000000      1.41700000                  #
# cp-cp Class2 Bond Coeff:         1.41700000    470.83610000   -627.61790000    1327.63450000               #
# cp-cp Morse  Bond Coeff:         29.63670000   2.90000000      1.41700000                                  #
##############################################################################################################
# bond radius to search          start,  end,  inc    decimal place round (Angstrom)
#radius = list(np.around(np.arange(0.8,   4,   0.05), 2))
radius = list(np.around(np.arange(0.8,    10,  0.1), 2))
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
r2 = 1.41700000
a1 = -0.1185
a2 = 6.3204
a3 = 0.0
def compute_middlebondtorsion(r, rot):
    p = np.array([math.radians(i) for i in rot])
    e = (r-r2)*( a1*np.cos(p) + a2*np.cos(2*p) + a3*np.cos(3*p) )
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
middlebondtorsion = np.array(compute_middlebondtorsion(np.ravel(X), np.ravel(Y))).reshape(X.shape)


# Plot surfaces
t = 0.5
ax.plot_surface(X, Y, class2bond, color='tab:orange', alpha=t, label='class2 bond')
ax.plot_surface(X, Y, morsebond, color='tab:green', alpha=t, label='morse bond')
ax.plot_surface(X, Y, middlebondtorsion, color='tab:blue', alpha=t, label='middlebondtorsion')


# Adjust axis and set label axis
ax.set_zlim3d(-1*d, 2*d)
ax.set_xlabel('Bond Radius (A)')
ax.set_ylabel('Torsion (Degrees)')
ax.set_zlabel('Energy (Kcal/mol)')



        