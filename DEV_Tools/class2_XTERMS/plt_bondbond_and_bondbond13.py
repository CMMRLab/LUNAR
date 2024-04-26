# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
April 13th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

##############################
# Import Necessary Libraries #
##############################
import matplotlib.pyplot as plt
import numpy as np
import math


##################################################################################################
# cp-cp-cp     bondbond Coeff:     68.28560000   1.41700000      1.41700000                      #
# cp-cp-cp-cp  bondbond13 Coeff:   53.00000000   1.41700000      1.41700000                      #   
# cp-cp Class2 Bond Coeff:         1.41700000    470.83610000   -627.61790000   1327.63450000    #
# cp-cp Morse  Bond Coeff:         29.63670000   2.90000000      1.41700000                      #
##################################################################################################
# bond radius to search          start,  end,  inc     decimal place round (Angstrom)
radius = list(np.around(np.arange(0.001,   4,  0.001), 3))


# Generate class2_bond_coeff data from radius
r0 = 1.41700000 
k2 = 470.83610000 
k3 = -627.61790000
k4 = 1327.63450000
class2_bond = []
for r in radius:
    e = k2*(r - r0)**2 + k3*(r - r0)**3 + k4*(r - r0)**4
    class2_bond.append(e)
    
    
# Generate class2_bond_coeff data from radius
r0 = 1.41700000 
a = 2.90000000 
d = 29.63670000
morse_bond = []
for r in radius:
    e = d*( 1 - math.exp(-a*(r-r0) ) )**2
    morse_bond.append(e)
    
# Generate bondbond_coeff data from radius
r1 = 1.41700000
r2 = 1.41700000
m = 68.28560000
bondbond = []
for r in radius:
    e = m*(r-r1)*(r-r2)
    bondbond.append(e)
    
# Generate bondbond_coeff data from radius
r1 = 1.41700000
r2 = 1.41700000
n = 53.00000000
bondbond13 = []
for r in radius:
    e = n*(r-r1)*(r-r2)
    bondbond13.append(e)
    

################
### Plotting ###
################
# Create new plot
fig, ax = plt.subplots()

# Turn on the grid
#ax.grid()

# Plot data and set up data format
lw = 2 # set line width
ax.plot(radius, class2_bond, linewidth=4*lw, color='tab:orange', label='class2 bond')
ax.plot(radius, morse_bond, linewidth=4*lw, color='tab:blue', label='morse bond')
ax.plot(radius, bondbond, linewidth=lw, color='tab:green', label='bondbond')
ax.plot(radius, bondbond13, linewidth=lw, color='tab:red', label='bondbond13')

# Adjust axis and set label axis
plt.xlim((0, max(radius)+0.1))
plt.ylim((0, 4*d))  
plt.xlabel('Bond Radius (A)')
plt.ylabel('Energy (Kcal/mol)')


# Set legend
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, ncol=4, fontsize=10)