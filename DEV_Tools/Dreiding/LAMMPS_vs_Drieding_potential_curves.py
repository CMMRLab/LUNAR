# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 18:08:24 2023

@author: jdkem
"""

##############################
# Import Necessary Libraries #
##############################
import matplotlib.pyplot as plt
import numpy as np
from math import e
import math


######################
### Math functions ###
######################
# Find cos in degrees
def cosd(angle):
    return math.cos(math.radians(angle))

# Find sin in degrees
def sind(angle):
    return math.sin(math.radians(angle))


#######################
### MORSE BOND TEST ###
#######################
# bond radius to search          start,   end,  inc      decimal place round (Angstrom)
radius = list(np.around(np.arange(0.0001,   5,  0.0001), 4))

# Set some arbitary morse parms
alpha = 2.0
r0 = 1.42
De = 85

# Create potentials
LAMMPS_morse = [] # LAMMPS morse bond formula
Drieding_morse = [] # from paper
for r in radius:
    LAMMPS_morse.append(   De*(1 - e**(  -alpha*(r - r0)  ) )**2 )
    Drieding_morse.append( De*(e**(  -alpha*(r - r0)  ) - 1)**2  )


# Create new plot
fig, ax = plt.subplots()

# Plot data and set up data format
plt.plot(radius, LAMMPS_morse, linewidth=6.0, label='LAMMPS Morse bond style')
plt.plot(radius, Drieding_morse, linewidth=2.0, label='DRIEDING Morse bond formula')


plt.xlim((0, max(radius)+0.1))
plt.ylim((-De, 3*De)) 

plt.xlabel('Bond Radius (A)')
plt.ylabel('Energy (Kcal/mol)')

# Set legend
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, ncol=3, fontsize=8.0)

# Set figure size (width, height, forward to GUI)
fig.set_size_inches(6, 4, forward=True)

# otherwise the legend gets clipped off of saved image
fig.tight_layout()



####################
### TOSRION TEST ###
####################
# bond radius to search      start,   end,  inc      decimal place round (Angstrom)
phi = list(np.around(np.arange(0.0,   360,  0.01), 2))

# Set some arbitary Dreiding torsion parms
Vjk = 45.0 # Kcal/mol
Njk = 2    # int
phi0 = 180 # degrees

Vjk = 2    # Kcal/mol
Njk = 3    # int
phi0 = 180 # degrees

Vjk = 25   # Kcal/mol
Njk = 2    # int
phi0 = 0   # degrees

Vjk = 1    # Kcal/mol
Njk = 6    # int
phi0 = 0   # degrees


# Create potentials
LAMMPS_tor_harmonic = [] # LAMMPS torsion harmonic formula
Drieding_torsion = [] # from paper
for angle in phi:
    
    # Set LAMMPS d-var for harmonic (Dreiding is opposite then LAMMPS, thus 180=direction:1 and 0=directions:-1)
    if phi0 == 180: 
        d = 1 # set d as 1
    elif phi0 == 0:
        d = 1 # set d as 1 (use to be neg 1 without re-adjusting phi0)
        phi0 = 180 # Adjust phi0
    
    
    # Set LAMMPS k-var for harmonic (Multiply by 0.5)
    k = 0.5*Vjk
        
    LAMMPS_tor_harmonic.append( k*(1 + d*cosd(Njk*angle)) )
    Drieding_torsion.append( 0.5*Vjk*(1 - cosd(Njk*angle - phi0) ) )
    
# Create new plot
fig, ax = plt.subplots()

# Plot data and set up data format
plt.plot(phi, LAMMPS_tor_harmonic, linewidth=6.0, label='LAMMPS harmonic torsion formula')
plt.plot(phi, Drieding_torsion, linewidth=2.0, label='DRIEDING torsion formula')


plt.xlabel('Dihedral angle (deg)')
plt.ylabel('Energy (Kcal/mol)')

# Set legend
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, ncol=3, fontsize=8.0)

# Set figure size (width, height, forward to GUI)
fig.set_size_inches(6, 4, forward=True)

# otherwise the legend gets clipped off of saved image
fig.tight_layout()
    