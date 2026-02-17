# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 18:08:24 2023

@author: jdkem
%matplotlib qt  
"""

##############################
# Import Necessary Libraries #
##############################
import matplotlib.pyplot as plt
import numpy as np
from math import e
import math
plt.close('all')


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
phi = list(np.around(np.arange(-180,   180,  0.01), 2))

# Set some arbitary Dreiding torsion parms
Vjk = 5.6250 # Kcal/mol
Njk = 2      # int
phi0 = 180   # degrees

# Vjk = 0.1111 # Kcal/mol
# Njk = 3      # int
# phi0 = 0     # degrees




# Create potentials
LAMMPS_tor_harmonic = [] # LAMMPS torsion harmonic formula
Drieding_torsion = [] # from paper
for angle in phi:
    
    # Set LAMMPS d-var for harmonic (Dreiding is opposite then LAMMPS, thus 180=direction:1 and 0=directions:-1)
    d = 0.5
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
    


######################
### VDW - LJ vs X6 ###
######################
def buck_mixing(iparms, jparms):
    r0_ii, d0_ii, x_ii = iparms
    r0_jj, d0_jj, x_jj = jparms
    
    A_ii = (6*d0_ii/(x_ii - 6))*np.exp(x_ii)
    B_ii = x_ii*d0_ii/(x_ii - 6)*r0_ii**6
    C_ii = x_ii/r0_ii
    print('A_ii   = ', A_ii)
    print('B_ii   = ', B_ii)
    print('C_ii   = ', C_ii)
    
    A_jj = (6*d0_jj/(x_jj - 6))*np.exp(x_jj)
    B_jj = x_jj*d0_jj/(x_jj - 6)*r0_jj**6
    C_jj = x_jj/r0_jj
    print('A_jj   = ', A_jj)
    print('B_jj   = ', B_jj)
    print('C_jj   = ', C_jj)
    
    A_ij   =  (A_ii*A_jj)**(1/2)
    B_ij   =  (B_ii*B_jj)**(1/2)
    C_ij   = 0.5*C_ii + 0.5*C_jj
    rho_ij = 1/C_ij
    return float(A_ij), float(B_ij), float(C_ij), float(rho_ij)


# VDW for H_
d0 = 0.01520  # kcal/mol
r0 = 3.1950   # Angstroms
xi = 12.382   # unitless

# # VDW for C_R
# d0 = 0.09510  # kcal/mol
# r0 = 3.8983   # Angstroms
# xi = 14.034   # unitless



# Compute LJ stuff
r       = np.linspace(1.0, 10.0, 1000)
c       = 2**(1.0/6.0)
epsilon = d0 
sigma   = r0/c
A       = 4*epsilon*sigma**12
B       = 4*epsilon*sigma**6
LJ      = 4*epsilon*((sigma/r)**12 - (sigma/r)**6)
print('epsilon = ', epsilon)
print('sigma   = ', sigma)

# Compute DREIDNG x6
rho = r/r0
x6_dreiding = d0*( (6/(xi - 6))*np.exp(xi*(1 - rho)) - (xi/(xi - 6))*rho**-6 )


A   = (6*d0/(xi - 6))*np.exp(xi)
B   = xi*d0/(xi - 6)*r0**6
C   = xi/r0
#rho = r0/xi
print('\n\nDREIDING parameters')
print('A   = ', A)
print('B   = ', B)
print('C   = ', C)
iparms = [r0, d0, xi]
print(buck_mixing(iparms, iparms))
x6_dreiding = A*np.exp(-C*r) - B*r**-6
x6_dreiding = A*np.exp(-C*r) - B/r**6

# Convert to LAMMPS parameters
rho = 1/C # Convert using DREIDING C first
C = B     # Swap DREIDNG C for LAMMPS C

# Manual over-ride for testing params, these
# params will always be commented out unless
# testing with them (H_ from Jacob)
# A = 2473.87282957 
# rho = 0.26625000
# C = 32.33692762
x6_lammps   = A*np.exp(-r/rho) - C/r**6
print('\n\nLAMMPS parameters')
print('A_LAMMPS   =', A)
print('rho_LAMMPS =', rho)
print('C_LAMMPS   =', C)




fig, ax = plt.subplots(1, 1, figsize=(1*6, 1*4))

ax.plot(r, LJ,          lw=4, ls='-', label='LJ - 12/6')
ax.plot(r, x6_dreiding, lw=4, ls='-', label='x6 - DREIDING')
ax.plot(r, x6_lammps,   lw=4, ls='--', label='buck - LAMMPS')


ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
ax.set_xlabel('Radius (Å)')
ax.set_ylabel('Potential Energy')
ax.legend()
ax.set_ylim( (-2*epsilon, 12*epsilon) )
fig.tight_layout()