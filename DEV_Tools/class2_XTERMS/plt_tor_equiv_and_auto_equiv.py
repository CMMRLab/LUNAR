# -*- coding: utf-8 -*-
"""
#torsion_3            cff91
> E = SUM(n=1,3) { V(n) * [ 1 - cos(n*Phi - Phi0(n)) ] }

#torsion_1            cff91_auto
> E = Kphi * [ 1 + cos(n*Phi - Phi0) ]

 2.0  2   *     cp_   s_    *         1.5000    2   180.0000
"""
##############################
# Import Necessary Libraries #
##############################
import matplotlib.pyplot as plt
import numpy as np
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


####################################################
# Function to convert from autoequiv to equivalent #
####################################################
def dihedral_autoequiv2equiv(kphi, n, phi0):
    # Set coefff list as zeros and update later on with new info
    coeff_list = 6*[0.0]
    
    # Build coeff list
    if int(n) != 0:

        # Adjust phi0 accordingly and set as phi
        if phi0 == 180: phi = 0.0; # reset phi if phi0 is 180 back to 0
        elif phi0 == 0.0: phi = 180; # reset phi if phi0 is 0 back to 180
        else: phi = phi0; # else keep phi as is
    
        # Rebuild coeffs list
        if int(n) == 1:
            coeff_list = [kphi, phi, 0.0, 0.0, 0.0, 0.0]
        elif int(n) == 2:
            coeff_list = [0.0, 0.0, kphi, phi, 0.0, 0.0]
        elif int(n) == 3:
            coeff_list = [0.0, 0.0, 0.0, 0.0, kphi, phi]
    return coeff_list


#######################################################
# Function for fitting autoequivalent dihedral coeff  #
#######################################################
def dihedral_autoequiv(kphi, n, phi0):
    # Set coefff list as zeros and update later on with new info
    coeff_list = 3*[0.0]
    
    if phi0 == 0:
        d = 1;
        coeff_list = [kphi, d, n]
    elif phi0 == 180:
        d = -1;
        coeff_list = [kphi, d, n]
    else: # For Drieding support
        d = 0;
        coeff_list = [kphi, d, n]
    return coeff_list



####################
### TOSRION TEST ###
####################
def create_and_plt_torsion(kphi, n, phi0):
    # bond radius to search      start,   end,  inc      decimal place round (Angstrom)
    phi = list(np.around(np.arange(-10,   370,  0.01), 2))
    
    # Create potentials
    auto_equiv = []; equivalent = [];
    for angle in phi:
        
        # Get kphi, d, n for auto equivalent (ae) and create potential
        # E = Kphi * [ 1 + cos(n*Phi - Phi0) ]
        kphi_ae, d_ae, n_ae = dihedral_autoequiv(kphi, n, phi0)
        auto_equiv.append( kphi_ae*(1 + d_ae*cosd(n_ae*angle)) )
        
        # v1, phi1, v2, phi2, v3, phi3 for equivalent (eq) and create potential
        # E = SUM(n=1,3) { V(n) * [ 1 - cos(n*Phi - Phi0(n)) ] }
        v1, phi1, v2, phi2, v3, phi3 = dihedral_autoequiv2equiv(kphi, n, phi0)
        n1 = 1; n2 = 2; n3 = 3;
        eq = v1*(1 - cosd(n1*angle - phi1)) + v2*(1 - cosd(n2*angle - phi2)) + v3*(1 - cosd(n3*angle - phi3))
        equivalent.append( eq )
        
    # Create new plot
    fig, ax = plt.subplots()
    
    # Plot data and set up data format
    plt.plot(phi, auto_equiv, linewidth=6.0, label='auto-equiv torsion_1 formula')
    plt.plot(phi, equivalent, linewidth=2.0, label='equivalent torsion_3 formula')
    
    # Set title and axis
    plt.title(f'Testing: kphi={kphi} n={n} phi0={phi0};   equiv={v1, phi1, v2, phi2, v3, phi3}', y=1.0, pad=30)
    plt.xlabel('Dihedral angle (deg)')
    plt.ylabel('Energy (Kcal/mol)')
    
    # Set legend
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), fancybox=True, ncol=3, fontsize=8.0)
    
    # Set figure size (width, height, forward to GUI)
    fig.set_size_inches(6, 4, forward=True)
    
    # otherwise the legend gets clipped off of saved image
    fig.tight_layout()
    return

# Testing:  2.0  2   *     cp_   s_    *         1.5000    2   180.0000
create_and_plt_torsion(kphi=1.5, n=2, phi0=180)

# Testing:  2.0  2   *     c=_   na_   *         0.0000    0     0.0000
create_and_plt_torsion(kphi=0.0, n=0, phi0=0.0)

# Testing:  2.0  2   *     c=_   si_   *         0.2110    3     0.0000
create_and_plt_torsion(kphi=0.2110, n=3, phi0=0.0)

# Testing:  2.0  2   *     c=_   n_    *         1.2500    2   180.0000
create_and_plt_torsion(kphi=1.2500, n=2, phi0=180)

    