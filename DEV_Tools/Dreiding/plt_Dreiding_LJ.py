# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 16:55:03 2023

@author: jdkem
"""


def compute_DREIDING_LJ_sigma(r0, d0):
    # Set intial r and r_inc to iterate through until E > 0
    r = 0; r_inc = 0.0000001;
    #r = 0; r_inc = 0.001; # comment/uncomment for quicker testing
    
    # Compute intial DREIDING LJ-Energy
    r += r_inc # Set intial incrementation to avoid zero errors
    rho = r/r0;
    E = d0*(rho**-12 - 2*rho**-6)

    while E > 0:
        
        # increment r by r_inc and recompute E
        r += r_inc
        rho = r/r0
        E = d0*(rho**-12 - 2*rho**-6)

    return r

def compute_DREIDING_LJ_sigma1(r0):
    # Set intial r and r_inc to iterate through
    r = 0; r_inc = 0.000001;
    
    # Compute intial r and rho
    r += r_inc # Set intial incrementation to avoid zero errors
    rho = r/r0;
    
    # Initalize rho_12 and rho_6 (rho_12 should always be larger near 0)
    rho_12 = rho**-12; rho_6 = 2*rho**-6;
    
    # While through until rho_6 gets larger
    while rho_12 > rho_6:
        r += r_inc
        rho = r/r0
        
        rho_12 = rho**-12
        rho_6 = 2*rho**-6

    return r

sigma = compute_DREIDING_LJ_sigma1(r0=3.4046)
print('sigma O_3', sigma)

# Pair Coeffs

#   1      0.095700     3.033154   # O_3
#   2      0.095100     3.472990   # C_3
#   3      0.095100     3.472990   # C_R
#   4      0.015200     2.846421   # H_

sigma = compute_DREIDING_LJ_sigma1(r0=3.4046)
print('sigma O_3', sigma)

sigma = compute_DREIDING_LJ_sigma1(r0=3.8983)
print('sigma C_3/C_R', sigma)

sigma = compute_DREIDING_LJ_sigma1(r0=3.195)
print('sigma H_', sigma)


