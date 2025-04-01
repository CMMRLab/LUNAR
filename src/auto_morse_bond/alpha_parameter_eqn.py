# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.5
November 13th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

Will find the best fit alpha paramerter for bond coeffs that will
be convertd to morse.
"""

##############################
# Import Necessary Libraries #
##############################
from math import e



# Function to find when bond is dissociated
def find_d_rcut(morse_lst, radius, D, r0):
    poverunder = 2/100 # 2% error
    Dlo = D - poverunder*D
    Dhi = D + poverunder*D
    rs = []
    for d, r in zip(morse_lst, radius):
        if Dlo < d < Dhi and r > r0: rs.append(r)
    try: rd = min(rs)
    except: rd = 2*r0
    return rd
    

############################################
# Fitting the alpha parameter              #
# https://docs.lammps.org/bond_class2.html #
# https://docs.lammps.org/bond_morse.html  #
############################################
class fitting:
    def __init__(self, m, bond_info, radius, alpha_specs, alpha_scale, ff_class, bondbreak_scale, include_rcut, log):
        self.e_harmonic = {} # {bond coeff id: potential for plotting }
        self.e_morse = {} # {bond coeff id: potential for plotting }
        self.morse_harmonic = {}  # {bond coeff id: lst[parms] }
        self.dpoint = {} # { bond coeff ID : dissociation break point }
        self.rcut = {i:0 for i in m.bond_coeffs} # { bond coeff ID : rcut }

        
        ##################################
        ### Finding the harmonic curve ###
        ##################################
        for i in m.bond_coeffs:
            # get coeffs after converting to quadratic (if desired)
            if ff_class in [1, '1']: 
                k, r0 = m.bond_coeffs[i].coeffs
            if ff_class in [2, '2']:
                r0, k2, k3, k4 = m.bond_coeffs[i].coeffs 
            self.e_harmonic[i] = [];
            self.e_morse[i] = [];
            self.morse_harmonic[i] = [];
            for r in radius:
                if ff_class in [1, '1']:
                    k_i = k*(r - r0)**2
                if ff_class in [2, '2']:
                    k2e = k2*(r - r0)**2
                    k3e = k3*(r - r0)**3
                    k4e = k4*(r - r0)**4
                    k_i = k2e + k3e + k4e
                self.e_harmonic[i].append(k_i)
                
                
        ###############################
        ### Finding the morse curve ###
        ###############################
        log.out('Fitting alpha parameters ....')
        for i in bond_info.data:
            bond = bond_info.data[i]
            morse = False
            
            ###############################
            ## Finding morse parameters ###
            ###############################
            if len(bond) == 2 and bond_info.updated[i]:
                morse = True
                if ff_class in [1, '1']:
                    k2, r0 = m.bond_coeffs[i].coeffs
                if ff_class in [2, '2']:
                    r0, k2, k3, k4 = m.bond_coeffs[i].coeffs 
                r0 = bond[0]
                D = round(bond[1]*0.239006, 4) # conver to Kcal/mol
                
                # Finding alpha via equation (use 2*K2 instead of k2 as LAMMPS parameters inlcude the 1/2 factor)
                alpha = round( (2*k2/(2*D))**(0.5), 1)
                
                # scale alpha
                alpha = round(alpha_scale*alpha, 1)
                
                # add new parms to self.morse_harmonic
                self.morse_harmonic[i] = ['morse', D, alpha, r0]
                
                # Creating the morse potential curve of the minimized alpha value
                for r in radius:                
                    k_i = D*(1 - e**(  -alpha*(r - r0)  ) )**2                         
                    self.e_morse[i].append(k_i)


            ######################################################
            ## Finding user defined over ride morse parameters ###
            ######################################################
            elif len(bond) == 3:
                morse = True
                # add new parms to self.morse_harmonic
                r0 = bond[0]
                D = round(bond[1]*0.239006, 4) # Kcal/mol
                alpha = bond[2]
                self.morse_harmonic[i] = ['morse', D, alpha, r0]
                if include_rcut:
                    rcut = float('{:.4f}'.format(bondbreak_scale*r0))
                    self.morse_harmonic[i].append(rcut)
                
                # Creating the morse potential curve of the minimized alpha value
                for r in radius:                
                    k_i = D*(1 - e**(  -alpha*(r - r0)  ) )**2                         
                    self.e_morse[i].append(k_i)
                    
                
            ####################################
            # IF morse = True add some options #
            ####################################
            if morse:
                self.dpoint[i] = find_d_rcut(self.e_morse[i], radius, D, r0)
                if bondbreak_scale > 0:
                    rcut = float('{:.4f}'.format(bondbreak_scale*r0))
                else:
                    rcut = float('{:.4f}'.format(self.dpoint[i]))
                self.rcut[i] = rcut
                if include_rcut: # New rcut and offset options
                    self.morse_harmonic[i].append(rcut)
                    if rcut < r0: log.warn( 'WARNING bond coeff {} rcut < r0; rcut = {}, r0 = {}, bondbreak_scale = {}'.format(i, rcut, r0, bondbreak_scale) )

                    
                    
            ###################################################################################################
            ### Else self morse_harmonic dict as class2 and save morse curve as harmonic curve for plotting ###
            ###################################################################################################
            else:
                if ff_class in [1, '1']: 
                    k, r0 = m.bond_coeffs[i].coeffs
                    self.morse_harmonic[i] = ['harmonic', k, r0]
                if ff_class in [2, '2']:
                    r0, k2, k3, k4 = m.bond_coeffs[i].coeffs 
                    self.morse_harmonic[i] = ['class2', r0, k2, k3, k4]
                self.e_morse[i] = self.e_harmonic[i]
