# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.4
July 14th, 2023
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
import math


##################################################
### Functions used for fitting alpha parameter ###
##################################################
# Function to compute distance in 2D
def distance(x1, y1, x2, y2):
    dx = x1 - x2; dy = y1 - y2;
    return math.sqrt(dx*dx + dy*dy)

# Function to reduce data set
def reduce_data(xdata, ydata, xlo, xhi):
    reducedx = []; reducedy = [];  
    for x, y in zip(xdata, ydata):
        if x >= xlo and x < xhi:
            reducedx.append(x); reducedy.append(y);
    return reducedx, reducedy

# function to find closest value in list
def closest(lst, value):
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-value))]

# Function to find width of harmonic potential at dissociation energy
def width_at_D(e_harmonic, D, radius, r0):
    lo_index = 0; hi_index = 0; width = 0;
    
    # find left and right data centered around r0
    r0_lo, e_lo = reduce_data(radius, e_harmonic, xlo=min(radius), xhi=r0)
    r0_hi, e_hi = reduce_data(radius, e_harmonic, xlo=r0, xhi=max(radius))
    
    # find closest harmonic energy point to dissociation energy on left and right
    closest_e_lo = closest(e_lo, D); closest_r_lo = r0_lo[e_lo.index(closest_e_lo)];
    closest_e_hi = closest(e_hi, D); closest_r_hi = r0_hi[e_hi.index(closest_e_hi)];
    
    # Finding lo and hi index in raduis  
    lo_index = radius.index(closest_r_lo); hi_index = radius.index(closest_r_hi);
    
    # Find width at dissociation energy
    width = abs(radius[lo_index] - radius[hi_index])
    return lo_index, hi_index, width


# Function to find when bond is dissociated
def find_d_rcut(morse_lst, radius, D, r0):
    poverunder = 3/100 # 3% error
    Dlo = D - poverunder*D
    Dhi = D + poverunder*D
    rs = []
    for d, r in zip(morse_lst, radius):
        if Dlo < d < Dhi and r > r0: rs.append(r)
    try: rd = min(rs)
    except: rd = 3*r0
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
        
        # Generate alpha_testing
        alpha_round = str(alpha_specs['increment'])[::-1].find('.')
        alpha_testing = [round(n*alpha_specs['increment']+alpha_specs['start'], alpha_round) for n in range( int((alpha_specs['end']-alpha_specs['start'])/alpha_specs['increment'])+1) ]
        
        ##################################
        ### Finding the harmonic curve ###
        ##################################
        for i in m.bond_coeffs:            
            # get coeffs after converting to quadratic (if desired)
            if ff_class == 1: 
                k, r0 = m.bond_coeffs[i].coeffs
            if ff_class == 2:
                r0, k2, k3, k4 = m.bond_coeffs[i].coeffs 
            self.e_harmonic[i] = [];
            self.e_morse[i] = [];
            self.morse_harmonic[i] = [];
            for r in radius:
                if ff_class == 1:
                    k_i = k*(r - r0)**2
                if ff_class == 2:
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
                r0 = bond[0]; D = round(bond[1]*0.239006, 4) # Kcal/mol
                morse = True
                
                # Putting morse curve data into holder dictionary for analysis
                holder = {alpha: [D*(1 - e**(  -alpha*(r - r0)  ) )**2 for r in radius] for alpha in alpha_testing}
                        
                # Finding width of harmonic curve at D
                lo, hi, width = width_at_D(self.e_harmonic[i], D, radius, r0)
                    
                ##################    
                ### Minimizing ###
                ##################
                sum_distances = [];
                for j in holder:
                    # Finding minimum sum of distances
                    dist = []; alphas = holder[j];
                    for val in range(lo, hi, 1):
                        p1 = [radius[val], self.e_harmonic[i][val]]
                        p2 = [radius[val], alphas[val]]
                        dist.append(distance(p1[0], p1[1], p2[0], p2[1]))   
                    sum_distances.append(sum(dist))
                
                # Finding index of minimum sum of distances
                minimized_index = sum_distances.index(min(sum_distances))
                
                # Finding minimized alpha and add 0.1 to minimized value to be consistent with manual method
                alpha = alpha_testing[minimized_index] + 0.1
                
                # scale alpha to allow user to adjust alpha via alpha_scale
                alpha = round(alpha_scale*alpha, alpha_round)
                
                # add new parms to self.morse_harmonic
                self.morse_harmonic[i] = ['morse', D, alpha, r0]
                
                # Creating the morse potential curve of the minimized alpha value
                for r in radius:                
                    k_i = D*(1 - e**(  -alpha*(r - r0)  ) )**2                         
                    self.e_morse[i].append(k_i)
                    
                # Warn if found alpha parameter is the same as the max or min aplha's in aplha_test
                if alpha <= min(alpha_testing):
                    log.out(f"\n\nBond Coeff {i} has found an alpha parameter near the minimum in alpha_testing.")
                    log.out(f"Please lower the start value in alpha_testing to something less than {alpha} and run the code once more.")
                elif alpha >= max(alpha_testing):
                    log.out(f"\n\nBond Coeff {i} has found an alpha parameter near the maximum in alpha_testing.")
                    log.out(f"Please raise the end value in alpha_testing to something greater than {alpha+0.1} and run the code once more.")


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
                
                # Creating the morse potential curve of the minimized alpha value
                for r in radius:                
                    k_i = D*(1 - e**(  -alpha*(r - r0)  ) )**2                         
                    self.e_morse[i].append(k_i)
                    
            ####################################
            # IF morse = True add some options #
            ####################################
            if morse:
                self.dpoint[i] = find_d_rcut(self.e_morse[i], radius, D, r0)
                if include_rcut: # New rcut and offset options
                    rcut = float('{:.4f}'.format(bondbreak_scale*r0))
                    self.morse_harmonic[i].append(rcut)
                    
                    
            ###############################################################################################################
            ### Else self morse_harmonic dict as harmonic or class2 and save morse curve as harmonic curve for plotting ###
            ###############################################################################################################
            else:
                if ff_class == 1: 
                    k, r0 = m.bond_coeffs[i].coeffs
                    self.morse_harmonic[i] = ['harmonic', k, r0]
                if ff_class == 2:
                    r0, k2, k3, k4 = m.bond_coeffs[i].coeffs 
                    self.morse_harmonic[i] = ['class2', r0, k2, k3, k4]
                self.e_morse[i] = self.e_harmonic[i]