# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 18th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

Atom typing rules for CVFF-IFF with user-defined assumed types
found in the assumed dictionary. These are stored inside each
force field atom-typing script because not all forcefields 
need to have the exact same assumed atom types and some
force fields may not even contain basic elements that others
may.
"""


##############################
# Import Necessary Libraries #
##############################
import src.atom_typing.typing.typing_functions as tf
from collections import OrderedDict
import os


######################################################################
# Function to perform atom-typing for CVFF-IFF forcefield atom types #
######################################################################
def nta(mm, basename, ff_name):
    tally = {'found':0, 'assumed':0, 'failed':0} # To tally findings
    
    ################################################
    # Set CVFF-IFF specific informations and flags #
    ################################################
    # Log supported types
    supported_types = {'Carbon':      ['ct', 'c+', 'cr', 'c-', 'c*', " c' ", 'c5', 'cs', 'cp', 'ci', 'c=',
                                       'c=1', 'c=2', 'c3h', 'c3m', 'c4h', 'c4m', 'ca', 'cg', 'co', 'coh',
                                       'cn', 'c1', 'c2', 'c3', 'c',
                                       'ce1 (Q)'], # NEW CVFF-IFF additions
                       
                       'Hydrogen':    ['hi', 'hc', 'hw', 'ho', 'hn', 'hs', 'hp',
                                       'he1 (Q)', 'ha1 (Q)'], # NEW CVFF-IFF additions
                       
                       'Oxygen':      [" o' ", 'o-', 'oz', 'o*', 'oe', 'oc', 'o', 'o3e', 'o4e', 'op',
                                       'oh', 'of',
                                       'oe1 (Q)', 'oa1 (Q)'], # NEW CVFF-IFF additions
                       
                       'Nitrogen':    ['nt', 'n=', 'n=1', 'n=2', 'n1', 'n2', 'ni', 'nho', 'npc', 'nh',
                                       'nh+', 'n3n', 'n3m', 'n4n', 'n4m', 'n', 'nn', 'na', 'no', 'n3',
                                       'nh+', 'n4'],
                       
                       'Sulfer':      ['s-', 's3e', 's4e', 'sp', " s' ", 'sc', 'sh', 's1'],
                        
                       'Silicone':    ['sz', 'si'],
                       
                       'Fluorine':    ['f'],
                       'Helium':      ['he'],
                       'Deuterium':   ['dw', 'd'],
                       'Chlorine':    ['cl'],
                       'Calcium':     ['ca+'],
                       'Bromine':     ['br', 'Br'],
                       'Argon':       ['ar'],
                       'Phosphorous': ['pz', 'p']
                       }
    
    # Set lst of heavy elements (NOTE: MUST BE UPDATED IF MORE CVFF-IFF ATOM TYPES GET CODED IN)
    heavies = ['C', 'O', 'N', 'S', 'F', 'Si', 'Xe', 'Ne', 'Kr', 'Cl', 'Br', 'Ar', 'P']
    
    # CVFF-IFF flags to help set correct atom-types (True to use False not to use)
    use_peo_poly_types = False # Option to use ce1/oe1/he1/ha1 atom-types else c2/oc/hc/ho2 will be used (charge will also be reset)
    
    # Update supported types lists if flags are used and set flag status in parentheses
    supported_types = tf.update_supported_types(supported_types, element='Carbon', atomtype='ce1 (Q)', flag=use_peo_poly_types)
    supported_types = tf.update_supported_types(supported_types, element='Hydrogen', atomtype='he1 (Q)', flag=use_peo_poly_types)
    supported_types = tf.update_supported_types(supported_types, element='Hydrogen', atomtype='ha1 (Q)', flag=use_peo_poly_types)
    supported_types = tf.update_supported_types(supported_types, element='Oxygen', atomtype='oe1 (Q)', flag=use_peo_poly_types)
    supported_types = tf.update_supported_types(supported_types, element='Oxygen', atomtype='oa1 (Q)', flag=use_peo_poly_types)

    
    
    ##########################################################################
    # Set Assumed CVFF-IFF atomtypes if code fails to find correct atom-type #
    ##########################################################################
    assumed = {'Carbon-2-connect-no-ring':    'cp',     # 2-connects or less: Sp1 carbon no ring
               'Carbon-2-connect-in-ring':    'ct',     # 2-connects or less: Sp1 carbon in ring
               'Carbon-3-connect-no-ring':    'c=2',    # 3-connects: Sp2 carbon no ring
               'Carbon-3-connect-in-ring':    'cp',     # 3-connects: Sp2 carbon in ring
               'Carbon-4-connect-no-ring':    'c2',     # 4-connects: Sp3 carbon no ring
               'Carbon-4-connect-in-ring':    'c2',     # 4-connects: Sp3 carbon in ring
               'Hydrogen-1-connect':          'hc',     # 1-connect : Hydrogen
               'Oxygen-1-connect':            " o' ",   # 1-connect : Oxygen
               'Oxygen-2-connect-no-ring':    'o',      # 2-connects: Oxygen no ring
               'Oxygen-2-connect-in-ring':    'op',     # 2-connects: Oxygen in ring
               'Nitrogen-1-connect':          'nt',     # 1-connect : Nitrogen
               'Nitrogen-2-connect-no-ring':  'n=',     # 2-connects: Nitrogen no ring
               'Nitrogen-2-connect-in-ring':  'np',     # 2-connects: Nitrogen in ring
               'Nitrogen-3-connect-no-ring':  'n3',     # 3-connects: Nitrogen no ring
               'Nitrogen-3-connect-in-ring':  'nh',     # 3-connects: Nitrogen in ring
               'Nitrogen-4-connect-no-ring':  'n4',     # 4-connects: Nitrogen no ring
               'Nitrogen-4-connect-in-ring':  'nh+',    # 4-connects: Nitrogen in ring
               'Sulfur-1-connect':            's-',     # 1-connect : Sulfur
               'Sulfur-2-connect-no-ring':    's',      # 2-connects: Sulfur no ring
               'Sulfur-2-connect-in-ring':    'sp',     # 2-connects: Sulfur in ring
              }


    ########################################################
    # Opening files to dump assumed/failed atom types into #
    ########################################################
    a = open(basename+'_assumed.txt', 'w')
    f = open(basename+'_failed.txt', 'w')

    ########################################################    
    # Loop through mm.atoms and start assigning atom-types #
    ########################################################
    for i in mm.atoms:
        # Set atom instance and pull out desired info
        atom = mm.atoms[i]; formula = atom.molecule.formula;
        nb = int(atom.nb); ring = int(atom.ring); element = atom.element;
        
        # Find 1st-neighbor info (denoted by nameN; where N=neighbor depth)
        elements1 = tf.neigh_extract(atom, depth=1, info='element') # Example: ['C', 'C', 'N']
        rings1 = tf.neigh_extract(atom, depth=1, info='ring') # Example: [6, 6, 0]
        nbs1 = tf.neigh_extract(atom, depth=1, info='nb') # Example: [3, 3, 3]
        
        # Find 2nd-neighbor info (denoted by nameN; where N=neighbor depth)
        elements2 = tf.neigh_extract(atom, depth=2, info='element') # Example: ['C', 'C', 'N']
        #rings2 = tf.neigh_extract(atom, depth=2, info='ring') # Example: [6, 6, 0]
        nbs2 = tf.neigh_extract(atom, depth=2, info='nb') # Example: [3, 3, 3]
        
        # Find 3rd-neighbor info (denoted by nameN; where N=neighbor depth)
        #rings3 = tf.neigh_extract(atom, depth=3, info='ring') # Example: [3, 3, 3]
        
        
        # Set intial .nta and .nta_comments and update later on if found
        atom.nta_type = '{}-type-yourself'.format(element)
        atom.nta_info = 'FAILED TO BE TYPED:  element: {}, ring: {}, nb: {}'.format(element, ring, nb)

        
        
        ######################################################################################
        # Sp1 Carbon atom typing (ordering of nested if/elif/else statements set precedence) #
        ######################################################################################
        if element == 'C' and nb <= 2:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#  
            # For 5 member-ring open valence
            if ring == 5:
                atom.nta_type = 'c5'; tally['found'] += 1;
                atom.nta_info = 'User-defined for Open Valence ReaxFF polymerization (element:C, ring:5, nb:2)'

            # For 6 member-ring open valence
            elif ring >= 6:
                atom.nta_type = 'cp'; tally['found'] += 1;
                atom.nta_info = 'User-defined for Open Valence ReaxFF polymerization (element:C, ring:>=6, nb:2)'
                
            #----------------------------------------------------------------------------#
            # Strict IFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------#
            # ct      12.01115      C          2        sp carbon involved in a triple bond
            elif ring == 0:
                tally['found'] += 1
                atom.nta_type = 'ct'
                atom.nta_info = 'Correctly found'
                
            #-----------------------------------------------------------------------------#
            # If user defined attempt fails and Strict IFF atom-typing fails used assumed #
            #-----------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Carbon-2-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Carbon-2-connect-no-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                    
        ######################################################################################
        # Sp2 Carbon atom typing (ordering of nested if/elif/else statements set precedence) #
        ######################################################################################
        elif element == 'C' and nb == 3:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # cp       12.01115      C          3        sp2 aromatic carbon
            # For 7 and 8 member rings on carbonized structures
            if ring > 6:
                atom.nta_type = 'cp'; tally['found'] += 1;
                atom.nta_info = 'User-defined for 7/8-ring carbonization (element:C, ring:>6, nb:3)'
                
            #----------------------------------------------------------------------------#
            # Strict IFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------#
            # c+       12.01115      C          3        C in guanidinium group
            elif ring == 0 and elements1.count('N') == 3 and formula == 'C1-H5-N3':              
                atom.nta_type = 'c+'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # cr        12.01115      C          3        C in neutral arginine 
            elif ring == 0 and elements1.count('C') == 2 and elements1.count('N') == 1:
                atom.nta_type = 'cr'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            elif ring == 0 and elements1.count('N') == 3:
                atom.nta_type = 'cr'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c*    12.011150     C           3        Carbon in carbonyl  group,   non_amides
            # c"    12.011150     C           3        Carbon in carbonyl  group,   non_amides
            # (O-index=2 in nbs1)
            elif ring == 0 and elements1.count('C') == 2 and elements1.count('O') == 1 and nbs1[2] == 0:
                atom.nta_type = " c* "; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c'     12.011150     C            3        Sp2 carbon in carbonyl (C=O) group  
            # (O-index=2 in nbs1)
            elif ring == 0 and elements1.count('C') >= 1 and elements1.count('O') == 1 and nbs1[2] == 0:
                atom.nta_type = " c' "; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c-        12.01115      C          3        C in charged carboxylate  
            # (N atom would be at index 1 and 2 in rings1 so check that rings[1 and 2] are zero)
            elif ring == 0 and elements1.count('C') == 1 and elements1.count('N') == 2 and rings1[1] == 0  and rings1[2] == 0:
                atom.nta_type = 'c-'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c5       12.01115      C          3        sp2 aromatic carbon in 5-membered ring  
            elif ring == 5 and elements1.count('S') == 0:
                atom.nta_type = 'c5'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # cs        12.01115      C          3        sp2 aromatic carbon in 5 membered ring next to S 
            elif ring == 5 and elements1.count('S') > 0:
                atom.nta_type = 'cs'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # cp       12.01115      C          3        sp2 aromatic carbon
            elif ring == 6:
                atom.nta_type = 'cp'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # ci        12.01115      C          3        sp2 aromatic carbon in charged imidazole ring (His+)
            # (N-index=2 for ring and element)
            elif ring >= 5 and rings1[2] >= 5 and elements1[2] =='N' and formula == 'C3-H4-N2':
                atom.nta_type = 'ci'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c=       12.01115      C          3        nonaromatic end doubly bonded carbon
            elif ring == 0 and nbs1.count(1) == 2:
                atom.nta_type = 'c='; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                                
            # c=1     12.01115      C          3        nonaromatic, next to end doubly bonded carbon 
            elif ring == 0 and nbs2.count(1) >= 2:
                atom.nta_type = 'c=1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c=2     12.01115      C          3        nonaromatic doubly bonded carbon  
            elif ring == 0:
                atom.nta_type = 'c=2'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict CVFF-IFF atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Carbon-3-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Carbon-3-connect-no-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                    
                    
        ######################################################################################
        # Sp3 Carbon atom typing (ordering of nested if/elif/else statements set precedence) #
        ######################################################################################
        elif element == 'C' and nb == 4:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #----------------------------------------------------------------------------#
            # Strict IFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------#
            # ce1     12.01115       C          4        Carbon atom in backbone or terminal group in PEO
            # Will check for topology in backbone and terminal unit for element carbon
            if ring == 0  and tf.check_peo_topo(atom, element_type='C') and use_peo_poly_types:              
                atom.nta_type = 'ce1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                atom.charge = 0.094000 # reset C-element charge for IFF if PEO atom-type
                
            # c3h      12.01115      C          4        sp3 carbon in 3-membered ring with hydrogens
            elif ring == 3 and 'H' in elements1:              
                atom.nta_type = 'c3h'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c3m     12.01115      C          4        sp3 carbon in 3-membered ring
            elif ring == 3:              
                atom.nta_type = 'c3m'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c4h      12.01115      C          4        sp3 carbon in 4-membered ring with hydrogens
            elif ring == 4 and 'H' in elements1:              
                atom.nta_type = 'c4h'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'

            # c4m     12.01115      C          4        sp3 carbon in 4-membered ring
            elif ring == 4:              
                atom.nta_type = 'c4m'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # ca      12.011150    C           4        General amino acid alpha carbon (sp3)
            elif ring == 0 and elements1.count('C') == 2 and elements1.count('H') == 1 and elements1.count('N') == 1:              
                atom.nta_type = 'ca'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # cg        12.01115      C          4        sp3 alpha carbon in glycine
            elif ring == 0 and formula == 'C2-H5-N1-O2':              
                atom.nta_type = 'cg'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # co        12.01115      C          4        sp3 carbon in acetals
            elif ring == 0 and elements1.count('C') == 2 and elements1.count('O') == 2:              
                atom.nta_type = 'co'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # coh     12.01115      C          4        sp3 carbon in acetals with hydrogen
            elif ring == 0 and elements1.count('C') == 1 and elements1.count('H') == 1 and elements1.count('O') == 2:              
                atom.nta_type = 'coh'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # cn      12.011150    C           4        Sp3 Carbon bonded to N 
            elif elements1.count('N') >= 1:              
                atom.nta_type = 'cn'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                

            #############################################################################
            # attempt more unqiue types 1st before appling general c1, c2, and c3 types #
            #############################################################################
            # c1        12.01115      C          4        sp3 carbon with 1 H 3 heavies
            elif elements1.count('H') == 1 and tf.count_heavies(elements1, heavies) == 3:              
                atom.nta_type = 'c1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c2        12.01115      C          4        sp3 carbon with 2 H's, 2 Heavy’s
            elif elements1.count('H') == 2 and tf.count_heavies(elements1, heavies) == 2:              
                atom.nta_type = 'c2'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c3        12.01115      C          4        sp3 carbon with 3 Hs 1 heavy
            elif elements1.count('H') == 3 and tf.count_heavies(elements1, heavies) == 1:              
                atom.nta_type = 'c3'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # General Sp3 C atom types
            # c         12.01115      C          4        generic SP3 carbon
            elif element == 'C':
                atom.nta_type = 'c'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict CVFF-IFF atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Carbon-4-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Carbon-4-connect-no-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                    
                    
        ####################################################################################
        # Hydrogen atom typing (ordering of nested if/elif/else statements set precedence) #
        ####################################################################################
        elif element == 'H' and nb == 1:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #---------------------------------------------------------------------------------#
            # Strict CVFF-IFF atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # he1     1.00797      H         1        Hydrogen atom in backbone or terminal group in PEO
            # Will check for topology in backbone and terminal unit for element hydrogen
            if ring == 0  and tf.check_peo_topo(atom, element_type='H') and use_peo_poly_types:              
                atom.nta_type = 'he1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                atom.charge = 0.053000 # reset H-element charge for IFF if PEO atom-type
                
            # ha1     1.00797      H         1        Hydrogen atom in terminal OH group in PEO (+0.4)
            # Will check for topology in terminal -OH unit for element hydrogen
            elif ring == 0  and tf.check_peo_topo(atom, element_type='H') and use_peo_poly_types:              
                atom.nta_type = 'ha1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                atom.charge = 0.400000 # reset H-element charge for IFF if PEO atom-type
                
            # hi       1.00797      H          1        Hydrogen in charged imidazole ring
            elif rings1.count(5) > 0 and tf.count_neigh(atom.neighbor_info[2], element='N', ring=5, nb=False) > 0 and formula == 'C3-H4-N2':
                atom.nta_type = 'hi'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # hc      1.00797      H          1        hydrogen bonded to carbon
            elif elements1[0] == 'C':
                atom.nta_type = 'hc'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # hw       1.00797     H         1        hydrogen in water (+0.41 in this model !)
            # h*       1.007970    H         1        Hydrogen in water molecule
            # hspc     1.007970    H         1        Hydrogen in SPC water molecule   
            # htip     1.007970    H         1        Hydrogen in TIP water molecule    
            elif elements1[0] == 'O' and formula == 'H2-O1':
                atom.nta_type = 'hw'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # ho      1.00797      H          1        hydrogen bonded to oxygen
            elif elements1[0] == 'O':
                atom.nta_type = 'ho'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # hn      1.00797      H          1        hydrogen bonded to nitrogen
            elif elements1[0] == 'N':
                atom.nta_type = 'hn'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # hs       1.00797      H         1        hydrogen bonded to sulfur
            elif elements1[0] == 'S':
                atom.nta_type = 'hs'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # hp       1.007970    H           1        Hydrogen bonded to P  
            elif elements1[0] == 'P':
                atom.nta_type = 'P'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
        ###########################################################################################
        # 1-bonded Oxygen atom typing (ordering of nested if/elif/else statements set precedence) #
        ###########################################################################################
        elif element == 'O' and nb == 1:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #----------------------------------------------------------------------------------#
            # Strict CVFF-IFFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------------#
            # o'    15.999400     O           1        Oxygen in carbonyl (C=O) group
            if ring == 0 and elements1.count('C') == 1 and 'C' in elements2 and 'O' not in elements2:              
                atom.nta_type = " o' "; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o-    15.999400    O           1        Oxygen in charged carboxylate (COO-) group
            elif ring == 0 and elements1.count('C') == 1 and 'O' in elements2 and 'C' not in elements2:              
                atom.nta_type = 'o-'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # oz    15.999400    O           1        Oxygen in Zeolite
            elif ring == 0 and elements1.count('Si') >= 1:              
                atom.nta_type = 'oz'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict CVFF-IFF atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                tag = 'Oxygen-1-connect'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                

        ###########################################################################################
        # 2-bonded Oxygen atom typing (ordering of nested if/elif/else statements set precedence) #
        ###########################################################################################
        elif element == 'O' and nb == 2:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #---------------------------------------------------------------------------------#
            # Strict CVFF-IFF atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # o*     15.999400    O           2        Oxygen in water molecule 
            # ospc   15.999400    O           2        Oxygen in SPC water molecule     
            # otip   15.999400    O           2        Oxygen in TIP3P water molecule 
            if formula == 'H2-O1':
                atom.nta_type = 'o*'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # oe1     15.99940      O          2        Oxygen atom in PEO backbone
            # oet5    14.99940      O          2        Oxygen atom in PEO backbone with lone pairs et5 attached (et5=-0.12)
            # Will check for topology in backbone and terminal unit for element oxygen
            elif ring == 0  and tf.check_peo_topo(atom, element_type='O') and use_peo_poly_types:              
                atom.nta_type = 'oe1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # oa1      15.99940      O          2       Oxygen atom in terminal OH group in PEO (-0.6)
            # oat5     14.999400     O          2       Oxygen atom in terminal OH group in PEO with lone pairs et5 attached (et5=-0.18)
            # Will check for topology in backbone and terminal unit for element oxygen
            elif ring == 0  and tf.check_peo_topo(atom, element_type='O') and use_peo_poly_types:              
                atom.nta_type = 'oa1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                atom.charge = -0.600000 # reset O-element charge for IFF if PEO atom-type
                
            # oe     15.999400    O           2        Sp3 oxygen in ester
            elif ring == 0  and len(rings1) == rings1.count(0) and elements1.count('C') == 2 and tf.count_neigh(atom.neighbor_info[2], element='O', ring=0, nb=1) == 1 and elements2.count('C') == 2:              
                atom.nta_type = 'oe'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # oc     15.999400     O           2        Sp3 oxygen in ether or acetals
            elif ring == 0 and elements1.count('C') == 2 and elements2.count('O') == 1:              
                atom.nta_type = 'oc'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o        15.999400    O           2        Sp3 oxygen in ether or ester groups
            elif ring == 0 and elements1.count('C') == 2:              
                atom.nta_type = 'o'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o3e     15.99940      O          2        sp3 oxygen  in three membered ring
            elif ring == 3:              
                atom.nta_type = 'o3e'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o4e     15.99940      O          2        sp3 oxygen  in  four  membered ring
            elif ring == 4:              
                atom.nta_type = 'o4e'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # op       15.99940      O          2        sp2 aromatic in 5 membered ring
            elif ring >= 5:              
                atom.nta_type = 'op'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # oh      15.999400    O           2        Oxygen in hydroxyl (OH) group
            elif ring == 0 and elements1.count('C') == 1 and elements1.count('H') == 1:              
                atom.nta_type = 'oh'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # of      15.999400    O           2        Oxygen 
            elif element == 'O':              
                atom.nta_type = 'of'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict CVFF-IFF atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Oxygen-2-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Oxygen-2-connect-no-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                    
                    
        #############################################################################################
        # 1-bonded Nitrogen atom typing (ordering of nested if/elif/else statements set precedence) #
        #############################################################################################
        elif element == 'N' and nb == 1:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #---------------------------------------------------------------------------------#
            # Strict CVFF-IFF atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # nt      14.00670      N          1        sp nitrogen involved in a triple bond 
            # nz      14.00670      N          1        sp3 nitrogen bonded to two atoms
            if element == 'N':
                atom.nta_type = 'nt'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict CVFF-IFF atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                tag = 'Nitrogen-1-connect'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                
                
        #############################################################################################
        # 2-bonded Nitrogen atom typing (ordering of nested if/elif/else statements set precedence) #
        #############################################################################################
        elif element == 'N' and nb == 2:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #----------------------------------------------------------------------------#
            # Strict IFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------#                
            # np       14.00670      N          2        sp2 nitrogen in 5- or 6- membered ring
            if ring >= 5:
                atom.nta_type = 'np'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n=       14.00670      N          2        non aromatic end doubly bonded nitrogen
            elif ring == 0 and nbs1.count(1) >= 1:
                atom.nta_type = 'n='; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n=1     14.00670      N          2        non aromatic, next to end doubly bonded carbon
            elif ring == 0 and nbs2.count(1) >= 1:
                atom.nta_type = 'n=1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n=2     14.00670      N          2        non aromatic doubly bonded nitrogen  
            elif ring == 0:
                atom.nta_type = 'n=2'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #-----------------------------------------------------------------------------#
            # If user defined attempt fails and Strict IFF atom-typing fails used assumed #
            #-----------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Nitrogen-2-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Nitrogen-2-connect-no-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                    
        #############################################################################################
        # 3-bonded Nitrogen atom typing (ordering of nested if/elif/else statements set precedence) #
        #############################################################################################
        elif element == 'N' and nb == 3:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #---------------------------------------------------------------------------------#
            # Strict CVFF-IFF atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # n1        14.00670      N          3        sp2 nitrogen in charged arginine
            if ring == 0 and formula == 'C6-H14-N4-O2':
                atom.nta_type = 'n1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n2        14.00670      N          3        sp2 nitrogen (NH2) in guanidinium group (HN=C(NH2)2)
            # nr        14.00670      N          3        sp2 nitrogen (NH2) in guanidinium group (HN=C(NH2)2)
            elif ring == 0 and formula == 'C1-H5-N3':
                atom.nta_type = 'n1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                  
            # n         14.006700    N          3        Sp2 nitrogen with 1 H, 2 heavy atoms (amide group)
            elif ring == 0 and elements1.count('H') == 1 and elements1.count('C') == 2:
                atom.nta_type = 'n'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # nn       14.00670       N          3        sp2 nitrogen in aromatic amines
            # nb       14.00670       N          3        sp2 nitrogen in aromatic amines
            elif len([i for i in rings1 if i > 0]) > 0 and 3 in nbs1:
                atom.nta_type = 'nn'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # na        14.00670      N          3        sp3 nitrogen in amines
            elif ring == 0 and elements1.count('H') == 2 and elements1.count('C') == 1:
                atom.nta_type = 'na'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            elif ring == 0 and elements1.count('C') == 3:
                atom.nta_type = 'na'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                  
            # ni         14.00670      N           3        nitrogen in charged imidazole ring
            elif ring == 5 and tf.count_neigh(atom.neighbor_info[1], element='C', ring=5, nb=3) >= 2 and formula == 'C3-H4-N2':
                atom.nta_type = 'ni'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # nho     14.00670      N           3        sp2 nitrogen in 6 membered ring next to a carbonyl
            elif ring == 6 and tf.count_neigh(atom.neighbor_info[2], element='O', ring=0, nb=1) >= 1 and 'C' in elements1:
                atom.nta_type = 'nho'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # npc     14.00670       N          3        sp2 nitrogen in 5- or 6- membered ring and with a heavy atom
            elif ring >= 5 and tf.count_heavies(elements1, heavies) > 0 or tf.count_heavies(elements2, heavies) > 0:  
                atom.nta_type = 'npc'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'   
                
            # nh       14.00670      N           3        sp2 nitrogen in 5 or 6 membered ring
            elif ring >= 5:  
                atom.nta_type = 'nh'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n3n      14.00670      N          3        sp2 nitrogen in 3- membered ring (ASSUME n means w/ N like H for c3h)
            elif ring == 3 and 'N' in elements1:
                atom.nta_type = 'n3n'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n3m     14.00670      N          3        sp3 nitrogen in 3- membered ring (ASSUME m means member like m for c3m)
            elif ring == 3:
                atom.nta_type = 'n3m'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n4n      14.00670      N          3        sp2 nitrogen in 4- membered ring (ASSUME n means w/ N like H for c4h)
            elif ring == 4 and 'N' in elements1 or 'N' in elements2:
                atom.nta_type = 'n4n'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n4m     14.00670      N          3        sp3 nitrogen in 4- membered ring (ASSUME m means member like m for c4m)
            elif ring == 4 and 'N' in elements1 or 'N' in elements2:
                atom.nta_type = 'n4m'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # no       14.006700    N          3        Sp2 nitrogen in nitro group 
            elif ring == 0 and elements1.count('O') == 2 and elements1.count('C') == 1:
                atom.nta_type = 'no'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n3       14.006700    N          3        Sp3 nitrogen with three substituents
            elif element == 'n':
                atom.nta_type = 'n3'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict CVFF-IFF atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Nitrogen-3-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Nitrogen-3-connect-no-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                    
        #############################################################################################
        # 4-bonded Nitrogen atom typing (ordering of nested if/elif/else statements set precedence) #
        #############################################################################################
        elif element == 'N' and nb == 4:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #---------------------------------------------------------------------------------#
            # Strict CVFF-IFF atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # nh+     14.00670      N           3        protonated nitrogen in 6 membered ring
            if ring == 6:
                atom.nta_type = 'nh+'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n4      14.00670      N          4        sp3 nitrogen in protonated amines
            # n+      14.00670      N          4        sp3 nitrogen in protonated amines 
            elif element == 'N':
                atom.nta_type = 'n4'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict CVFF-IFF atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Nitrogen-4-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Nitrogen-4-connect-no-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                    
                    
        ###########################################################################################
        # 1-bonded Sulfur atom typing (ordering of nested if/elif/else statements set precedence) #
        ###########################################################################################
        elif element == 'S' and nb == 1:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #---------------------------------------------------------------------------------#
            # Strict CVFF-IFF atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#                
            # s-      32.06400      S          1        partial double sulfur
            if element == 'S':
                atom.nta_type = "s-"; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict CVFF-IFF atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                tag = 'Sulfur-1-connect'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                
                
        ###########################################################################################
        # 2-bonded Sulfur atom typing (ordering of nested if/elif/else statements set precedence) #
        ###########################################################################################
        elif element == 'S' and nb == 2:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #---------------------------------------------------------------------------------#
            # Strict CVFF-IFF atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # s3e    32.06400      S          2        sulfur  in three membered ring
            if ring == 3:
                atom.nta_type = 's3e'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # s4e    32.06400      S          2        sulfur  in four membered ring
            elif ring == 4:
                atom.nta_type = 's4e'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                   
            # sp      32.06400      S          2        sulfur in an aromatic ring (e.g. thiophene)
            elif ring > 4:
                atom.nta_type = 'sp'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # s'      32.064000    S           2        Sulfur in thioketone (>C=S) group
            elif ring == 0 and elements1.count('C') == 1:
                atom.nta_type = 'sc'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # sc      32.06400      S          2        sp3 sulfur in methionines (C-S-C) group
            # s       32.064000    S           2        Sulfur in methionine (C-S-C) group
            elif ring == 0 and elements1.count('C') == 2:
                atom.nta_type = 'sc'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # sh      32.06400      S          2        sp3 sulfur in sulfhydryl (-SH) group (e.g. cysteine)
            elif ring == 0 and elements1.count('H') >= 1 and elements1.count('C') >= 1:
                atom.nta_type = 'sh'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # s1      32.06400      S          2        sp3 sulfur involved in (S-S) group of disulfides
            elif elements1.count('S') >= 1:
                atom.nta_type = 's1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict CVFF-IFF atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Sulfur-2-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Sulfur-2-connect-no-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                    
                    
        ############################################################################################
        # 4-bonded Silicon atom typing (ordering of nested if/elif/else statements set precedence) #
        ############################################################################################
        elif element == 'Si':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #----------------------------------------------------------------------------#
            # Strict IFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------#
            # sz    28.086000    Si          1        Silicon atom in a Zeolite or Silicate
            if nb == 1:          
                atom.nta_type = 'sz'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # si    28.086000    Si          4        Silicon 
            elif nb == 4:          
                atom.nta_type = 'si'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'

                    
        ###########################################################################
        # Phosphorous (ordering of nested if/elif/else statements set precedence) #
        ###########################################################################
        elif element == 'P':
            # pz    30.973800     P           1        Phosphorous atom in ALPO type structure 
            if nb == 1:
                atom.nta_type = 'pz'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # p     30.973800      P           4        General phosphorous atom
            elif nb == 4:
                atom.nta_type = 'p'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
        ####################################################################################
        # Fluorine atom typing (ordering of nested if/elif/else statements set precedence) #
        ####################################################################################
        elif element == 'F':
            # f     18.998400    F           1        Fluorine bonded to a carbon
            atom.nta_type = 'f'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ######################################################################
        # Helium (ordering of nested if/elif/else statements set precedence) #
        ######################################################################
        elif element == 'He':
            # he     4.002600    He          0        Helium atom
            atom.nta_type = 'he'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        #########################################################################
        # Deuterium (ordering of nested if/elif/else statements set precedence) #
        #########################################################################
        elif element == 'D':
            # dw     2.014000    D           1        Deuterium in heavy water
            if formula == 'D2-O1':
                atom.nta_type = 'dw'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            
            # d        2.014000    D           1        General Deuterium Atom
            elif element == 'D':
                atom.nta_type = 'd'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
        ########################################################################
        # Chlorine (ordering of nested if/elif/else statements set precedence) #
        ########################################################################
        elif element == 'Cl':
            # Cl    35.453000    Cl          1        Chloride ion  Cl-
            atom.nta_type = 'cl'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ########################################################################
        # Calcium (ordering of nested if/elif/else statements set precedence) #
        ########################################################################
        elif element == 'Ca':
            # ca+     40.08000     Ca          1        calcium ion 
            atom.nta_type = 'ca+'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        #######################################################################
        # Bromine (ordering of nested if/elif/else statements set precedence) #
        #######################################################################
        elif element == 'Br':
            # br    79.909000    Br          1        Bromine bonded to a carbon
            if 'C' in elements1:
                atom.nta_type = 'br'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'

            # Br    79.904000    Br          1        Bromide ion   Br-
            elif element == 'Br':
                atom.nta_type = 'Br'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
        #####################################################################
        # Argon (ordering of nested if/elif/else statements set precedence) #
        #####################################################################
        elif element == 'Ar':
            # ar    39.948       Ar          0        Argon atom
            atom.nta_type = 'ar'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'

                

                
                    
                    
                    
                    
                    
                    
                    
        ############################################################################
        # else all fails print to screen to tell user they must type by themselves #
        ############################################################################
        else:
            # print(f'atom-id {i} failed to be characterized');
            tf.write2failed(f, atom, ff_name, i); tally['failed'] += 1;
            
        ##################
        # DEBUG PRINTING #
        ##################
        #print(i, element, ring, elements1, elements2, rings1, rings1.count(0), atom.nta_type)
                
        
    ####################################################################################
    # Sort final mm.atoms[ID] to ensure written files have atomids in asscending order #
    ####################################################################################
    mm.atoms = dict(OrderedDict(sorted(mm.atoms.items())))
    
    
    ##################################################################
    # Add to the following instances to mm for data logging purposes #
    ##################################################################
    mm.supported_types = supported_types; mm.tally = tally; mm.assumed = assumed;
    
    
    ##########################################################################
    # Close the open files and check if they are zero in size; if so delete #
    #########################################################################
    a.close(); f.close(); asize = os.path.getsize(a.name); fsize = os.path.getsize(f.name);
    try:
        if asize == 0: os.remove(a.name);
        if fsize == 0: os.remove(f.name);
    except: a.close(); f.close();
    
    return mm