# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
December 15, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

Atom typing rules for PCFF with user-defined assumed types
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


##################################################################
# Function to perform atom-typing for PCFF forcefield atom types #
##################################################################
def nta(mm, basename, ff_name):
    tally = {'found':0, 'assumed':0, 'failed':0} # To tally findings
    
    ############################################
    # Set PCFF specific informations and flags #
    ############################################
    # Log supported types
    supported_types = {'Carbon':      ['ct', 'c+', 'cr', 'c-', 'c5', 'cs', 'cp', 'c_0', 'c_1', 'c_2',
                                       'cz', 'ci', 'c=', 'c=1', 'c=2', 'c3h', 'c3m', 'c4h', 'c4m', 'c_a',
                                       'cg', 'coh', 'c1', 'c2', 'c3'],
                       
                       'Hydrogen':    ['hi', 'hc', 'hw', 'hos', 'ho2', 'ho', 'hn2', 'hn', 'h*', 'hsi',
                                       'hs', 'h'],
                       
                       'Oxygen':      ['o_1', 'oo', 'o=', 'o-', 'o*', 'oz', 'o_2', 'oc', 'o3e', 'o4e',
                                       'op', 'osh', 'osi', 'oh', 'ob'],
                       
                       'Nitrogen':    ['nt', 'np', 'n=', 'n=1', 'n=2', 'n1', 'n2', 'n_2', 'nn', 'na',
                                       'nho', 'ni', 'npc', 'nh', 'n3n', 'n3m', 'n4n', 'n4m', 'n', 'nh+',
                                       'n4'],
                       
                       'Sulfer':      ["s'", 's-', 's3e', 's4e', 'sp', 'sc', 'sh', 's1', 's', 'sf'],
                        
                       'Silicone':    ['sio', 'si'],
                       
                       'Fluorine':    ['f'],
                       'Xenon':       ['xe'],
                       'Neon':        ['ne'],
                       'Krypton':     ['kr'],
                       'Helium':      ['he'],
                       'Deuterium':   ['dw'],
                       'Chlorine':    ['cl'],
                       'Calcium':     ['ca'],
                       'Bromine':     ['br'],
                       'Argon':       ['ar'],
                       'Phosphorous': ['p', 'p='],
                       'Molybdenum':  ['Mo']
                       }
    
    # Set lst of heavy elements (NOTE: MUST BE UPDATED IF MORE IFF ATOM TYPES GET CODED IN)
    heavies = ['C', 'O', 'N', 'S', 'F', 'Si', 'Xe', 'Ne', 'Kr', 'Cl', 'Br', 'Ar', 'P']
    
    # PCFF flags to help set correct atom-types (True to use False not to use)

    
    
    ######################################################################
    # Set Assumed PCFF atomtypes if code fails to find correct atom-type #
    ######################################################################
    assumed = {'Carbon-2-connect-no-ring':    'ct',   # 2-connects or less: Sp1 carbon no ring
               'Carbon-2-connect-in-ring':    'cp',   # 2-connects or less: Sp1 carbon in ring
               'Carbon-3-connect-no-ring':    'c=2',  # 3-connects: Sp2 carbon no ring
               'Carbon-3-connect-in-ring':    'cp',   # 3-connects: Sp2 carbon in ring
               'Carbon-4-connect-no-ring':    'c2',   # 4-connects: Sp3 carbon no ring
               'Carbon-4-connect-in-ring':    'c2',   # 4-connects: Sp3 carbon in ring
               'Hydrogen-1-connect':          'h',    # 1-connect : Hydrogen
               'Oxygen-1-connect':            'o=',   # 1-connect : Oxygen
               'Oxygen-2-connect-no-ring':    'o',    # 2-connects: Oxygen no ring
               'Oxygen-2-connect-in-ring':    'op',   # 2-connects: Oxygen in ring
               'Oxygen-3-connect-no-ring':    'ob',   # 3-connects: Oxygen no ring
               'Oxygen-3-connect-in-ring':    'ob',   # 3-connects: Oxygen in ring
               'Nitrogen-1-connect':          'nt',   # 1-connect : Nitrogen
               'Nitrogen-2-connect-no-ring':  'n=',   # 2-connects: Nitrogen no ring
               'Nitrogen-2-connect-in-ring':  'np',   # 2-connects: Nitrogen in ring
               'Nitrogen-3-connect-no-ring':  'n',    # 3-connects: Nitrogen no ring
               'Nitrogen-3-connect-in-ring':  'nh',   # 3-connects: Nitrogen in ring
               'Nitrogen-4-connect-no-ring':  'n2',   # 4-connects: Nitrogen no ring
               'Nitrogen-4-connect-in-ring':  'nh+',  # 4-connects: Nitrogen in ring
               'Sulfur-1-connect':            's-',   # 1-connect : Sulfur
               'Sulfur-2-connect-no-ring':    's',    # 2-connects: Sulfur no ring
               'Sulfur-2-connect-in-ring':    'sp',   # 2-connects: Sulfur in ring
               'Sulfur-4-connect-no-ring':    's_m',  # 4-connects: Sulfur no ring
               'Sulfur-4-connect-in-ring':    's_m',  # 4-connects: Sulfur in ring
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
                
            #-----------------------------------------------------------------------------#
            # Strict PCFF atom-typing that occurs after User defined atom-typing attempts #
            #-----------------------------------------------------------------------------#
            # ct      12.01115      C          2        sp carbon involved in a triple bond
            elif ring == 0:
                tally['found'] += 1
                atom.nta_type = 'ct'
                atom.nta_info = 'Correctly found'
                
            #------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict PCFF atom-typing fails used assumed #
            #------------------------------------------------------------------------------#
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
                
            #-----------------------------------------------------------------------------#
            # Strict PCFF atom-typing that occurs after User defined atom-typing attempts #
            #-----------------------------------------------------------------------------#
            # c+       12.01115      C          3        C in guanidinium group
            elif ring == 0 and elements1.count('N') == 3 and formula == 'C1-H5-N3':              
                atom.nta_type = 'c+'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # cr        12.01115      C          3        C in neutral arginine 
            elif ring == 0 and elements1.count('C') == 2 and elements1.count('N') == 1 and formula == 'C6-H14-N4-O2':
                atom.nta_type = 'cr'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            elif ring == 0 and elements1.count('N') == 3 and formula == 'C6-H14-N4-O2':
                atom.nta_type = 'cr'; tally['found'] += 1;
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
                
            # c_0     12.01115      C          3        carbonyl carbon of aldehydes, ketones  
            elif ring == 0 and elements1.count('C') == 2 and elements1.count('O') == 1:
                atom.nta_type = 'c_0'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            # Ketone (O-index=2 in nbs1)
            elif ring == 0 and elements1.count('C') == 2 and elements1.count('O') == 1 and nbs1[2] == 0:
                atom.nta_type = 'c_0'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            # Aldehyde (O-index=2 in nbs1)
            elif ring == 0 and elements1.count('C') == 1 and elements1.count('O') == 1 and elements1.count('H') == 1 and nbs1[2] == 0:
                atom.nta_type = 'c_0'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c_1     12.01115      C          3        carbonyl carbon of acid, ester, amide 
            # Amide (C-index=0, N-index=1, O-index=2 in elements1)
            elif ring == 0 and elements1[0] == 'C' and elements1[1] == 'N' and elements1[2] == 'O':# and rings1.count(0) == 3:
                atom.nta_type = 'c_1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            # Ester, Carboxylic acid (C-index=0, O-index=1, O-index=2 in elements1)
            elif ring == 0 and elements1[0] == 'C' and elements1[1] == 'O' and elements1[2] == 'O':# and rings1.count(0) == 3:
                atom.nta_type = 'c_1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c_2     12.01100      C          3        carbonyl carbon of carbamate, urea
            # Carbamate (N-index=0, O-index=1, O-index=2 in elements1)
            elif ring == 0 and elements1[0] == 'N' and elements1[1] == 'O' and elements1[2] == 'O':# and rings1.count(0) == 3:
                atom.nta_type = 'c_2'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            # Urea (N-index=0, N-index=1, O-index=2 in elements1)
            elif ring == 0 and elements1[0] == 'N' and elements1[1] == 'N' and elements1[2] == 'O':# and rings1.count(0) == 3:
                atom.nta_type = 'c_2'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # cz        12.01100      C          3        carbonyl carbon of carbonate
            elif ring == 0 and elements1.count('O') == 3:
                atom.nta_type = 'cz'; tally['found'] += 1;
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
            elif ring == 0 and nbs2.count(1) == 2:
                atom.nta_type = 'c=1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c=2     12.01115      C          3        nonaromatic doubly bonded carbon  
            elif ring == 0:
                atom.nta_type = 'c=2'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict PCFF atom-typing fails used assumed #
            #------------------------------------------------------------------------------#
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
            # # For 6 member-ring open valence
            # if ring >= 6 and elements1.count('C') == 4 and rings1.count(6) >= 3:
            #     atom.nta_type = 'cp'; tally['found'] += 1;
            #     atom.nta_info = 'User-defined for Graphite xlinks (element:C, ring:>=6, nb:4)'
            
            #-----------------------------------------------------------------------------#
            # Strict PCFF atom-typing that occurs after User defined atom-typing attempts #
            #-----------------------------------------------------------------------------#
            # c3h      12.01115      C          4        sp3 carbon in 3-membered ring with hydrogens
            if ring == 3 and 'H' in elements1:              
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
                
            # c_a      12.01115      C          4        general amino acid alpha carbon (sp3)
            elif ring == 0 and elements1.count('C') == 2 and elements1.count('H') == 1 and elements1.count('N') == 1:              
                atom.nta_type = 'c_a'; tally['found'] += 1;
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
                
            # c4o     12.01115       C          4        Carbon atom, sp3, bonded to oxygen (+0.054, compass)
            # Setting near end so all other types can be tried 1st since it is not very unique
            elif 'O' in elements1:              
                atom.nta_type = 'c4o'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # General Sp3 C atom types
            # c         12.01115      C          4        generic SP3 carbon
            # c4        12.01115      C          4        Carbon atom, sp3, generic, 4 bonds (compass)
            elif element == 'C':
                atom.nta_type = 'c'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict PCFF atom-typing fails used assumed #
            #------------------------------------------------------------------------------#
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
            
            #-----------------------------------------------------------------------------#
            # Strict PCFF atom-typing that occurs after User defined atom-typing attempts #
            #-----------------------------------------------------------------------------#
            # hi       1.00797      H          1        Hydrogen in charged imidazole ring
            if rings1.count(5) > 0 and tf.count_neigh(atom.neighbor_info[2], element='N', ring=5, nb=False) > 0 and formula == 'C3-H4-N2':
                atom.nta_type = 'hi'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # hc      1.00797      H          1        hydrogen bonded to carbon
            elif elements1[0] == 'C':
                atom.nta_type = 'hc'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'

                
            # hw      1.00797      H         1        hydrogen in water (+0.41 in this model !)
            elif elements1[0] == 'O' and formula == 'H2-O1':
                atom.nta_type = 'hw'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # hos     1.00782      H          1        hydrogen atom in terminal hydroxyl group on silicon
            elif elements1[0] == 'O' and 'Si' in elements2 and 'Mg' not in elements2 and 'Al' not in elements2:
                atom.nta_type = 'hos'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # ho2    1.00800      H          1        hydroxyl hydrogen
            elif elements1[0] == 'O' and 'C' in elements2 or 'S' in elements2:
                atom.nta_type = 'ho2'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # ho      1.00797      H          1        hydrogen bonded to oxygen
            elif elements1[0] == 'O':
                atom.nta_type = 'ho'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # # hn2    1.00800      H          1        amino hydrogen
            # # An amino group is a nitrogen atom bonded to two hydrogen atoms
            # elif elements1[0] == 'N' and len(elements2) == 2 and elements2[0] == 'C' and elements2[1] == 'H':
            #     atom.nta_type = 'hn2'; tally['found'] += 1;
            #     atom.nta_info = 'Correctly found'

            # hn      1.00797      H          1        hydrogen bonded to nitrogen
            elif elements1[0] == 'N':
                atom.nta_type = 'hn'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # h*     1.00797      H          1        hydrogen bonded to nitrogen, Oxygen
            elif elements1[0] == 'N' or elements1[0] == 'O':
                atom.nta_type = 'h*'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # hsi      1.00800      H         1        silane hydrogen
            elif elements1[0] == 'S' and formula == 'H4-S1':
                atom.nta_type = 'hsi'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # hs       1.00797      H         1        hydrogen bonded to sulfur
            elif elements1[0] == 'S':
                atom.nta_type = 'hs'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # h         1.00797      H         1        generic hydrogen bound to C, Si,or H
            elif elements1[0] == 'C' or elements1[0] == 'Si' or elements1[0] == 'H':
                atom.nta_type = 'h'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict PCFF atom-typing fails used assumed #
            #------------------------------------------------------------------------------#
            else:
                tag = 'Hydrogen-1-connect'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                
                
        ###########################################################################################
        # 1-bonded Oxygen atom typing (ordering of nested if/elif/else statements set precedence) #
        ###########################################################################################
        elif element == 'O' and nb == 1:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #-----------------------------------------------------------------------------#
            # Strict PCFF atom-typing that occurs after User defined atom-typing attempts #
            #-----------------------------------------------------------------------------#
            # o_1     15.99940      O          1        oxygen in carbonyl group
            if ring == 0 and elements1.count('C') == 1 and 'C' in elements2:# and 'O' not in elements2:              
                atom.nta_type = 'o_1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # oo       15.99940      O          1        oxygen in carbonyl group, carbonate only
            elif ring == 0 and elements1.count('C') == 1 and 'O' in elements2 and 'C' not in elements2:              
                atom.nta_type = 'oo'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o=       15.99940      O          1        oxygen double bonded to O, C, S, N, P 
            elif ring == 0 and elements1.count('O') > 0 or elements1.count('C') > 0 or elements1.count('S') > 0  or elements1.count('N') > 0  or elements1.count('P') > 0:              
                atom.nta_type = 'o='; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o-        15.99940       O          1        partial double oxygen
            elif element == 'O':              
                atom.nta_type = 'o-'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict PCFF atom-typing fails used assumed #
            #------------------------------------------------------------------------------#
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
            
            #----------_------------------------------------------------------------------#
            # Strict PCFF atom-typing that occurs after User defined atom-typing attempts #
            #-----------_-----------------------------------------------------------------#
            # o*       15.99940      O          2        oxygen in water (-0.82 in this moddel !)
            if formula == 'H2-O1':
                atom.nta_type = 'o*'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # oz       15.99940      O          2        ester oxygen in carbonate
            elif ring == 0  and len(rings1) == rings1.count(0) and elements1.count('C') == 2 and tf.count_neigh(atom.neighbor_info[2], element='O', ring=0, nb=1) == 1 and tf.count_neigh(atom.neighbor_info[2], element='O', ring=0, nb=2) == 1:              
                atom.nta_type = 'oz'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o_2     15.99940      O          2        ester oxygen
            # oe       15.99940      O          2        sp3 oxygen  in ester
            elif ring == 0  and len(rings1) == rings1.count(0) and elements1.count('C') == 2 and tf.count_neigh(atom.neighbor_info[2], element='O', ring=0, nb=1) == 1 and elements2.count('C') >= 1:              
                atom.nta_type = 'o_2'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # oc       15.99940      O          2        sp3 oxygen  in ether or acetals
            elif ring == 0 and elements1.count('C') == 2:              
                atom.nta_type = 'oc'; tally['found'] += 1;
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
                
            # osh     15.99491      O          2        oxygen atom in terminal hydroxyl group on silicon
            elif ring == 0 and elements1.count('C') == 1 and elements1.count('Si') == 1:              
                atom.nta_type = 'osh'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # osi     16.00000      O          2        siloxane oxygen
            # oss     15.99491      O          2        oxygen atom betweem two silicons
            elif ring == 0 and 'Si' in elements1 and 'O' in elements2:              
                atom.nta_type = 'osi'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # oh       15.99940      O          2        oxygen bonded to hydrogen 
            elif 'H' in elements1:
                atom.nta_type = 'oh'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o         15.99940      O          2        generic SP3 oxygen 
            elif element == 'O':
                atom.nta_type = 'o'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'

                
            #------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict PCFF atom-typing fails used assumed #
            #------------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Oxygen-2-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Oxygen-2-connect-no-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)

        
        ###########################################################################################
        # 3-bonded Oxygen atom typing (ordering of nested if/elif/else statements set precedence) #
        ###########################################################################################
        elif element == 'O' and nb == 3:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #-----------------------------------------------------------------------------#
            # Strict PCFF atom-typing that occurs after User defined atom-typing attempts #
            #-----------------------------------------------------------------------------#
            # ob      15.99491      O          3        oxygen atom in bridging hydroxyl group
            if element == 'O':
                atom.nta_type = 'ob'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict PCFF atom-typing fails used assumed #
            #------------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Oxygen-3-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Oxygen-3-connect-no-ring'; tally['assumed'] += 1;
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
            
            #-----------------------------------------------------------------------------#
            # Strict PCFF atom-typing that occurs after User defined atom-typing attempts #
            #-----------------------------------------------------------------------------#
            # nt      14.00670      N          1        sp nitrogen involved in a triple bond 
            # nz      14.00670      N          1        sp3 nitrogen bonded to two atoms
            if element == 'N':
                atom.nta_type = 'nt'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict PCFF atom-typing fails used assumed #
            #------------------------------------------------------------------------------#
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
            
            #-----------------------------------------------------------------------------#
            # Strict PCFF atom-typing that occurs after User defined atom-typing attempts #
            #-----------------------------------------------------------------------------#
            # np       14.00670      N          2        sp2 nitrogen in 5- or 6- membered ring
            if ring >= 5:
                atom.nta_type = 'np'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n=       14.00670      N          2        non aromatic end doubly bonded nitrogen
            elif ring == 0 and nbs1.count(1) >= 1:
                atom.nta_type = 'n='; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n=1     14.00670      N          2        non aromatic, next to end doubly bonded carbon
            elif ring == 0 and nbs2.count(1) == 1:
                atom.nta_type = 'n=1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n=2     14.00670      N          2        non aromatic doubly bonded nitrogen  
            elif ring == 0:
                atom.nta_type = 'n=2'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict PCFF atom-typing fails used assumed #
            #------------------------------------------------------------------------------#
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
            
            #-----------------------------------------------------------------------------#
            # Strict PCFF atom-typing that occurs after User defined atom-typing attempts #
            #-----------------------------------------------------------------------------#
            # n1        14.00670      N          3        sp2 nitrogen in charged arginine
            if ring == 0 and formula == 'C6-H14-N4-O2':
                atom.nta_type = 'n1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n2        14.00670      N          3        sp2 nitrogen (NH2) in guanidinium group (HN=C(NH2)2)
            # nr        14.00670      N          3        sp2 nitrogen (NH2) in guanidinium group (HN=C(NH2)2)
            elif ring == 0 and formula == 'C1-H5-N3':
                atom.nta_type = 'n1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n_2     14.01000       N          3        nitrogen of urethane
            elif ring == 0 and formula == 'C3-H7-N1-O2':
                atom.nta_type = 'n_2'; tally['found'] += 1;
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
            elif ring == 0 and elements1.count('H') == 1 and elements1.count('C') == 2:
                atom.nta_type = 'na'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            elif ring == 0 and elements1.count('C') == 3:
                atom.nta_type = 'na'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # nho     14.00670      N           3        sp2 nitrogen in 6 membered ring next to a carbonyl
            elif ring == 6 and tf.count_neigh(atom.neighbor_info[2], element='O', ring=0, nb=1) >= 1 and 'C' in elements1:
                atom.nta_type = 'nho'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # ni         14.00670      N           3        nitrogen in charged imidazole ring
            elif ring == 5 and tf.count_neigh(atom.neighbor_info[1], element='C', ring=5, nb=3) >= 2 and formula == 'C3-H4-N2':
                atom.nta_type = 'ni'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # npc     14.00670       N          3        sp2 nitrogen in 5- or 6- membered ring and with a heavy atom
            elif ring >= 5 and tf.count_heavies(elements1, heavies) > 0 or ring >= 5 and tf.count_heavies(elements2, heavies) > 0:  
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
                
            # n          14.00670      N          3        generic sp2 nitrogen (in amids))
            elif element == 'n':
                atom.nta_type = 'n'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict PCFF atom-typing fails used assumed #
            #------------------------------------------------------------------------------#
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
            
            #----------------------------------------------------------------------------#
            # Strict IFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------#
            # nh+     14.00670      N           3        protonated nitrogen in 6 membered ring
            if ring == 6:
                atom.nta_type = 'nh+'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n4      14.00670      N          4        sp3 nitrogen in protonated amines
            # n+      14.00670      N          4        sp3 nitrogen in protonated amines 
            elif element == 'N':
                atom.nta_type = 'n4'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #-----------------------------------------------------------------------------#
            # If user defined attempt fails and Strict IFF atom-typing fails used assumed #
            #-----------------------------------------------------------------------------#
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
            
            #----------------------------------------------------------------------------#
            # Strict IFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------#
            # s'      32.06400      S          1        S in thioketone group
            if 'C' in elements1:
                atom.nta_type = "s'"; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # s-      32.06400      S          1        partial double sulfur
            elif element == 'S':
                atom.nta_type = "s-"; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #-----------------------------------------------------------------------------#
            # If user defined attempt fails and Strict IFF atom-typing fails used assumed #
            #-----------------------------------------------------------------------------#
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
            
            #----------------------------------------------------------------------------#
            # Strict IFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------#
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
                
            # sc      32.06400      S          2        sp3 sulfur in methionines (C-S-C) group
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
                
            # s       32.06400       S          2        sp3 sulfur
            elif element == 'S':
                atom.nta_type = 's'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #-----------------------------------------------------------------------------#
            # If user defined attempt fails and Strict IFF atom-typing fails used assumed #
            #-----------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Sulfur-2-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Sulfur-2-connect-no-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                    
                    
        ###########################################################################################
        # 4-bonded Sulfur atom typing (ordering of nested if/elif/else statements set precedence) #
        ###########################################################################################
        elif element == 'S' and nb == 4:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #----------------------------------------------------------------------------#
            # Strict IFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------#
            # sf      32.06400      S          1        S in sulfonate group
            if elements1.count('O') == 3:
                atom.nta_type = 'sf'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #-----------------------------------------------------------------------------#
            # If user defined attempt fails and Strict IFF atom-typing fails used assumed #
            #-----------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Sulfur-4-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Sulfur-4-connect-no-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                    
            
        ############################################################################################
        # 4-bonded Silicon atom typing (ordering of nested if/elif/else statements set precedence) #
        ############################################################################################
        elif element == 'Si' and nb == 4:
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #----------------------------------------------------------------------------#
            # Strict IFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------#
            # sio     28.08600     Si          4         siloxane silicon
            if ring == 0 and 'O' in elements1 and 'Si' in elements2:          
                atom.nta_type = 's_m'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # si        28.08600     Si          4         silicon atom
            elif element == 'Si':
                atom.nta_type = 'si'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
                
        ####################################################################################
        # Fluorine atom typing (ordering of nested if/elif/else statements set precedence) #
        ####################################################################################
        elif element == 'F':
            # f       18.99840      F          1        fluorine  atom
            atom.nta_type = 'f'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        #####################################################################
        # Xenon (ordering of nested if/elif/else statements set precedence) #
        #####################################################################
        elif element == 'Xe':
            # xe     131.30000     Xe          0        Xenon
            atom.nta_type = 'xe'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ####################################################################
        # Neon (ordering of nested if/elif/else statements set precedence) #
        ####################################################################
        elif element == 'Ne':
            # ne      20.18300     Ne          0        Neon
            atom.nta_type = 'ne'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        #######################################################################
        # Krypton (ordering of nested if/elif/else statements set precedence) #
        #######################################################################
        elif element == 'Kr':
            # kr      83.80000     Kr          0        Krypton
            atom.nta_type = 'kr'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ######################################################################
        # Helium (ordering of nested if/elif/else statements set precedence) #
        ######################################################################
        elif element == 'He':
            # he       4.00300     He          0        Helium
            atom.nta_type = 'he'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        #########################################################################
        # Deuterium (ordering of nested if/elif/else statements set precedence) #
        #########################################################################
        elif element == 'D':
            # dw       2.01400      D          1        deuterium in heivy water
            atom.nta_type = 'dw'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ########################################################################
        # Chlorine (ordering of nested if/elif/else statements set precedence) #
        ########################################################################
        elif element == 'Cl':
            # cl      35.45300       Cl          1        chlorine atom
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
            # br      79.90900     Br          1        bromine atom
            atom.nta_type = 'br'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        #####################################################################
        # Argon (ordering of nested if/elif/else statements set precedence) #
        #####################################################################
        elif element == 'Ar':
            # ar      39.94400     Ar          0        Argon
            atom.nta_type = 'ar'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ###########################################################################
        # Phosphorous (ordering of nested if/elif/else statements set precedence) #
        ###########################################################################
        elif element == 'P':
            # p       30.97380      P          4        general phosphorous atom
            if nb == 4:
                atom.nta_type = 'p'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # p=      30.97380      P          5        phosphazene phosphorous atom
            elif nb == 5:
                atom.nta_type = 'p='; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                    
        ##########################################################################
        # Molybdenum (ordering of nested if/elif/else statements set precedence) #
        ##########################################################################
        elif element == 'Mo':
            # Mo      95.94000     Mo          0        Molybdenum metal
            atom.nta_type = 'Mo'; tally['found'] += 1;
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