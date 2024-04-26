# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 18th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

Atom typing rules for compass with user-defined assumed types
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
    supported_types = {'Carbon':      ['c1o', 'c2=', 'c3a', "c3'", 'c41o', 'c43o', 'c4z', 'c43',
                                       'c44', 'c4o', 'c4'],
                       
                       'Hydrogen':    ['h1h', 'h1o', 'h1'],
                       
                       'Oxygen':      ['o1=', 'o1=*', 'o1c', 'o2n', 'o1n', 'o1o', 'o12', 'o2s', 'o2h',
                                       'o2e', 'o2z', 'o2'],
                       
                       'Nitrogen':    ['n1z', 'n1n', 'n1o', 'n2t', 'n2z', 'n2o', 'n2=', 'n3o', 'n3m'],
                       
                       'Sulfer':      ['s1=', 's2='],
                        
                       'Silicone':    ['si4c', 'si4'],
                       
                       'Xenon':       ['xe'],
                       'Neon':        ['ne'],
                       'Krypton':     ['kr'],
                       'Helium':      ['he'],
                       'Argon':       ['ar'],
                       'Phosphorous': ['p4=']
                       }
    
    # Set lst of heavy elements (NOTE: MUST BE UPDATED IF MORE IFF ATOM TYPES GET CODED IN)
    heavies = ['C', 'O', 'N', 'S', 'F', 'Si', 'Xe', 'Ne', 'Kr', 'Cl', 'Br', 'Ar', 'P']
    
    # Set lst of elements that are considered as functional groups to set alpha carbonds and maybe others ...
    functional_groups = {1: ['O', 'N', 'Br'],      # 1st neighs
                         2: ['O', 'N', 'C', 'H'],  # 2nd neighs
                         }
    
    # PCFF-IFF flags to help set correct atom-types (True to use False not to use)
    use_c4o_alpha_carbon = False # Option to use c4o alpha carbon atom-type (this one is difficult to get correct other atom types could be used in its place, hence the option)
    
    # compass flags to help set correct atom-types (True to use False not to use)
    supported_types = tf.update_supported_types(supported_types, element='Carbon', atomtype='c4o', flag=use_c4o_alpha_carbon)

    
    
    #########################################################################
    # Set Assumed compass atomtypes if code fails to find correct atom-type #
    #########################################################################
    assumed = {'Carbon-2-connect-no-ring':    'c1o',   # 2-connects or less: Sp1 carbon no ring
               'Carbon-2-connect-in-ring':    'c3a',   # 2-connects or less: Sp1 carbon in ring
               'Carbon-3-connect-no-ring':    "c3'",   # 3-connects: Sp2 carbon no ring
               'Carbon-3-connect-in-ring':    'c3a',   # 3-connects: Sp2 carbon in ring
               'Carbon-4-connect-no-ring':    'c4',    # 4-connects: Sp3 carbon no ring
               'Carbon-4-connect-in-ring':    'c4',    # 4-connects: Sp3 carbon in ring
               'Hydrogen-1-connect':          'h1',    # 1-connect : Hydrogen
               'Oxygen-1-connect':            'o1c',   # 1-connect : Oxygen
               'Oxygen-2-connect-no-ring':    'o2',    # 2-connects: Oxygen no ring
               'Oxygen-2-connect-in-ring':    'o2',    # 2-connects: Oxygen in ring
               'Nitrogen-1-connect':          'n1n',   # 1-connect : Nitrogen
               'Nitrogen-2-connect-no-ring':  'n2=',   # 2-connects: Nitrogen no ring
               'Nitrogen-2-connect-in-ring':  'n2=',   # 2-connects: Nitrogen in ring
               'Nitrogen-3-connect-no-ring':  'n3m',   # 3-connects: Nitrogen no ring
               'Nitrogen-3-connect-in-ring':  'n3m',   # 3-connects: Nitrogen in ring
               'Sulfur-1-connect':            's1=',   # 1-connect : Sulfur
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
        elements3 = tf.neigh_extract(atom, depth=3, info='element') # Example: ['C', 'C', 'N']
        #rings3 = tf.neigh_extract(atom, depth=3, info='ring') # Example: [3, 3, 3]
        
        # Build functional group tests for alpha carbon test (not the best way, but not sure of any other)
        functional_group_tests = {} # { neigh depth : [lst of elements that meet criteria in functional_groups] }
        alpha_carbon_functional_flag = True # Update if any depth is empty
        for depth in functional_groups:
            functional_group_tests[depth] = []
            elementsdepth = tf.neigh_extract(atom, depth=depth, info='element')
            for j in elementsdepth:
                if j in functional_groups[depth]:
                    functional_group_tests[depth].append(j)
        for j in functional_group_tests:
            if not functional_group_tests[j]:
                alpha_carbon_functional_flag = False
        #print(i, element, functional_group_tests, alpha_carbon_functional_flag)
            
        
        
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
            # For 5/6 member-ring open valence
            if ring >= 5:
                atom.nta_type = 'c3a'; tally['found'] += 1;
                atom.nta_info = 'User-defined for Open Valence ReaxFF polymerization (element:C, ring:>=5, nb:2)'
                
            #--------------------------------------------------------------------------------#
            # Strict compass atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # c1o  12.01115     C          carbon in CO
            elif formula == 'C1-O1':
                tally['found'] += 1
                atom.nta_type = 'c1o'
                atom.nta_info = 'Correctly found'
                
            # c2=  12.01115     C          carbon in CO2 and CS2
            elif formula == 'C1-O2':
                tally['found'] += 1
                atom.nta_type = 'c2='
                atom.nta_info = 'Correctly found'
                
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict compass atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
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
            # c3a  12.01115     C          aromatic carbon
            # For 7 and 8 member rings on carbonized structures
            if ring > 6:
                atom.nta_type = 'c3a'; tally['found'] += 1;
                atom.nta_info = 'User-defined for 7/8-ring carbonization (element:C, ring:>6, nb:3)'
                
            #--------------------------------------------------------------------------------#
            # Strict compass atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # c3a  12.01115     C          aromatic carbon
            elif ring >= 5:
                tally['found'] += 1
                atom.nta_type = 'c3a'
                atom.nta_info = 'Correctly found'
                
            # c3’   12.01115     C          carbonyl carbon [one polar substituent]
            elif ring == 0 and elements1.count('C') == 2 and elements1.count('O') == 1:
                atom.nta_type = "c3'"; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict compass atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
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
            
            #--------------------------------------------------------------------------------#
            # Strict compass atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # c41o  12.01115     C          carbon, sp3, in methanol 
            if formula == 'C1-H4-O1':
                atom.nta_type = 'c41o'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c43o  12.01115     C          carbon, sp3 in secondary alcohols
            elif ring == 0 and elements1.count('C') == 2 and elements1.count('O') == 1 and elements2.count('H') >= 1:
                atom.nta_type = 'c43o'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c4z    12.01115     C          carbon, sp3, bonded to -N3 (azides)
            elif ring == 0 and elements1.count('N') == 1 and elements2.count('N') >= 1 and elements3.count('N') >= 1:
                atom.nta_type = 'c4z'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c43    12.01115     C          sp3 carbon with three heavy atoms attached
            elif tf.count_heavies(elements1, heavies) == 3:
                atom.nta_type = 'c43'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c44  12.01115     C          sp3 carbon with four heavy atoms attached
            elif tf.count_heavies(elements1, heavies) == 4:
                atom.nta_type = 'c44'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c4o    12.01115     C          alpha carbon (e.g. alpha to oxygen in ethers and alcohols)
            elif ring == 0 and functional_group_tests and use_c4o_alpha_carbon:
                atom.nta_type = 'c4o'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # c4      12.01115     C          generic sp3 carbon
            elif element == 'C':
                atom.nta_type = 'c4'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict compass atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
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
        elif element == 'H':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #--------------------------------------------------------------------------------#
            # Strict compass atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # h1h   1.00797     H          hydrogen in H2      
            if formula == 'H2':
                atom.nta_type = 'h1h'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # h1o   1.00797     H          strongly polar hydrogen, bonded to O,F   
            elif elements1.count('O') >= 1 or elements1.count('F') >= 1:
                atom.nta_type = 'h1o'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
              
            # h1     1.00797     H          nonpolar hydrogen   
            elif element == 'H':
                atom.nta_type = 'h1'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict compass atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
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
            
            #--------------------------------------------------------------------------------#
            # Strict compass atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # o1=   15.99940     O          oxygen in NO2 and SO2 [and carbonyl]
            # Carbonyl
            if ring == 0 and elements1.count('C') == 1 and 'C' in elements2:               
                atom.nta_type = 'o1='; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            # NO2
            elif ring == 0 and formula == 'N1-O2':               
                atom.nta_type = 'o1='; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            # SO2
            elif ring == 0 and formula == 'O2-S1':               
                atom.nta_type = 'o1='; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o1=* 15.99940     O          oxygen in CO2
            elif ring == 0 and formula == 'C1-O2':               
                atom.nta_type = 'o1=*'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o1c   15.99940     O          oxygen in CO
            elif ring == 0 and formula == 'C1-O2':               
                atom.nta_type = 'o1c'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            
            # o1n   15.99940     O          oxygen in NO
            elif ring == 0 and formula == 'N1-O1':               
                atom.nta_type = 'o1n'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o1o   15.99940     O          oxygen in O2
            elif ring == 0 and formula == 'O2':               
                atom.nta_type = 'o1o'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o2n  15.99940     O          oxygen in nitrates
            elif ring == 0 and formula == 'N1-O3':               
                atom.nta_type = 'o2n'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o12   15.99940     O          oxygen in nitro group (-NO2)
            elif elements1.count('N') == 1 and elements2.count('O') == 1 and elements2.count('C') == 1:               
                atom.nta_type = 'o12'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict compass atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
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
            
            #--------------------------------------------------------------------------------#
            # Strict compass atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # o2s  15.99940     O          ester oxygen
            if ring == 0  and len(rings1) == rings1.count(0) and elements1.count('C') == 2 and tf.count_neigh(atom.neighbor_info[2], element='O', ring=0, nb=1) == 1 and elements2.count('C') == 2:              
                atom.nta_type = 'o2s'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o2h      15.99940      O          2        Oxygen atom in hydroxyl group (-0.57, compass) 
            elif ring == 0 and elements1.count('C') == 1 and elements1.count('H') == 1:              
                atom.nta_type = 'o2h'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            
            # o2e  15.99940     O          ether oxygen
            elif ring == 0 and elements1.count('C') == 2:              
                atom.nta_type = 'o2e'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o2z  15.99940     O          oxygen, in siloxanes and zeolites
            elif ring == 0 and 'Si' in elements1 and 'O' in elements2:              
                atom.nta_type = 'o2z'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # o2    15.99940     O          generic oxygen with two bonds attached
            elif element == 'O':
                atom.nta_type = 'o2'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict compass atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
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
            
            #----------------------------------------------------------------------------#
            # Strict IFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------#
            # n1z  14.00670     N          nitrogen, terminal atom in -N3 (azides)
            if ring == 0 and elements1.count('N') == 1 and elements2.count('N') >= 1 and elements3.count('C') == 1:
                atom.nta_type = 'n1z'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n1n  14.00670     N          nitrogen in N2
            elif ring == 0 and formula == 'N2':               
                atom.nta_type = 'n1n'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                     
            # n1o  14.00670     N          nitrogen in NO
            elif ring == 0 and formula == 'N1-O1':               
                atom.nta_type = 'n1o'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict compass atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
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
            
            #--------------------------------------------------------------------------------#
            # Strict compass atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # n2t    14.00670     N          nitrogen, central atom in -N3 (azides)
            if ring == 0 and elements1.count('N') == 2 and elements2.count('C') >= 1:
                atom.nta_type = 'n2t'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n2z    14.00670     N          nitrogen, first atom in -N3 (azides)
            elif ring == 0 and elements1.count('N') == 1 and elements1.count('C') >= 1 and elements2.count('N') >= 1:
                atom.nta_type = 'n2z'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n2o   14.00670     N          nitrogen in NO2
            elif ring == 0 and formula == 'N1-O2':
                atom.nta_type = 'n2o'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n2=   14.00670     N          nitrogen
            elif element == 'N':
                atom.nta_type = 'n2='; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict compass atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
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
            
            #----------------------------------------------------------------------------#
            # Strict IFF atom-typing that occurs after User defined atom-typing attempts #
            #----------------------------------------------------------------------------#
            # n3o    14.00670     N          nitrogen in nitro group
            if ring == 0 and elements1.count('O') == 2 and elements1.count('C') == 1:
                atom.nta_type = 'n3o'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # n3m   14.00670     N          sp3 nitrogen in amides without hydrogen
            elif ring == 0 and elements1.count('C') == 3:
                atom.nta_type = 'n3m'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict compass atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
            else:
                if ring > 0:
                    tag = 'Nitrogen-3-connect-in-ring'; tally['assumed'] += 1;
                    atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                    tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                else:
                    tag = 'Nitrogen-3-connect-no-ring'; tally['assumed'] += 1;
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
            
            #--------------------------------------------------------------------------------#
            # Strict compass atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # s1=  32.06400     S          sulfur in CS2
            if ring == 0 and formula == 'C1-S2':
                atom.nta_type = 's1='; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
            
            # s2=  32.06400     S          sulfur in SO2
            elif ring == 0 and formula == 'O2-S1':
                atom.nta_type = 's2='; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict compass atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
            else:
                tag = 'Sulfur-1-connect'; tally['assumed'] += 1;
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
            # si4c 28.08600     Si         a subset of si4, non-hydrogen atom attached [siloxanes]
            if ring == 0 and elements1.count('H') == 0:
                atom.nta_type = 'si4c'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
            # si4   28.08600     Si         generic silicon with four bonds attached
            else:
                atom.nta_type = 'si4c'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                
        #####################################################################
        # Xenon (ordering of nested if/elif/else statements set precedence) #
        #####################################################################
        elif element == 'Xe':
            # xe   1 31.30000     Xe         xenon
            atom.nta_type = 'xe'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ####################################################################
        # Neon (ordering of nested if/elif/else statements set precedence) #
        ####################################################################
        elif element == 'Ne':
            # ne    20.18300      Ne         neon
            atom.nta_type = 'ne'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        #######################################################################
        # Krypton (ordering of nested if/elif/else statements set precedence) #
        #######################################################################
        elif element == 'Kr':
            # kr     83.80000      Kr          krypton
            atom.nta_type = 'kr'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ######################################################################
        # Helium (ordering of nested if/elif/else statements set precedence) #
        ######################################################################
        elif element == 'He':
            # he    4.00300        He         helium
            atom.nta_type = 'he'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        #####################################################################
        # Argon (ordering of nested if/elif/else statements set precedence) #
        #####################################################################
        elif element == 'Ar':
            # ar    39.94400      Ar         argon
            atom.nta_type = 'ar'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ###########################################################################
        # Phosphorous (ordering of nested if/elif/else statements set precedence) #
        ###########################################################################
        elif element == 'P':
            # p4=  30.97380     P           phosphorous 
            atom.nta_type = 'p4='; tally['found'] += 1;
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