# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
May 24th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

Atom typing rules for DRIEDING with user-defined assumed types
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
    
    ################################################
    # Set DREIDING specific informations and flags #
    ################################################
    # Log supported types
    supported_types = {'Carbon':      ['C_11', 'C_1', 'C_R', 'C_RB', 'C_21Mid', 'C_2Mid', 'C_21', 'C_2',
                                       'C_33', 'C_32', 'C_31', 'C_3'],
                       
                       'Hydrogen':    ['H_'],
                       
                       'Oxygen':      ['O_1', 'O_R', 'O_R1', 'O_2', 'O_31', 'O_3'],
                       
                       'Nitrogen':    ['N_1', 'N_R', 'N_R1', 'N_21Mid', 'N_2Mid', 'N_21', 'N_2',
                                       'N_33', 'N_32', 'N_31', 'N_3'],
                       
                       'Sulfer':      ['S_3'],
                        
                       'Silicone':    ['Si3'],
                       
                       'Fluorine':    ['F_'],
                       'Chlorine':    ['Cl'],
                       'Calcium':     ['Ca'],
                       'Bromine':     ['Br'],
                       'Phosphorous': ['P_3']
                       }

    
    # DREIDING flags to help set correct atom-types (True to use False not to use)
    use_implicit_Hs_on_Carbon = True
    use_implicit_Hs_on_Oxygen = True
    use_implicit_Hs_on_Nitrogen = True
    
    
    # Update supported types lists if flags are used and set flag status in parentheses
    supported_types = tf.update_supported_types(supported_types, element='Carbon', atomtype='C_11', flag=use_implicit_Hs_on_Carbon)
    supported_types = tf.update_supported_types(supported_types, element='Carbon', atomtype='C_21Mid', flag=use_implicit_Hs_on_Carbon)
    supported_types = tf.update_supported_types(supported_types, element='Carbon', atomtype='C_21', flag=use_implicit_Hs_on_Carbon)
    supported_types = tf.update_supported_types(supported_types, element='Carbon', atomtype='C_33', flag=use_implicit_Hs_on_Carbon)
    supported_types = tf.update_supported_types(supported_types, element='Carbon', atomtype='C_32', flag=use_implicit_Hs_on_Carbon)
    supported_types = tf.update_supported_types(supported_types, element='Carbon', atomtype='C_31', flag=use_implicit_Hs_on_Carbon)
    supported_types = tf.update_supported_types(supported_types, element='Oxygen', atomtype='O_R1', flag=use_implicit_Hs_on_Oxygen)
    supported_types = tf.update_supported_types(supported_types, element='Oxygen', atomtype='O_31', flag=use_implicit_Hs_on_Oxygen)
    supported_types = tf.update_supported_types(supported_types, element='Nitrogen', atomtype='N_R1', flag=use_implicit_Hs_on_Nitrogen)
    supported_types = tf.update_supported_types(supported_types, element='Nitrogen', atomtype='N_21Mid', flag=use_implicit_Hs_on_Nitrogen)
    supported_types = tf.update_supported_types(supported_types, element='Nitrogen', atomtype='N_21', flag=use_implicit_Hs_on_Nitrogen)
    supported_types = tf.update_supported_types(supported_types, element='Nitrogen', atomtype='N_33', flag=use_implicit_Hs_on_Nitrogen)
    supported_types = tf.update_supported_types(supported_types, element='Nitrogen', atomtype='N_32', flag=use_implicit_Hs_on_Nitrogen)
    supported_types = tf.update_supported_types(supported_types, element='Nitrogen', atomtype='N_31', flag=use_implicit_Hs_on_Nitrogen)

    
    
    ##########################################################################
    # Set Assumed DREIDING atomtypes if code fails to find correct atom-type #
    ##########################################################################
    assumed = {'Sp1-Carbon':   'C_1',
               'Sp2-Carbon':   'C_2',
               'Sp3-Carbon':   'C_3',
               'Sp1-Oxygen':   'O_1',
               'Sp2-Oxygen':   'O_2',
               'Sp3-Oxygen':   'O_3',
               'Sp1-Nitrogen': 'N_1',
               'Sp2-Nitrogen': 'N_2',}


    ########################################################
    # Opening files to dump assumed/failed atom types into #
    ########################################################
    a = open(basename+'_assumed.txt', 'w')
    f = open(basename+'_failed.txt', 'w')
    
    # Function to count number of difference of IDs (meant for ringIDs and fusedringID to DRIEDING's _RB types)
    def count_diff_of_IDs(ID, IDs):
        count = 0
        if ID in IDs:
            for i in IDs:
                if i != ID: count += 1
        return count

    ########################################################    
    # Loop through mm.atoms and start assigning atom-types #
    ########################################################
    for i in mm.atoms:
        # Set atom instance and pull out desired info
        atom = mm.atoms[i]; formula = atom.molecule.formula;
        nb = int(atom.nb); ring = int(atom.ring); element = atom.element;
        ringID = atom.ringID; ringformula = atom.ringformula;
        hybridization = atom.hybridization; fusedringID = atom.fusedringID;
        
        
        # Find 1st-neighbor info (denoted by nameN; where N=neighbor depth)
        elements1 = tf.neigh_extract(atom, depth=1, info='element') # Example: ['C', 'C', 'N']
        rings1 = tf.neigh_extract(atom, depth=1, info='ring') # Example: [6, 6, 0]
        nbs1 = tf.neigh_extract(atom, depth=1, info='nb') # Example: [3, 3, 3]
        ringIDs1 = tf.get_neigh_ringIDs(i, mm, depth=1, criteria={'minringsize':1}) # Example: [1, 2, 2] # 1st neighs belong to ringIDs 1, 2, 2 
        fusedringIDs1 = [mm.atoms[j].fusedringID for j in atom.neighbor_ids[1] if mm.atoms[j].ring != 0] # first neigh fusedRingIDs if neigh is in ring
        hybrid1 = [mm.atoms[j].hybridization for j in mm.atoms[i].neighbor_ids[1]] # Example: ['Sp2', 'Sp2', 'Sp3']
        
        # Find 2nd-neighbor info (denoted by nameN; where N=neighbor depth)
        #elements2 = tf.neigh_extract(atom, depth=2, info='element') # Example: ['C', 'C', 'N'], if all H's molecule edge found
        #rings2 = tf.neigh_extract(atom, depth=2, info='ring') # Example: [6, 6, 0]
        #nbs2 = tf.neigh_extract(atom, depth=2, info='nb') # Example: [3, 3, 3]
        
        # Find 3rd-neighbor info (denoted by nameN; where N=neighbor depth)
        elements3 = tf.neigh_extract(atom, depth=3, info='element') # Example: ['H', 'H', 'H'], if all H's molecule edge found
        #rings3 = tf.neigh_extract(atom, depth=3, info='ring') # Example: [3, 3, 3]

        
        
        # Set intial .nta and .nta_comments and update later on if found
        atom.nta_type = '{}-type-yourself'.format(element)
        atom.nta_info = 'FAILED TO BE TYPED:  element: {}, ring: {}, nb: {}'.format(element, ring, nb)
        #print(i, ringID, ringIDs1, 'fused:  ', fusedringID, fusedringIDs1)


        ######################################################################################
        # Sp1 Carbon atom typing (ordering of nested if/elif/else statements set precedence) #
        ######################################################################################
        if element == 'C' and hybridization == 'Sp1':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#           
                
            #---------------------------------------------------------------------------------#
            # Strict DREIDING atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # C_11        13.019000   C         1            Sp1 hybridized Carbon with 1-implicit hydrogens
            if nb == 1 and use_implicit_Hs_on_Carbon:
                tally['found'] += 1
                atom.nta_type = 'C_11'
                atom.nta_info = 'Correctly found'
                
            # C_1         12.011000   C         2            Sp1 hybridized Carbon with 0-implicit hydrogens
            elif nb <= 2:
                tally['found'] += 1
                atom.nta_type = 'C_1'
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict DREIDING atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                tag = 'Sp1-Carbon'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                
                
        ######################################################################################
        # Sp2 Carbon atom typing (ordering of nested if/elif/else statements set precedence) #
        ######################################################################################
        elif element == 'C' and hybridization == 'Sp2':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#           
                
            #---------------------------------------------------------------------------------#
            # Strict DREIDING atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # C_RB        12.011000   C         3            Sp2 aromatic Carbon with 0-implicit hydrogens "B"ridging bond between two aromatic rings like in propene (will assign torsions differently)
            if nb >= 2 and ring >= 5 and ringID != 0 and count_diff_of_IDs(ringID, ringIDs1) == 1 and count_diff_of_IDs(fusedringID, fusedringIDs1) == 1:
                tally['found'] += 1
                atom.nta_type = 'C_RB'
                atom.nta_info = 'Correctly found'
                
            # C_R         12.011000   C         3            Sp2 aromatic Carbon with 0-implicit hydrogens
            elif nb >= 2 and ring >= 5:
                tally['found'] += 1
                atom.nta_type = 'C_R'
                atom.nta_info = 'Correctly found'
                
            # C_21Mid     13.019000   C         2            Sp2 non-aromatic Carbon with 1-implicit hydrogens in the "middle" of a molecule bonded to an Sp3 atom (uniquely named by Josh to help all2lmp set correct torsions)
            elif nb <= 2 and use_implicit_Hs_on_Carbon and 'Sp3' in hybrid1 and len(elements3) == 0: # if no 3rd-neighbors it must be at the center like in propene (CH2-CH1-CH3) or acetate anion (CH3-C-OO)
                tally['found'] += 1
                atom.nta_type = 'C_21Mid'
                atom.nta_info = 'Correctly found'
                
            # C_2Mid      12.011000   C         3            Sp2 non-aromatic Carbon with 0-implicit hydrogens in the "middle" of a molecule bonded to an Sp3 atom (uniquely named by Josh to help all2lmp set correct torsions)
            elif nb == 3 and 'Sp3' in hybrid1 and len(elements3) == 0: # if no 3rd-neighbors it must be at the center like in propene (CH2-CH1-CH3) or acetate anion (CH3-C-OO)
                tally['found'] += 1
                atom.nta_type = 'C_2Mid'
                atom.nta_info = 'Correctly found'
            
            # C_21        13.019000   C         2            Sp2 non-aromatic Carbon with 1-implicit hydrogens
            elif nb <= 2 and use_implicit_Hs_on_Carbon:
                tally['found'] += 1
                atom.nta_type = 'C_21'
                atom.nta_info = 'Correctly found'
                
            # C_2         12.011000   C         3            Sp2 non-aromatic Carbon with 0-implicit hydrogens
            elif nb == 3:
                tally['found'] += 1
                atom.nta_type = 'C_2'
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict DREIDING atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                tag = 'Sp2-Carbon'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                
                
        ######################################################################################
        # Sp3 Carbon atom typing (ordering of nested if/elif/else statements set precedence) #
        ######################################################################################
        elif element == 'C' and hybridization == 'Sp3':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#           
                
            #---------------------------------------------------------------------------------#
            # Strict DREIDING atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # C_33        15.035000   C         1            Sp3 hybridized Carbon with 3-implicit hydrogens
            if nb == 1 and use_implicit_Hs_on_Carbon:
                tally['found'] += 1
                atom.nta_type = 'C_33'
                atom.nta_info = 'Correctly found'
                
            # C_32        14.027000   C         2            Sp3 hybridized Carbon with 2-implicit hydrogens
            elif nb == 2 and use_implicit_Hs_on_Carbon:
                tally['found'] += 1
                atom.nta_type = 'C_32'
                atom.nta_info = 'Correctly found'
                
            # C_31        13.019000   C         3            Sp3 hybridized Carbon with 1-implicit hydrogens
            elif nb == 3 and use_implicit_Hs_on_Carbon:
                tally['found'] += 1
                atom.nta_type = 'C_31'
                atom.nta_info = 'Correctly found'
                
            # C_3         12.011000   C         4            Sp3 hybridized Carbon with 0-implicit hydrogens
            elif nb == 4:
                tally['found'] += 1
                atom.nta_type = 'C_3'
                atom.nta_info = 'Correctly found'
                

            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict DREIDING atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                tag = 'Sp3-Carbon'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                
                
        ####################################################################################
        # Hydrogen atom typing (ordering of nested if/elif/else statements set precedence) #
        ####################################################################################
        elif element == 'H':
            # H_           1.008000   H         1            Generic hydrogen
            atom.nta_type = 'H_'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
            
        ######################################################################################
        # Sp1 Oxygen atom typing (ordering of nested if/elif/else statements set precedence) #
        ######################################################################################
        elif element == 'O' and hybridization == 'Sp1':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #---------------------------------------------------------------------------------#
            # Strict DREIDING atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # Sp1 hybridized Oxygen with 0-implicit hydrogens
            if nb == 1:
                tally['found'] += 1
                atom.nta_type = 'O_1'
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict DREIDING atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                tag = 'Sp1-Oxygen'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                
                
        ######################################################################################
        # Sp2 Oxygen atom typing (ordering of nested if/elif/else statements set precedence) #
        ######################################################################################
        elif element == 'O' and hybridization == 'Sp2':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #---------------------------------------------------------------------------------#
            # Strict DREIDING atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # O_R         15.999000   0         2            Sp2 aromatic Oxygen with 0-implicit hydrogens
            if ring >= 5 and nb == 2:
                tally['found'] += 1
                atom.nta_type = 'O_R'
                atom.nta_info = 'Correctly found'
                
            # O_R1        17.007000   0         1            Sp2 aromatic Oxygen with 1-implicit hydrogens
            elif ring >= 5 and nb == 1 and use_implicit_Hs_on_Oxygen:
                tally['found'] += 1
                atom.nta_type = 'O_R1'
                atom.nta_info = 'Correctly found'
                
            # O_2         15.999000   0         2            Sp2 non-aromatic Oxygen with 0-implicit hydrogens
            elif nb == 2:
                tally['found'] += 1
                atom.nta_type = 'O_2'
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict DREIDING atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                tag = 'Sp2-Oxygen'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                
                
        ######################################################################################
        # Sp3 Oxygen atom typing (ordering of nested if/elif/else statements set precedence) #
        ######################################################################################
        elif element == 'O' and hybridization == 'Sp3':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #---------------------------------------------------------------------------------#
            # Strict DREIDING atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # O_31        17.007000   0         1            Sp3 hybridized Oxygen with 1-implicit hydrogens
            if nb == 1 and use_implicit_Hs_on_Oxygen:
                tally['found'] += 1
                atom.nta_type = 'O_31'
                atom.nta_info = 'Correctly found'
                
            # O_3         15.999000   0         2            Sp3 hybridized Oxygen with 0-implicit hydrogens
            elif nb == 2:
                tally['found'] += 1
                atom.nta_type = 'O_3'
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict DREIDING atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                tag = 'Sp3-Oxygen'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                
                
        ########################################################################################
        # Sp1 Nitrogen atom typing (ordering of nested if/elif/else statements set precedence) #
        ########################################################################################
        elif element == 'N' and hybridization == 'Sp1':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY
            
            #---------------------------------------------------------------------------------#
            # Strict DREIDING atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # N_1         14.007000   N         2            Sp1 hybridized Nitrogen with 0-implicit hydrogens
            if nb == 2:
                tally['found'] += 1
                atom.nta_type = 'N_1'
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict DREIDING atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                tag = 'Sp1-Nitrogen'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                
                
        ########################################################################################
        # Sp2 Nitrogen atom typing (ordering of nested if/elif/else statements set precedence) #
        ########################################################################################
        elif element == 'N' and hybridization == 'Sp2':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#           
                
            #---------------------------------------------------------------------------------#
            # Strict DREIDING atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#                
            # N_R         14.007000   N         3            Sp2 aromatic Nitrogen with 0-implicit hydrogens
            if nb == 2 and ring >= 5:
                tally['found'] += 1
                atom.nta_type = 'N_R'
                atom.nta_info = 'Correctly found'
                
            # N_R1        15.015000   N         2            Sp2 aromatic Nitrogen with 1-implicit hydrogens
            elif nb == 1 and ring >= 5 and use_implicit_Hs_on_Nitrogen:
                tally['found'] += 1
                atom.nta_type = 'N_R'
                atom.nta_info = 'Correctly found'
                
            # N_21Mid     15.015000   N         2            Sp2 non-aromatic Nitrogen with 1-implicit hydrogens in the "middle" of a molecule bonded to an Sp3 atom (uniquely named by Josh to help all2lmp set correct torsions)
            elif nb <= 1 and use_implicit_Hs_on_Nitrogen and 'Sp3' in hybrid1 and len(elements3) == 0: # if no 3rd-neighbors it must be at the center like in propene (CH2-CH1-CH3) or acetate anion (CH3-C-OO)
                tally['found'] += 1
                atom.nta_type = 'N_21Mid'
                atom.nta_info = 'Correctly found'
                
            # N_2Mid      14.007000   N         3            Sp2 non-aromatic Nitrogen with 0-implicit hydrogens in the "middle" of a molecule bonded to an Sp3 atom (uniquely named by Josh to help all2lmp set correct torsions)
            elif nb == 3 and 'Sp3' in hybrid1 and len(elements3) == 0: # if no 3rd-neighbors it must be at the center like in propene (CH2-CH1-CH3) or acetate anion (CH3-C-OO)
                tally['found'] += 1
                atom.nta_type = 'N_2Mid'
                atom.nta_info = 'Correctly found'
            
            # N_21        15.015000   N         2            Sp2 non-aromatic Nitrogen with 1-implicit hydrogens
            elif nb <= 1 and use_implicit_Hs_on_Nitrogen:
                tally['found'] += 1
                atom.nta_type = 'N_21'
                atom.nta_info = 'Correctly found'
                
            # N_2         14.007000   N         3            Sp2 non-aromatic Nitrogen with 0-implicit hydrogens
            elif nb == 2:
                tally['found'] += 1
                atom.nta_type = 'N_2'
                atom.nta_info = 'Correctly found'
                
            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict DREIDING atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                tag = 'Sp2-Nitrogen'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                
                
        ########################################################################################
        # Sp3 Nitrogen atom typing (ordering of nested if/elif/else statements set precedence) #
        ########################################################################################
        elif element == 'N' and hybridization == 'Sp3':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#           
                
            #---------------------------------------------------------------------------------#
            # Strict DREIDING atom-typing that occurs after User defined atom-typing attempts #
            #---------------------------------------------------------------------------------#
            # N_33        17.031000   N         1            Sp3 hybridized Nitrogen with 3-implicit hydrogens
            if nb == 0 and use_implicit_Hs_on_Nitrogen:
                tally['found'] += 1
                atom.nta_type = 'N_33'
                atom.nta_info = 'Correctly found'
                
            # N_32        16.023000   N         2            Sp3 hybridized Nitrogen with 2-implicit hydrogens
            elif nb == 1 and use_implicit_Hs_on_Nitrogen:
                tally['found'] += 1
                atom.nta_type = 'N_32'
                atom.nta_info = 'Correctly found'
                
            # N_31        15.015000   N         3            Sp3 hybridized Nitrogen with 1-implicit hydrogens
            elif nb == 2 and use_implicit_Hs_on_Nitrogen:
                tally['found'] += 1
                atom.nta_type = 'N_31'
                atom.nta_info = 'Correctly found'
                
            # N_3         14.007000   N         4            Sp3 hybridized Nitrogen with 0-implicit hydrogens
            elif nb == 3:
                tally['found'] += 1
                atom.nta_type = 'N_3'
                atom.nta_info = 'Correctly found'
                

            #----------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict DREIDING atom-typing fails used assumed #
            #----------------------------------------------------------------------------------#
            else:
                tag = 'Sp3-Nitrogen'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                
        ##################################################################################
        # Sulfur atom typing (ordering of nested if/elif/else statements set precedence) #
        ##################################################################################
        elif element == 'S':
            # S_3         32.060000   S         2            Sulfur, e.g. in thiols, thioethers
            atom.nta_type = 'S_3'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ####################################################################################
        # Silicone atom typing (ordering of nested if/elif/else statements set precedence) #
        ####################################################################################
        elif element == 'Si':
            # Si3         28.085000   Si        4            Sp3 Silicon, e.g. in silanes
            atom.nta_type = 'Si3'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ####################################################################################
        # Fluorine atom typing (ordering of nested if/elif/else statements set precedence) #
        ####################################################################################
        elif element == 'F':
            # F_          18.998400   F         1            Fluorine
            atom.nta_type = 'F_'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ####################################################################################
        # Chlorine atom typing (ordering of nested if/elif/else statements set precedence) #
        ####################################################################################
        elif element == 'Cl':
            # Cl          35.450000   Cl        1            Chlorine
            atom.nta_type = 'Cl'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ###################################################################################
        # Calcium atom typing (ordering of nested if/elif/else statements set precedence) #
        ###################################################################################
        elif element == 'Cl':
            # Ca          40.078000   Ca        1            Calcium(++) ion
            atom.nta_type = 'Ca'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        ###################################################################################
        # Bromine atom typing (ordering of nested if/elif/else statements set precedence) #
        ###################################################################################
        elif element == 'Cl':
            # Br          79.904000   Br        1            Bromine
            atom.nta_type = 'Br'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            
        #######################################################################################
        # Phosphorous atom typing (ordering of nested if/elif/else statements set precedence) #
        #######################################################################################
        elif element == 'Cl':
            # P_3         30.973700   P         3            Sp3 Phosphorous, e.g. in phosphanes
            atom.nta_type = 'P_3'; tally['found'] += 1;
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