# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 18th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

Atom typing rules for Clay-FF with user-defined assumed types
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


#####################################################################
# Function to perform atom-typing for Clay-FF forcefield atom types #
#####################################################################
def nta(mm, basename, ff_name):
    tally = {'found':0, 'assumed':0, 'failed':0} # To tally findings
    
    ###############################################
    # Set Clay-FF specific informations and flags #
    ###############################################
    # Log supported types
    supported_types = {'Hydrogen':    ['h* (Q)', 'ho (Q)'],
                       
                       'Oxygen':      ['o* (Q)', 'oh- (Q)', 'oh (Q)', 'obts (Q)', 'obos (Q)',
                                       'ob (Q)', 'obss (Q)'],
                       
                       'Silicone':    ['st (Q)'],
                       'Chlorine':    ['Cl (Q)'],
                       'Calcium':     ['cah (Q)', 'cao (Q)', 'Ca (Q)'],
                       'Aluminum':    ['ao (Q)', 'at (Q)'],
                       'Magnesium':   ['mgh (Q)', 'mgo (Q)', 'Mg (Q)'],
                       'Lithium ':    ['lio (Q)'],
                       'Iron':        ['feo (Q)'],
                       'Sodium':      ['Na (Q)'],
                       'Potassium':   ['K (Q)'],
                       'Cesium':      ['Cs (Q)'],
                       'Barium':      ['Ba (Q)'],
                       'Strontium':   ['Sr (Q)'],
                       'Lead':        ['Pb (Q)'],
                       }
    
    # Clay-FF flags to help set correct atom-types (True to use False not to use)

    
    # Update supported types lists if flags are used and set flag status in parentheses
    reset_charge_w_FF_specifc_charge = True # Will reset all charges to FF specific charges (True or False)


    
    
    #########################################################################
    # Set Assumed Clay-FF atomtypes if code fails to find correct atom-type #
    #########################################################################
    #          Key to access assumed type     (type, charge)
    assumed = {'Hydrogen':                     ('ho', 0.425),    # 1-connect : Hydrogen
               'Oxygen':                       ('ob', -0.105),   # Any-connect : Oxygen
               'Aluminium':                    ('ao', 1.575),    # Any-connect: Aluminium
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
        rings3 = tf.neigh_extract(atom, depth=3, info='ring') # Example: [3, 3, 3]
        
        
        # Set intial .nta and .nta_comments and update later on if found
        atom.nta_type = '{}-type-yourself'.format(element)
        atom.nta_info = 'FAILED TO BE TYPED:  element: {}, ring: {}, nb: {}'.format(element, ring, nb)

        
        
        ####################################################################################
        # Hydrogen atom typing (ordering of nested if/elif/else statements set precedence) #
        ####################################################################################
        if element == 'H':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY

                
            #--------------------------------------------------------------------------------#
            # Strict Clay-FF atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # h*     1.007970     H          1        spc water H          0.41
            if elements1[0] == 'O' and formula == 'H2-O1':
                atom.nta_type = 'h*'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = 0.41 # reset H-element charge
                
            # ho     1.007970     H          1        hydroxyl H           0.425
            elif 'O' in elements1:
                atom.nta_type = 'ho'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = 0.425 # reset H-element charge
                
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict Clay-FF atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
            else:
                tag = 'Hydrogen'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag][0]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = assumed[tag][1] # reset H-element charge


        ##################################################################################
        # Oxygen atom typing (ordering of nested if/elif/else statements set precedence) #
        ##################################################################################
        if element == 'O':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY

                
            #--------------------------------------------------------------------------------#
            # Strict Clay-FF atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # o*         15.99940     O          2        spc water O               -0.82
            if formula == 'H2-O1':
                atom.nta_type = 'o*'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = -0.82 # reset O-element charge
                    
                
            # oh-        15.99940     O          1        hydroxide O 
            elif nb == 1:
                atom.nta_type = 'oh-'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = -0.95 # reset O-element charge (using oh charge)
                    
            # oh         15.99940     O          2        hydroxyl O              -0.95 
            # ohs       15.99940      O          2        hydroxyl O sub.         -1.0808
            elif ring == 0 and nb == 2 and elements1.count('H') == 1:              
                atom.nta_type = 'oh'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = -0.95 # reset O-element charge
              
            # obts     15.99940     O          2        oxygen tet. sub.        -1.1688
            elif nb == 2 and 4 in nbs1:
                atom.nta_type = 'obts'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = -1.1688 # reset O-element charge
                    
            # obos     15.99940    O          2        oxygen oct. sub.        -1.1808
            elif nb == 2 and 6 in nbs1:
                atom.nta_type = 'obos'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = -1.1808 # reset O-element charge
                    
            # ob        15.99940     O          2        bridging oxygen         -1.05
            elif nb == 2:
                atom.nta_type = 'ob'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = -1.05 # reset O-element charge
                    
            # obss     15.99940     O          3        oxygen double sub.  -1.2996
            elif nb == 2:
                atom.nta_type = 'obss'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = -1.2996 # reset O-element charge
                    
                    
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict Clay-FF atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
            else:
                tag = 'Oxygen'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag][0]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = assumed[tag][1] # reset O-element charge
                    
                    
        ####################################################################################
        # Silicone atom typing (ordering of nested if/elif/else statements set precedence) #
        ####################################################################################
        elif element == 'Si':
            # st     28.08550    Si          4        tetrahedral silicon  2.1
            atom.nta_type = 'st'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            if reset_charge_w_FF_specifc_charge:
                atom.charge = 2.1 # reset Si-element charge
                    
                
        ########################################################################
        # Chlorine (ordering of nested if/elif/else statements set precedence) #
        ########################################################################
        elif element == 'Cl':
            # Cl     35.45300    Cl          0        chloride ion        -1.0
            atom.nta_type = 'Cl'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            if reset_charge_w_FF_specifc_charge:
                atom.charge = -1.0 # reset Cl-element charge
                
        #######################################################################
        # Calcium (ordering of nested if/elif/else statements set precedence) #
        #######################################################################
        elif element == 'Ca':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY

                
            #--------------------------------------------------------------------------------#
            # Strict Clay-FF atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # cah     40.08000    Ca          6        hydroxide calcium    1.05
            # for non-terminal
            if nb == 6 and elements1.count('O') == 6 and elements2.count('H') == 6:
                atom.nta_type = 'cah'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = 1.05 # reset Ca-element charge
            # terminal
            elif nb == 2 and elements1.count('O') == 2 and elements2.count('H') == 2:
                atom.nta_type = 'cah'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = 1.05 # reset Ca-element charge
                    
            # cao     40.08000    Ca          6        octahedral calcium
            elif nb == 6:
                atom.nta_type = 'cao'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = 1.05 # reset Ca-element charge (Using cah charge)
                    
            # Ca     40.07980     Ca           0        calcium ion                2.0 
            elif nb == 0:
                atom.nta_type = 'Ca'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = 2.0 # reset Ca-element charge
                    
                    
        #########################################################################
        # Aluminium (ordering of nested if/elif/else statements set precedence) #
        #########################################################################
        elif element == 'Al':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY

                
            #--------------------------------------------------------------------------------#
            # Strict Clay-FF atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # ao    26.98154    Al          6        octahedral aluminum  1.575  
            if nb == 6:
                atom.nta_type = 'ao'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = 1.575 # reset Al-element charge
                    
            # at     26.98154    Al          4        tetrahedral aluminum 1.575    
            elif nb == 4:
                atom.nta_type = 'at'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = 1.575 # reset Al-element charge
                    
            #---------------------------------------------------------------------------------#
            # If user defined attempt fails and Strict Clay-FF atom-typing fails used assumed #
            #---------------------------------------------------------------------------------#
            else:
                tag = 'Aluminium'; tally['assumed'] += 1;
                atom.nta_type = assumed[tag][0]; atom.nta_info = 'ASSUMED: {}'.format(assumed[tag]);
                tf.write2assumed(a, tag, assumed[tag], atom, ff_name, i)
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = assumed[tag][1] # reset O-element charge
                    
                    
        #########################################################################
        # Magnesium (ordering of nested if/elif/else statements set precedence) #
        #########################################################################
        elif element == 'Mg':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#
            # PASS NO USER DEFINED TYPES CURRENTLY

                
            #--------------------------------------------------------------------------------#
            # Strict Clay-FF atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # mgh     24.30500    Mg          6        hydroxide magnesium  1.05
            # for non-terminal
            if nb == 6 and elements1.count('O') == 6 and elements2.count('H') == 6:
                atom.nta_type = 'mgh'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = 1.05 # reset Ca-element charge
            # terminal
            elif nb == 2 and elements1.count('O') == 2 and elements2.count('H') == 2:
                atom.nta_type = 'mgh'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = 1.05 # reset Ca-element charge
                    
            # mgo     24.30500    Mg          6        octahedral magnesium 1.36 
            elif nb == 6:
                atom.nta_type = 'mgo'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = 1.36 # reset Ca-element charge
                    
            # Mg       24.3050      Mg          0        magnesium ion               2.0   
            elif nb == 0:
                atom.nta_type = 'mgo'; tally['found'] += 1;
                atom.nta_info = 'Correctly found'
                if reset_charge_w_FF_specifc_charge:
                    atom.charge = 2.0 # reset Ca-element charge
                    
                    
        ####################################################################
        # Iron (ordering of nested if/elif/else statements set precedence) #
        ####################################################################
        elif element == 'Fe':
            # feo     55.84700    Fe          6        octahedral iron      1.575
            atom.nta_type = 'feo'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            if reset_charge_w_FF_specifc_charge:
                atom.charge = 1.575 # reset Fe-element charge
                    
                    
        #######################################################################
        # lithium (ordering of nested if/elif/else statements set precedence) #
        #######################################################################
        elif element == 'Li':
            # lio     6.941000   Li          6        octahedral lithium   0.525
            atom.nta_type = 'lio'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            if reset_charge_w_FF_specifc_charge:
                atom.charge = 0.525 # reset Li-element charge


        ######################################################################
        # Sodium (ordering of nested if/elif/else statements set precedence) #
        ######################################################################
        elif element == 'Na':
            # Na     22.99000    Na          0        sodium ion           1.0
            atom.nta_type = 'Na'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            if reset_charge_w_FF_specifc_charge:
                atom.charge = 1.0 # reset Na-element charge
                
                
        #########################################################################
        # Potassium (ordering of nested if/elif/else statements set precedence) #
        #########################################################################
        elif element == 'K':
            # K     39.10        K          0        potassium ion        1.0
            atom.nta_type = 'K'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            if reset_charge_w_FF_specifc_charge:
                atom.charge = 1.0 # reset K-element charge
                

        ######################################################################
        # Cesium (ordering of nested if/elif/else statements set precedence) #
        ######################################################################
        elif element == 'Cs':
            # Cs     132.9100    Cs          0        cesium ion           1.0
            atom.nta_type = 'Cs'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            if reset_charge_w_FF_specifc_charge:
                atom.charge = 1.0 # reset Cs-element charge
                

        ######################################################################
        # Barium (ordering of nested if/elif/else statements set precedence) #
        ######################################################################
        elif element == 'Ba':
            # Ba     137.3300    Ba          0        barium ion           2.0
            atom.nta_type = 'Ba'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            if reset_charge_w_FF_specifc_charge:
                atom.charge = 2.0 # reset Ba-element charge
                

        #########################################################################
        # Strontium (ordering of nested if/elif/else statements set precedence) #
        #########################################################################
        elif element == 'Sr':
            # Sr      87.6200    Sr          0        strontium ion        2.0
            atom.nta_type = 'Sr'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            if reset_charge_w_FF_specifc_charge:
                atom.charge = 2.0 # reset Sr-element charge
                

        ####################################################################
        # Lead (ordering of nested if/elif/else statements set precedence) #
        ####################################################################
        elif element == 'Pb':
            # Pb     207.2000    Pb          0        lead ion             2.0
            atom.nta_type = 'Pb'; tally['found'] += 1;
            atom.nta_info = 'Correctly found'
            if reset_charge_w_FF_specifc_charge:
                atom.charge = 2.0 # reset Sr-element charge
                
                    
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