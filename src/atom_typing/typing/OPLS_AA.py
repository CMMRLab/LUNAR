# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 18th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

Atom typing rules for OPLS-AA with user-defined assumed types
found in the assumed dictionary. These are stored inside each
force field atom-typing script because not all forcefields 
need to have the exact same assumed atom types and some
force fields may not even contain basic elements that others
may.

PCFF-IFF also has PCFF-IFF flags to override the usage of standard
PCFF atom-types for better defined atom-types. Using these flags
will also automatically set the charges based on PCFF-IFF charges
found in the IFF DATABASE folder.
"""


##############################
# Import Necessary Libraries #
##############################
import src.atom_typing.typing.typing_functions as tf
from collections import OrderedDict
import os


#####################################################################
# Function to perform atom-typing for OPLS-AA forcefield atom types #
#####################################################################
def nta(mm, basename, ff_name):
    tally = {'found':0, 'assumed':0, 'failed':0} # To tally findings
    
    # print warning about OPLS-AA ff file from msi2lmp
    print('\n\n')
    print('Using OPLS-AA atom typing module. The atom-types set by this module are consistant')
    print('with the expeirmental and not fully supported OPLS-AA .frc file provided with msi2lmp.')
    print('The OPLS-AA extension for atom_typing.py and all2lmp.py have the same limited atom')
    print('types and functionality that the current version of msi2lmp (v3.9.9 - 2018-11-05) has.')
    print()
    print('If a full OPLS-AA .frc file where to be developed the addition of new atom-types will')
    print('have to be coded into OPLS_AA.py automatic typing code. The OPLS_AA.py automatic typing')
    print('code serves as a rudimentary attempt for OPLS-AA atom typing and should be treated as')
    print('such.')
    
    ###############################################
    # Set OPLS-AA specific informations and flags #
    ###############################################
    # Log supported types
    supported_types = {'Carbon':      ['CT', 'CG', 'CC'],
                       
                       'Hydrogen':    ['HM', 'HC', 'HO', 'HS'],
                       
                       'Oxygen':      ['OH', 'OC'],
                       
                       
                       'Sulfer':      ['SH', 'S'],
                       
                       'Chlorine':    ['CL'],

                       }
    

    # OPLS-AA flags to help set correct atom-types (True to use False not to use)

    
    
    #########################################################################
    # Set Assumed OPLS-AA atomtypes if code fails to find correct atom-type #
    #########################################################################
    assumed = {}


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

        
        
        ##################################################################################
        # Carbon atom typing (ordering of nested if/elif/else statements set precedence) #
        ##################################################################################
        if element == 'C':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#

                
            #--------------------------------------------------------------------------------#
            # Strict OPLS-AA atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # CT    12.011150    C           4        Aliphatic C
            if nb == 4:
                tally['found'] += 1
                atom.nta_type = 'CT'
                atom.nta_info = 'Correctly found'
                
            # CG    12.011150    C           1        Graphite C
            elif ring >= 5:
                tally['found'] += 1
                atom.nta_type = 'CG'
                atom.nta_info = 'Correctly found'
                
            # CC    12.011150    C           3        Carbonate ion C, AMBER
            elif ring == 0 and formula == 'C1-O3':
                tally['found'] += 1
                atom.nta_type = 'CC'
                atom.nta_info = 'Correctly found'
                
                
        ####################################################################################
        # Hydrogen atom typing (ordering of nested if/elif/else statements set precedence) #
        ####################################################################################
        if element == 'H':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#

                
            #--------------------------------------------------------------------------------#
            # Strict OPLS-AA atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # HM     1.007970    H           1        H(C), CH3OH
            if 'C' in elements1 and 'H' in elements2:
                tally['found'] += 1
                atom.nta_type = 'HM'
                atom.nta_info = 'Correctly found'
                
            # HC     1.007970    H           1        H, RH, alkanes
            elif 'C' in elements1:
                tally['found'] += 1
                atom.nta_type = 'HC'
                atom.nta_info = 'Correctly found'
                
            # HO     1.007970    H           1        H(O), ROH
            elif 'O' in elements1:
                tally['found'] += 1
                atom.nta_type = 'HO'
                atom.nta_info = 'Correctly found'

            # HS     1.007970    H           1        H(S), RSH
            elif 'S' in elements1:
                tally['found'] += 1
                atom.nta_type = 'HS'
                atom.nta_info = 'Correctly found'
                
                
        ##################################################################################
        # Oxygen atom typing (ordering of nested if/elif/else statements set precedence) #
        ##################################################################################
        if element == 'O':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#

                
            #--------------------------------------------------------------------------------#
            # Strict OPLS-AA atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # OH    15.999400    O           2        O, ROH
            if nb == 2 and 'C' in elements1 and 'H' in elements1:
                tally['found'] += 1
                atom.nta_type = 'OH'
                atom.nta_info = 'Correctly found'
                
            # OC    15.035060    O           1        Carbonate ion O, AMBER
            elif nb == 1 and formula == 'C1-O3':
                tally['found'] += 1
                atom.nta_type = 'OC'
                atom.nta_info = 'Correctly found'
                

        ##################################################################################
        # Sulfur atom typing (ordering of nested if/elif/else statements set precedence) #
        ##################################################################################
        if element == 'S':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#

                
            #--------------------------------------------------------------------------------#
            # Strict OPLS-AA atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # SH    32.064000    S           2        S, RSH
            if nb == 2 and 'C' in elements1 and 'H' in elements1:
                tally['found'] += 1
                atom.nta_type = 'SH'
                atom.nta_info = 'Correctly found'
                
            # S     32.064000    S           2        S, RSR
            elif nb == 2 and elements1.count('C') == 2:
                tally['found'] += 1
                atom.nta_type = 'S'
                atom.nta_info = 'Correctly found'


        ####################################################################################
        # Chlorine atom typing (ordering of nested if/elif/else statements set precedence) #
        ####################################################################################
        if element == 'Cl':
            #----------------------------------------------------------------------------#
            # User defined intial attempts (For ReaxFF with open valence polymerization) #
            #----------------------------------------------------------------------------#

                
            #--------------------------------------------------------------------------------#
            # Strict OPLS-AA atom-typing that occurs after User defined atom-typing attempts #
            #--------------------------------------------------------------------------------#
            # CL    35.453000    Cl          1        Cl, RCl
            if nb == 1 and 'C' in elements1:
                tally['found'] += 1
                atom.nta_type = 'CL'
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