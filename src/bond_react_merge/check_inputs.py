#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
February 1st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to check user inputs and let them know what input is incorrect
def safety(files, parent_dir, atom_style, generate_map_file, write_rxn_mol2files,
           write_rxn_datafiles, write_moleculefiles, print_options, map_near_edge_rxn_charges,
           commandline_inputs, log):
    
    # Intialize as True and update as False as needed
    return_boolean = True
    
    # Check if files is a dictionary and is not empty
    if not type(files) is dict:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR files variable: ', files, ' is not a python dictionary'))
    if len(files) == 0:
        return_boolean = False
        log.error('{}'.format('ERROR files dictionary is empty'))
        
    # Check if tags are supported in files dictionary
    supported_tags = ['pre', 'post', 'data', 'infile']
    for file in files:
        # strip last number of filename. Assumes tags either:
        # - data1, data2, ... dataN  ->  data
        # - pre1, pre2, ... preN     ->  pre
        # - post1, post2, ... postN  ->  post
        tmpname = ''.join([i for i in file if i.isalpha()])
        
        if tmpname not in supported_tags:
            return_boolean = False
            log.error('{} {} {}'.format('ERROR filetag: ', file, ' is not a supported filetag'))
            
    # Check if atom style is supported
    atomstyles = ['full', 'charge', 'molecular', 'angle', 'bond', 'atomic', 'dipole', 'dpd', 'line']
    if atom_style not in atomstyles:
        return_boolean = False
        log.error('{} {} {} {}'.format('ERROR atom_style: ', atom_style, 'not supported. Currently supported atom styles: ', '  '.join(atomstyles)))
        
    # Check if generate_map_file is True or False
    if generate_map_file not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR generate_map_file: ', generate_map_file, 'is not True or False, update to either'))
        
    # Check if map_near_edge_rxn_charges is False or an acceptable int
    acceptable_ints = [i for i in range(21)]; acceptable_bools = [False];
    acceptable_edge_inputs = acceptable_ints + acceptable_bools
    if map_near_edge_rxn_charges not in acceptable_edge_inputs:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR map_near_edge_rxn_charges: ', map_near_edge_rxn_charges, 'is not False or an integer value (likely string), update to either'))
        
    # Check if write_rxn_mol2files is True or False
    if write_rxn_mol2files not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR write_rxn_mol2files: ', write_rxn_mol2files, 'is not True or False, update to either'))
        
    # Check if write_rxn_datafiles is True or False
    if write_rxn_datafiles not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR write_rxn_datafiles: ', write_rxn_datafiles, 'is not True or False, update to either'))
        
    # Check if write_moleculefiles is True or False
    if write_moleculefiles not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR write_moleculefiles: ', write_moleculefiles, 'is not True or False, update to either'))
        
    # Check if print_options is True or False
    if print_options not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR print_options: ', print_options, 'is not True or False, update to either'))
    return return_boolean
