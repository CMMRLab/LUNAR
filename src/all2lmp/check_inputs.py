#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
October 4th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to check user inputs and let them know what input is incorrect
def safety(topofile, nta_file, frc_file, assumed, parent_dir, newfile, atom_style, ff_class, use_auto_equivalence, use_assumed_auto_fill, reset_molids,
           reset_charges, write_txt_comments, write_bond_react, print_options, use_morse_bonds, include_type_labels, log):
    
    # Intialize as True and update as False as needed
    return_boolean = True
    
    # Check if topofile is valid. If passed is empty switch return_boolean and print error
    topo_exts = ['.mol', '.data', '.dat', '.mdf', 'sdf', '.mol2', '.pdb'] # valid extensions
    passed = [i for i in topo_exts if i in topofile]
    if not passed:
        return_boolean = False
        log.error('{} {} {} {}'.format('ERROR topofile: ', topofile, 'not supported. Currently supported file extensions: ', '  '.join(topo_exts)))
        

    # Check if nta_file is valid. If passed is empty switch return_boolean and print error
    nta_exts = ['.nta', '.car', 'topofile', 'PCFF-IFF', 'PCFF', 'compass', 'CVFF-IFF', 'CVFF', 'Clay-FF', 'DREIDING', 'OPLS-AA'] # valid extensions
    passed = [i for i in nta_exts if i in nta_file]
    if not passed:
        return_boolean = False
        log.error('{} {} {} {}'.format('ERROR nta_file: ', nta_file, 'not supported. Currently supported file extensions: ', '  '.join(nta_exts)))
        
        
    # Check if frc_file is valid. If passed is empty switch return_boolean and print error
    frc_exts = ['.frc'] # valid extensions
    passed = [i for i in frc_exts if i in frc_file]
    if not passed:
        return_boolean = False
        log.error('{} {} {} {}'.format('ERROR frc_file: ', frc_file, 'not supported. Currently supported file extensions: ', '  '.join(frc_exts)))


    # Check if atom style is supported
    atomstyles = ['full', 'charge', 'molecular', 'angle', 'bond', 'atomic', 'dipole', 'dpd', 'line']
    if atom_style not in atomstyles:
        return_boolean = False
        log.error('{} {} {} {}'.format('ERROR atom_style: ', atom_style, 'not supported. Currently supported atom styles: ', '  '.join(atomstyles)))
        
    # Check if FF class is supported
    ff_classes = [0, 1, 2, 'i', 'ilmp', 'd', 's1', 's2']
    if ff_class not in ff_classes:
        return_boolean = False
        log.error('{} {} {} {}'.format('ERROR ff_class: ', ff_class, 'not supported. Currently supported FF classes: ', '  '.join([str(i) for i in ff_classes])))
    
    # Check if use_auto_equivalence is True or False
    if use_auto_equivalence not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR use_auto_equivalence: ', use_auto_equivalence, 'is not True or False, update to either'))
        
    # Check if use_assumed_auto_fill is True or False
    if use_assumed_auto_fill not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR use_assumed_auto_fill: ', use_assumed_auto_fill, 'is not True or False, update to either'))
    
    # Check if reset_molids is True or False
    if reset_molids not in [True, False]:
        integer_input = False
        try:
            if reset_molids.isdigit():
                integer_input = True
        except: pass
        if not integer_input:
            return_boolean = False
            log.error('{} {} {}'.format('ERROR reset_molids: ', reset_molids, 'is not True or False or an integer, update to either'))

    # Check if reset_charges is True or False
    if reset_charges not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR reset_charges: ', reset_charges, 'is not True or False, update to either'))
        
    # Check if write_bond_react is True or False
    if write_bond_react not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR write_bond_react: ', write_bond_react, 'is not True or False, update to either'))
        
    # Check if print_options is True or False
    if print_options not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR print_options: ', print_options, 'is not True or False, update to either'))
        
    # Check if use_auto_equivalence is True or False
    if use_morse_bonds not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR use_morse_bonds: ', use_morse_bonds, 'is not True or False, update to either'))
        
    # Check if include_type_labels is True or False
    if include_type_labels not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR include_type_labels: ', include_type_labels, 'is not True or False, update to either'))
        
    # Check if write_txt_comments is True or False
    if write_txt_comments not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR write_txt_comments: ', write_txt_comments, 'is not True or False, update to either'))
    return return_boolean
