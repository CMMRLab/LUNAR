
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
September 19th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to check user inputs and let them know what input is incorrect
def safety(topofile, bondfile, parent_dir, newfile, ff_name, delete_atoms, mass_map, bondorder,
           maxbonded, boundary, vdw_radius_scale, reset_charges, print_options, log):
    
    # Intialize as True and update as False as needed
    return_boolean = True
    
    # Check if topofile is valid. If passed is empty switch return_boolean and print error
    topo_exts = ['.mol', '.data', '.dat', 'sdf', '.mol2', '.smiles', '.pdb'] # valid extensions
    passed = [i for i in topo_exts if i in topofile]
    if not passed:
        return_boolean = False
        log.error('{} {} {} {}'.format('ERROR topofile: ', topofile, 'not supported. Currently supported file extensions: ', '  '.join(topo_exts)))
                
    # Check if frc_file is valid. If passed is empty switch return_boolean and print error
    ff_names = ['PCFF-IFF', 'PCFF', 'compass', 'CVFF-IFF', 'CVFF', 'Clay-FF', 'DREIDING', 'OPLS-AA', 'general:0', 'general:1', 'general:2', 'general:3', 'general:4'] # valid ff_names
    passed = [i for i in ff_names if i in ff_name]
    if not passed:
        return_boolean = False
        log.error('{} {} {} {}'.format('ERROR ff_name: ', ff_name, 'not supported. Currently supported force field names: ', '  '.join(ff_names)))
        
    # Check if reset_charges is True or False
    if reset_charges not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR reset_charges: ', reset_charges, 'is not True or False, update to either'))
        
    # Check if print_options is True or False
    if print_options not in [True, False]:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR print_options: ', print_options, 'is not True or False, update to either'))
        
    # Check delete_atoms has proper keys and values
    if 'method' not in delete_atoms:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR delete_atoms: ', delete_atoms, 'is missing "method" key'))
    if 'criteria' not in delete_atoms:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR delete_atoms: ', delete_atoms, 'is missing "criteria" key'))
    if delete_atoms['method'] not in ['mass', 'size']:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR delete_atoms: ', delete_atoms, '"method" is not "mass" or "size"'))
    try: 
        float(delete_atoms['criteria'])
    except:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR delete_atoms: ', delete_atoms, '"criteria" value is not a float or int value'))
        
    # Check vdw_radius_scale
    try: 
        float(vdw_radius_scale)
    except:
        return_boolean = False
        log.error('{} {} {}'.format('ERROR vdw_radius_scale: ', vdw_radius_scale, 'value is not a float or int value'))
        
    # check boundary conditions
    pflags = boundary.split() # split boundary
    count = pflags.count('f') + pflags.count('p')
    if len(pflags) != 3 and count != 3:
        return_boolean = False
        log.error('ERROR boundary does not contain 3-dimensions or boundary flags are not f or p')
    return return_boolean
