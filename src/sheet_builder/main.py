# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
June 4th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.sheet_builder.add_pi_electrons as add_pi_electrons
import src.sheet_builder.misc_functions as misc_functions
import src.sheet_builder.command_line as command_line
import src.sheet_builder.build_sheets as build_sheets
import src.sheet_builder.build_tubes as build_tubes
import src.sheet_builder.atoms2lmp as atoms2lmp
import src.io_functions as io_functions
import src.write_lmp as write_lmp
import time
import sys
import os



#####################################
### Main function to build sheets ###
#####################################
def main(sheet_basename, symmetric_tube_basename, chiral_tube_basename, run_mode, parent_directory, length_in_perpendicular, length_in_edgetype, sheet_edgetype, types,
         bond_length, sheet_layer_spacing, sheet_nlayers, stacking, plane, tube_edgetype, tube_layer_spacing, symmetric_ntubes, symmetric_length, diameter, n, m,
         chiral_length, symmetric_tube_axis, chiral_tube_axis, find_bonds, periodic_bonds, charges, masses, commandline_inputs=[], log=io_functions.LUNAR_logger()):
    start_time = time.time()
    
    # Configure log (default is level='production', switch to 'debug' if debuging)
    log.configure(level='production')
    #log.configure(level='debug')
    
    #########################
    # Command Line Override #
    #########################
    # if -opt or -man option is in commandline_inputs print options and stop code execution
    if '-opt' in commandline_inputs or  '-man' in commandline_inputs:
        # call man page and exit if '-opt' or '-man' is provided at the command line
        command_line.print_man_page(sheet_basename, symmetric_tube_basename, chiral_tube_basename, run_mode, parent_directory, length_in_perpendicular, length_in_edgetype, sheet_edgetype, types,
                                    bond_length, sheet_layer_spacing, sheet_nlayers, stacking, plane, tube_edgetype, tube_layer_spacing, symmetric_ntubes, symmetric_length, diameter, n, m,
                                    chiral_length, symmetric_tube_axis, chiral_tube_axis, find_bonds, periodic_bonds)
        sys.exit()
        
    ###################################################################################
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters #
    ###################################################################################
    if '-opt' not in commandline_inputs and commandline_inputs:    
        # call inputs for commandline over rides
        over_rides = command_line.inputs(commandline_inputs, sheet_basename, symmetric_tube_basename, chiral_tube_basename, run_mode, parent_directory, length_in_perpendicular, length_in_edgetype,
                                         sheet_edgetype, types, bond_length, sheet_layer_spacing, sheet_nlayers, stacking, plane, tube_edgetype, tube_layer_spacing, symmetric_ntubes, symmetric_length,
                                         diameter, n, m, chiral_length, symmetric_tube_axis, chiral_tube_axis, find_bonds, periodic_bonds)
        
        # Set new inputs from over_rides class
        sheet_basename = over_rides.sheet_basename
        symmetric_tube_basename = over_rides.symmetric_tube_basename
        chiral_tube_basename = over_rides.chiral_tube_basename
        run_mode = over_rides.run_mode
        parent_directory = over_rides.parent_directory
        length_in_perpendicular = over_rides.length_in_perpendicular
        length_in_edgetype = over_rides.length_in_edgetype
        sheet_edgetype = over_rides.sheet_edgetype
        types = over_rides.types
        bond_length = over_rides.bond_length
        sheet_layer_spacing = over_rides.sheet_layer_spacing
        sheet_nlayers = over_rides.sheet_nlayers
        stacking = over_rides.stacking
        plane = over_rides.plane
        tube_edgetype = over_rides.tube_edgetype
        tube_layer_spacing = over_rides.tube_layer_spacing
        symmetric_ntubes = over_rides.symmetric_ntubes
        symmetric_length = over_rides.symmetric_length
        diameter = over_rides.diameter
        n = over_rides.n
        m = over_rides.m
        chiral_length = over_rides.chiral_length
        symmetric_tube_axis = over_rides.symmetric_tube_axis
        chiral_tube_axis = over_rides.chiral_tube_axis
        find_bonds = over_rides.find_bonds
        periodic_bonds = over_rides.periodic_bonds
    
    
    ###########################################
    # Initialize some preliminary information #
    ###########################################
    # set version and print starting information to screen
    version = 'v1.0 / 4 June 2024'
    log.out(f'\n\nRunning sheet_builder {version} in mode {run_mode}')
    log.out(f'Using Python version {sys.version}')
    
    
    #######################################################################################
    # Determine atom types and/or pi-electrons naming = 'atomtype|pi-electron', where the #
    # '|' character seperates the two. If no '|', then no pi-electrons will be added.     #
    #######################################################################################
    atom_types = {}; pi_electrons = {};
    for i in types:
        string = types[i]
        if '|' in string:
            split_string = string.split('|')
            atom_type = split_string[0]
            pi_electron = split_string[1]
        else:
            atom_type = string
            pi_electron = ''
        atom_types[i] = atom_type
        pi_electrons[i] = pi_electron
    
    
    ################################
    # Generate sheets or nanotubes #
    ################################
    if run_mode == 'sheet':
        if plane == 'xy': boundary = 'p p f'
        elif plane == 'xz': boundary = 'p f p'
        elif plane == 'yz': boundary = 'f p p'
        else: log.error(f"ERROR plane {plane} is not supported. Supported planes are 'xy' or 'xz' or 'yz'")
        basename = sheet_basename
        pflag = True
        atoms, box = build_sheets.generate(length_in_perpendicular, length_in_edgetype, bond_length, sheet_edgetype, atom_types, sheet_layer_spacing, sheet_nlayers, stacking, plane, periodic_bonds, pflag, log)
    elif run_mode in ['chiral-tube', 'symmetric-tube']:
        if run_mode == 'symmetric-tube':
            if symmetric_tube_axis == 'x': boundary = 'p f f'
            elif symmetric_tube_axis == 'y': boundary = 'f p f'
            elif symmetric_tube_axis == 'z': boundary = 'f f p'
            basename = symmetric_tube_basename
            atoms, box = build_tubes.generate_MWCNT(symmetric_length, diameter, bond_length, atom_types, tube_edgetype, tube_layer_spacing, symmetric_ntubes, symmetric_tube_axis, periodic_bonds, log)
        if run_mode == 'chiral-tube':
            if chiral_tube_axis == 'x': boundary = 'p f f'
            elif chiral_tube_axis == 'y': boundary = 'f p f'
            elif chiral_tube_axis == 'z': boundary = 'f f p'
            basename = chiral_tube_basename
            atoms, box = build_tubes.generate_chiral(n, m, chiral_length, bond_length, atom_types, chiral_tube_axis, periodic_bonds, log)
    else: log.error(f"ERROR run_mode {run_mode} is not supported. Supported modes are 'sheet' or 'symmetric-tube' or 'chiral-tube'")    
    log.out(f'  Created: {len(atoms)} atoms')
    
    
    ##############################
    # find bonds if user desires #
    ##############################
    if find_bonds:
        if not periodic_bonds: boundary = 'f f f'
        pflag = True # print bonding results
        max_bonds_per_atom = 3 # only use the first 3 closest atoms to generate bonds
        tolerance = bond_length/2 # max_distance = r0*tolerance, which tolerance can be large since code will use 3 closets atoms to create bonds
        domain_size = 5*bond_length # generate a moderate domain size for domain decomposition (best of all worlds for performance)
        log.out('\n\nFinding bonds ...')
        bonds = misc_functions.find_bonds(atoms, box, boundary, bond_length, tolerance, max_bonds_per_atom, domain_size, pflag, log)
        log.out(f'  Created: {len(bonds)} bonds')
    else: bonds = {}
    
    
    #####################################################################
    # Add pi-electrons if any of the four types is not any empty string #
    #####################################################################
    electron_types = list(pi_electrons.values())
    if electron_types.count('') != 4:
        if find_bonds:
            log.out(f'\n\nAdding pi-electrons of types = {str(pi_electrons)} ...')
            atoms, bonds = add_pi_electrons.add(atoms, bonds, box, pi_electrons, log)
        else:
            log.warn(f'WARNING pi-electrons of types = {str(pi_electrons)}, but find_bonds is False. Set find_bonds as True to add pi-electrons.')
          
            
    ###############################################################
    # Attempting adding charges to the system, if they dont exist #
    # in the charges dictionary charges will remain as zeros      #
    ###############################################################
    atoms = misc_functions.charge_system(atoms, charges, log)
    
    
    #######################################
    # Setting up directories and where to #
    # write final files and results to    #
    #######################################
    # Find present working directory
    pwd = os.getcwd()
    
    # Find/create paths to store code results
    path = os.path.join(pwd, parent_directory)
        
    # Check if path exits. IF not create
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
            
    # Change the current working directory to path2
    # so all files get written to that directory
    os.chdir(path)
    
    
    ##############################
    # Write new lammps data file #
    ##############################
    header = 'HEADER, created with sheet_builder {}'.format(version)
    if atoms and box:
        # Create m-object like read_lmp would generate
        m = atoms2lmp.Molecule_File(basename, atoms, bonds, box, masses, header)
        
        # Write LAMMPS datafile
        atom_style = 'full'
        include_type_labels = False
        write_lmp.file(m, basename+'.data', header, atom_style, include_type_labels, log)
        
        # Write .nta file only if find_bonds is True
        if find_bonds:
            style = 'type' # 'type' or 'id'
            misc_functions.write_nta(m, style, version, charges, basename+'.nta')
    
    ###########
    # Wrap up #
    ###########
    # Print file locations
    log.out(f'\n\nAll outputs can be found in {path} directory')
    
    # Print completion of code
    log.out('\n\nNormal program termination\n\n')
    
    # Script run time
    execution_time = (time.time() - start_time)
    log.out('Execution time in seconds: ' + str(execution_time))  
    
    # Show number of warnings and errors
    log.out_warnings_and_errors()
    
    # write log
    log.write_logged(basename+'.log.lunar')
    
    # Change back to the intial directory to keep directory free for deletion
    os.chdir(pwd)
    return