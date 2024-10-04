# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.10
October 4th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.bond_react_merge.merge_coeffs as merge_coeffs
import src.cell_builder.write_lmp_dump as write_lmp_dump
import src.cell_builder.misc_functions as misc_functions
import src.cell_builder.lmp_inscript as lmp_inscript
import src.cell_builder.command_line as command_line
import src.cell_builder.group_atoms as group_atoms
import src.cell_builder.subcell as subcell
import src.io_functions as io_functions
import src.write_lmp as write_lmp
import src.read_lmp as read_lmp
import time
import sys
import os


                
    
##############################################
### Main function to perform cell building ###
##############################################
def main(topofiles, force_field_joining, duplicate, distance_scale, newfile, atom_style, parent_directory, max_rotations,
         reset_molids, unwrap_atoms_via_image_flags, include_type_labels, group_monomers_locally,
         seed, domain, maxtry, tolerance, mixing_rule, boundary, commandline_inputs=[], log=io_functions.LUNAR_logger()):
    start_time = time.time()
    
    # Configure log (default is level='production', switch to 'debug' if debuging)
    log.configure(level='production')
    #log.configure(level='debug', print2console=False)
    
    
    ########################################################
    # set version and print starting information to screen #
    ########################################################
    version = 'v1.10 / 4 October 2024'
    log.out(f'\n\nRunning cell_builder {version}')
    log.out(f'Using Python version {sys.version}')
    
    
    #########################
    # Command Line Override #
    #########################
    # if -opt or -man option is in commandline_inputs print options and stop code execution
    if '-opt' in commandline_inputs or  '-man' in commandline_inputs:
        # call man page and exit if '-opt' or '-man' is provided at the command line
        command_line.print_man_page(topofiles, force_field_joining, duplicate, distance_scale, newfile, atom_style, parent_directory, max_rotations,
                                    reset_molids, unwrap_atoms_via_image_flags, include_type_labels, group_monomers_locally, seed, domain, maxtry,
                                    tolerance, mixing_rule, boundary)
        sys.exit()
        
        
    ###################################################################################
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters #
    ###################################################################################
    if '-opt' not in commandline_inputs and commandline_inputs:    
        # call inputs for commandline over rides
        over_rides = command_line.inputs(topofiles, force_field_joining, duplicate, distance_scale, newfile, atom_style, parent_directory, max_rotations,
                                          reset_molids, unwrap_atoms_via_image_flags, include_type_labels, group_monomers_locally,
                                          seed, domain, maxtry, tolerance, mixing_rule, boundary, commandline_inputs)
        
        # Set new inputs from over_rides class
        topofiles = over_rides.topofiles
        force_field_joining = over_rides.force_field_joining
        parent_directory = over_rides.parent_directory
        atom_style = over_rides.atom_style
        include_type_labels = over_rides.include_type_labels
        newfile = over_rides.newfile
        duplicate = over_rides.duplicate
        distance_scale = over_rides.distance_scale
        reset_molids = over_rides.reset_molids
        max_rotations = over_rides.max_rotations
        unwrap_atoms_via_image_flags = over_rides.unwrap_atoms_via_image_flags
        group_monomers_locally = over_rides.group_monomers_locally
        seed = over_rides.seed
        domain = over_rides.domain
        maxtry = over_rides.maxtry
        tolerance = over_rides.tolerance 
        mixing_rule = over_rides.mixing_rule
        boundary = over_rides.boundary
                
    
    
    ###############################################################
    # Read in files and set fileID to qty map. FileIDs will be    #
    # set when reading files and will start at 1 and go to nfiles #
    ###############################################################
    if len(topofiles) == 0: log.error('ERROR files dictionary is empty')
    # Unique style hints logger
    unique_style_hints = {'Masses': [], 'Pair_Coeffs': [], 'Bond_Coeffs': [],
                          'Angle_Coeffs': [], 'Dihedral_Coeffs': [],
                          'Improper_Coeffs': [], 'BondBond_Coeffs': [],
                          'BondAngle_Coeffs': [], 'AngleAngleTorsion_Coeffs': [],
                          'EndBondTorsion_Coeffs': [], 'MiddleBondTorsion_Coeffs': [],
                          'BondBond13_Coeffs': [], 'AngleTorsion_Coeffs': [], 'AngleAngle_Coeffs': []}
    files = {} # { fileID : file-object }
    qtys = {} # { fileID : n-quantity } 
    molecule_insertion = {} # { fileID : [inserted, desired] }
    moleculespans = [] # [lst to record max molecule spans of]
    ntypes = {'atom':set(),'bond':set(),'angle':set(),'dihedral':set(),'improper':set()}
    n = 0; grouping_file = ''; defined_groups = set(); # defaults for grouping file
    for topofile in topofiles:
        qty = topofiles[topofile]
        if topofile.endswith('data'):
            if os.path.isfile(topofile):
                # read file and recenter about (0, 0, 0 and find maxspan)
                m = read_lmp.Molecule_File(topofile, method='forward', sections=['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Velocities'])
                
                # Find unique style hints 
                unique_style_hints['Masses'].append(m.mass_coeffs_style_hint)
                unique_style_hints['Pair_Coeffs'].append(m.pair_coeffs_style_hint)
                unique_style_hints['Bond_Coeffs'].append(m.bond_coeffs_style_hint)
                unique_style_hints['Angle_Coeffs'].append(m.angle_coeffs_style_hint)
                unique_style_hints['Dihedral_Coeffs'].append(m.dihedral_coeffs_style_hint)
                unique_style_hints['Improper_Coeffs'].append(m.improper_coeffs_style_hint)
                unique_style_hints['BondBond_Coeffs'].append(m.bondbond_coeffs_style_hint)
                unique_style_hints['BondAngle_Coeffs'].append(m.bondangle_coeffs_style_hint)
                unique_style_hints['AngleAngleTorsion_Coeffs'].append(m.angleangletorsion_coeffs_style_hint)
                unique_style_hints['EndBondTorsion_Coeffs'].append(m.endbondtorsion_coeffs_style_hint)
                unique_style_hints['MiddleBondTorsion_Coeffs'].append(m.middlebondtorsion_coeffs_style_hint)
                unique_style_hints['BondBond13_Coeffs'].append(m.bondbond13_coeffs_style_hint)
                unique_style_hints['AngleTorsion_Coeffs'].append(m.angletorsion_coeffs_style_hint)
                unique_style_hints['AngleAngle_Coeffs'].append(m.angleangle_coeffs_style_hint)
                

                # Unwrap atoms and geometrically center molecules that will be added to the system
                if qty > 0:
                    if unwrap_atoms_via_image_flags:
                        m = misc_functions.unwrap_atoms_via_iflags(m)
                    m = misc_functions.center_about_000(m)
                    moleculespans.append(m.maxspan)
                
                    # Get header info such as Regions or Rotations
                    m.Regions, m.Rotations, m.Shfits = misc_functions.get_dict(m.header, log)
                    
                # Check for Shifts from files with qty = 0 (ONLY USE THE Shifts dict as others dont
                # make sense to use for files that are meant to define a system being built.)
                if qty == 0:
                    Regions, Rotations, Shift = misc_functions.get_dict(m.header, log)
                    if Shift: m = misc_functions.shift_system(m, Shift, log)
                
                # Save and and set fileID to qty map
                n += 1
                files[n] = m; qtys[n] = topofiles[topofile]; molecule_insertion[n] = [0, qty*duplicate]
                
                # Check ntypes
                ntypes['atom'].add(len(m.masses))
                ntypes['bond'].add(len(m.bond_coeffs))
                ntypes['angle'].add(len(m.angle_coeffs))
                ntypes['dihedral'].add(len(m.dihedral_coeffs))
                ntypes['improper'].add(len(m.improper_coeffs))
            else: log.error(f'ERROR lammps datafile: {topofile} does not exist')
        else:
            if os.path.isfile(topofile):
                log.out(f'The qty associated with filename {topofile} is ZERO. Using file to group atoms.')
                grouping_file = topofile
            else: log.error(f'ERROR file: {topofile} does not exist')
    
    
    #########################################
    # Set up groups if grouping_file exists #
    #########################################
    if grouping_file:
        files, defined_groups = group_atoms.add_group_to_atoms(files, grouping_file, log)
        log.debug(defined_groups)
    
    
    ###############################################
    # Set up how to deal with force field joining #
    ###############################################
    if force_field_joining == 'none':
        offset_coeff_types = False
    elif force_field_joining == 'merge':
        # warn if system already seems merged
        for i in ntypes:
            if ntypes[i] and len(ntypes[i]) == 1:
                log.warn(f'WARNING consistent number of {i} types between files. Perhaps  force_field_joining method "merge" was not meant to be used.')
        
        # merge coeffs and update ntypes
        offset_coeff_types = False
        new = merge_coeffs.merged(files, log)
        merge_coeffs.print_merged(new, log, skip_printing_cross_terms=True)
        ntypes = {'atom':set(),'bond':set(),'angle':set(),'dihedral':set(),'improper':set()}
        for i in files:
            m = files[i]
            m = merge_coeffs.update_TypeIDs(m, new, log) # Update coeffs from merge
            files[i] = m
            ntypes['atom'].add(len(m.masses))
            ntypes['bond'].add(len(m.bond_coeffs))
            ntypes['angle'].add(len(m.angle_coeffs))
            ntypes['dihedral'].add(len(m.dihedral_coeffs))
            ntypes['improper'].add(len(m.improper_coeffs))
    elif force_field_joining == 'offset':
        # warn if system already seems merged
        for i in ntypes:
            if ntypes[i] and len(ntypes[i]) == 1:
                log.warn(f'WARNING consistent number of {i} types between files. Perhaps  force_field_joining method "offset" was not meant to be used.')
        offset_coeff_types = True
    else: log.error(f'ERROR force_field_joining {force_field_joining} option not supported')
    
        
    #######################################################################################
    # Loop through style hints to make sure they are all consistant across all files.     #
    # Additionally if not using offset_coeff_types, check for consistent number of coeffs #
    #######################################################################################
    for coeff in unique_style_hints:
        hints = unique_style_hints[coeff]
        no_dups = set(hints)
        if 'N/A' not in hints and len(no_dups) > 1:
            log.warn(f'WARNING style hints between some files are different for {coeff} Coeff. Style Hints {str(no_dups)}')
    #if not offset_coeff_types:
    if not offset_coeff_types:
        for i in ntypes:
            if ntypes[i] and len(ntypes[i]) != 1:
                log.error(f'ERROR inconsistent number of {i} types between files. Either use force_field_joining method "merge" or "offset".')
                
    
    ##################################################################################################
    # Start building new simulation cell out of subcells w/ randomly orientated and placed molecules #
    ##################################################################################################
    # Execute local grouping first if user desires
    occurrences = 0; molids_based_on_files = {}; molids_based_on_offsets = {}; coeff_offsets = {};  subcells = []; # to keep track of if locally grouping
    if group_monomers_locally:
        log.out('Grouping monomers locally ...')
        if hasattr(m,'Regions'): m.Regions = {}; log.warn('WARNING Regions header info is not compatable with group_monomers_locally. Ignoring Regions.')
        if hasattr(m,'Rotations'): m.Rotations = {}; log.warn('WARNING Rotations header info is not compatable with group_monomers_locally. Ignoring Rotations.')
        if domain.count('A') == 3:
            log.error(f'ERROR can not used group_monomers_locally and random insertion setup with domain {domain}.')
        tmp_files = {}; tmp_qtys = {}; tmp_moleculespans = []; tmp_duplicate = duplicate; tmp_molecule_insertion = {} # { fileID : [inserted, desired] }
        for n in range(duplicate):
            locally_grouped = subcell.constructor(files, qtys, offset_coeff_types, 1, max_rotations, distance_scale, moleculespans, seed, reset_molids, domain, molecule_insertion,
                                                  occurrences, subcells, maxtry, tolerance, mixing_rule, boundary, log, pflag=False, grouping=True)
            locally_grouped = misc_functions.center_about_000(locally_grouped)
            tmp_moleculespans.append(locally_grouped.maxspan)
            locally_grouped.Regions = {}; locally_grouped.Rotations = {};
            tmp_files[n+1] = locally_grouped; tmp_qtys[n+1] = 1; subcells.append(locally_grouped)
            tmp_molecule_insertion[n+1] = [1, 1]
            if n == 0:
                coeff_offsets = locally_grouped.coeff_offsets
                molids_based_on_files = locally_grouped.molids_based_on_files
                molids_based_on_offsets = locally_grouped.molids_based_on_offsets
        duplicate = 1; files = tmp_files; qtys = tmp_qtys; molecule_insertion = tmp_molecule_insertion
        moleculespans.extend(tmp_moleculespans); occurrences += 1

    # Build larger simulation cell
    m = subcell.constructor(files, qtys, offset_coeff_types, duplicate, max_rotations, distance_scale, moleculespans, seed, reset_molids, domain, molecule_insertion,
                            occurrences, subcells, maxtry, tolerance, mixing_rule, boundary, log, pflag=True, grouping=False)
    occurrences += 1
    
    # if occurance equals 2 over ride current dicts based on local grouping dicts
    if occurrences == 2:
        duplicate = tmp_duplicate
        m.coeff_offsets = coeff_offsets
        m.molids_based_on_files = molids_based_on_files
        m.molids_based_on_offsets = molids_based_on_offsets


    #########################
    # Write parameters used #
    #########################
    log.out('\n\n\n-----------------------------------------------------')
    log.out('   Cell building and randomization parameters used')
    log.out('-----------------------------------------------------')
    log.out('{:<26} {:>18}'.format('unwrap atoms by iflags: ', str(unwrap_atoms_via_image_flags)))
    log.out('{:<26} {:>18}'.format('locally group atoms: ', str(group_monomers_locally)))
    log.out('{:<26} {:>18}'.format('force field joining: ', str(force_field_joining)))
    log.out('{:<26} {:>18}'.format('distance scale: ', str(distance_scale)))
    log.out('{:<26} {:>18}'.format('duplicate: ', str(duplicate)))
    log.out('{:<26} {:>18}'.format('max rotations: ', str(max_rotations)))
    log.out('{:<26} {:>18}'.format('seed: ', str(seed)))
    log.out('{:<26} {:>18}'.format('domain: ', str(domain)))
    log.out('{:<26}'.format('System files used: '))
    for file in topofiles:
        log.out('   {} qty: {}'.format(file, topofiles[file]))
        try:
            for coeff in m.coeff_offsets[file]:
                offset = m.coeff_offsets[file][coeff]
                coeff_type = '{} TypeIDs offset'.format(coeff)
                log.out('       {:<24} {}'.format(coeff_type, offset))
        except: pass
        
    # write molIDs based on files if reset_molids was files
    if reset_molids == 'files':
        log.out('{:<18}'.format("reset_molids = 'files': "))
        for filename in m.molids_based_on_files:
            log.out('    all atoms in {} have their molIDs set as {}'.format(filename, m.molids_based_on_files[filename]))
            
    # write molIDs based on offset if reset_molids was files
    if reset_molids == 'offset':
        log.out('{:<18}\n'.format("reset_molids = 'offset': "))
        for filename in m.molids_based_on_offsets:
            log.out('    all atoms in {} have their molIDs offset by {}'.format(filename, m.molids_based_on_offsets[filename]))    

    ################################
    # Reset molids is user desires #
    ################################
    if reset_molids == 'clusters' and len(m.bonds) != 0:
        m = misc_functions.update_molids(m, log)
        
        
    ################################
    # Print new system information #
    ################################
    log.out('\n\n\n--------------------------------------')
    log.out('Atoms/Bonds/Angles/Dihedrals/Impropers')
    log.out('--------------------------------------')
    log.out('{:<18} {:^5}'.format('natoms: ', m.natoms))
    log.out('{:<18} {:^5}'.format('nbonds: ', m.nbonds))
    log.out('{:<18} {:^5}'.format('nangles: ', m.nangles))
    log.out('{:<18} {:^5}'.format('ndihedrals: ', m.ndihedrals))
    log.out('{:<18} {:^5}'.format('nimpropers: ', m.nimpropers))
    lx, ly, lz, volume, mass, density = misc_functions.compute_mass_volume_density(m)
    log.out('\n')
    log.out('---------------------------------------')
    log.out('System Unit cell, volume, mass, density')
    log.out('---------------------------------------')
    log.out('{:<10} {:^10.5f} {:<10}'.format('Lx:', lx, 'angstrom'))
    log.out('{:<10} {:^10.5f} {:<10}'.format('Ly:', ly, 'angstrom'))
    log.out('{:<10} {:^10.5f} {:<10}'.format('Lz:', lz, 'angstrom'))
    log.out('{:<10} {:10.4E} {:<10}'.format('volume:', volume, 'cm^3'))
    log.out('{:<10} {:10.4E} {:<10}'.format('mass:', mass, 'grams'))
    log.out('{:<10} {:^10.5f} {:<10}' .format('density:', density, 'g/cc'))
    
    
    ########################################################################
    # Setting up directories and where to write final files and results to #
    ########################################################################
    # Find present working directory
    pwd = os.getcwd()
    
    # Find/create paths to store code results
    path = os.path.join(pwd, parent_directory)
    
    # If parent_directory == 'topofile' use topofile path as path
    if 'topofile' in parent_directory:
        log.out('Using path from 1st file in files dictionary to set parent_directory ...')
        topofile = list(topofiles.keys())[0]
        path = io_functions.get_dir_from_topofile(topofile, parent_directory)
    
    # Check if path exits. IF not create
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
        
    # Change the current working directory to path1
    # so all files get written to that directory
    os.chdir(path)
        
    
    ################################################################################
    # Write new LAMMPS datafile, new mol2 file, and a generalized lmp input script #
    ################################################################################
    #newfile = '{}_seed_{}'.format(newfile, seed)
    header = '{} > cell_builder version: {}'.format(m.header, version)
    write_lmp.file(m, newfile+'.data', header, atom_style, include_type_labels, log)
    lmp_inscript.write('in.'+newfile, version, atom_style, unique_style_hints, log)
    
    # If grouping file, write out lmp groups file
    if grouping_file: group_atoms.write(newfile, defined_groups, m, log, write_group_lmp_file=True)
    
    # Internal flag to write a dump file of all subcells
    write_grp_local_subcells = True
    if write_grp_local_subcells and subcells:
        write_lmp_dump.subcells2dump(newfile, subcells)
    
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
    log.write_logged(newfile+'.log.lunar')
    
    # Change back to the intial directory to keep directory free for deletion
    os.chdir(pwd)
    return