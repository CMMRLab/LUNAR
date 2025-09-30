# -*- coding: utf-8 -*-
"""
Revision 1.12
April 14, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

 
    *********************************************************
    * Requirements:                                         *
    *   python 3.7+                                         *
    *                                                       *   
    * Dependencies:                                         *
    *   python tqdm module:                                 *
    *    - pip3 install tqdm (if pip manager is installed)  *
    *                                                       *
    *   python numpy module:                                *
    *    - pip3 install numpy (if pip manager is installed) *
    *                                                       *
    *********************************************************
"""


##############################
# Import Necessary Libraries #
##############################
import src.auto_morse_bond.compute_bonds_distances as compute_bonds_distances
import src.auto_morse_bond.class2xe_conversion as class2xe_conversion
import src.auto_morse_bond.potential_plotting as potential_plotting
import src.auto_morse_bond.alpha_parameter_min as alpha_parameter
import src.auto_morse_bond.merge_input_files as merge_input_files
import src.auto_morse_bond.cluster_analysis as cluster_analysis
import src.auto_morse_bond.zero_xterm_r0s as zero_xterm_r0s
import src.auto_morse_bond.ring_analysis as ring_analysis
import src.auto_morse_bond.coordination as coordination
import src.auto_morse_bond.command_line as command_line
import src.auto_morse_bond.out2console as out2console
import src.auto_morse_bond.atom_info as atom_info
import src.auto_morse_bond.bond_info as bond_info
import src.io_functions as io_functions
import src.write_lmp as write_lmp
import glob
import time
import sys
import os


######################################################
### Main function to perform all atom-typing tasks ###
######################################################
def main(topofile, morsefile, parent_directory, newfile, mass_map, min_bond_length, coeffs2skip, radius_specs, alpha_specs, alpha_scale, files2write, atom_style,
         zero_effected_xterms, bondbreak_scale, ff_class, include_type_labels, class2xe_update, include_rcut, commandline_inputs=[], log=None):
    
    #########################
    # Command Line Override #
    #########################
    # if -opt or -man option is in commandline_inputs print options and stop code execution
    if '-opt' in commandline_inputs or  '-man' in commandline_inputs:
        # call man page and exit if '-opt' or '-man' is provided at the command line
        command_line.print_man_page(topofile, morsefile, parent_directory, newfile, mass_map, min_bond_length, coeffs2skip,
                                    radius_specs, alpha_specs, alpha_scale, files2write, atom_style,
                                    zero_effected_xterms, bondbreak_scale, ff_class, include_type_labels,
                                    class2xe_update, include_rcut)
        sys.exit()
        
    ###################################################################################
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters #
    ###################################################################################
    if '-opt' not in commandline_inputs and commandline_inputs:    
        # call inputs for commandline over rides
        over_rides = command_line.inputs(topofile, morsefile, parent_directory, newfile, min_bond_length, coeffs2skip,
                                         alpha_scale, atom_style, zero_effected_xterms, bondbreak_scale,
                                         ff_class, include_type_labels, class2xe_update, include_rcut, commandline_inputs)
        
        # Set new inputs from over_rides class
        topofile = over_rides.topofile
        morsefile = over_rides.morsefile
        parent_directory = over_rides.parent_directory
        newfile = over_rides.newfile
        atom_style = over_rides.atom_style
        ff_class = over_rides.ff_class
        include_type_labels = over_rides.include_type_labels
        min_bond_length = over_rides.min_bond_length
        coeffs2skip = over_rides.coeffs2skip
        alpha_scale = over_rides.alpha_scale
        zero_effected_xterms = over_rides.zero_effected_xterms
        bondbreak_scale = over_rides.bondbreak_scale
        class2xe_update = over_rides.class2xe_update
        include_rcut = over_rides.include_rcut
    
    #####################################################
    # Set up Tristan's "array" analysis using recursion #
    #####################################################
    if not os.path.isfile(str(topofile)):
        if log is None: log = io_functions.LUNAR_logger()
        log.configure(level='production', print2console=True, write2log=True)
        files = glob.glob(topofile); array_time = time.time()
        if files:
            for n, file in enumerate(files, 1):
                log.clear_all()
                log.out('\n\nUsing array input option:')
                log.out(' - topofile     : {}'.format(topofile))
                log.out(' - matched file : {}'.format(file))
                log.out(' - elapsed time : {:.2f} (seconds)'.format(time.time() - array_time))
                log.out(' - progress     : {} of {} ({:.2f}%)'.format(n, len(files), 100*(n/len(files))))
                root = file[:file.rfind('.')]
                if io_functions.check_outfile_existance(file, newfile, parent_directory, filetype='topofile'):
                    log.warn(f' - WARNING matched file {file} already has been processed and was skipped')
                    continue
                if newfile.endswith(':') and root.startswith(newfile[:-1]): # prefix
                    log.warn(f' - WARNING matched file {file} already has newfile {newfile} extension and was skipped')
                    continue
                elif newfile.startswith(':') and root.endswith(newfile[1:]): # suffix
                    log.warn(f' - WARNING matched file {file} already has newfile {newfile} extension and was skipped')
                    continue
                try: # we dont want crashes to exit this loop
                    main(file, morsefile, parent_directory, newfile, mass_map, min_bond_length, coeffs2skip, radius_specs, alpha_specs, alpha_scale, files2write, atom_style,
                         zero_effected_xterms, bondbreak_scale, ff_class, include_type_labels, class2xe_update, include_rcut, commandline_inputs, log)
                except: pass
            print('\a') # Alert
        else: log.error(f'ERROR topofile: {topofile} unwrapped zero files or does not exist')
        return
    else:
        # To remove recursion, un-indedent this all the way to the return and delete the lines above
        start_time = time.time()
        
        # Configure log (default is level='production', switch to 'debug' if debuging)
        if log is None:
            log = io_functions.LUNAR_logger()
        log.configure(level='production')
        #log.configure(level='debug', print2console=False)
    
    
        ###########################################
        # Initialize some preliminary information #
        ###########################################
        # set version and print starting information to screen
        version = 'v1.12 / 14 April 2025'
        log.out(f'\n\nRunning auto_morse_bond: {version}')
        log.out(f'Using Python version: {sys.version}')
        log.out(f'Using Python executable: {sys.executable}')
        
        
        ######################################################################################################
        # Read .data file w/bonds:                                                                           #
        #      - m.atoms[ID] instance .element = element; .comment = element                                 #
        ######################################################################################################
        m = merge_input_files.merge(topofile, mass_map, ff_class, log)   
        
        
        ###################################################################################
        # Find clusters and add cluster info to m. Will add the following instances to m: #
        #   m.molecules -> .data .atoms .clusters .formula .molids                        # 
        #   m.atoms[ID] -> .molecule -> .mass .size .formula                              #
        ###################################################################################
        m = cluster_analysis.add_molecule_data2m(m, log)
    
            
        #####################################################################################################
        # Find and add neighbor_ids information to mm.atoms as mm.atoms[ID].neighbor_ids; were neighbor_ids #
        # is a dictionary as:                                                                               #
        #    {1st-neighs,   2nd-neighs,   Nth-neighs,}                                                      #
        #    {1: [1, 2, 3], 2: [4, 5, 6], ....}                                                             #
        #    Use a modified BFS algorithm to find Nth neighbors away from each atom to find neighboring IDs #
        #    set the Nth_neigh_depth as 4 and find up to 4 neighs away from atomid (adjust as needed)       #
        #####################################################################################################
        m = coordination.neighbors(m, Nth_neigh_depth=4)
        
        
        ########################################################################################################
        # Find and add ring data to m class. ring data will be found via information in find_rings dictionary  #
        # below; where 'elements2walk' sets the element types that can belong in the ring and 'rings2check are #
        # sizes of rings to look for. Code provided by Jake Gissinger with modifications and additions by Josh #
        ########################################################################################################
        find_rings = {'elements2walk': list(mass_map.keys()), # Elements to walk (set by elements found in mass_map)
                      'rings2check': [3, 4, 5, 6],            # Ring sizes to check
                      'fused-rings': False,                   # True or False to run ring connectivty analysis ('perform' key must be True).
                      'fused2check': [5, 6]}
        m = ring_analysis.add_ring_data2m(m, find_rings, log)
        
        
        ##########################################################################################################
        # Find final atom info needed to perform atom-typing and add tom m.atoms[ID] instance for each atom.     #
        # Full verbose set of instances that mm.atoms[ID] will contain for determing atom type:                  #
        #     .molecule.formula = formula atomID belongs to (sorted by element with '-' delimiter. I.E. C1-H4)   #
        #     .element = element type of atomID                                                                  #
        #     .neighbor_ids = {1: [1, 2, 3], 2: [4, 5, 6], ....} dict or neighbors at certain depths             #  
        #     .nb = number of bonded 1st neighbors                                                               #
        #     .info = [element, ringsize, nb]                                                                    #
        #     .neighbor_info = {1: [], 2: []}; lst = [['C', 6, 3], ['H', 0, 1]]; sorted by nb -> ring -> element #
        ##########################################################################################################
        m = atom_info.add(m)
        
        
        #################################################
        # Find bond information and fit alpha parameter #
        #################################################
        radius = [round(n*radius_specs['increment']+radius_specs['start'], str(radius_specs['increment'])[::-1].find('.')) for n in range( int((radius_specs['end']-radius_specs['start'])/radius_specs['increment'])+1) ]
        m.bonds_lengths = compute_bonds_distances.compute(m)
        m.bond_info = bond_info.topology(m, min_bond_length, coeffs2skip, ff_class, morsefile, class2xe_update, log)
        m.alpha_parameter = alpha_parameter.fitting(m, m.bond_info, radius, alpha_specs, alpha_scale, ff_class, bondbreak_scale, include_rcut, log)
        
        
        ############################################################
        # Create .harmonic attribute and update .coeffs with morse #
        ############################################################
        bond2style = {} # {bondTypeID:'style'}
        for i in m.bond_coeffs:
            m.bond_coeffs[i].harmonic = m.bond_coeffs[i].coeffs
            m.bond_coeffs[i].coeffs = m.alpha_parameter.morse_harmonic[i]
            bond2style[i] = m.bond_coeffs[i].coeffs[0]
        
        
        ############################
        # log file current results #
        ############################
        out2console.out(m, log, version, min_bond_length, coeffs2skip, zero_effected_xterms, alpha_specs, alpha_scale, include_type_labels)
        
        
        ###########################
        # Class2 crossterm issues #
        ###########################
        # Log any changes in potential styles to log to user
        potential_styles = [] # [(potential, style), ...]
        
        # Check that zero_effected_xterms and class2_morse are not both True
        if zero_effected_xterms and class2xe_update:
            log.error('ERROR zero_effected_xterms and class2xe_update are both True. Only One can be used.')
        
        # Zero any effect crossterms if users desire
        if zero_effected_xterms and ff_class in [2, '2']: 
            m = zero_xterm_r0s.zero(m, log)
            
        # Convert class2 crossterms to class2xe variant
        if class2xe_update and ff_class in [2, '2']:
            m = class2xe_conversion.update(m, morsefile, potential_styles, include_rcut, log)
            
        # Update bond_coeffs.coeffs if only a single bond style is used
        bond_styles = {m.bond_coeffs[i].coeffs[0] for i in m.bond_coeffs} # will be 'harmonic' or 'class2' or 'morse'
        if len(bond_styles) == 1:
            for i in m.bond_coeffs:
                coeffs = m.bond_coeffs[i].coeffs
                coeffs.pop(0) # remove style
                m.bond_coeffs[i].coeffs = coeffs
            bond_style = list(bond_styles)[0]
            m.bond_coeffs_style_hint = bond_style
            potential_styles.insert(0, ('bond', bond_style))
        else: 
            bond_style = ' '.join(bond_styles)
            m.bond_coeffs_style_hint = bond_style
            potential_styles.insert(0, ('bond', 'hybrid '+bond_style))
            
        # Write to user the change on potential_styles
        log.out('\n\nNew potential styles to use in LAMMPS:')
        for potential, style in potential_styles:
            log.out( '  {: <8} style {}'.format( potential, style) )
        
        # Write new combination of parameters
        log.out('\n\nNew combinations of harmonic and morse bond coeffs:')   
        for i in m.bond_coeffs:
            coeffs = m.bond_coeffs[i].coeffs
            string = '  '.join([str(i) for i in coeffs])
            log.out('  {} {}'.format(i, string))
        
        
        #######################################
        # Setting up directories and where to #
        # write final files and results to    #
        #######################################
        # Find present working directory
        pwd = os.getcwd()
        
        # Find/create paths to store code results
        path = os.path.join(pwd, parent_directory)
        
        # If parent_directory == 'topofile' use topofile path as path
        if 'topofile' in parent_directory:
            log.out('Using path from topofile to set parent_directory ...')
            path = io_functions.get_dir_from_topofile(topofile, parent_directory)
        
        # Check if path exists. IF not create.
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
            
        # Change the current working directory to path1
        # to write all new files to outputs directory
        os.chdir(path)
    
        # Plotting Morse and Harmonic potentials
        basename = io_functions.get_basename(topofile, newfile=newfile, character=':', pflag=True)
        m.plotting = potential_plotting.figure(m, radius, basename, files2write, version, bondbreak_scale, ff_class, bond2style, class2xe_update)
    
        # Write lammps datafile 
        if files2write['write_datafile']:
            header = '{} > auto_morse_bond version: {}'.format(m.header, version)
            write_lmp.file(m, basename+'.data', header, atom_style, include_type_labels, log)
            
        # Write lammps datafile w/ force field parameters only
        if files2write['write_forcefield']:
            header = '{} > auto_morse_bond version: {}'.format(m.header, version)
            write_lmp.file(m, basename+'.FF.data', header, atom_style, include_type_labels, log, force_field_only=True)
        
        
        # Print file locations
        log.out(f'\n\nAll outputs can be found in {path} directory')
        
        # Print shifted morse bond note
        if not include_rcut:
            log.out('\n\n*NOTE: To use the shifted morse bond potential, modification to the LAMMPS source code')
            log.out('is required! You must subtract out d0[type] from ebond in src/MOLECULE/bond_morse.cpp *')
        else:
            log.out('\n\n*NOTE: To use the shifted morse bond potential, it is required to use the morse_bond.cpp')
            log.out('file edited by Pieter in t Veld from BASF where an optional rcut is provided by this code')
            log.out('that will then be used to automatically shift the morse potential by the energy at the')
            log.out('rcut value. The rcut value can be adjusted by the bondbreak_scale varaible.')
            
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