# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.7
April 14, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

 
    ***************************************************
    * Requirements:                                   *
    *   python 3.7+                                   *
    *                                                 *
    ***************************************************
"""


##############################
# Import Necessary Libraries #
##############################
import src.add_pi_electrons.convert_coeffs as convert_coeffs
import src.add_pi_electrons.pi_electrons as pi_electrons
import src.add_pi_electrons.command_line as command_line
import src.io_functions as io_functions
import src.write_lmp as write_lmp
import src.read_lmp as read_lmp
import glob
import time
import sys
import os



#############################################################
### Main function to perform all adding pi-electrons task ###
#############################################################
def main(topofile, types2convert, atom_style, reset_charges, net_zero_charge, convert2cg1, add_pi_electrons, parent_directory, newfile, 
         include_type_labels, neighbor_charge_constraint, reset_simulation_cell, commandline_inputs=[], log=None):
    
    #########################
    # Command Line Override #
    #########################
    # if -opt or -man option is in commandline_inputs print options and stop code execution
    if '-opt' in commandline_inputs or  '-man' in commandline_inputs:
        # call man page and exit if '-opt' or '-man' is provided at the command line
        command_line.print_man_page(topofile, types2convert, atom_style, reset_charges, net_zero_charge, convert2cg1, add_pi_electrons,
                                    parent_directory, newfile, include_type_labels, neighbor_charge_constraint, reset_simulation_cell)
        sys.exit()
        
    ###################################################################################
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters #
    ###################################################################################
    if '-opt' not in commandline_inputs and commandline_inputs:    
        # call inputs for commandline over rides
        over_rides = command_line.inputs(commandline_inputs, topofile, types2convert, atom_style, reset_charges, net_zero_charge, convert2cg1, add_pi_electrons,
                                         parent_directory, newfile, include_type_labels, neighbor_charge_constraint, reset_simulation_cell)
        
        # Set new inputs from over_rides class
        topofile = over_rides.topofile
        parent_directory = over_rides.parent_directory
        newfile = over_rides.newfile
        atom_style = over_rides.atom_style
        types2convert = over_rides.types2convert
        reset_charges = over_rides.reset_charges
        net_zero_charge = over_rides.net_zero_charge
        convert2cg1 = over_rides.convert2cg1
        add_pi_electrons = over_rides.add_pi_electrons
        include_type_labels = over_rides.include_type_labels
        neighbor_charge_constraint = over_rides.neighbor_charge_constraint
        reset_simulation_cell = over_rides.reset_simulation_cell
    
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
                if newfile.endswith(':') and root.startswith(newfile[:-1]): # prefix
                    log.warn(f' - WARNING matched file {file} already has newfile {newfile} extension and was skipped')
                    continue
                elif newfile.startswith(':') and root.endswith(newfile[1:]): # suffix
                    log.warn(f' - WARNING matched file {file} already has newfile {newfile} extension and was skipped')
                    continue
                try: # we dont want crashes to exit this loop
                   main(file, types2convert, atom_style, reset_charges, net_zero_charge, convert2cg1, add_pi_electrons, parent_directory,
                        newfile, include_type_labels, neighbor_charge_constraint, reset_simulation_cell, commandline_inputs, log)
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
        version = 'v1.7 / 14 April 2025'
        log.out(f'\n\nRunning add_pi_electrons {version}')
        log.out(f'Using Python version {sys.version}')
        
        
        ##########################################################
        # Read LAMMPS datafile into memory as class "m" applying #
        # forward mapping of type labels (if applicable)         #
        ##########################################################
        if os.path.isfile(topofile):
            m = read_lmp.Molecule_File(topofile, method='forward')
            log.out(f'Read in {m.filename} LAMMPS datafile')
        else:
            log.error(f'ERROR lammps datafile: {topofile} does not exist')
        basename = io_functions.get_basename(topofile, newfile=newfile, character=':', pflag=True)
        
        
        ########################################
        # Charges to set new graphene atoms to #
        ########################################
        charges = {'cg1': 0.250, 'cge': -0.125}
        
        
        ##############################################
        # Generate graph to find 1st-neighs later-on #
        ##############################################
        graph = {i:[] for i in m.atoms}
        for i in m.bonds:
            id1, id2 = m.bonds[i].atomids
            graph[id1].append(id2)
            graph[id2].append(id1)
            
            
        ###########################################################
        # Check for string type in types2convert and switch it to #
        # atomTypeID (if possible - requires comments to exist)   #
        ###########################################################
        types2remove = []; types2add = [];
        for i in types2convert:
            if type(i) is str:
                TypeIDs = []; j = str(i);
                j = j.strip()
                j = j.strip("'")
                j = j.strip('"')
                j = j.strip()
                for k in m.masses:
                    if m.masses[k].type == j:
                        TypeIDs.append(k)
                if TypeIDs: 
                    types2remove.append(i)
                    for j in TypeIDs:
                        types2add.append(j)
                        log.out('{} in types2convert was found to be a string, and maps to atomTypeID {} in {}.data'.format(i, j, basename))
                else: log.out('{} in types2convert was found to be a string, but no atomTypeIDs currently have this comment in {}.data'.format(i, basename))
        types2convert.extend(types2add)
        types2convert = list(filter(lambda i: i not in types2remove, types2convert))
            
            
        #########################################################
        # Reset charges and switch atom.comment if user desires #
        #########################################################
        if add_pi_electrons:
            log.out('add_pi_electrons and reset_charge were both True. Setting reset_charge = False, since add_pi_electrons will adjust charge.')
            reset_charges = False
            
        if reset_charges:
            log.out('Resetting cg1 charges.')
        # switch atom type comment and charge based on inputs
        for i in m.atoms:
            atom = m.atoms[i]
            if atom.type in types2convert:
                if reset_charges: atom.charge = charges['cg1'] # update charge
                if convert2cg1: atom.comment = 'cg1/C'         # update comment w/bond_react_merge comptable comment
        
        
        #########################
        # Convert coeffs to cg1 #
        #########################
        if convert2cg1:
            log.out('Converting coeffs to cg1.')
            m = convert_coeffs.cg1(m, types2convert, graph)
        
        
        ###############################################
        # Add pi-electrons and add atoms/bonds/angles #
        ###############################################
        if add_pi_electrons:
            if net_zero_charge and neighbor_charge_constraint in ['check-neighbors', 'accumulate-carbon', 'accumulate-pi-electron']:
                log.warn(f'WARNING net_zero_charge is {net_zero_charge} and neighbor_charge_constraint is {neighbor_charge_constraint}. The net_zero_charge should likely be set to False.')
            log.out('Adding pi electrons and updating coeffs/charges. Created:')
            m = pi_electrons.add(m, types2convert, charges, graph, log, neighbor_charge_constraint, reset_simulation_cell)
            
            
        #####################################################################################################
        # Make system charge net zero if user desires by adjusting charges not set by IFF cg1 and cge types #
        #####################################################################################################
        if net_zero_charge:
            natoms = 0; system_charge = 0; # Tallys to find how much charges needs to be added where
            if add_pi_electrons: types2skip = sorted([m.cg1_type, m.cge_type] + types2convert) # types to leave charge as is
            else: types2skip = sorted(types2convert) # types to leave charge as is
            log.out( 'Zeroing system net charge by adding a fixed value to atoms that DO NOT belong to the following types: {}.'.format(' '.join([str(i) for i in types2skip])) )
            for i in m.atoms:
                if m.atoms[i].type not in types2skip: natoms += 1
                system_charge += m.atoms[i].charge
            if natoms > 0:
                fixed_charge = -system_charge/natoms # Charge required to neutralize system
                log.out('    Old system net charge: {:>4.4f} that will be neutralized by adding {:>4.4f} to {} atoms.'.format(system_charge, fixed_charge, natoms))
                system_charge = 0; # Reset to tally charge again
                typesscaled = sorted([str(i) for i in m.masses if i not in types2skip])
                for i in m.atoms:
                    if m.atoms[i].type not in types2skip:
                        m.atoms[i].charge += fixed_charge
                    system_charge += m.atoms[i].charge
                log.out('    New system net charge: {:>4.4f} by adding fixed charge to the following types: {}.'.format(system_charge, ' '.join(typesscaled) ))
            else:
                log.warn('WARNING net_zero_charge could not be enforced since all atoms are defined by the atom types in types2convert = {}. System charge = {}'.format( ' '.join([str(i) for i in types2skip]), system_charge) )
            
            
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
            
        # Check if path exits. IF not create
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
                
        # Change the current working directory to path2
        # so all files get written to that directory
        os.chdir(path)
        
        
        ##############################
        # Write new lammps data file #
        ##############################
        header = '{} > add_pi_electrons version: {} w/types2convert={}'.format(m.header, version, str(types2convert))
        write_lmp.file(m, basename+'.data', header, atom_style, include_type_labels, log)
        
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