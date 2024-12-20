# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
April 1st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

 
    *********************************************************
    * Requirements:                                         *
    *   python 3.7+                                         *
    *                                                       *   
    * Dependencies:                                         *
    *********************************************************
"""


##############################
# Import Necessary Libraries #
##############################
import src.atom_removal.command_line as command_line
import src.atom_removal.rm_atoms as rm_atoms
import src.io_functions as io_functions
import src.write_lmp as write_lmp
import src.read_lmp as read_lmp
import time
import sys
import os


######################################################
### Main function to perform all atom-typing tasks ###
######################################################
def main(topofile, parent_directory, newfile, atom_style, atoms2remove, include_type_labels, method='atomIDs', commandline_inputs=[], log=None):
    start_time = time.time()
    
    # Configure log (default is level='production', switch to 'debug' if debuging)
    if log is None:
        log = io_functions.LUNAR_logger()
    log.configure(level='production')
    #log.configure(level='debug')
    
    #########################
    # Command Line Override #
    #########################
    # if -opt or -man option is in commandline_inputs print options and stop code execution
    if '-opt' in commandline_inputs or  '-man' in commandline_inputs:
        # call man page and exit if '-opt' or '-man' is provided at the command line
        command_line.print_man_page(topofile, parent_directory, newfile, atom_style, atoms2remove, include_type_labels, method)
        sys.exit()
        
    ###################################################################################
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters #
    ###################################################################################
    if '-opt' not in commandline_inputs and commandline_inputs:    
        # call inputs for commandline over rides
        over_rides = command_line.inputs(topofile, parent_directory, newfile, atom_style, atoms2remove, include_type_labels, method, commandline_inputs)
        
        # Set new inputs from over_rides class
        topofile = over_rides.topofile
        parent_directory = over_rides.parent_directory
        atom_style = over_rides.atom_style
        include_type_labels = over_rides.include_type_labels
        newfile = over_rides.newfile
        atoms2remove = over_rides.atoms2remove
        include_type_labels = over_rides.include_type_labels 
        method = over_rides.method


    ###########################################
    # Initialize some preliminary information #
    ###########################################
    # set version and print starting information to screen
    version = 'v1.3 / 1 April 2023'
    log.out(f'\n\nRunning atom_removal {version}')
    log.out(f'Using Python version {sys.version}')
    
    
    ############################################
    # read topofile with or without typelabels #
    ############################################
    if os.path.isfile(topofile):
        m = read_lmp.Molecule_File(topofile, method='forward', sections=['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Velocities'])
        log.out(f'Read in {m.filename} LAMMPS datafile')
    else:
        log.error(f'ERROR lammps datafile: {topofile} does not exist')
        
    ####################################
    # Remove atoms and rebuild m class #
    ####################################
    mm = rm_atoms.constructor(m, atoms2remove, log, method=method, pflag=True)
    

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

    # Write lammps datafile file
    basename = io_functions.get_basename(topofile, newfile=newfile, character=':', pflag=True, log=log)
    header = '{} > atom_removal: {} atoms2remove={}'.format(m.header, version, str(atoms2remove))
    write_lmp.file(mm, basename+'.data', header, atom_style, include_type_labels, log)
    
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