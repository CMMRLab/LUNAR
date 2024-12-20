# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.4
April 2nd, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
 
    ********************************************************
    * Requirements:                                        *
    *   python 3.7+                                        *
    *                                                      *   
    * Dependencies:                                        *
    *   correct atom-typing                                *
    ********************************************************
"""

##############################
# Import Necessary Libraries #
##############################
import src.bond_react_merge_prep.add_comments2m as add_comments2m
import src.bond_react_merge_prep.command_line as command_line
import src.bond_react_merge_prep.read_cta as read_cta
import src.io_functions as io_functions
import src.write_lmp as write_lmp
import src.read_lmp as read_lmp
import time
import sys
import os




####################################################################
### Main function to perform all topology and FF assigning tasks ###
####################################################################
def main(topofile, cta_file, newfile, atom_style, parent_directory, rm_unused_coeffs, commandline_inputs=[], log=None):
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
        command_line.print_man_page(topofile, cta_file, newfile, atom_style, parent_directory, rm_unused_coeffs)
        sys.exit()
        
    ###################################################################################
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters #
    ###################################################################################
    if '-opt' not in commandline_inputs and commandline_inputs:    
        # call inputs for commandline over rides
        over_rides = command_line.inputs(commandline_inputs, topofile, cta_file, newfile, atom_style, parent_directory, rm_unused_coeffs)
        
        # Set new inputs from over_rides class
        topofile = over_rides.topofile
        cta_file = over_rides.cta_file
        parent_directory = over_rides.parent_dir
        newfile = over_rides.newfile
        atom_style = over_rides.atom_style
        rm_unused_coeffs = over_rides.rm_unused_coeffs



    ###########################################
    # Initialize some preliminary information #
    ###########################################
    # set version and print starting information to screen
    version = 'v1.4 / 2 April 2023'
    log.out(f'\n\nRunning bond_react_merge_prep {version}')
    log.out(f'Using Python version {sys.version}')
    
    
    ####################
    # Read input files #
    ####################
    # Read lammps data file
    if os.path.isfile(topofile):
        m = read_lmp.Molecule_File(topofile, method='forward', sections=['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Velocities'])
        log.out(f'Read in {m.filename} LAMMPS datafile')
    else: log.error(f'ERROR lammps datafile: {topofile} does not exist')
        
    # Read .cta file for .data
    if cta_file == 'topofile':
        base = topofile.split('.')[0]; ext = 'cta'
        cta_file = '{}.{}'.format(base, ext)
        log.out(f'Using path and topfile name from topofile to set nta_file = {cta_file}')
    if os.path.isfile(cta_file):
        try:
            cta = read_cta.atoms(cta_file, m, log)
            log.out(f'Read in {cta_file} new-type-assignment file')
        except: log.error(f'ERROR something wrong with nta file: {cta_file}')
    else: log.error(f'ERROR .cta file: {cta_file} does not exist')
        
    
    #####################################
    # Add nessacary comments to m class #
    #####################################
    log.out('Adding necessary comments to datafile class')
    m = add_comments2m.comment(m, cta, rm_unused_coeffs, log)
        
    
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
    # Set new file name and write new datafile and version number
    basename = io_functions.get_basename(topofile, newfile=newfile, character=':', pflag=True)
       
    # write lammps datafile
    log.out('Coeff Style Hints should be transferred')
    log.out('Writing LAMMPS datafile')
    header = '{} > bond_react_merge_prep: {}'.format(m.header, version)
    include_type_labels = False
    write_lmp.file(m, basename+'.data', header, atom_style, include_type_labels, log)
    #write_lmp.datafile(basename, atom_style, m, version)
    
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
    return m