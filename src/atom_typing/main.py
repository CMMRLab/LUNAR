# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.9
January 5th, 2024
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
import src.atom_typing.Gasteiger.Gasteiger_algorithm as Gasteiger_algorithm
import src.atom_typing.merge_input_files as merge_input_files
import src.atom_typing.cluster_analysis as cluster_analysis
import src.atom_typing.ring_analysis as ring_analysis
import src.atom_typing.hybridization as hybridization
import src.atom_typing.coordination as coordination
import src.atom_typing.check_inputs as check_inputs
import src.atom_typing.command_line as command_line
import src.atom_typing.out2console as out2console
import src.atom_typing.atom_info as atom_info
import src.atom_typing.write_lmp as write_lmp
import src.atom_typing.write_nta as write_nta
import src.atom_typing.write_pdb as write_pdb
import src.io_functions as io_functions
import time
import sys
import os



######################################################
### Main function to perform all atom-typing tasks ###
######################################################
def main(topofile, bondfile, parent_directory, newfile, ff_name, delete_atoms, mass_map, bondorder, maxbonded, boundary,
         vdw_radius_scale, reset_charges, print_options, commandline_inputs, bonds_via_distance_override, pdb_file, chargefile,
         log=io_functions.LUNAR_logger()):  
    start_time = time.time()
    
    # Configure log (default is level='production', switch to 'debug' if debuging)
    log.configure(level='production')
    #log.configure(level='debug')
    
    #########################
    # Command Line Override #
    #########################
    # if -opt or -man option is in commandline_inputs print options and stop code execution
    if '-opt' in commandline_inputs or '-man' in commandline_inputs or print_options:
        # call man page and exit if '-opt' or '-man' is provided at the command line
        command_line.print_man_page(topofile, bondfile, parent_directory, newfile, ff_name, delete_atoms,
                                    mass_map, bondorder, maxbonded, boundary, vdw_radius_scale, reset_charges,
                                    print_options, pdb_file, chargefile)
        sys.exit()
    
    ###################################################################################
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters #
    ###################################################################################
    if '-man' not in commandline_inputs and commandline_inputs:    
        # call inputs for commandline over rides
        over_rides = command_line.inputs(topofile, bondfile, parent_directory, newfile, ff_name, delete_atoms, mass_map, bondorder, maxbonded,
                                         boundary, vdw_radius_scale, reset_charges, print_options, pdb_file, chargefile, commandline_inputs)
        
        # Set new inputs from over_rides class
        topofile = over_rides.topofile
        bondfile = over_rides.bondfile
        parent_directory = over_rides.parent_directory
        newfile = over_rides.newfile
        ff_name = over_rides.ff_name
        delete_atoms = over_rides.delete_atoms
        mass_map = over_rides.mass_map
        bondorder = over_rides.bondorder
        maxbonded = over_rides.maxbonded
        boundary = over_rides.boundary
        vdw_radius_scale = over_rides.vdw_radius_scale
        reset_charges = over_rides.reset_charges
        print_options = over_rides.print_options
        pdb_file = over_rides.pdb_file
        chargefile = over_rides.chargefile
    
    
    ######################################################
    # Check for valid inputs; if not valid exit the code #
    ######################################################
    if not check_inputs.safety(topofile, bondfile, parent_directory, newfile, ff_name, delete_atoms, mass_map, bondorder,
                               maxbonded, boundary, vdw_radius_scale, reset_charges, print_options, log):
        log.error('ERROR safety check failed')


    ###########################################
    # Initialize some preliminary information #
    ###########################################
    # set version and print starting information to screen
    version = 'v1.9 / 5 January 2024'
    log.out(f'\n\nRunning atom_typing {version}')
    log.out(f'Using Python version {sys.version}')
    log.out(f'Assigning {ff_name} atom-types')
    
    
    ######################################################################################################
    # Read .data file or .mol or .mol2 or .sdf file based on extension and add following information:    #
    #   for .mol or .sdf or .mol2:                                                                       #
    #      - m.masses instance: {atomtype : coeffs class}; coeffs class .type .coeffs                    #
    #      - m.atoms[ID] instance: .type = atomtype set by element; .comment set by element              #
    #                                                                                                    #
    #   for .data w/bonds:                                                                               #
    #      - m.atoms[ID] instance .element = element; .comment = element                                 #
    #                                                                                                    #
    #   for .data w/o bonds (LAMMPS ReaxFF, use bondfile):                                               #
    #      - m.atoms[ID] instance .element = element; .comment = element                                 #
    #      - m.bonds[ID] instance .type = 1 (arbitrary) .atomids = [id1, id2] (from BO file and cutoffs) #
    #      - m.nbonds = number of bonds found from BO file and bond-order cutoffs applied                #
    #      - m.nonbdtypes = 1 (arbitrarly set every bond type as 1 since this is just a place holder)    #
    ######################################################################################################
    m = merge_input_files.merge(topofile, bondfile, mass_map, bondorder, maxbonded, boundary, vdw_radius_scale, reset_charges, bonds_via_distance_override, log)
    
    
    ###################################################################################
    # Find clusters and add cluster info to m. Will add the following instances to m: #
    #   m.molecules -> .data .atoms .clusters .formula .molids                        # 
    #   m.atoms[ID] -> .molecule -> .mass .size .formula                              #
    ###################################################################################
    m = cluster_analysis.add_molecule_data2m(m, log)
    
    
    ###################################################################################################
    # Extract clusters and rebuild m class as mm with all necessary information and contigous atomids #
    # with for .atoms and .bonds, transfering box info, header, filename, reaxff info, molecule info  #
    # .... method='extract' to find large clusters, method='by-product' to find by-products           #
    ###################################################################################################
    mm = cluster_analysis.cluster_removal(m, log, delete_atoms, method='extract') # mm=large kept molecules
    bp = cluster_analysis.cluster_removal(m, log, delete_atoms, method='by-product') # bp=removed by-products
    
    #####################################################################################################
    # FROM HERE ON OUT ONLY ANALYZING mm CLASS SINCE IT CONTAINS THE ATOMS WERE CARE ABOUT              #
    # Find and add neighbor_ids information to mm.atoms as mm.atoms[ID].neighbor_ids; were neighbor_ids #
    # is a dictionary as:                                                                               #
    #    {1st-neighs,   2nd-neighs,   Nth-neighs,}                                                      #
    #    {1: [1, 2, 3], 2: [4, 5, 6], ....}                                                             #
    #    Use a modified BFS algorithm to find Nth neighbors away from each atom to find neighboring IDs #
    #    set the Nth_neigh_depth as 4 and find up to 4 neighs away from atomid (adjust as needed)       #
    #####################################################################################################
    mm = coordination.neighbors(mm, Nth_neigh_depth=4)
    
    
    ########################################################################################################
    # Find and add ring data to mm class. ring data will be found via information in find_rings dictionary #
    # below; where 'elements2walk' sets the element types that can belong in the ring and 'rings2check are #
    # sizes of rings to look for. Code provided by Jake Gissinger with modifications and additions by Josh #
    # Then add .rings .ring .ringID .ringformula to mm.atoms[ID]                                           # 
    # .rings = [lst or ringsizes atom belongs too]                                                         #
    # .ring = integer value of partioned ring                                                              #
    # .ringID = ringID (zero if not logical)                                                               #
    # .ringformula = ring formula (blank if not logical)                                                   #
    ########################################################################################################
    #supported_elements = ['C', 'H', 'O', 'N', 'S', 'F', 'Si','Xe', 'Ne', 'Kr', 'He', 'D', 'Cl', 'Ca', 'Br', 'Ar', 'P']
    find_rings = {'elements2walk': mm.elements,    # List of elements to walk (recommended to walk along all elements in system).
                  'rings2check': [3, 4, 5, 6, 7],  # List of ring sizes to check for in main ring analysis code.
                  'fused-rings': True,             # True or False to run ring connectivty analysis ('perform' key must be True).
                  'fused2check': [5, 6]}           # List of ring sizes to check for in fused ring analysis code.
    # Shut off fused ring if not used with certain FFs, b/c it is very time consuming to run
    no_fusedrings = ['CVFF', 'CVFF-IFF', 'Clay-FF', 'compass', 'general:0', 'general:1', 'general:2', 'general:3', 'general:4', 'OPLS-AA', 'PCFF-IFF', 'PCFF']
    if ff_name in no_fusedrings and len(mm.atoms) > 1000:
        find_rings['fused-rings'] = False
    elif ff_name not in no_fusedrings and len(mm.atoms) >= 1000:
        log.out(f'   WARNING very large system and FF: {ff_name} requires a fused ring analysis which may significantly slow down the code ...')
    mm = ring_analysis.add_ring_data2m(mm, find_rings, log)
    
    
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
    mm = atom_info.add(mm)

    
    ####################################################################################################################
    # Find hybrdization and add to mm.atoms[ID].hybridization. hybridization will be intialized as 'Unknown' and if    #
    # found will be set as 'Sp0' or 'Sp1' or 'Sp2' or 'Sp3'. Where the meanings are as follows:                        #
    #    'Sp1' is a linear geometry set by searching for angles around 180 degrees                                     # 
    #    'Sp2' is a Trigonal Planar geometry set in the order 1.) ring size 2.) search for angles around 120 degrees   #
    #    'Sp3' is a Tetrahedral geometry set in the order of 1.) C/N nb==4 2.) search for angles around 109.5          #
    #    'Terminal' is an atom that is terminating and does not fall into the Sp1/Sp2/Sp3 categories domain (H, F, ..) #
    #'   'unknown' is the intialized hybridization state and is the default if the hybridization could not be found    #
    #                                                                                                                  #
    # Also add mm.atoms[ID].avg_angle to atoms dict, since this maybe useful outside of the hybridization              #
    # characterization. If the atom is not in the center of an angle this will be set to zero.                         #
    #                                                                                                                  #
    # The method above is used because it can not be assumed that charge is set and correct nor can it be assumed      #
    # that we have access to lone pairs, thus the hybridization will be set via by topological and geometric           #
    # considerations. It will 1st be attempted by topological considerations, because it is independant of atomic      #
    # positions and I usually catch individuals performing MD simulations starting out with poor geometries (like      #
    # having a molecule that has Sp3 hybridized atoms, but only having the molecule with 2D coordinates, even though   #
    # modern chemical drawing tools have the ability to intialize in 3D ....). The main purpose of finding the         #
    # hybridization state before atom-typing is for force fields like DREIDING that rely on the known hybridization    #
    # state to make atom-typing easier. DREIDING also posses another issue and that is standard topolical              #
    # considerations are not enough when implicit Hydrogens are used. IF any atom-typing using this hybridization      #
    # info to perform atom-typing the USER WILL BE WARNED that imporper atomic positions may result in incorrect       #
    # atom types being used. It will also serve as a useful starting point for a hand built Gasteiger Charge Method    #
    # since this info is needed to set the intial Gasteiger parameters.                                                #
    ####################################################################################################################
    mm = hybridization.partially_topological_partially_geometric(mm, log) # Find per atom hybridization
    mm = hybridization.hybridization_data(mm) # Tally and find global hybridization data table

    
    
    #####################################################################
    # Find Gasteirger charges and update mm.atoms[atomID].charge object #
    #####################################################################
    if reset_charges:
        mm = Gasteiger_algorithm.compute_charges(mm, chargefile, log)
    
    
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
        
    # Change the current working directory to path
    # to write all new files to outputs directory
    os.chdir(path)
    
    
    #######################################################################################
    # Perform atom-typing based of ff_name. Will add the following info to mm.atoms[ID]:  #
    #     .nta_type = atom-type that was determined (nta=new-type-assignment)             #
    #     .nta_info = comment for how .nta type was set or found                          #
    # Will also add .supported_types to mm class; were supported_types = {'element': [lst #
    # of atom-types that belong to that element]}; IE. {'C': ['cp', 'c'], 'H': ['hc']}    #
    #######################################################################################
    if 'PCFF-IFF' in ff_name:
        import src.atom_typing.typing.PCFF_IFF as typing
    elif 'PCFF' in ff_name:
        import src.atom_typing.typing.PCFF as typing
    elif 'DREIDING' in ff_name:
        import src.atom_typing.typing.DREIDING as typing
    elif 'compass' in ff_name:
        import src.atom_typing.typing.compass as typing
    elif 'CVFF-IFF' in ff_name:
        import src.atom_typing.typing.CVFF_IFF as typing
    elif 'CVFF' in ff_name:
        import src.atom_typing.typing.CVFF as typing
    elif 'OPLS-AA' in ff_name:
        import src.atom_typing.typing.OPLS_AA as typing
    elif 'Clay-FF' in ff_name:
        import src.atom_typing.typing.Clay_FF as typing
    elif 'general' in ff_name:
        import src.atom_typing.typing.general as typing
    else: log.error(f'ERROR requesting unsupported Force Feild name: {ff_name}')
    
    # Start performing atom-typing based on force field typing module imported
    basename = io_functions.get_basename(topofile, newfile=newfile, character=':', pflag=True)
    mm = typing.nta(mm, basename, ff_name)

               
    ####################
    # Print to console #
    ####################
    out2console.out(mm, bp, log, version, ff_name)

    
    ##################
    # Write lmp file #
    ##################
    filename  = '{}.data'.format(basename)
    write_lmp.file(mm, filename, version, ff_name) # Write lmp .data file
    
    
    ##################
    # Write nta file #
    ##################
    filename  = '{}.nta'.format(basename)
    write_nta.file(mm, filename, version, ff_name) # Write lmp .data file
    
    
    ##################
    # Write pdb file #
    ##################
    if pdb_file in ['types', 'typeIDs']:
        if len(mm.atoms) <= 99_999:
            filename  = '{}_packmol'.format(basename)
            write_pdb.file(mm, filename, version, ff_name, pdb_file, log) # Write .pdb file
        else: log.warn('WARNING could not write .pdb file since number of atoms > 99,999.')
    elif pdb_file == 'skip':
        log.out(f'pdb_file = {pdb_file}, so skipping writing of optional .pdb file')
    else: log.error(f'ERROR unsupported pdb_file input {pdb_file}')
    
        
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
    return mm