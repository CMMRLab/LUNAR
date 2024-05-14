# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.15
April 3rd, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.bond_react_merge.molecule_file_options_builder as molecule_file_options_builder
import src.bond_react_merge.auto_gen_map_file as auto_gen_map_file
import src.bond_react_merge.near_edge_charges as near_edge_charges
import src.bond_react_merge.read_input_file as read_input_file
import src.bond_react_merge.reduce_topology as reduce_topology
import src.bond_react_merge.write_molecule as write_molecule
import src.bond_react_merge.command_line as command_line
import src.bond_react_merge.check_inputs as check_inputs
import src.bond_react_merge.lmp_inscript as lmp_inscript
import src.bond_react_merge.merge_coeffs as merge_coeffs
import src.bond_react_merge.write_mol2 as write_mol2
import src.io_functions as io_functions
import src.write_lmp as write_lmp
import src.read_lmp as read_lmp
import collections
import time
import sys
import os

       
##################################################
### Main function to perform all merging tasks ###
##################################################
def main(files, parent_dir, newfile, atom_style, generate_map_file, write_rxn_mol2files, write_rxn_datafiles,
         write_moleculefiles, print_options, commandline_inputs, map_near_edge_rxn_charges, molecule_file_options,
         include_type_labels, log=io_functions.LUNAR_logger()):
    start_time = time.time()
    
    # Configure log (default is level='production', switch to 'debug' if debuging)
    log.configure(level='production')
    #log.configure(level='debug')
    
    # set version and print starting information to screen
    version = 'v1.15 / 3 April 2024'
    log.out(f'\n\nRunning bond_react_merge {version}')
    log.out(f'Using Python version {sys.version}')


    #########################
    # Command Line Override #
    #########################
    # if -opt or -man option is in commandline_inputs print options and stop code execution
    if '-opt' in commandline_inputs or  '-man' in commandline_inputs or print_options:
        # call man page and exit if '-opt' or '-man' is provided at the command line
        command_line.print_man_page(files, parent_dir, newfile, atom_style, generate_map_file, write_rxn_mol2files,
                                    write_rxn_datafiles, write_moleculefiles, print_options, map_near_edge_rxn_charges,
                                    molecule_file_options, include_type_labels)
        sys.exit()
        
    ###################################################################################
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters #
    ###################################################################################
    if '-opt' not in commandline_inputs and commandline_inputs:    
        # call inputs for commandline over rides
        over_rides = command_line.inputs(files, parent_dir, newfile, atom_style, generate_map_file, write_rxn_mol2files,
                                         write_rxn_datafiles, write_moleculefiles, print_options, map_near_edge_rxn_charges,
                                         molecule_file_options, include_type_labels, commandline_inputs)
        
        # Set new inputs from over_rides class
        files = over_rides.files
        parent_dir = over_rides.parent_dir
        newfile = over_rides.newfile
        atom_style = over_rides.atom_style
        generate_map_file = over_rides.generate_map_file
        write_rxn_mol2files = over_rides.write_rxn_mol2files
        write_rxn_datafiles = over_rides.write_rxn_datafiles
        write_moleculefiles = over_rides.write_moleculefiles
        print_options = over_rides.print_options
        map_near_edge_rxn_charges = over_rides.map_near_edge_rxn_charges
        molecule_file_options = over_rides.molecule_file_options
        include_type_labels = over_rides.include_type_labels

        
    ######################################################
    # Check for valid inputs; if not valid exit the code #
    ######################################################
    if not check_inputs.safety(files, parent_dir, atom_style, generate_map_file, write_rxn_mol2files,
                               write_rxn_datafiles, write_moleculefiles, print_options, map_near_edge_rxn_charges,
                               commandline_inputs, log):
        log.error('ERROR safety check failed')
        
        
    #####################################################################################
    # Check if inputfile is in files, if so regenerate files dictionary from input file #
    #####################################################################################
    if 'infile' in files:
        files, parent_directory = read_input_file.file_inputs(files['infile'], log)
        if parent_directory != '': parent_dir = parent_directory # update parenet directory if applicable

    
    ##############################
    # Read in all files in files #
    ##############################
    merge = {}; filenames = [];
    for i in files:
        # create dictionary
        filenames.append(files[i]);
        
        # find file and add to merge dictionary
        if os.path.isfile(files[i]):
            merge[i] = read_lmp.Molecule_File(files[i], method='forward', sections=['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Velocities'])
        else: log.error(f'ERROR lammps datafile: {files[i]} does not exist')
        
    # Check for duplicate filenames
    duplicates = [item for item, count in collections.Counter(filenames).items() if count > 1]
    if duplicates:
        for i in duplicates:
            log.error(f'ERROR filename: {i} is used more then once in')
    
    # Print start of table of file tag to file
    log.out('\n\n---------------------------------------------------------------------------------------------')
    log.out('|                read in files key to value dictionary map to organize files                |')
    log.out('---------------------------------------------------------------------------------------------')
    log.out('| {:^16} | {:^70} |'.format('key to', 'value to dictionary'))
    log.out('| {:^16} | {:^70} |'.format('dictionary', 'filename'))
    log.out('---------------------------------------------------------------------------------------------')
    for i in merge:
        log.out('| {:^16} | {:^70} |'.format(i, os.path.basename(merge[i].filename)[0:70])) # 70 character limit and then reject the rest
    log.out('---------------------------------------------------------------------------------------------\n\n')
    
    
    ##################################
    # Reduce templates if applicable #
    ##################################
    template_pairs = reduce_topology.find_rxn_pairs(merge) # Find template pairs
    pairids = sorted(list(template_pairs.keys())) # sort such that pairID 1 is found before pairID 2 ...
    for i in pairids:
        pair = sorted(template_pairs[i], reverse=True) # sort such that list will be order as ['preN', 'postN']
        if len(pair) == 2:
            pre = merge[pair[0]]; post = merge[pair[1]];
            BondingIDs, CreateIDs, Reduce, Remove, Keep = auto_gen_map_file.get_lmp_header_info(pre, log)
            if len(Reduce) in [2, 5] or Remove or Keep:
                pre_reduced, post_reduced = reduce_topology.template(pre, post, BondingIDs, CreateIDs, Reduce, Remove, Keep, log)
                if pre_reduced != '' and post_reduced != '' and pre_reduced.natoms == post_reduced.natoms:
                    merge[pair[0]] = pre_reduced; merge[pair[1]] = post_reduced; print();
                if pre_reduced == '' or post_reduced == '': log.error('ERROR template reduce failed due to unknown cause')
            if len(Reduce) not in [0, 2, 5]: log.error('ERROR reduce option used, but length of list was not 2 or 5')
    
    
    ############################
    # Call merged class as new #
    ############################
    new = merge_coeffs.merged(merge, log)
    
    # Print out new coeff types in debug mode
    merge_coeffs.print_merged(new, log, skip_printing_cross_terms=True)

    
    #######################################
    # Setting up directories and where to #
    # write final files and results to    #
    #######################################
    # Find present working directory
    pwd = os.getcwd()
    
    # Find/create paths to store code results
    path = os.path.join(pwd, parent_dir)
    
    # If parent_directory == 'topofile' use topofile path as path
    if 'topofile' in parent_dir:
        log.out('Using path from 1st file in files dictionary to set parent_directory ...')
        topofile = list(files.values())[0]
        path = io_functions.get_dir_from_topofile(topofile, parent_dir)

    
    # Check if path exits. IF not create
    if not os.path.isdir(path): 
        os.makedirs(path, exist_ok=True)
        
    # Change the current working directory to path so all files get written to that directory
    os.chdir(path)


    ###################################################################
    # Try auto generating the map file based on a cost fitting method #
    ###################################################################
    # Function to test is var is False or an Integer for map_near_edge_rxn_charges var
    def check_var_for_int(var):
        string_var = str(var); return_boolean = False
        try: 
            int(string_var)
            return_boolean = True
        except: pass
        return return_boolean
    if generate_map_file:       
        log.out('\n\n************************************************************************************************************')
        log.out('************************************************************************************************************')
        log.out('************************************************************************************************************')
        log.out('***  Using generate_map_file, which auto-generates the map file based on tag-ID for preN/postN pairings. ***')
        log.out('***  PLEASE CONFIRM preN tagged file is consistent with postN tagged file. Each section is set by the    ***')
        log.out('***  following methods:                                                                                  ***')
        log.out('***      InitiatorIDs                                                                                    ***')
        log.out('***          Will be initialized as "ID1" and "ID2", then if two molecues are found in the preN tagged   ***')
        log.out('***          file and edgeIDs were found, a 1st attempt at setting the InitiatorIDs is to search 4-6     ***')
        log.out('***          neighbors away from all edgeIDs and select to most unique atom-type in each molecule. If    ***')
        log.out('***          the code can find all equivalences (which are assumed to be correct), a 2nd attempt is then ***')
        log.out('***          tried where the code applies a backwards map of all atomIDs in the postN tagged file to the ***')
        log.out('***          preN tagged file. After the backwards mapping the code tries to find the bond responsible   ***')
        log.out('***          that links the preN tagged-file two molecules together in the postN tagged-file. If this    ***')
        log.out('***          "linking" bond can be found, the InitiatorIDs are updated with to the new "linking" atomIDs.***')   
        log.out('***      EdgeIDs                                                                                         ***')
        log.out('***          Set by searching for the number of bonded atoms to each atom type in the dataN tagged files ***')
        log.out('***          and comparing them to the preN/postN tagged files, if the number of bonded atoms for each   ***')
        log.out('***          atom type is found to be different IT MUST BE AN EDGE ATOM (assumes correct atom-typing,    ***')
        log.out('***          thus if no edge atoms are found it is an indication that atom types are not consistent      ***')
        log.out('***          between the dataN tagged files and the preN tagged files).                                  ***')
        log.out('***      DeleteIDs                                                                                       ***')
        log.out('***           Will automatically be found by the code if the code can find all equivalences (which are   ***')
        log.out('***           assumed to be correct). The code will warn the user that deleteIDs were detected which are ***')
        log.out('***           determined by if the postN tagged files has more then one molecular fragment, the largest  ***')
        log.out('***           frament is kept and the smaller ones are are assumed to be deleteIDs. Once the smaller     ***')
        log.out('***           fragments are found the atomIDs are mapped backwards onto the preN tagged atomIDs. Thus if ***')
        log.out('***           the Equivalences are found to be wrong its possible that the DeleteIDs are wrong, however  ***')
        log.out('***           in the *_commented.txt map file there will be a comment on the atomID it is in the postN   ***')
        log.out('***           tagged file which can be used to help reset the DeleteIDs if needed during checking of the ***')
        log.out('***           generated map file.                                                                        ***')
        log.out('***      Equivalences                                                                                    ***')
        log.out('***           Set by a pseudo-graph-theory/cost-based algorithm were each preN/postN tagged pair are     ***')
        log.out('***           converted to graphs. The atom-types, number of bonded neighbors (nb), elements, atomIDs,   ***')
        log.out('***           molID (reset by code), molecule formula, and relative location of atom in the template is  ***')
        log.out('***           used to assign a map cost. Then an incremented approach of slowly increase an acceptable   ***')
        log.out('***           cost is then used where the lower cost mappings are assigned first (best matches) and the  ***')
        log.out('***           atomIDs that have the worst map are set last (due to creation of new topologies).          ***')
        log.out('***  In general, it is up to the user to check all information in the generated map files is correct.    ***')
        log.out('***  NOTE: terminating atoms like H do NOT require perfect mapping since nothing is bonded "through"     ***')
        log.out('***  them. The code is designed to try getting everything correct, but there are scenerios where the     ***')
        log.out('***  code may not be able to due to the complexity present in the generation of the map files.           ***')
        log.out('************************************************************************************************************')
        log.out('************************************************************************************************************')
        log.out('************************************************************************************************************')

        
        # Find pre/post pairing to auto the generate map file
        template_pairs = {} # {1:[pre1, post1], 2:[pre2, post2], ... N:[] }
        dataN_tagged = {} # { file-tag : molecule class }
        postN_tagged = {} # { file-tag : molecule class }
        preN_tagged = {} # { file-tag : molecule class }
        dataN_types2nb_lst = {} # { atom-type: [lst of nb atoms to type]}
        dataN_types = set(); postN_types = set(); preN_types = set(); # 2 check if any preN tagged files has new types compared to dataN tagged file
        for file in merge:    
            # strip last number of filename. Assumes tags either:
            # - data1, data2, ... dataN  ->  data
            # - pre1, pre2, ... preN     ->  pre
            # - post1, post2, ... postN  ->  post
            tmpname = ''.join([i for i in file if i.isalpha()])
            
            # String begining tag call from file. Assumes tags either:
            # - data1, data2, ... dataN  ->  1, 2, ... N
            # - pre1, pre2, ... preN     ->  1, 2, ... N
            # - post1, post2, ... postN  ->  1, 2, ... N
            tmpid = ''.join([i for i in file if i.isdigit()])
            
            # Log pairs if tmpid exists, else create key/valued list if tmpname is not data
            if tmpname != 'data':
                if tmpid in template_pairs:
                    template_pairs[tmpid].append(file)
                else: template_pairs[tmpid] = [file]
                
            # Find dataN_types2nb_lst from datafiles
            if tmpname == 'data':
                dataN_types2nb_lst = auto_gen_map_file.dataN_type2nb_map(merge[file], dataN_types2nb_lst)
                dataN_tagged[file] = merge[file]
                
                # Get dataN tagged types
                for i in merge[file].masses:
                    dataN_types.add(merge[file].masses[i].type)
                
            # Find postN file and add to postN_tagged
            if tmpname == 'post':
                postN_tagged[file] = merge[file]
                
                # Get postN tagged types
                for i in merge[file].masses:
                    postN_types.add(merge[file].masses[i].type)
                
            # Find preN file and get atom types
            if tmpname == 'pre':
                preN_tagged[file] = merge[file]
                
                # Get preN tagged types
                for i in merge[file].masses:
                    preN_types.add(merge[file].masses[i].type)
        
        # No preN file should have newer types then dataN or postN tagged file
        if preN_types:
            for i in preN_types:
                if i not in dataN_types and i not in postN_types:
                    log.warn(f'WARNING atom type {i} in a preN tagged file does not exist in any dataN or postN tagged file.')
                    for filetag in preN_tagged:
                        for j in preN_tagged[filetag].masses:
                            if preN_tagged[filetag].masses[j].type == i:
                                log.warn('  offending filetag {} -> {}'.format(filetag, preN_tagged[filetag].filename))
                    sys.exit()
                if i not in dataN_types and i in postN_types:
                    log.warn(f'WARNING atom type {i} in a preN tagged file does not exist in any dataN tagged file, but exists in a postN tagged file.')
                
        # Find most frequent number of bonded atoms to each atom type in dataN_types2nb_lst to use for
        # detecting edge atoms in molecule file (Terminating elements like H, F, ... will also be used)
        dataN_types2nb_map = {} # { atom-type: most-freqeunty-nb (from all read in datafiles) }
        for atomtype in dataN_types2nb_lst:
            nb_lst = dataN_types2nb_lst[atomtype]
            dataN_types2nb_map[atomtype] = set(nb_lst)
            
        # Write generalized LAMMPS input file for fix bond/react
        if template_pairs:
            lmp_inscript.bond_react('in.fix_bond_react.script', newfile, version, merge, atom_style, new, pairids, template_pairs, log)
                    
        # Loop through template_pairs and try finding map
        pairids = sorted(list(template_pairs.keys())) # sort such that pairID 1 is found before pairID 2 ...
        for i in pairids:
            # sort such that list will be order as ['preN', 'postN']
            pair = sorted(template_pairs[i], reverse=True)
            
            # Only try if length of pair is two
            if len(pair) == 2:
                pre = merge[pair[0]]; post = merge[pair[1]];
                pre_rxnid = ''.join([i for i in pair[0] if i.isdigit()])
                post_rxnid = ''.join([i for i in pair[1] if i.isdigit()])
                if pre_rxnid == post_rxnid:
                    filename = '{}-{}_rxn-map'.format(pair[0], pair[1])
                    
                    # Print header for rxn-map
                    log.out('\n\n')
                    log.out('****************************************************')
                    log.out('* Finding rxn map/super-impose file for rxnid: {:^3} *'.format(pre_rxnid))
                    log.out('****************************************************')
                    
                    # Find pre2post class
                    pre2post = auto_gen_map_file.find(pre, post, filename, pair, dataN_types2nb_map, log)
                    log.out('  {} -> {} rxn map equivalences converged in {} iterations with a max cost of {}\n'.format(pair[0], pair[1], pre2post.iterations, pre2post.max_cost))
                    
                    # Find new title and write file
                    title = 'Map file: bond_react_merge: {} auto-generated mapfile for {}  ->  {}    rxnid: {}'.format(version, os.path.basename(pre.filename), os.path.basename(post.filename), pre_rxnid)
                    auto_gen_map_file.write_map(filename+'_commented.txt', pre2post, title, version, comment_flag=True)
                    auto_gen_map_file.write_map(filename+'_uncommented.txt', pre2post, title, version, comment_flag=False)
                    
                    
                    # Apply any molecule_file_options if lst is not empty
                    if molecule_file_options:
                        try: pre, post = molecule_file_options_builder.add2molecules(molecule_file_options, pre, post, pre2post, newfile, log)
                        except: log.warn('    WARNING something went wrong with molecule file option usage.')
                
                    # Try finding edge atom mapping to update charges, by determining if var is an intger and set ndepth_update by in value
                    if check_var_for_int(map_near_edge_rxn_charges): 
                        try:
                            # Only attempt it there are at least 1-edgeID found and at least 1-equivalence found
                            if len(pre2post.pre_edge) > 0 and len(pre2post.map) > 0:
                                # Find edge atom charges and create new attributes in atoms dict class:
                                #    - mapped_charge = x.xxxx
                                #    - mapped_comment = coment
                                log.out('  Attempting to find edge atom equivalences from dataN or postN tagged files to update edge atom charge(s) or near edge atom charge(s)')
                                fileN_tagged = {**dataN_tagged, **postN_tagged}  # Join dataN_tagged and postN_tagged classes such that dataN will be attempted 1st. This allows for reactions to be built on top of one another
                                near_edge_charges.find(pre, post, pre2post, fileN_tagged, log, ndepth_update=map_near_edge_rxn_charges)
                            else:
                                log.warn('    WARNING map_near_edge_rxn_charges was used to try updating edge atom (or near edge atom) charges, but')
                                log.out('    NO edge atomIDs or NO equivalences could be found from the generate_map_file option. Skipping usage')
                                log.out('    of map_near_edge_rxn_charges, since the code would not have anything to map since it does understand')
                                log.out(f'    the edge atoms or equivs in your {os.path.basename(pre.filename)} -> {os.path.basename(post.filename)}')
                                log.out('    rxn files.')
                        except:
                            log.out('')
                            log.warn('    WARNING map_near_edge_rxn_charges failed for an unknown reason. Most likely reason is that a dataN tagged')
                            log.out('    file had multiple molecules in it and the code assumes each dataN tagged file only has a single molecule in')
                            log.out('    in. Edge atom charges will have to be updated by the user manually for these set of files.') 
                else: log.warn('    WARNING generate_map_file was used but code messed up the logic somewhere')
            else: log.warn(f'WARNING generate_map_file was used but N-tag: {i} does NOT have exactly two pairs: {str(pair)}')
    
    
    ###################################################################################################################
    # Write outputted files (change coeff numbering on the fly when writing each file - using read in comments from   #
    # read_lmp (string types) and map to new numeric types using *_types_map; *=atom, bond, angle, dihedral, and      #
    # improper to set new type). All coeff types will be inserted into every datafile                                 #
    ###################################################################################################################
    # Write generalized LAMMPS input file for create_atoms
    if write_moleculefiles:
        lmp_inscript.create_atoms('in.create_atoms.script', newfile, version, merge, atom_style, new, log)
    filenames = []
    for file in merge:
        m = merge[file]
        m = merge_coeffs.update_TypeIDs(m, new, log) # Update coeffs before writing
        header = '{} > bond_react_merge: {} datafile (filetag: {})'.format(m.header, version, file)
        
        # Find old name and make new name
        oldname = os.path.basename(m.filename)
        filenames.append(oldname)
        if ':' in newfile or newfile == '':
            newname = io_functions.get_basename(oldname, newfile=newfile, character=':', pflag=True)
        else: newname = ('{}{}'.format(file, newfile))
        
        # strip last number of filename. Assumes tags either:
        # - data1, data2, ... dataN  ->  data
        # - pre1, pre2, ... preN     ->  pre
        # - post1, post2, ... postN  ->  post
        tmpname = ''.join([i for i in file if i.isalpha()])
    
        # Write remaing datafile info (w/o or w/ coeffs depending on coeff_flag)
        if tmpname == 'data' or tmpname == 'Data':
                write_lmp.file(m, newname+'.data', header, atom_style, include_type_labels, log)
                if write_moleculefiles:
                    write_molecule.file(newname, m, new, file, version, include_type_labels)
                    write_lmp.file(m, 'force_field.data', header, atom_style, include_type_labels, log, force_field_only=True)
                    
                
        # Write molecule files
        elif tmpname == 'pre' or tmpname == 'Pre' or tmpname == 'post' or tmpname == 'Post':
            write_molecule.file(newname, m, new, file, version, include_type_labels)
            if write_rxn_datafiles: write_lmp.file(m, newname+'.data', header, atom_style, include_type_labels, log)
            if write_rxn_mol2files: write_mol2.file(newname, m, file, version)
    
    
    # Print warning about charges on edge atoms
    if check_var_for_int(map_near_edge_rxn_charges):
        if generate_map_file:
            log.out('\n\n****************************************************************************************************************************** ')
            log.out('*** map_near_edge_rxn_charges was used. Please check all molecule templates to confirm that charges were mapped properly!  *** ')
            log.out('****************************************************************************************************************************** ')
        else:
            log.out('\n\n**************************************************************************************************************************************************** ')
            log.out('*** map_near_edge_rxn_charges was used, but generate_map_file was not. Thus the molecule template may need to have the edge atom charges updated *** ') 
            log.out('**************************************************************************************************************************************************** ')
    else:
        log.out('\n\n************************************************************************************************  ')
        log.out('*** Only modification to molecule template is to update charges on edge atoms if not updated *** ')
        log.out('************************************************************************************************ ')

    
    # Print file locations
    log.out(f'\n\nAll outputs can be found in {path} directory')
    
    
    # Print completion of code
    log.out('\n\nNormal program termination\n\n')
    
    # Script run time
    execution_time = (time.time() - start_time)
    log.out('Execution time in seconds: ' + str(execution_time))
    
    # Show number of warnings and errors
    log.out_warnings_and_errors()
    
    # write log (set basename based on first and last file, if only one file set as single file basename)
    if len(filenames) >= 2:
        file1 = filenames[0]
        file2 = filenames[-1]
        root1 = file1[:file1.rfind('.')]
        root2 = file2[:file2.rfind('.')]
        basename = '{}-{}-nfiles={}'.format(root1, root2, len(filenames))
    else:
        file1 = filenames[0]
        root1 = file1[:file1.rfind('.')]
        basename = '{}-nfiles={}'.format(root1, len(filenames))
    log.write_logged(basename+'.log.lunar')
    
    # Change back to the intial directory to keep directory free for deletion
    os.chdir(pwd)
    return