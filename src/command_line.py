# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
September 7th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

File contains command line options/man page and
command line over riding class to reset hard coded
inputs.
"""

##############################
# Import Necessary Libraries #
##############################
import sys


#################################
# For single file clusters code #
#################################
def clusters_single(topofile, N0, txtfile, fav, commandline_inputs):
    # if -opt or -man option is in commandline_inputs print options and stop code execution
    if '-opt' in commandline_inputs or  '-man' in commandline_inputs:
        print('\n\ncluster_analysis has been run with -opt or -man option to show the optional command line overrides available. Command line')
        print('option summary [-tag <tag-input>]:')
        print('python3 cluster_analysis.py [-topo <topo-filename>] [-txt <T|F>] [-n0 <int>] [-fav <float>] <-gui> <-opt>|<-man>')


        print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the cluster_analysis.py')
        print('file will be enforced. The command line inputs option gives more flexibility to the code depending on IDE usage or')
        print('command line usage. The user may specify as many or as few command line options as desired, but if the command line')
        print('option is used and not all options are provided by the user, the hard coded inputs will be enforced, and the user')
        print('will be warned.*')
        
        print('\nTo use command line override options there must be alteration between tags started with the - character and the')
        print('input value. Ordering of tag/input pairs do not matter. Example:')
        print('    python3 cluster_analysis.py -tag1 tag1-input  -tag2 tag2-input  -tag3 tag3-input  ...N-arbitrary tagsN/tagN-input pairs')
        
        # print options
        print('\n\nOptions:')
        
        # print -topo option
        print(f'\n -topo or -t <topo-filename>  cluster_analysis variable: topofile    hard coded: {topofile}')
        print('    Command line option to set topofile name of topo/coord file to be read into cluster_analysis. Currently')
        print('    supported topofile formats are:')
        print('         .data = .data or dat file from LAMMPS')
        print('    Example usage:')
        print('        python3 cluster_analysis.py -topo EXAMPLE_TOPO-FILE.data')
        
        # print -n0 option
        print(f'\n -n0 <int>  cluster_analysis variable: N0    hard coded: {N0}')
        print('    Command line option to set the number of intial molecules before polymerization for computing p. Example usage:')
        print('        python3 cluster_analysis.py -n0 100')
        
        # print -fav option
        print(f'\n -fav <float>  cluster_analysis variable: fav    hard coded: {fav}')
        print('    Command line option to set the average number of functional groups present per monomer unit for computing p, pg, Xn. Example usage:')
        print('        python3 cluster_analysis.py -fav 1.5')
        
        # print -txt option
        print(f'\n -txt <T|F>  cluster_analysis variable: txtfile    hard coded: {txtfile}')
        print('    Command line option to write topofile BASENAME.txt file of information found (T or F). Example usage:')
        print('        python3 cluster_analysis.py -txt T')
        
        # print -opt or -man option
        print('\n -opt or -man')
        print('    Will tell cluster_analysis to only print out avaiable command line options known as tagN and tagN-inputs. This is the')
        print('    only tag that doesnt require a tagN-input following the tag since the code will only look for if -opt is in')
        print('    command line input list. Example usage:')
        print('        python3 cluster_analysis.py -man')
            
        # print -gui option
        print('\n -gui <no addition info required> cluster_analysis variable: use_GUI')
        print('    Command line option flag to load GUI. Example usage:')
        print('        python3 cluster_analysis.py -gui')
        
        # print requirements and dependencies
        print('\n\nRequirements and Dependencies:')
        print('    Requirements: python 3.7+ to guarantee ordering of dictionaries which stores most of the data')
        
        # print buffer and exit code
        print('\n\n')
        sys.exit()
        
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters
    if '-opt' not in commandline_inputs and commandline_inputs:   
        # Check that the given command line inputs are even for alternating tags/tag-inputs
        if (len(commandline_inputs) % 2) != 0:
            print('\nERROR command line option over ride used but odd number of command line arguments provided (missing tag or tag-input value).\n')
            sys.exit()
        
        # Check that tags are alternating between tag and tag-input
        nonalternating_tags = [commandline_inputs[i] for i in range(0, len(commandline_inputs), 2) if commandline_inputs[i][0] != '-']
        if nonalternating_tags:
            print(f'\nERROR tags are not alernating between tag and tag-input. Incorrect tag(s) list: {str(nonalternating_tags)} in {str(commandline_inputs)}\n')
            sys.exit()
        
        # Check that tag is supported and log if tag from the command line set supported tags
        supported_tags = ['-topo', '-n0', '-fav', '-txt']
        shortcut_tags = {'-t':'-topo', '-n0':'-n0', '-fav':'-fav', '-txt':'-txt'}
        default_variables ={'-topo': topofile, '-n0':N0, '-fav': fav, '-txt': txtfile}
        tags = {i:'' for i in supported_tags}
        
        # Update tags and raise exception if user provided unsupported tag
        for i in range(0, len(commandline_inputs), 2):
            user_tag = commandline_inputs[i]
            tag_input = commandline_inputs[i+1]
            if user_tag in shortcut_tags: 
                user_tag = shortcut_tags[user_tag]
            if user_tag not in supported_tags:
                print(f'\nERROR requesting unsupported command line tag   {str(user_tag)}\n')
                sys.exit()
            tags[user_tag] = tag_input
            
        # Loop through tag_checks and warn user hard coded variables will be enforced
        print('\n\nCommand line run option override checks (will warn if command line run option is used but not all options are provided at the command line):')
        for i in tags:
            if not tags[i]:
                print('WARNING override option   {:<18} not provided at the command line. Hard coded input being enforced: {}'.format(i, default_variables[i]))
        
        ##########################################
        # Start applying command line over rides #
        ##########################################
        # Start changing python variables based on command line over rides
        print('\n\nCommand line run option override inputs confirmation (will confirm to user that the command line override does override the hard coded inputs):')
        def T_F_string2boolean(tag, string):
            if string == 'F': return False
            elif string == 'T': return True
            else: print(f'\nERROR Provided tag: {tag} string is {string} which is not T or F. Use T or F string\n'); sys.exit()
            
        # set new -topo option and print confirmation
        if tags['-topo']:
            topofile = tags['-topo']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-topo', topofile))
            
        # set new -n0 option and print confirmation
        if tags['-n0']:
            try: N0 = int(tags['-n0']) # try getting int
            except: 
                print('ERROR -n0 input could not be converted to an int value')
                sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-n0', N0))    
            
        # set new -fav option and print confirmation
        if tags['-fav']:
            try: fav = float(tags['-fav']) # try getting float
            except: 
                print('ERROR -fav input could not be converted to an float value')
                sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-fav', fav))  
            
        # set new -txt option and print confirmation
        if tags['-txt']:
            txtfile = T_F_string2boolean('-txt', (tags['-txt']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-txt', txtfile))
            
        # print buffer
        print('\n\n')
    
    return topofile, N0, txtfile, fav

###############################
# For auto file clusters code #
###############################
def clusters_auto(files_directory, N0, txtfile, fav, newfile, commandline_inputs):
    # if -opt or -man option is in commandline_inputs print options and stop code execution
    if '-opt' in commandline_inputs or  '-man' in commandline_inputs:
        print('\n\nauto_cluster_analysis has been run with -opt or -man option to show the optional command line overrides available. Command line')
        print('option summary [-tag <tag-input>]:')
        print('python3 auto_cluster_analysis.py [-dir <dir w/ .data files>] [-txt <T|F>] [-n0 <int>] [-fav <float>] [-new-file <new filename>] <-gui> <-opt>|<-man>')


        print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the auto_cluster_analysis.py')
        print('file will be enforced. The command line inputs option gives more flexibility to the code depending on IDE usage or')
        print('command line usage. The user may specify as many or as few command line options as desired, but if the command line')
        print('option is used and not all options are provided by the user, the hard coded inputs will be enforced, and the user')
        print('will be warned.*')
        
        print('\nTo use command line override options there must be alteration between tags started with the - character and the')
        print('input value. Ordering of tag/input pairs do not matter. Example:')
        print('    python3 auto_cluster_analysis.py -tag1 tag1-input  -tag2 tag2-input  -tag3 tag3-input  ...N-arbitrary tagsN/tagN-input pairs')
        
        # print options
        print('\n\nOptions:')
        
        # print -dir option
        print(f'\n -dir or -d <new directory name>   auto_cluster_analysis variable: files_directory    hard coded: {files_directory}')
        print('    Command line option to to set directory where LAMMPS datafiles are stored. The code will change to this directory')
        print('    to read ALL files in this directory with the *.data extension. This means ONLY have LAMMPS *.data files in this')
        print('    directory you want to read into this code to automate the molecule logging process. Example usage:')
        print('        python3 auto_cluster_analysis.py -dir EXAMPLE_DIRECTORY')
        
        # print -new-file option
        print(f'\n -new-file or -nf <new filename>   auto_cluster_analysis variable: newfile    hard coded: {newfile}')
        print('    Command line option to set new filename. Example usage:')
        print('        python3 auto_cluster_analysis.py -new-file csv_log_file')
        
        # print -n0 option
        print(f'\n -n0 <int>  auto_cluster_analysis variable: N0    hard coded: {N0}')
        print('    Command line option to set the number of intial molecules before polymerization for computing p. Example usage:')
        print('        python3 auto_cluster_analysis.py -n0 100')
        
        # print -fav option
        print(f'\n -fav <float>  auto_cluster_analysis variable: fav    hard coded: {fav}')
        print('    Command line option to set the average number of functional groups present per monomer unit for computing p, pg, Xn. Example usage:')
        print('        python3 auto_cluster_analysis.py -fav 1.5')
        
        # print -txt option
        print(f'\n -txt <T|F>  auto_cluster_analysis variable: txtfile    hard coded: {txtfile}')
        print('    Command line option to write topofile BASENAME.txt file of information found (T or F). Example usage:')
        print('        python3 auto_cluster_analysis.py -txt T')
        
        # print -opt or -man option
        print('\n -opt or -man')
        print('    Will tell auto_cluster_analysis to only print out avaiable command line options known as tagN and tagN-inputs. This is the')
        print('    only tag that doesnt require a tagN-input following the tag since the code will only look for if -opt is in')
        print('    command line input list. Example usage:')
        print('        python3 auto_cluster_analysis.py -man')
            
        # print -gui option
        print('\n -gui <no addition info required> auto_cluster_analysis variable: use_GUI')
        print('    Command line option flag to load GUI. Example usage:')
        print('        python3 auto_cluster_analysis.py -gui')
        
        # print requirements and dependencies
        print('\n\nRequirements and Dependencies:')
        print('    Requirements: python 3.7+ to guarantee ordering of dictionaries which stores most of the data')
        
        # print buffer and exit code
        print('\n\n')
        sys.exit()
        
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters
    if '-opt' not in commandline_inputs and commandline_inputs:   
        # Check that the given command line inputs are even for alternating tags/tag-inputs
        if (len(commandline_inputs) % 2) != 0:
            print('\nERROR command line option over ride used but odd number of command line arguments provided (missing tag or tag-input value).\n')
            sys.exit()
        
        # Check that tags are alternating between tag and tag-input
        nonalternating_tags = [commandline_inputs[i] for i in range(0, len(commandline_inputs), 2) if commandline_inputs[i][0] != '-']
        if nonalternating_tags:
            print(f'\nERROR tags are not alernating between tag and tag-input. Incorrect tag(s) list: {str(nonalternating_tags)} in {str(commandline_inputs)}\n')
            sys.exit()
        
        # Check that tag is supported and log if tag from the command line set supported tags
        supported_tags = ['-dir', '-n0', '-fav', '-txt', '-new-file']
        shortcut_tags = {'-d':'-dir', '-n0':'-n0', '-fav':'-fav', '-txt':'-txt', '-nf':'-new-file'}
        default_variables ={'-dir': files_directory, '-n0':N0, '-fav': fav, '-txt': txtfile, '-new-file': newfile}
        tags = {i:'' for i in supported_tags}
        
        # Update tags and raise exception if user provided unsupported tag
        for i in range(0, len(commandline_inputs), 2):
            user_tag = commandline_inputs[i]
            tag_input = commandline_inputs[i+1]
            if user_tag in shortcut_tags: 
                user_tag = shortcut_tags[user_tag]
            if user_tag not in supported_tags:
                print(f'\nERROR requesting unsupported command line tag   {str(user_tag)}\n')
                sys.exit()
            tags[user_tag] = tag_input
            
        # Loop through tag_checks and warn user hard coded variables will be enforced
        print('\n\nCommand line run option override checks (will warn if command line run option is used but not all options are provided at the command line):')
        for i in tags:
            if not tags[i]:
                print('WARNING override option   {:<18} not provided at the command line. Hard coded input being enforced: {}'.format(i, default_variables[i]))
        
        ##########################################
        # Start applying command line over rides #
        ##########################################
        # Start changing python variables based on command line over rides
        print('\n\nCommand line run option override inputs confirmation (will confirm to user that the command line override does override the hard coded inputs):')
        def T_F_string2boolean(tag, string):
            if string == 'F': return False
            elif string == 'T': return True
            else: print(f'\nERROR Provided tag: {tag} string is {string} which is not T or F. Use T or F string\n'); sys.exit()
            
        # set new -dir option and print confirmation
        if tags['-dir']:
            if tags['-dir'].isspace() or tags['-dir'] == '':
                print('-dir flag was used, with an empty directory string. The directory string must contain alpha-numeric characters')
            else:
                files_directory = tags['-dir']
                print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-dir', files_directory))
            
        # set new -n0 option and print confirmation
        if tags['-n0']:
            try: N0 = int(tags['-n0']) # try getting int
            except: 
                print('ERROR -n0 input could not be converted to an int value')
                sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-n0', N0))    
            
        # set new -fav option and print confirmation
        if tags['-fav']:
            try: fav = float(tags['-fav']) # try getting float
            except: 
                print('ERROR -fav input could not be converted to an float value')
                sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-fav', fav))  
            
        # set new -txt option and print confirmation
        if tags['-txt']:
            txtfile = T_F_string2boolean('-txt', (tags['-txt']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-txt', txtfile))
            
        # set -new-file option and print confirmation
        if tags['-new-file']:
            newfile = tags['-new-file']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-new-file', newfile))
            
        # print buffer
        print('\n\n')
    
    return files_directory, N0, txtfile, fav, newfile

##########################
# For lmp2SYBYLmol2 code #
##########################
def lmp2SYBYLmol2_interface(topofile, parent_directory, remove_PBC_bonds, commandline_inputs):
    # if -opt or -man option is in commandline_inputs print options and stop code execution
    if '-opt' in commandline_inputs or  '-man' in commandline_inputs:
        print('\n\nlmp2SYBYLmol2 has been run with -opt or -man option to show the optional command line overrides available. Command line')
        print('option summary [-tag <tag-input>]:')
        print('python3 lmp2SYBYLmol2.py [-topo <topo-filename>] [-dir <new directory name>] [-rm-pbc-bonds <T|F>] <-gui> <-opt>|<-man>')


        print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the lmp2SYBYLmol2.py')
        print('file will be enforced. The command line inputs option gives more flexibility to the code depending on IDE usage or')
        print('command line usage. The user may specify as many or as few command line options as desired, but if the command line')
        print('option is used and not all options are provided by the user, the hard coded inputs will be enforced, and the user')
        print('will be warned.*')
        
        print('\nTo use command line override options there must be alteration between tags started with the - character and the')
        print('input value. Ordering of tag/input pairs do not matter. Example:')
        print('    python3 lmp2SYBYLmol2.py -tag1 tag1-input  -tag2 tag2-input  -tag3 tag3-input  ...N-arbitrary tagsN/tagN-input pairs')
        
        # print options
        print('\n\nOptions:')
        
        # print -topo option
        print(f'\n -topo or -t <topo-filename>  lmp2SYBYLmol2 variable: topofile    hard coded: {topofile}')
        print('    Command line option to set topofile name of topo/coord file to be read into cluster_analysis. Currently')
        print('    supported topofile formats are:')
        print('         .data = .data or dat file from LAMMPS')
        print('    Example usage:')
        print('        python3 lmp2SYBYLmol2.py -topo EXAMPLE_TOPO-FILE.data')
        
        # print -dir option
        print(f'\n -dir or -d <new directory name>   lmp2SYBYLmol2 variable: parent_directory    hard coded: {parent_directory}')
        print('    Command line option to set new directory name to store all files bond_react_template_merge can write. The file')
        print('    management is based on the present working directory (pwd) and will store files one directory deep inside the')
        print('    pwd. The -dir option will specify the directory that is inside of pwd. This will allow for easier file ')
        print('    management and allows the user to merge multiple files for a certain system or research project and keep them')
        print('    all organized into a single directory.')
        print('    python variables*. Example usage:')
        print('        python3 lmp2SYBYLmol2.py -dir EXAMPLE_TEST_SYSTEM_DIRECTORY')
        
        # print -rm-pbc-bonds option
        print(f'\n -rm-pbc-bonds or -rm <T|F>   lmp2SYBYLmol2 variable: remove_PBC_bonds    hard coded: {remove_PBC_bonds}')
        print('    Command line option to remove periodic bonds from the written .mol2 file. Typically you would only use this option to')
        print('    visualize the converted files in ChemDraw, VMD, Avogadro, Avogadro2, ... since periodic bonds will be visualized as')
        print('    spanning the imaginary simulation cell. If you are using lmp2SYBYLmol2 for adding atoms for more MD simulations DO NOT')
        print('    USE this option (T or F). Example usage:')
        print('        python3 lmp2SYBYLmol2.py -rm-pbc-bonds T')
        
        # print -opt or -man option
        print('\n -opt or -man')
        print('    Will tell lmp2SYBYLmol2 to only print out avaiable command line options known as tagN and tagN-inputs. This is the')
        print('    only tag that doesnt require a tagN-input following the tag since the code will only look for if -opt is in')
        print('    command line input list. Example usage:')
        print('        python3 lmp2SYBYLmol2.py -man')
            
        # print -gui option
        print('\n -gui <no addition info required> lmp2SYBYLmol2 variable: use_GUI')
        print('    Command line option flag to load GUI. Example usage:')
        print('        python3 lmp2SYBYLmol2.py -gui')
        
        # print requirements and dependencies
        print('\n\nRequirements and Dependencies:')
        print('    Requirements: python 3.7+ to guarantee ordering of dictionaries which stores most of the data')
        
        # print buffer and exit code
        print('\n\n')
        sys.exit()
        
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters
    if '-opt' not in commandline_inputs and commandline_inputs:   
        # Check that the given command line inputs are even for alternating tags/tag-inputs
        if (len(commandline_inputs) % 2) != 0:
            print('\nERROR command line option over ride used but odd number of command line arguments provided (missing tag or tag-input value).\n')
            sys.exit()
        
        # Check that tags are alternating between tag and tag-input
        nonalternating_tags = [commandline_inputs[i] for i in range(0, len(commandline_inputs), 2) if commandline_inputs[i][0] != '-']
        if nonalternating_tags:
            print(f'\nERROR tags are not alernating between tag and tag-input. Incorrect tag(s) list: {str(nonalternating_tags)} in {str(commandline_inputs)}\n')
            sys.exit()
        
        # Check that tag is supported and log if tag from the command line set supported tags
        supported_tags = ['-topo', '-dir', '-rm-pbc-bonds']
        shortcut_tags = {'-t':'-topo', '-d':'-dir', '-rm':'-rm-pbc-bonds'}
        default_variables ={'-topo': topofile, '-dir':parent_directory, '-rm-pbc-bonds':remove_PBC_bonds}
        tags = {i:'' for i in supported_tags}
        
        # Update tags and raise exception if user provided unsupported tag
        for i in range(0, len(commandline_inputs), 2):
            user_tag = commandline_inputs[i]
            tag_input = commandline_inputs[i+1]
            if user_tag in shortcut_tags: 
                user_tag = shortcut_tags[user_tag]
            if user_tag not in supported_tags:
                print(f'\nERROR requesting unsupported command line tag   {str(user_tag)}\n')
                sys.exit()
            tags[user_tag] = tag_input
            
        # Loop through tag_checks and warn user hard coded variables will be enforced
        print('\n\nCommand line run option override checks (will warn if command line run option is used but not all options are provided at the command line):')
        for i in tags:
            if not tags[i]:
                print('WARNING override option   {:<18} not provided at the command line. Hard coded input being enforced: {}'.format(i, default_variables[i]))
        
        ##########################################
        # Start applying command line over rides #
        ##########################################
        # Start changing python variables based on command line over rides
        print('\n\nCommand line run option override inputs confirmation (will confirm to user that the command line override does override the hard coded inputs):')
        def T_F_string2boolean(tag, string):
            if string == 'F': return False
            elif string == 'T': return True
            else: print(f'\nERROR Provided tag: {tag} string is {string} which is not T or F. Use T or F string\n'); sys.exit()
            
        # set new -topo option and print confirmation
        if tags['-topo']:
            topofile = tags['-topo']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-topo', topofile))
            
        # set new -dir option and print confirmation
        if tags['-dir']:
            if tags['-dir'].isspace() or tags['-dir'] == '':
                print('-dir flag was used, with an empty directory string. The directory string must contain alpha-numeric characters')
            else:
                parent_directory = tags['-dir']
                print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-dir', parent_directory))
                
        # set new -rm-pbc-bonds option and print confirmation
        if tags['-rm-pbc-bonds']:
            remove_PBC_bonds = T_F_string2boolean('-rm-pbc-bonds', (tags['-rm-pbc-bonds']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-rm-pbc-bonds', remove_PBC_bonds))
            
            
        # print buffer
        print('\n\n')
    
    return topofile, parent_directory, remove_PBC_bonds