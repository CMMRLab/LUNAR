#!/usr/bin/env python
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


def print_man_page(topofile, parent_directory, newfile, atom_style, atoms2remove, include_type_labels, method):

    
    # print general command line options
    print('\n\natom_removal has been run with -opt or -man option to show the optional command line overrides available. Command line')
    print('option summary [-tag <tag-input>]:')
    print('python3 atom_removal.py [-topo <topo-filename>] [-dir <new directory name>] [-newfile <string>] [-atomstyle <atomstyle>]')
    print('                        [-type-labels<T|F>] [-types <string of IDs>] [-method <atomIDs or TypeIDs or ...>] <-gui> <-opt>|<-man>')


    print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the atom_removal.py')
    print('file will be enforced. The command line inputs option gives more flexibility to the code depending on IDE usage or')
    print('command line usage. The user may specify as many or as few command line options as desired, but if the command line')
    print('option is used and not all options are provided by the user, the hard coded inputs will be enforced, and the user')
    print('will be warned.*')
    
    print('\nTo use command line override options there must be alteration between tags started with the - character and the')
    print('input value. Ordering of tag/input pairs do not matter. Example:')
    print('    python3 atom_removal.py -tag1 tag1-input  -tag2 tag2-input  -tag3 tag3-input  ...N-arbitrary tagsN/tagN-input pairs')
    
    # print options
    print('\n\nOptions:')
    
    # print -topo option
    print(f'\n -topo or -t <topo-filename>  atom_removal variable: topofile    hard coded: {topofile}')
    print('    Command line option to set topofile name of topo/coord file to be read into convert2graphite. Currently')
    print('    supported topofile formats are:')
    print('         .data = .data or dat file from LAMMPS')
    print('    Example usage:')
    print('        python3 atom_removal.py -topo EXAMPLE_TOPO-FILE.data')

    # print -dir option
    print(f'\n -dir or -d <new directory name>   atom_removal variable: parent_directory    hard coded: {parent_directory}')
    print('    Command line option to set new directory name to store all files bond_react_template_merge can write. The file')
    print('    management is based on the present working directory (pwd) and will store files one directory deep inside the')
    print('    pwd. The -dir option will specify the directory that is inside of pwd. This will allow for easier file ')
    print('    management and allows the user to merge multiple files for a certain system or research project and keep them')
    print('    all organized into a single directory.')
    print('    python variables*. Example usage:')
    print('        python3 atom_removal.py -dir EXAMPLE_TEST_SYSTEM_DIRECTORY')
    
    # print -newfile option
    print(f'\n -newfile or -nf <string>   atom_removal variable: newfile    hard coded: {newfile}')
    print('    Command line option to set the new output filename(s). The following options exist for using the newfile string')
    print('    for setting the output file basenames: ')
    print()
    print("      if newfile starts with ':' or ends with ':'")
    print("        The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix added to")
    print("        the file basename. The following are examples:")
    print("          Suffix (newfile = ':_rm_atoms'  and topofile = 'detda.data')")
    print("            basename = 'detda_rm_atoms', where the ':' character acts as a placeholder for the topofile basename.")
    print("          Prefix (newfile = 'rm_atoms-:'  and topofile = 'detda.data')")
    print("            basename = 'rm_atoms-detda', where the ':' character acts as a placeholder for the topofile basename.")
    print('        Recommended usage: common and safe as this method is safe and output filename(s) carry similar names to')
    print('        input filename(s).')
    print()
    print("      if newfile == 'ANYTEXT'")
    print("        The output filename(s) will be set as the file tag of each file and 'ANYTEXT' will be a suffix appended")
    print("        to the file tag. For example:")
    print("          newfile = 'detda_renamed' and topofile = 'detda.data'")
    print("            basename = 'detda_renamed'")
    print('        Recommended usage: occasional and safe as output filename(s) may no longer carry similar names as output')
    print('        filename(s), but is safe as output filename(s) will not overwrite input filename(s)')
    print()
    print("      if newfile == ''")
    print("        The output filename(s) will be the same as the input filename(s). This can be a dangerous option as you")
    print("        may inadvertently overwrite a file and then must assume the file to contain certain information, but it")
    print("        contains other information.")
    print("        Recommended usage: rare and dangerous as this could lead to using incorrect files if not paying attention")
    print("        very carefully.")
    print('    Example usage:')
    print('        python3 atom_removal.py -newfile  :_rm_atoms')
    print('        python3 atom_removal.py -newfile ": rm_atoms"  # note the use of the " " bounding characters to pass a whitespace')
    
    # print -atomstyle option
    print(f'\n -atomstyle or -as <atomstyle>   atom_removal variable: atom_style    hard coded: {atom_style}')
    print('    Command line option to set atom style in written LAMMPS .data file. Currently supported atom styles are full,')
    print('    charge, molecular. *NOTE: If another atomstyle is desired it can be added in the write_lmp.py script.* Example')
    print('    usage:')
    print('        python3 atom_removal.py -atomstyle full')
    
    # print -atoms2remove option
    print(f'\n -atoms2remove or -a2r <string of IDs> atom_removal variable: atoms2remove    hard coded: {atoms2remove}')
    print('    Command line option to set the list of atomIDs or atomTypeIDs or mass cutoff or size cutoff to remove from the')
    print('    system. NOTE for cluster-size or cluster-mass only a single value maybe supplied. Example usage:')
    print('        python3 atom_removal.py -atoms2remove 1,2')
    
    # print -method option
    print(f'\n -method or -m <string>   atom_removal variable: include_type_labels    hard coded: {method}')
    print('    Command line option to set method of identifying which atoms to remove from the system. The following methods are supported:')
    print("       'atomIDs' will identify atoms based on atomIDs")
    print("       'typeIDs' will identify atoms based on atom typeIDs")
    print("       'cluster-mass' will perform cluster analysis and identify atoms based on a cutoff value set in")
    print('                      atoms2remove, where all cluster mass less than or equal to the cutoff value will')
    print('                      be removed')
    print("       'cluster-size' will perform cluster analysis and identify atoms based on a cutoff value set in")
    print('                      atoms2remove, where all cluster sizes (number of atoms) less than or equal to the')
    print('                      cutoff value will be removed')
    print('    Example usage:')
    print('        python3 atom_removal.py -method atomIDs')
    print('        python3 atom_removal.py -method cluster-mass -atoms2remove 25.55')
    
    # print -type-labels option
    print(f'\n -type-labels or -tl <T|F>   atom_removal variable: include_type_labels    hard coded: {include_type_labels}')
    print('    Command line option to inlcude LAMMPS new type labels in written datafile to make bond/angle/dihedrals/improper')
    print('    types if applicable consistent between datafiles or molecule files (primary use for bond/react file generation)')
    print('    *NOTE: this option will also change how -write-bond-react molecule files are written to be consistent between')
    print('    the written .data file and the written molecule .lmpmol file.* Example usage:')
    print('        python3 atom_removal.py -type-labels T')
    
    # print -opt or -man option
    print('\n -opt or -man')
    print('    Will tell atom_removal to only print out avaiable command line options known as tagN and tagN-inputs. This is the')
    print('    only tag that doesnt require a tagN-input following the tag since the code will only look for if -opt is in')
    print('    command line input list. Example usage:')
    print('        python3 atom_removal.py -man')
        
    # print -gui option
    print('\n -gui <no addition info required> convert2graphite variable: use_GUI')
    print('    Command line option flag to load GUI. Example usage:')
    print('        python3 atom_removal.py -gui')
    
    # print requirements and dependencies
    print('\n\nRequirements and Dependencies:')
    print('    Requirements: python 3.7+ to guarantee ordering of dictionaries which stores most of the data')
    
    # print buffer and exit code
    print('\n\n')
    return

# Class to update inputs
class inputs:
    def __init__(self, topofile, parent_directory, newfile, atom_style, atoms2remove, include_type_labels, method, commandline_inputs):
        
        # Give access to inputs (update later on if command line over ride is given)
        self.commandline_inputs = commandline_inputs
        self.topofile = topofile
        self.parent_directory = parent_directory
        self.newfile = newfile
        self.atom_style = atom_style
        self.atoms2remove = atoms2remove
        self.include_type_labels = include_type_labels 
        self.method = method
        
        
        # Check that the given command line inputs are even for alternating tags/tag-inputs
        if (len(commandline_inputs) % 2) != 0:
            print('\nERROR command line option over ride used but odd number of command line arguments provided (missing tag or tag-input value).\n')
            sys.exit()
        
        # Check that tags are alternating between tag and tag-input
        nonalternating_tags = [commandline_inputs[i] for i in range(0, len(commandline_inputs), 2) if commandline_inputs[i][0] != '-']
        if nonalternating_tags:
            print(f'\nERROR tags are not alernating between tag and tag-input. Incorrect tag(s) list: {str(nonalternating_tags)} in {str(commandline_inputs)}\n')
            sys.exit()
        
        
        # Check that tag is supported and log if tag from the command line
        # set supported tags
        supported_tags = ['-topo', '-dir', '-newfile', '-atomstyle', '-type-labels', '-atoms2remove', '-method']
        
        # set shortcut_tags mapping
        shortcut_tags = {'-t':'-topo', '-d':'-dir', '-nf':'-newfile', '-as':'-atomstyle', '-tl':'-type-labels', '-a2r':'-atoms2remove', '-m':'-method'}
        
        # set default variables
        default_variables ={'-topo': self.topofile, '-dir':self.parent_directory, '-newfile': self.newfile, '-atomstyle': self.atom_style,
                            '-atoms2remove':self.atoms2remove, '-type-labels': self.include_type_labels, '-method':self.method}
        
        # set tag/tag-input pair as empty string and update
        tags = {i:'' for i in supported_tags}
        
        # Update tags and raise exception if user provided unsupported tag
        for i in range(0, len(commandline_inputs), 2):
            user_tag = commandline_inputs[i]
            tag_input = commandline_inputs[i+1]
            
            # Find if shortcut was used and update to full tag if in shortcuts
            if user_tag in shortcut_tags: 
                user_tag = shortcut_tags[user_tag]
            
            # 1st check make sure it is a tag
            if user_tag not in supported_tags:
                print(f'\nERROR requesting unsupported command line tag   {str(user_tag)}\n')
                sys.exit()
                
            # If 1st check clears add input to tags
            tags[user_tag] = tag_input
            
        # Loop through tag_checks and warn user hard coded variables will be enforced
        print('\n\nCommand line run option override checks (will warn if command line run option is used but not all options are provided at the command line):')
        for i in tags:
            if not tags[i]:
                print('WARNING override option   {:<18} not provided at the command line. Hard coded input being enforced: {}'.format(i, default_variables[i]))
        
                
        # Start changing python variables based on command line over rides
        print('\n\nCommand line run option override inputs confirmation (will confirm to user that the command line override does override the hard coded inputs):')
        def T_F_string2boolean(tag, string):
            if string == 'F': return False
            elif string == 'T': return True
            else: print(f'\nERROR Provided tag: {tag} string is {string} which is not T or F. Use T or F string\n'); sys.exit()
        
        
        ###############################################
        # set new -topo option and print confirmation #
        ###############################################
        if tags['-topo']:
            self.topofile = tags['-topo']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-topo', self.topofile))
            
        # set new -dir option and print confirmation
        if tags['-dir']:
            if tags['-dir'].isspace() or tags['-dir'] == '':
                print('-dir flag was used, with an empty directory string. The directory string must contain alpha-numeric characters')
            else:
                self.parent_directory = tags['-dir']
                print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-dir', self.parent_directory))
            
        # set new -newfile option and print confirmation
        if tags['-newfile']:
            self.newfile = tags['-newfile']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-newfile', self.newfile))
        
        # set new -atomstyle option and print confirmation
        if tags['-atomstyle']:
            self.atom_style = tags['-atomstyle']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-atomstyle', self.atom_style))
            
        # set new -atoms2remove option and print confirmation
        if tags['-atoms2remove']:
            try: self.atoms2remove = [eval(i) for i in tags['-atoms2remove'].split(',')]
            except:
                print(f'ERROR tag -atoms2remove atomID list could not be converted to ints or floats {tags["-atoms2remove"]}'); sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-atoms2remove', self.atoms2remove))
        
        # set new -type-labels option and print confirmation
        if tags['-type-labels']:
            self.include_type_labels = T_F_string2boolean('-type-labels', (tags['-type-labels']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-type-labels', self.include_type_labels))
        
        # set new -method option and print confirmation
        if tags['-method']:
            self.method = tags['-method']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-method', self.method))
        
        # print buffer
        print('\n\n')