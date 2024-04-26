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

python3 bond_react_merge_prep.py -topo smp20230330_IFF_rlxd0_testv3_8proc_1ts.data -cta smp20230330_IFF_rlxd0_testv3_8proc_1ts.cta -dir test -ext test_cta -as charge
"""

##############################
# Import Necessary Libraries #
##############################
import sys


def print_man_page(topofile, cta_file, newfile, atom_style, parent_dir, rm_unused_coeffs):

    
    # print general command line options
    print('\n\nbond_react_merge_prep has been run with -opt or -man option to show the optional command line overrides available. Command line')
    print('option summary [-tag <tag-input>]:')
    print('python3 bond_react_merge_prep.py [-topo <topo-filename>] [-cta <cta-filename>] [-dir <new directory name>] [-newfile <string>]')
    print('                                 [-atomstyle <atomstyle>] <-gui> <-opt>|<-man>')


    print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the bond_react_merge_prep.py')
    print('file will be enforced. The command line inputs option gives more flexibility to the code depending on IDE usage or')
    print('command line usage. The user may specify as many or as few command line options as desired, but if the command line')
    print('option is used and not all options are provided by the user, the hard coded inputs will be enforced, and the user')
    print('will be warned.*')
    
    print('\nTo use command line override options there must be alteration between tags started with the - character and the')
    print('input value. Ordering of tag/input pairs do not matter. Example:')
    print('    python3 bond_react_merge_prep.py -tag1 tag1-input  -tag2 tag2-input  -tag3 tag3-input  ...N-arbitrary tagsN/tagN-input pairs')
    
    # print options
    print('\n\nOptions:')
    
    # print -topo option
    print(f'\n -topo or -t <topo-filename>  bond_react_merge_prep variable: topofile    hard coded: {topofile}')
    print('    Command line option to set topofile name of topo/coord file to be read into bond_react_merge_prep. Currently')
    print('    supported topofile formats are:')
    print('         .data = .data or dat file from LAMMPS')
    print('    Example usage:')
    print('        python3 bond_react_merge_prep.py -topo EXAMPLE_TOPO-FILE.data')
    
    # print -cta option
    print(f'\n -cta or -c <cta-filename>   bond_react_merge_prep variable: cta_file    hard coded: {cta_file}')
    print('    Command line option to set nta filename of file containing new type assignment (nta) file to be read into')
    print('    all2lmp for conversion. This file contains atom-ids and atom-types in columns, the atom-type in each row will')
    print('    be mapped onto the corresponding atom-id. The mapping between atom-ids and atom-types allows for the changing')
    print('    of atom-types or setting of atom-types. *NOTE: This code assumes correct atom-typing since the atom-types will')
    print('    set the atom, bond, angle, dihedral, improper coeff types and coeff parameters. See Materials Studio software')
    print('    below for msi2lmp conversion. Example usage:')
    print('        python3 bond_react_merge_prep.py -cta EXAMPLE_CTA-FILE.cta')

    # print -dir option
    print(f'\n -dir or -d <new directory name>   bond_react_merge_prep variable: parent_directory    hard coded: {parent_dir}')
    print('    Command line option to set new directory name to store all files bond_react_template_merge can write. The file')
    print('    management is based on the present working directory (pwd) and will store files one directory deep inside the')
    print('    pwd. The -dir option will specify the directory that is inside of pwd. This will allow for easier file ')
    print('    management and allows the user to merge multiple files for a certain system or research project and keep them')
    print('    all organized into a single directory.')
    print('    python variables*. Example usage:')
    print('        python3 bond_react_merge_prep.py -dir EXAMPLE_TEST_SYSTEM_DIRECTORY')
    
    # print -newfile option
    print(f'\n -newfile or -nf <string>   bond_react_merge_prep variable: newfile    hard coded: {newfile}')
    print('    Command line option to set the new output filename(s). The following options exist for using the newfile string')
    print('    for setting the output file basenames: ')
    print()
    print("      if newfile starts with ':' or ends with ':'")
    print("        The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix added to")
    print("        the file basename. The following are examples:")
    print("          Suffix (newfile = ':_cta'  and topofile = 'detda.data')")
    print("            basename = 'detda_cta', where the ':' character acts as a placeholder for the topofile basename.")
    print("          Prefix (newfile = 'cta-:'  and topofile = 'detda.data')")
    print("            basename = 'cta-detda', where the ':' character acts as a placeholder for the topofile basename.")
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
    print('        python3 bond_react_merge_prep.py -newfile  :_cta')
    print('        python3 bond_react_merge_prep.py -newfile ": cta"  # note the use of the " " bounding characters to pass a whitespace')
    
    # print -atomstyle option
    print(f'\n -atomstyle or -as <atomstyle>   bond_react_merge_prep variable: atom_style    hard coded: {atom_style}')
    print('    Command line option to set atom style in written LAMMPS .data file. Currently supported atom styles are full,')
    print('    charge, molecular. *NOTE: If another atomstyle is desired it can be added in the write_lmp.py script.* Example')
    print('    usage:')
    print('        python3 bond_react_merge_prep.py -atomstyle full')
    
    # print -rm-na-coeffs option
    print(f'\n -rm-na-coeffs or -rm <T|F> bond_react_merge_prep variable: rm_unused_coeffs    hard coded: {rm_unused_coeffs}')
    print('    Command line option to remove any coeffIDs where no atomIDs do not currently use to coeffID. There are')
    print('    cases were LAMMPS datafiles may be built where not all coeffIDs are used by atomIDs in the current')
    print('    LAMMPS datafile. bond_react_merge_prep.py will be default assign "N/A" for atom CoeffIDs, "N/A  N/A" for')
    print('    bond coeff IDs, etc ... for any coeffIDs in which no atomIDs currently used the CoeffIDs. bond_react_merge.py')
    print('    is compatible with the "N/A" comments where if they exist in the file read by bond_react_merge.py will just')
    print('    not assign any atomID lists to these coeffs, however for cleanliness purposes it may be desired to remove')
    print('    these coeffIDs. This option will remove the unused CoeffIDs. T will use option and F will not. Example usage:')
    print('        python3 bond_react_merge_prep.py -rm-na-coeffs T')
    
    # print -opt or -man option
    print('\n -opt or -man')
    print('    Will tell bond_react_merge_prep to only print out avaiable command line options known as tagN and tagN-inputs.')
    print('    This is the only tag that doesnt require a tagN-input following the tag since the code will only look for if')
    print('     -opt is in command line input list. Example usage:')
    print('        python3 bond_react_merge_prep.py -man')
        
    # print -gui option
    print('\n -gui <no addition info required> bond_react_merge_prep variable: use_GUI')
    print('    Command line option flag to load GUI. Example usage:')
    print('        python3 bond_react_merge_prep.py -gui')
    
    # print requirements and dependencies
    print('\n\nRequirements and Dependencies:')
    print('    Requirements: python 3.7+ to guarantee ordering of dictionaries which stores most of the data')
    print('    Dependencies: CORRECT ATOM-TYPING in the .cta file and correct format of the .cta file')
    
    # print buffer and exit code
    print('\n\n')
    return

# Class to update inputs
class inputs:
    def __init__(self, commandline_inputs, topofile, cta_file, newfile, atom_style, parent_dir, rm_unused_coeffs):
        
        # Give access to inputs (update later on if command line over ride is given)
        self.commandline_inputs = commandline_inputs
        self.topofile = topofile
        self.cta_file = cta_file
        self.parent_dir = parent_dir
        self.newfile = newfile
        self.atom_style = atom_style
        self.rm_unused_coeffs = rm_unused_coeffs
        
        
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
        supported_tags = ['-topo', '-cta', '-dir', '-newfile', '-atomstyle', '-rm-na-coeffs']
        
        # set shortcut_tags mapping
        shortcut_tags = {'-t':'-topo', '-c':'-cta', '-d':'-dir', '-nf':'-newfile', '-as':'-atomstyle', '-rm':'-rm-na-coeffs'}
        
        # set default variables
        default_variables ={'-topo': self.topofile, '-cta': self.cta_file, '-dir':self.parent_dir, '-newfile': self.newfile,
                            '-atomstyle': self.atom_style, '-rm-na-coeffs':self.rm_unused_coeffs}
        
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
        
        # set new -cta option and print confirmation
        if tags['-cta']:
            self.cta_file = tags['-cta']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-cta', self.cta_file))
            
        # set new -dir option and print confirmation
        if tags['-dir']:
            if tags['-dir'].isspace() or tags['-dir'] == '':
                print('-dir flag was used, with an empty directory string. The directory string must contain alpha-numeric characters')
            else:
                self.parent_dir = tags['-dir']
                print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-dir', self.parent_dir))
            
        # set new -newfile option and print confirmation
        if tags['-newfile']:
            self.newfile = tags['-newfile']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-newfile', self.newfile))
        
        # set new -atomstyle option and print confirmation
        if tags['-atomstyle']:
            self.atom_style = tags['-atomstyle']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-atomstyle', self.atom_style))
            
        # set new -rm-na-coeffs option and print confirmation
        if tags['-rm-na-coeffs']:
            self.rm_unused_coeffs = T_F_string2boolean('-rm-na-coeffs', (tags['-rm-na-coeffs']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-rm-na-coeffs', self.rm_unused_coeffs))
        
        # print buffer
        print('\n\n')