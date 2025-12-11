#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
December 11, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

File contains command line options/man page and
command line over riding class to reset hard coded
inputs.

convert2cg1
add_pi_electrons

python3 add_pi_electrons.py -topo PFA_carb_2300K_2GPa_2000ps_typed_IFF_GT.data -dir pi -ext add_piE -as charge -rq T
python3 add_pi_electrons.py -topo PFA_carb_2300K_2GPa_2000ps_typed_IFF_GT.data -dir pi -ext add_piE -as charge -q0 T
"""

##############################
# Import Necessary Libraries #
##############################
import sys


def print_man_page(topofile, types2convert, atom_style, reset_charges, net_zero_charge, convert2cg1, add_pi_electrons,
                   parent_directory, newfile, include_type_labels, neighbor_charge_constraint, reset_simulation_cell):

    
    # print general command line options
    print('\n\nadd_pi_electrons has been run with -opt or -man option to show the optional command line overrides available. Command line')
    print('option summary [-tag <tag-input>]:')
    print('python3 add_pi_electrons.py [-topo <topo-filename>] [-cta <cta-filename>] [-dir <new directory name>] [-newfile <string>]')
    print('                            [-atomstyle <atomstyle>] [-types <string of IDs>] [-reset-charges <T|F>] [-charge0 <T|F>]')
    print('                            [-neigh_charge <none|check-neighbors|accumulate-carbon|accumulate-pi-electron>]')
    print('                            [-convert2cg1 <T|F>] [-pi-electrons <T|F>] [-type-labels<T|F>] <-gui> <-opt>|<-man>')


    print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the add_pi_electrons.py')
    print('file will be enforced. The command line inputs option gives more flexibility to the code depending on IDE usage or')
    print('command line usage. The user may specify as many or as few command line options as desired, but if the command line')
    print('option is used and not all options are provided by the user, the hard coded inputs will be enforced, and the user')
    print('will be warned.*')
    
    print('\nTo use command line override options there must be alteration between tags started with the - character and the')
    print('input value. Ordering of tag/input pairs do not matter. Example:')
    print('    python3 add_pi_electrons.py -tag1 tag1-input  -tag2 tag2-input  -tag3 tag3-input  ...N-arbitrary tagsN/tagN-input pairs')
    
    # print options
    print('\n\nOptions:')
    
    # print -topo option
    print(f'\n -topo or -t <topo-filename>  add_pi_electrons variable: topofile    hard coded: {topofile}')
    print('    Command line option to set topofile name of topo/coord file to be read into add_pi_electrons. Currently')
    print('    supported topofile formats are:')
    print('         .data = .data or dat file from LAMMPS')
    print('    Example usage:')
    print('        python3 add_pi_electrons.py -topo EXAMPLE_TOPO-FILE.data')

    # print -dir option
    print(f'\n -dir or -d <new directory name>   add_pi_electrons variable: parent_directory    hard coded: {parent_directory}')
    print('    Command line option to set new directory name to store all files bond_react_template_merge can write. The file')
    print('    management is based on the present working directory (pwd) and will store files one directory deep inside the')
    print('    pwd. The -dir option will specify the directory that is inside of pwd. This will allow for easier file ')
    print('    management and allows the user to merge multiple files for a certain system or research project and keep them')
    print('    all organized into a single directory.')
    print('    python variables*. Example usage:')
    print('        python3 add_pi_electrons.py -dir EXAMPLE_TEST_SYSTEM_DIRECTORY')
    
    # print -newfile option
    print(f'\n -newfile or -nf <string>   add_pi_electrons variable: newfile    hard coded: {newfile}')
    print('    Command line option to set the new output filename(s). The following options exist for using the newfile string')
    print('    for setting the output file basenames: ')
    print()
    print("      if newfile starts with ':' or ends with ':'")
    print("        The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix added to")
    print("        the file basename. The following are examples:")
    print("          Suffix (newfile = ':_pi_electrons'  and topofile = 'cnt.data')")
    print("            basename = 'cnt_pi_electrons', where the ':' character acts as a placeholder for the topofile basename.")
    print("          Prefix (newfile = 'pi_electrons-:'  and topofile = 'cnt.data')")
    print("            basename = 'pi_electrons-cnt', where the ':' character acts as a placeholder for the topofile basename.")
    print('        Recommended usage: common and safe as this method is safe and output filename(s) carry similar names to')
    print('        input filename(s).')
    print()
    print("      if newfile == 'ANYTEXT'")
    print("        The output filename(s) will be set as the file tag of each file and 'ANYTEXT' will be a suffix appended")
    print("        to the file tag. For example:")
    print("          newfile = 'cnt_renamed' and topofile = 'cnt.data'")
    print("            basename = 'cnt_renamed'")
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
    print('        python3 add_pi_electrons.py -newfile  :_pi_electrons')
    print('        python3 add_pi_electrons.py -newfile ": pi_electrons"  # note the use of the " " bounding characters to pass a whitespace')
    
    # print -atomstyle option
    print(f'\n -atomstyle or -as <atomstyle>   add_pi_electrons variable: atom_style    hard coded: {atom_style}')
    print('    Command line option to set atom style in written LAMMPS .data file. Currently supported atom styles are full,')
    print('    charge, molecular. *NOTE: If another atomstyle is desired it can be added in the write_lmp.py script.* Example')
    print('    usage:')
    print('        python3 add_pi_electrons.py -atomstyle full')
    
    # print -types option
    print(f'\n -types or -t2c <string of IDs> add_pi_electrons variable: types2convert    hard coded: {types2convert}')
    print('    Command line option to set the atom typeIDs to either convert to cg1 or to add pi electrons to. Note')
    print('    if the read-in topofile has comments of the atom types in the Masses section you may also specify a')
    print('    string type of the atom type, which will be converted to the corresponding atomTypeID interanlly.')
    print('    Example usage:')
    print('        python3 add_pi_electrons.py -types 1,2')
    print('        python3 add_pi_electrons.py -types 1,cp,c5')
    print('        python3 add_pi_electrons.py -types cp,c5')
    print('        python3 add_pi_electrons.py -types 1')
    
    # print -reset-charges option
    print(f'\n -reset-charges or -rq <T|F>   add_pi_electrons variable: reset_charges    hard coded: {reset_charges}')
    print('    Command line option to reset charges of all atom TypeIDs listed in types2convert to the cg1 charge. *NOTE convert2cg1')
    print('    will not automatically reset charge and add_pi_eletrons will automatically reset charge, so use accordingly.*')
    print('    T is for True and F is for False. Example usage:')
    print('        python3 add_pi_electrons.py -reset-charges T')
    
    # print -neigh_charge option
    print(f'\n -neigh_charge or -nq <none|check-neighbors|accumulate-carbon|accumulate-pi-electron>   add_pi_electrons variable: neighbor_charge_constraint    hard coded: {neighbor_charge_constraint}')
    print('    Command line option to set how to handle charge on "compounds" (materials that are not either pure graphite or CNT or')
    print('    fullerene, where there are other atom types bonded to the aromatic carbon atoms). The IFF charge model using the virtual')
    print('    pi-electrons was originally formulated only for pure graphitic systems and thus when trying to use the pi-electron charge')
    print('    with system that have different atom types other than aromatic carbons, the total system charge will become none-neutral.')
    print('    One method around this is to enforce the system charge to be neutral via the net_zero_charge option, however this method')
    print('    can have adverse effects such as completely removing charge from all none-aromatic atoms or arbitrarily scaling all none')
    print('    -aromatic atoms by the same charge value to achieve charge neutrality. This option provides a greater level of control to')
    print('    enforce charge neutrality. If your system is pure graphite or CNT or fullerene this option can be set')
    print('    to any supported options as none of the options effect this type of system. The following options exist:')
    print("      'none'                    which does not apply any constraint to how neighbor charges are handled. If there are any first")
    print('                                neighbors bonded to the aromatic carbon atoms, this method may result in a none-charge neutral system.')
    print()
    print("      'check-neighbors'         which checks that all neighbors are aromatic and only places the pi-electron if all neighbors are aromatic.")
    print()
    print("      'accumulate-carbon'       which will accumulate any residual charge into the carbon atom, to ensure the local grouping of atoms")
    print('                                stays charge neutral and thus the entire system will remain charge neutral.')
    print()
    print("      'accumulate-pi-electron' which will accumulate any residual charge into the pi-electron atoms, to ensure the local grouping of atoms")
    print('                               stays charge neutral and thus the entire system will remain charge neutral.')
    print()
    print("      'accumulate-neighbor'    which will accumulate any residual charge into the first neighboring atom, to ensure the local grouping of atoms")
    print('                               stays charge neutral and thus the entire system will remain charge neutral.')
    print('    Example usage:')
    print('        python3 add_pi_electrons.py -neigh_charge check-neighbors')
    
    # print -charge0 option
    print(f'\n -charge0 or -q0 <T|F>   add_pi_electrons variable: net_zero_charge    hard coded: {net_zero_charge}')
    print('    Command line option to make the system charge neutral after adding in pi-electrons via add_pi_electrons. If the system')
    print('    is a simple graphite sheet(s) or CNT then the system charge will already be charge neutral and this option can be left as')
    print('    False. However, if you have functional groups attached to graphite sheet(s) or CNTs and the system was charged via another')
    print('    method or the system is amorphous/glassy carbon from ReaxFF there will be a residual non-zero net charge of the system and')
    print('    this option shohuld be set as True. T is for True and F is for False. Example usage:')
    print('        python3 add_pi_electrons.py -charge0 T')
    
    # print -convert2cg1 option
    print(f'\n -convert2cg1 or -cg1 <T|F>   add_pi_electrons variable: convert2cg1    hard coded: {convert2cg1}')
    print('    Command line option to convert all bonds, angles, dihedrals, impropers, and crossterms that have TypeIDs')
    print('    in types2convert to "cg1" consistent parameters. T is for True and F is for False. Example usage:')
    print('        python3 add_pi_electrons.py -convert2cg1 T')
    
    # print -reset-box option
    print(f'\n -reset-box or -rb <T|F>   add_pi_electrons variable: reset_simulation_cell    hard coded: {reset_simulation_cell}')
    print('    Command line option to to reset the simulation cell size after adding pi-electrons or to not reset the simulation cell size.')
    print('    If the Boolean is True, the simulation cell size will be reset and if the Boolean is False, the simulation cell size will not')
    print('    be reset. This option is useful for when add_pi_electrons is True for some system types. When a pi-electron is added it inherits')
    print('    the image flag from the carbon atom that it is added to and LAMMPS will rewrap and atoms that are outside of the simulation cell')
    print('    when reading the file. However, some systems like a graphite system may have a simulation cell size that is too small for adequate')
    print('    re-wrapping of the pi-electron and it may cause errors. In these cases, it is beneficial to reset the simulation cell size. When')
    print('    the simulation cell size is reset the maximum distance between the any added pi-electron and the any orginal atom is used to')
    print('    increment the simulation cell to be larger in each direction.')

    print('    Example usage:')
    print('        python3 add_pi_electrons.py -reset-box T')
    
    # print -pi-electrons option
    print(f'\n -pi-electrons or -pie <T|F>   add_pi_electrons variable: convert2cg1    hard coded: {convert2cg1}')
    print('    Command line option to add pi electrons to atom TypeIDs in types2convert. Usage  will result in the')
    print('    addition of virtual pi electron atom, new bonds, and new angles. The atom types, bond types, and angle types will')
    print('    also be updated to be consistent with PCFF-IFF v1.5. Lastly, charges will also be updated and existing atoms and')
    print('    set to the new pi electrons. T is for True and F is for False. Example usage:')
    print('        python3 add_pi_electrons.py -pi-electrons T')
    
    # print -type-labels option
    print(f'\n -type-labels or -tl <T|F>   add_pi_electrons variable: include_type_labels    hard coded: {include_type_labels}')
    print('    Command line option to inlcude LAMMPS new type labels in written datafile to make bond/angle/dihedrals/improper')
    print('    types if applicable consistent between datafiles or molecule files (primary use for bond/react file generation)')
    print('    *NOTE: this option will also change how -write-bond-react molecule files are written to be consistent between')
    print('    the written .data file and the written molecule .lmpmol file.* Example usage:')
    print('        python3 add_pi_electrons.py -type-labels T')
    
    # print -opt or -man option
    print('\n -opt or -man')
    print('    Will tell add_pi_electrons to only print out avaiable command line options known as tagN and tagN-inputs. This is the')
    print('    only tag that doesnt require a tagN-input following the tag since the code will only look for if -opt is in')
    print('    command line input list. Example usage:')
    print('        python3 add_pi_electrons.py -man')
        
    # print -gui option
    print('\n -gui <no addition info required> add_pi_electrons variable: use_GUI')
    print('    Command line option flag to load GUI. Example usage:')
    print('        python3 add_pi_electrons.py -gui')
    
    # print requirements and dependencies
    print('\n\nRequirements and Dependencies:')
    print('    Requirements: python 3.7+ to guarantee ordering of dictionaries which stores most of the data')
    
    # print buffer and exit code
    print('\n\n')
    return

# Class to update inputs
class inputs:
    def __init__(self, commandline_inputs, topofile, types2convert, atom_style, reset_charges, net_zero_charge, convert2cg1, add_pi_electrons,
                 parent_directory, newfile, include_type_labels, neighbor_charge_constraint, reset_simulation_cell):
        
        # Give access to inputs (update later on if command line over ride is given)
        self.commandline_inputs = commandline_inputs
        self.topofile = topofile
        self.parent_directory = parent_directory
        self.newfile = newfile
        self.atom_style = atom_style
        self.types2convert = types2convert
        self.reset_charges = reset_charges
        self.net_zero_charge = net_zero_charge
        self.convert2cg1 = convert2cg1
        self.add_pi_electrons = add_pi_electrons
        self.include_type_labels = include_type_labels 
        self.neighbor_charge_constraint = neighbor_charge_constraint
        self.reset_simulation_cell = reset_simulation_cell
        
        
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
        supported_tags = ['-topo', '-dir', '-newfile', '-atomstyle', '-types', '-reset-charges', '-charge0', '-convert2cg1',
                          '-pi-electrons', '-type-labels', '-neigh_charge', '-reset-box']
        
        # set shortcut_tags mapping
        shortcut_tags = {'-t':'-topo', '-d':'-dir', '-nf':'-newfile', '-as':'-atomstyle', '-t2c':'-types', '-tl':'-type-labels', '-rb':'-reset-box',
                         '-rq':'-reset-charges', '-q0':'-charge0', '-cg1':'-convert2cg1', '-pie':'-pi-electrons', '-nq':'-neigh_charge'}
        
        # set default variables
        default_variables ={'-topo': self.topofile, '-dir':self.parent_directory, '-newfile': self.newfile, '-atomstyle': self.atom_style,
                            '-types':self.types2convert, '-reset-charges': self.reset_charges, '-charge0':self.net_zero_charge,
                            '-convert2cg1':self.convert2cg1, '-pi-electrons':self.add_pi_electrons, '-type-labels': self.include_type_labels,
                            '-neigh_charge':self.neighbor_charge_constraint, '-reset-box':self.reset_simulation_cell}
        
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
            
        # Function to convert strings2strings and ints2ints
        def s2s_i2i(inlst):
            outlst = []
            for i in inlst:
                try: j = int(i)
                except: j = str(i)
                outlst.append(j)
            return outlst
        
        
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
            
        # set new -neigh_charge option and print confirmation
        if tags['-neigh_charge']:
            self.neighbor_charge_constraint = tags['-neigh_charge']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-neigh_charge', self.neighbor_charge_constraint))
            
        # set new -types option and print confirmation
        if tags['-types']:
            try: self.types2convert = s2s_i2i(tags['-types'].split(','))
            except:
                print(f'ERROR tag -types atom TypeID list could not be converted to ints {tags["-types"]}'); sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-types', self.types2convert))
            
        # set new -reset-charges option and print confirmation
        if tags['-reset-charges']:
            self.reset_charges = T_F_string2boolean('-reset-charges', (tags['-reset-charges']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-reset-charges', self.reset_charges))
            
        # set new -charge0 option and print confirmation
        if tags['-charge0']:
            self.net_zero_charge = T_F_string2boolean('-charge0', (tags['-charge0']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-charge0', self.net_zero_charge))
            
        # set new -convert2cg1 option and print confirmation
        if tags['-convert2cg1']:
            self.convert2cg1 = T_F_string2boolean('-convert2cg1', (tags['-convert2cg1']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-convert2cg1', self.convert2cg1))
            
        # set new -pi-electrons option and print confirmation
        if tags['-pi-electrons']:
            self.add_pi_electrons = T_F_string2boolean('-pi-electrons', (tags['-pi-electrons']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-pi-electrons', self.add_pi_electrons))
            
        # set new -type-labels option and print confirmation
        if tags['-type-labels']:
            self.include_type_labels = T_F_string2boolean('-type-labels', (tags['-type-labels']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-type-labels', self.include_type_labels))
            
        # set new -reset-box option and print confirmation
        if tags['-reset-box']:
            self.reset_simulation_cell = T_F_string2boolean('-reset-box', (tags['-reset-box']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-reset-box', self.reset_simulation_cell))
        
        # print buffer
        print('\n\n')