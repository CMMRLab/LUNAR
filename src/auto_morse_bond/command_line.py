#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
November 13th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

File contains command line options/man page and
command line over riding class to reset hard coded
inputs.

python3 auto_morse_bond_update.py -topo detda_typed_IFF.data -min-r0 1.3 -skip 1,2,3,4,5,6
python3 auto_morse_bond_update.py -topo detda_typed_IFF.data -min-r0 1.3 -skip 1,2,3,4,5,6 -alpha 1.1
python3 auto_morse_bond_update.py -topo detda_typed_IFF.data -min-r0 1.3 -skip 1,2,3,4,5,6 -alpha 1.1 -zex F
python3 auto_morse_bond_update.py -topo detda_typed_IFF.data -min-r0 1.3 -skip 1,2,3,4,5,6 -alpha 1.1 -xterms F -bondbreak 1.75 -rcut T
"""

##############################
# Import Necessary Libraries #
##############################
import sys


def print_man_page(topofile, morsefile, parent_directory, newfile, mass_map, min_bond_length, coeffs2skip,
                   radius_specs, alpha_specs, alpha_scale, files2write, atom_style,
                   zero_effected_xterms, bondbreak_scale, ff_class, include_type_labels,
                   class2xe_update, include_rcut):

    
    # print general command line options
    print('\n\nauto_morse_bond has been run with -opt or -man option to show the optional command line overrides available. Command line option summary [-tag <tag-input>]:')
    print('python3 auto_morse_bond_update.py [-topo <topo-filename>] [-dir <new directory name>] [-newfile <string>] [-atomstyle <atomstyle>]')
    print('                                  [-class <|1|2>] [-type-labels <T|F>] [-min-r0 <float>] [-skip <string of IDs>] [-alpha <float>]')
    print('                                  [-rcut <T|F>] [-bondbreak <float>] [-xterms <T|F>] [-class2xe <T|F>] [-morse <morse parameter filename>]')
    print('                                  <-gui> <-opt>|<-man> [*NOTE: Not all options found in auto_morse_bond_update.py are currently supported')
    print('                                  via command line overrides.*] ')

    print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the auto_morse_bond_update.py')
    print('file will be enforced. The command line inputs option gives more flexibility to the code depending on IDE usage or command line')
    print('usage. The user may specify as many or as few command line options as desired, but if the command line option is used and not all')
    print('options are provided by the user, the hard coded inputs will be enforced, and the user will be warned.*')
    
    print('\nTo use command line override options there must be alteration between tags started with the - character and the')
    print('input value. Ordering of tag/input pairs do not matter. Example:')
    print('    python3 auto_morse_bond_update.py -tag1 tag1-input  -tag2 tag2-input  -tag3 tag3-input  ...N-arbitrary tagsN/tagN-input pairs')
    
    # print options
    print('\n\nOptions:')
    
    # print -topo option
    print(f'\n -topo or -t <topo-filename>  auto_morse_bond_update variable: topofile    hard coded: {topofile}')
    print('    Command line option to set topofile name of topo/coord file to be read into auto_morse_bond_update. Currently')
    print('    supported topofile formats are:')
    print('         .data = .data or dat file from LAMMPS (with bond coeffs either harmonic or class2)')
    print('    Example usage:')
    print('        python3 auto_morse_bond_update.py -topo EXAMPLE_TOPO-FILE.data')
    
    # print -morse option
    print(f'\n -morse or -m <morse parameter filename>  auto_morse_bond_update variable: morsefile    hard coded: {morsefile}')
    print('    Command line option to set morsefile name that contains the bond typing rules and the corresponding Morse bond')
    print('    dissociation energy coefficients. The morsefile currently has most of the bond typing rules set for all elemental')
    print('    and hybridization bonding configurations of the elements set in the mass_map dictionary. The morsefile may be added')
    print('    onto or adjusted as desired to users liking, however, the file is meant to be as comprehensive as possible and should')
    print('    rarely need to be modified. Example usage:')
    print('        python3 auto_morse_bond_update.py -morse frc_files/Morse_parameters.txt')
    
    # print -dir option
    print(f'\n -dir or -d <new directory name>   auto_morse_bond_update variable: parent_directory    hard coded: {parent_directory}')
    print('    Command line option to set new directory name to store all files bond_react_template_merge can write. The file')
    print('    management is based on the present working directory (pwd) and will store files one directory deep inside the')
    print('    pwd. The -dir option will specify the directory that is inside of pwd. This will allow for easier file ')
    print('    management and allows the user to merge multiple files for a certain system or research project and keep them')
    print('    all organized into a single directory.')
    print('    python variables*. Example usage:')
    print('        python3 auto_morse_bond_update.py -dir EXAMPLE_TEST_SYSTEM_DIRECTORY')
    
    # print -newfile option
    print(f'\n -newfile or -nf <string>   auto_morse_bond_update variable: newfile    hard coded: {newfile}')
    print('    Command line option to set the new output filename(s). The following options exist for using the newfile string')
    print('    for setting the output file basenames: ')
    print()
    print("      if newfile starts with ':' or ends with ':'")
    print("        The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix added to")
    print("        the file basename. The following are examples:")
    print("          Suffix (newfile = ':_morse_bond'  and topofile = 'detda.data')")
    print("            basename = 'detda_morse_bond', where the ':' character acts as a placeholder for the topofile basename.")
    print("          Prefix (newfile = 'morse_bond-:'  and topofile = 'detda.data')")
    print("            basename = 'morse_bond-detda', where the ':' character acts as a placeholder for the topofile basename.")
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
    print('        python3 auto_morse_bond_update.py -newfile  :_morse_bond')
    print('        python3 auto_morse_bond_update.py -newfile ": morse_bond"  # note the use of the " " bounding characters to pass a whitespace')
    
    # print -atomstyle option
    print(f'\n -atomstyle or -as <atomstyle>   auto_morse_bond_update variable: atom_style    hard coded: {atom_style}')
    print('    Command line option to set atom style in written LAMMPS .data file. Currently supported atom styles are full,')
    print('    charge, molecular. *NOTE: If another atomstyle is desired it can be added in the write_lmp.py script.* Example')
    print('    usage:')
    print('        python3 auto_morse_bond_update.py -atomstyle full')
    
    # print -class option
    print(f'\n -class or -c <1|2>   auto_morse_bond_update variable: ff_class    hard coded: {ff_class}')
    print('    Command line option to set the class of forcefield used. This will determine how the bonds/angles/diherdals/')
    print('    impropers are found, like wise how to coeffs are found and defined. There are 4 supported class types, 0, 1, 2')
    print('    and r, each of these MUST BE used with their correct .frc file. If the code detects the .frc file and class is')
    print('    inconsistent it will exit the code execution. Most cases have been tested, but it still is up to the user to')
    print('    make sure their inputs are consistent, and the outputs make sense. The class options are as follows with known')
    print('    FFs that use each class:')
    print('      1 = class1 (for cvff, clayff force field file in frc_files directory)')
    print('      2 = class2 (for PCFF-IFF, PCFF, compass force field file in frc_files directory)')
    print('    Example usage:')
    print('        python3 auto_morse_bond_update.py -class 2')
        
    # print -type-labels option
    print(f'\n -type-labels or -tl <T|F>   auto_morse_bond_update variable: include_type_labels    hard coded: {include_type_labels}')
    print('    Command line option to inlcude LAMMPS new type labels in written datafile to make bond/angle/dihedrals/improper')
    print('    types if applicable consistent between datafiles or molecule files (primary use for bond/react file generation)')
    print('    *NOTE: this option will also change how -write-bond-react molecule files are written to be consistent between')
    print('    the written .data file and the written molecule .moltemp file.* Example usage:')
    print('        python3 auto_morse_bond_update.py -type-labels T')
        
    # print -min-r0 option
    print(f'\n -min-r0 or -mr0 <float> auto_morse_bond_update variable: min_bond_length    hard coded: {min_bond_length}')
    print('    Command line option to set the minimum bond length r0 that a bond coeff must have to update from a')
    print('    harmonic bond to a morse bond. Example usage:')
    print('        python3 auto_morse_bond_update.py -min-r0 1.2')
        
    # print -skip option
    print(f'\n -skip or -c2s <string of IDs> auto_morse_bond_update variable: coeffs2skip    hard coded: {coeffs2skip}')
    print('    Command line option to set the any bond coeffID to skip updating from harmonic to morse bond. This info will')
    print('    be passed on through a list when running auto_morse_bond_update.py in IDE mode or a tag with the value of the')
    print('    tag seperated by commas with no whitespace. If the following list is empty it means that no coeffs2skip will')
    print('    be enforced. Example usage:')
    print('        python3 auto_morse_bond_update.py -skip 1,2,4')
        
    # print -alpha option
    print(f'\n -alpha or -a <float> auto_morse_bond_update variable: alpha_scale    hard coded: {alpha_scale}')
    print('    Command line option to set the alpha_scale which is a multiplier value to adjust the best fit alpha parameter.')
    print('    An alpha scale of 1.2 and 0.8 will increase the fit alpha by 20% and decrease the alpha by 20% respectively.')
    print('    A default value of 1.0 should be used, where adjusting this is only a special case scenerio. Example usage:')
    print('        python3 auto_morse_bond_update.py -alpha 1.1')
            
    # print -xterms option
    print(f'\n -xterms or -zex <T|F> auto_morse_bond_update variable: zero_effected_xterms    hard coded: {zero_effected_xterms}')
    print('    Command line option to set the boolean flag to zero the effected crossterms. This option is only applicable')
    print('    for class 2 FFs (ff_class = 2), since class 1 FFs do not contain any crossterms. Zeroing the crossterms when')
    print('    adding morse bonds to a class 2 force field helps keep the simulation from crashing since class 2 FFs assume')
    print('    that the harmonic bond can constrain the bond lengths to keep the crossterms from increasing dramatically.)')
    print('    T is for True and F is for False. Example usage:')
    print('        python3 auto_morse_bond_update.py -xterms T')
    
    # print -xterms option
    print(f'\n -class2xe or -2xe <T|F> auto_morse_bond_update variable: class2xe_update    hard coded: {class2xe_update}')
    print('    Command line option to convert class2 force field to a class2xe (x=crossterms, e=exponential), to allow for the')
    print('    crossterms to dissociate like a Morse bond. The following crossterms are updated to:')
    print('       E_bondbond = D*(1-e^(-alpha(r-r1)))*(1-e^(-alpha(r-r2)))')
    print('       E_bondangle = D1*(1-e^(-alpha(r-r1)))*(theta-theta0) + (1-e^(-alpha(r-r2)))*(theta-theta0)')
    print('       E_middlebondtorsion = (1-e^(-alpha(r-r2)))*[A1cos(phi) + A2cos(2*phi) + A3cos(3*phi)]')
    print('       E_endbondtorsion = (1-e^(-alpha(r-r1)))*[B1cos(phi) + B2cos(2*phi) + B3cos(3*phi)] +')
    print('                          (1-e^(-alpha(r-r3)))*[C1cos(phi) + C2cos(2*phi) + C3cos(3*phi)]')
    print('       E_bondbond13 = D*(1-e^(-alpha(r-r1)))*(1-e^(-alpha(r-r3)))')
    print()
    print('    If class2xe_update is T (True), the crossterms above will be reparmaterized and if F (False), the crossterms')
    print('    will be left alone. NOTE that when using class2xe_update the min_bond_length and coeffs2skip criteria')
    print('    will not be used and every bond coeff that is possible to update to a Morse bond will be updated, as')
    print('    the implementation of the class2xe potential requires all bonds to be Morse bonds. Example usage:')
    print('        python3 auto_morse_bond_update.py -class2xe_update T')
            
    # print -bondbreak option
    print(f'\n -bondbreak or -bb <float> auto_morse_bond_update variable: bondbreak_scale    hard coded: {bondbreak_scale}')
    print('    Command line option to set the bondbreak_scale which is a multiplier value to adjust the multiplier')
    print('    value of the bond coeff r0. The default bond break scale should be 1.75 (175% of r0). This value will')
    print('    be used to generate the LAMMPS include script with fix bond/break settings with the rcut value set as')
    print('    r0*bondbreak_scale. If include_rcut = True (used) the bond break scale will also be used to set the')
    print('    value of rcut in the bond coeffs. The rcut value in the bond coeffs section will be used to shift the')
    print('    morse potential. Alternatively if bondbreak_scale is set to 0, the Rmax will be set at the')
    print('    "dissociation point" of the Morse potential. Example usage:')
    print('        python3 auto_morse_bond_update.py -bondbreak 1.75')
    
    # print -rcut option
    print(f'\n -rcut or -irc <T|F> auto_morse_bond_update variable: include_rcut    hard coded: {include_rcut}')
    print('    Command line option to set the boolean flag to include the rcut value in the morse bond coeff. The')
    print('    rcut value in the bond coeffs section will be used to shift the morse potential. Example usage:')
    print('    T is for True and F is for False. Example usage:')
    print('        python3 auto_morse_bond_update.py -rcut T')
    
    # print -opt or -man option
    print('\n -opt or -man')
    print('    Will tell auto_morse_bond_update to only print out avaiable command line options known as tagN and tagN-inputs.')
    print('    This id the only tag that doesnt require a tagN-input following the tag since the code will only look for if')
    print('    -opt is in command line input list. Example usage:')
    print('        python3 auto_morse_bond_update.py -man')
    
    # print -gui option
    print('\n -gui <no addition info required> auto_morse_bond_update variable: use_GUI')
    print('    Command line option flag to load GUI. Example usage:')
    print('        python3 auto_morse_bond_update.py -gui')
    
    # print requirements and dependencies
    print('\n\nRequirements and Dependencies:')
    print('    Requirements: python 3.7+ to guarantee ordering of dictionaries which stores most of the data')
    print('    Dependencies:')
    print('      * python matplotlib module:')
    print('          - pip3 install matplotlib (if pip manager is installed - for plotting the overlay of the harmonic and morse bond)')

    
    # print buffer and exit code
    print('\n\n')
    return

# Class to update inputs
class inputs:
    def __init__(self, topofile, morsefile, parent_directory, newfile, min_bond_length, coeffs2skip,
                       alpha_scale, atom_style, zero_effected_xterms, bondbreak_scale,
                       ff_class, include_type_labels, class2xe_update, include_rcut, commandline_inputs):
        # Give access to inputs (update later on if command line over ride is given)
        self.commandline_inputs = commandline_inputs
        self.topofile = topofile
        self.morsefile = morsefile
        self.parent_directory = parent_directory
        self.newfile = newfile
        self.atom_style = atom_style
        self.ff_class = ff_class
        self.include_type_labels = include_type_labels
        self.min_bond_length = min_bond_length
        self.coeffs2skip = coeffs2skip
        self.alpha_scale = alpha_scale
        self.zero_effected_xterms = zero_effected_xterms
        self.bondbreak_scale = bondbreak_scale
        self.include_rcut = include_rcut
        self.class2xe_update = class2xe_update
        
        
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
        supported_tags = ['-topo', '-dir', '-newfile', '-atomstyle', '-class', '-type-labels', '-min-r0', 
                          '-skip', '-alpha', '-xterms', '-bondbreak', '-rcut', '-morse', '-class2xe']
        
        # set shortcut_tags mapping
        shortcut_tags = {'-t':'-topo', '-d':'-dir', '-nf':'-newfile', '-as':'-atomstyle', '-c':'-class',
                         '-tl': '-type-labels', '-mr0':'-min-r0', '-c2s':'-skip', '-a':'-alpha',
                         '-zex':'-xterms', '-bb':'-bondbreak', '-irc':'-rcut', '-m':'-morse', '-2xe':'-class2xe'}
        
        # set default variables
        default_variables = {'-topo': self.topofile, '-dir':parent_directory, '-newfile': self.newfile, '-atomstyle': self.atom_style,
                             '-class': self.ff_class, '-type-labels': self.include_type_labels, '-min-r0': self.min_bond_length,
                             '-skip': self.coeffs2skip, '-alpha': self.alpha_scale, '-xterms':self.zero_effected_xterms,
                             '-bondbreak':self.bondbreak_scale, '-rcut':self.include_rcut, '-morse':self.morsefile,'-class2xe':self.class2xe_update}
        
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
        
        # set new -morse option and print confirmation
        if tags['-morse']:
            self.morsefile = tags['-morse']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-morse', self.morsefile))
        
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
        
        # set new -class option and print confirmation
        if tags['-class']:
            self.ff_class = str(tags['-class'])
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-class', self.ff_class))
            
        # set new -type-labels option and print confirmation
        if tags['-type-labels']:
            self.include_type_labels = T_F_string2boolean('-type-labels', (tags['-type-labels']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-type-labels', self.include_type_labels))
            
        # set new -class2xe option and print confirmation
        if tags['-class2xe']:
            self.class2xe_update = T_F_string2boolean('-class2xe', (tags['-class2xe']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-class2xe', self.class2xe_update))
            
        # set new '-min-r0' option and print confirmation (ERROR checks will occur in the next step so only try float except set as input)
        if tags['-min-r0']:
            try: self.min_bond_length = float(tags['-min-r0'])
            except: 
                self.min_bond_length = tags['-min-r0']
                print('ERROR supplied -min-r0 or -mr0 could not be converted to a float'); sys.exit();
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-min-r0', self.min_bond_length))
            
        # set new -skip option and print confirmation
        if tags['-skip']:
            try: self.coeffs2skip = [int(i) for i in tags['-skip'].split(',')]
            except:
                print(f'ERROR tag -skip atomID TypeIDs list could not be converted to ints {tags["-skip"]}'); sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-skip', self.coeffs2skip))
            
        # set new '-alpha' option and print confirmation (ERROR checks will occur in the next step so only try float except set as input)
        if tags['-alpha']:
            try: self.alpha_scale = float(tags['-alpha'])
            except:
                print('ERROR supplied -alpha or -a could not be converted to a float'); sys.exit();
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-alpha', self.alpha_scale))
            
        # set new -xterms option and print confirmation
        if tags['-xterms']:
            self.zero_effected_xterms = T_F_string2boolean('-zex', (tags['-xterms']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-xterms', self.zero_effected_xterms))
            
        # set new '-bondbreak' option and print confirmation (ERROR checks will occur in the next step so only try float except set as input)
        if tags['-bondbreak']:
            try: self.bondbreak_scale = float(tags['-bondbreak'])
            except: 
                print('ERROR supplied -bondbreak or -bb could not be converted to a float'); sys.exit();
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-bondbreak', self.bondbreak_scale))
            
        # set new -rcut option and print confirmation
        if tags['-rcut']:
            self.include_rcut = T_F_string2boolean('-rcut', (tags['-rcut']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-rcut', self.include_rcut))
            
        # print buffer
        print('\n\n')