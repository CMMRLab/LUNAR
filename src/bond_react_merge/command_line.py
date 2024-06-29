# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
March 27th, 2024
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


def print_man_page(files, parent_dir, newfile, atom_style, generate_map_file, write_rxn_mol2files,
                   write_rxn_datafiles, write_moleculefiles, print_options, map_near_edge_rxn_charges,
                   molecule_file_options, include_type_labels):
    
    # function to string dict and once string get to long bump to next line
    def string_dict(dictionary):
        string = ''; max_length = 120; length_tally = 0; str_lst = []
        for tag in dictionary:
            file = dictionary[tag]; tag = '{}'.format(tag);
            string += '{}:{},'.format(tag, file)
            length_tally += len(string)
            if length_tally > max_length:
                length_tally = 0; str_lst.append(string); string = '';
            else: str_lst.append(string)
        
        # print to screen and remove last comma
        for n, i in enumerate(str_lst):
            if n+1 == len(str_lst): i = i[:-1];
            print(f'    {i}')
        return string[:-1], str_lst

    
    # print general command line options
    print('\n\nbond_react_merge has been run with -opt or -man option to show the optional command line overrides available.')
    print('Command line option summary [-tag <tag-input>]:\n')
    print('python3 bond_react_merge.py [-files <files-string>] [-dir <new directory name>] [-atomstyle <atomstyle>]')
    print('                            [-map <T|F>] [-write-rxn-datafiles <T|F>] [-type-labels<T|F>] [-newfile <string>]')
    print('                            [-write-rxn-mol2files <T|F>] [-edge <F|0|1|2|N>] [-write-moleculefiles <T|F>]')
    print('                            <-gui> <-opt>|<-man>')

    print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the')
    print('bond_react_template_merge.py file will be enforced. The command line inputs option gives more flexibility to the')
    print('code depending on IDE usage or command line usage. The user may specify as many or as few command line options')
    print('as desired, but if the command line option is used and not all options are provided by the user, the hard coded')
    print('inputs will be enforced, and the user will be warned.*')
    
    print('\nTo use command line override options there must be alteration between tags started with the - character and the')
    print('input value. Ordering of tag/input pairs do not matter. Example:')
    print('    python3 bond_react_merge.py -tag1 tag1-input  -tag2 tag2-input ...N-arbitrary tagsN/tagN-input pairs')
    
    # print options
    print('\n\nOptions:')

    # print -files option
    print('\n -files or -f <files string>   bond_react_merge variable: files    hard coded (converted from dict):')
    string_dict(files)
    print('\n    Command line option to pass in files information that will then be converted into a python dictionary. This')
    print('    will allow the user to provide the files/file-tags needed to run bond_react_template_merge at the command line.')
    print('    The <tag-input> MUST BE A SINGLE STRING with the following format (with comma delimiter and NO whitespaces):')
    print('        file-tag1:file1.data,file-tag2:file2.data,file-tagN:fileN.data   or  infile:merge_files.txt')
    print('    were the file-tag has to be preN or postN or dataN depending on associated and desired outputs. The N-value MUST')
    print('    BE CONSISTANT between each preN/postN reaction pair of files. Example usege:')
    print('        python3 bond_react_merge.py  -files  data1:mol1.data,data2:mol2.data,pre1:pre1.data,post1:post1.data')
    print('                                      or')
    print('        python3 bond_react_merge.py  -files  infile:merge_files.txt')
    
    # print -dir option
    print(f'\n -dir or -d <new directory name>   bond_react_merge variable: parent_directory    hard coded: {parent_dir}')
    print('    Command line option to set new directory name to store all files bond_react_template_merge can write. The file')
    print('    management is based on the present working directory (pwd) and will store files one directory deep inside the')
    print('    pwd. The -dir option will specify the directory that is inside of pwd. This will allow for easier file ')
    print('    management and allows the user to merge multiple files for a certain system or research project and keep them')
    print('    all organized into a single directory.')
    print('    python variables*. Example usage:')
    print('        python3 bond_react_merge.py -dir EXAMPLE_TEST_SYSTEM_DIRECTORY')
    
    # print -newfile option
    print(f'\n -newfile or -nf <string>   bond_react_merge variable: newfile    hard coded: {newfile}')
    print('    Command line option to set the new output filename(s). The following options exist for using the newfile string')
    print('    for setting the output file basenames: ')
    print()
    print("      if newfile starts with ':' or ends with ':'")
    print("        The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix added to")
    print("        the file basename. The following are examples:")
    print("          Suffix (newfile == ':_merged'  and topofile == 'detda.data')")
    print("            basename = 'detda_merged', where the ':' character acts as a placeholder for the topofile basename.")
    print("          Prefix (newfile == 'merged-:'  and topofile == 'detda.data')")
    print("            basename = 'merged-detda', where the ':' character acts as a placeholder for the topofile basename.")
    print('        Recommended usage: common and safe as this method is safe and output filename(s) carry similar names to')
    print('        input filename(s).')
    print()
    print("      if newfile == 'ANYTEXT'")
    print("        The output filename(s) will be set as the file tag of each file and 'ANYTEXT' will be a suffix appended")
    print("        to the file tag. For example:")
    print("          newfile == '-merged', topofile == 'detda.data', and tag == 'data1'")
    print("            basename = 'data1-merged'")
    print('        Recommended usage: occasional and safe as output filename(s) no longer carry similar names as output')
    print('        filename(s), but is safe as output filename(s) will not overwrite input filename(s)')
    print()
    print("      if newfile == ''")
    print("        The output filename(s) will be the same as the input filename(s). This can be a dangerous option as you")
    print("        may inadvertently overwrite a file and then must assume the file to contain certain information, but it")
    print("        contains other information.")
    print("        Recommended usage: rare and dangerous as this could lead to using incorrect files if not paying attention")
    print("        very carefully.")
    print('    Example usage:')
    print('        python3 bond_react_merge.py -newfile  :_merged')
    print('        python3 bond_react_merge.py -newfile ": merged"  # note the use of the " " bounding characters to pass a whitespace')
    
    # print -atomstyle option
    print(f'\n -atomstyle or -as <atomstyle>   bond_react_merge variable: atom_style    hard coded: {atom_style}')
    print('    Command line option to set atom style in written LAMMPS .data file. Currently supported atom styles are full,')
    print('    charge, molecular. *NOTE: If another atomstyle is desired it can be added in the write_lmp.py script.* Example')
    print('    usage:')
    print('        python3 bond_react_merge.py -atomstyle full')
    
    # print -map option
    print(f'\n -map or -m <T|F>   bond_react_merge variable: generate_map_file    hard coded: {generate_map_file}')
    print('    Command line option to auto generate the rxn map file based on a pseudo graph-theory/cost based algorithm.')
    print('    This option is ONLY AVAILABLE IF THE datafiles were generated using all2lmp v1.4+ since that atom sections')
    print('    has the styled atom comment of atom-type/element, because the elemental information and atom-type information')
    print('    will be used to generate the most accurate map file as possible. In general the generated map file should be')
    print('    ready for use in bond/react, IT IS STILL UP TO THE USER TO MAKE SURE THE MAP FILE IS CORRECT!. Example Usage:')
    print('        python3 bond_react_merge.py -map T')
    
    # print -edge option
    print(f'\n -edge or -e <F|0|1|2|N>   bond_react_merge variable: generate_map_file    hard coded: {map_near_edge_rxn_charges}')
    print('    Command line option to map charges near edge atoms from dataN or postN tagged files. *NOTE: It is assumed that')
    print('    there is a dataN or postN file that has an exact topological match to the preN tagged files. Once the preN tagged')
    print('    file is found to be a subset of the dataN or postN tagged file the charges will be mapped up to N-neighbors away')
    print('    and the postN tagged files will be updated via the the equivalences found via the generate_map_file option. THIS')
    print('    MEANS TO RUN THIS OPTION. THE generate_map_file OPTION MUST BE True* The following option are valid for input and')
    print('    their meanings:')
    print('        F = DO NOT attempt to match topologies for preN tagged files to dataN or postN tagged files and leave')
    print('            charges as is in preN/postN tagged files')
    print('        0 = match preN tagged files to dataN or postN tagged file and update only edge atom charges (int value')
    print('            sets depth from edge)')
    print('        1 = match preN tagged files to dataN or postN tagged file and update edge atom charges + 1 neighbor-deeper')
    print('            (int value sets depth from edge)')
    print('        2 = match preN tagged files to dataN or postN tagged file and update edge atom charges + 2 neighbor-deeper')
    print('            (int value sets depth from edge)')
    print('        N = match preN tagged files to dataN or postN tagged file and update edge atom charges + N neighbor-deeper')
    print('            (int value sets depth from edge)')
    print('    MAX INTGER VALUE CAN BE 20 currently since the topology matching looks for exact topological matches 20-neighbors')
    print('    away from any identified edge atoms. This also implies that if a topology is matched from any preN tagged file to')
    print('    any postN tagged file that the topology can be a guaranteed match of up to 16-neighbors deep from the edge atoms.')
    print('    Since the postN tagged file charges will also be updated via exact topological matching it means that the max-depth')
    print('    you should ever set with the map_near_edge_rxn_charges variable is one that is the max depth the the preN topology')
    print('    matches the postN topology, otherwise you will charges on the postN topology that should be left unchanged. A depth')
    print('    of 0 is the safest option since the topologies and connections of preN-to-postN tagged file must be identical. The')
    print('    matching seems to work well, but sometimes the code can get confused if there are symmetries on the preN tagged')
    print('    topology.')
    print()
    print('    The near edge mapping of charges will use either a dataN tagged file or a postN tagged file to try to match the')
    print('    preN tagged file tologogies to. The reason for this is if you are building a set of rxn-templates were the pre-rxn1')
    print('    topology generates a post-rxn1 topology, but then the post-rxn1 topology gets used as a molecule in a pre-rxn2')
    print('    topology then the only way to find the charges for the pre-rxn2 topology is to look at the post-rxn1 topology. The')
    print('    code handles this type of scenerio by updating the charges on the lowest N-tagged pre/post pair before moving onto')
    print('    the N+1-tagged pair. This method assumes that if the rxn-templates have a "sequential" progression, that the N-tagged')
    print('    value increases as the "sequential" progression increases. Therefore when setting up the file tags IT IS OF IMPORTANCE')
    print('    IF YOU PLAN TO USE THIS OPTION THAT ANY "SEQUENTIAL" TEMPLATE PAIRS HAVE N-TAGGED LABELS THAT ENCODE THE ORDERING OF')
    print('    THE EVOLUTION OF TOPOLOGIES OTHER WISE THE CODE WILL NOT UPDAE THE CHARGES PROPERPLY!')
    print()
    print('    The rational to this option is that when the partial topologies are made the edge atoms themselves are cut away')
    print('    from other parts of the orginal molecules found in the the dataN tagged files. Then depending on the charge method')
    print('    used to set the charges the edge atoms and possibly some N-depth deeper from the neighbors may have incorrect')
    print('    charging. The LAMMPS commnd bond/react currently super-imposes all attributes set in the molecule.moltemp files')
    print('    INCLUDING the edge atom charges. So these charges should be matched to the dataN tagged files manually or via some')
    print('    other methods. This option was built to help facilitate in automation of this issue in bond/react template creation.')
    print('    IT IS ULTIMATELY UP TO THE USER TO CHECK IF THE CHARGES WERE MAPPED PROPERLY. The code writes comments to each atomID')
    print('    in charge section of the molecule.moltemp file that was changed to due to the usage of this option.')
    print('        python3 bond_react_merge.py -edge 0')   
    
    # print -molecule option
    print(f'\n -molecule or -mol <string of options>   bond_react_merge variable: molecule_file_options    hard coded: {molecule_file_options}')
    print('    Options to add to molecule files (generate_map_file MUST BE True for these options to work):')
    print('    fragment/ID/OPTION (ordering matters and the delimiter is the / character.)')
    print('         index0 = fragment envokes an option to add fragment section to the molecule file')
    print('         index1 = ID is user defined ID of fragment')
    print('         index2 = OPTION is user defined option.')
    print('    Currently available options:')
    print('         custom_charges/depth_from_edge_N where custom_charges specifies the the code to build options for bond/reacts')
    print('                                          custom_charges option and depth_from_edge_N tells the code to build a fragmentID')
    print('                                          that only has atomic charges a certain depth N away from any edge atom. This is')
    print('                                          useful if the charge method used could not assign a proper charge to the edge atoms')
    print('                                          of the molecule template. bond_react_merge.py has another option:')
    print('                                          map_near_edge_rxn_charges that works fairly well to map charges from any dataN or')
    print('                                          postN tagged files onto edge atoms and a certain depth from edge atoms, but this is')
    print('                                          another work around to the charge issue of atoms near or edge atoms or edge atoms')
    print('                                          themselves. Much like map_near_edge_rxn_charges the meaning of depth for this option')
    print('                                          is as follows:')
    print('                                              0 = preN tagged files will have a .moltemp molecule template file with a fragmentID')
    print('                                                  with all atoms specified EXCEPT edge atoms at depth 0 from any edge.')
    print('                                              1 = preN tagged files will have a .moltemp molecule template file with a fragmentID')
    print('                                                  with all atoms specified EXCEPT edge atoms at depth 1 from any edge.')
    print('                                              N = preN tagged files will have a .moltemp molecule template file with a fragmentID')
    print('                                                  with all atoms specified EXCEPT edge atoms at depth N from any edge')
    print('                                          Example python string specifying this option to have a fragmentID in the preN/postN')
    print('                                          tagged molecule files that only has atoms in the fragmentID that are AT LEAST 2 deep')
    print('                                          from any edge atom:')
    print('                                                fragment/charge_edge_example/custom_charges/depth_from_edge_2')
    print()
    print('         custom_charges/equiv_cost_N where custom_charges specifies the the code to build options for bond/reacts custom_charges')
    print('                                     option and equiv_cost_N tells the code to build a fragmentID that only has atomic charges that')
    print('                                     have a total cost of equivalence fitting greater then N. This is useful to help gauge which')
    print('                                     charges should be updated based on the cost of the equivalence fit. If the cost is >0 it means')
    print('                                     that atom-type is new or its neighbors up to 2-neighbors away have changed and can be useful')
    print('                                     to help decide which charges should be updated. The meaning of N is as follows:')
    print('                                         0 = build a fragment whose total cost of equivalence fit for the map file is greater then 0.')
    print('                                         1 = build a fragment whose total cost of equivalence fit for the map file is greater then 1.')
    print('                                         N = build a fragment whose total cost of equivalence fit for the map file is greater then N.')
    print('                                     Example python string specifying this option to have a fragmentID in the preN/postN tagged molecule')
    print('                                     files that only has atoms in the fragmentID that were fit with a cost greater then 0:')
    print('                                     fragment/charge_cost_example/custom_charges/equiv_cost_0')
    print()
    print('    molecule will find moleculeIDs via a cluster analysis and order them largest-to-smallest and then adds a molecule section to')
    print('    the molecule files.')
    print()
    print('    This info will be passed on through a list when running bond_react_merge.py in IDE mode or a tag with the value of the tag')
    print('    seperated by commas with no whitespace. If the following list is empty it means that no molecule file options will be added')
    print('    to the written .moltemp files.')
    print('        python3 bond_react_merge.py -molecule fragment/charge_edge_example/custom_charges/depth_from_edge_0,molecule')   

    # print -write-rxn-datafiles option
    print(f'\n -write-rxn-datafiles or -wrd <T|F>    bond_react_merge variable: write_rxn_datafiles    hard coded: {write_rxn_datafiles}')
    print('    Command line option to write template datafiles for preN/postN tagged files. This option is useful to make sure')
    print('    all coeffs were mapped properly and allows the user to visualize the merge tempalte file in OVTIO of VMD. T is')
    print('    for True and use option and F is for False. Example usage:')
    print('        python3 bond_react_merge.py -write-rxn-datafiles T')

    # print -write-rxn-mol2files option
    print(f'\n -write-rxn-mol2files or -wrm <T|F>    bond_react_merge variable: write_rxn_mol2files    hard coded: {write_rxn_mol2files}')
    print('    Command line option to write template mol2 files for preN/postN tagged files. This option is useful to make sure')
    print('    for visualize the preN/postN reaction in VMD or ChemDraww. Morse specifically ChemDraw, because ChemDraw shows a')
    print('    table of atom-ids and allows multiple files to be open at once. Usage of this option and -map or -m option go')
    print('    hand and hand for quick checking to map sure the Equivalence section of the auto-generate map file makes sense.')
    print('    T is for True and use option and F is for False. Example usage:')
    print('        python3 bond_react_merge.py -write-rxn-mol2files T')
    
    # print -write-moleculefiles option
    print(f'\n -write-moleculefiles or -wmf <T|F>    bond_react_merge variable: write_moleculefiles    hard coded: {write_moleculefiles}')
    print('    Command line option to to write out a *_merged.lmpmol file for any file that has a dataN tag. This file will have all the')
    print('    coeff types and all topological types are mapped properly and can be used with the LAMMPS “create_atoms” command to build')
    print('    large random systems. T is for True and use option and F is for False. Example usage:')
    print('        python3 bond_react_merge.py -write-moleculefiles T')
    
    # print -type-labels option
    print(f'\n -type-labels or -tl <T|F>   bond_react_merge variable: include_type_labels    hard coded: {include_type_labels}')
    print('    Command line option to inlcude LAMMPS new type labels in written datafile to make bond/angle/dihedrals/improper')
    print('    types if applicable consistent between datafiles or molecule files (primary use for bond/react file generation)')
    print('    *NOTE: this option will also change how -write-bond-react molecule files are written to be consistent between')
    print('    the written .data file and the written molecule .moltemp file.* Example usage:')
    print('        python3 bond_react_merge.py -type-labels T')

    # print -opt or -man option
    print(f'\n -opt or -man   bond_react_merge variable: print_options    hard coded: {print_options}')
    print('    Will tell bond_react_merge to only print out avaiable command line options known as tagN and tagN-inputs. This is the')
    print('    only tag that doesnt require a tagN-input following the tag since the code will only look for if -opt is in')
    print('    command line input list. Example usage:')
    print('        python3 bond_react_merge.py -man')
    
    # print -gui option
    print('\n -gui <no addition info required> bond_react_merge variable: use_GUI')
    print('    Command line option flag to load GUI. Example usage:')
    print('        python3 bond_react_merge.py -gui')
    
    # print requirements and dependencies
    print('\n\nRequirements and Dependencies:')
    print('    Requirements: python 3.7+ to guarantee ordering of dictionaries for generate_map_file algorithm')
    
    # print buffer and exit code
    print('\n\n')
    return

# Class to update inputs
class inputs:
    def __init__(self, files, parent_dir, newfile, atom_style, generate_map_file, write_rxn_mol2files,
                 write_rxn_datafiles, write_moleculefiles, print_options, map_near_edge_rxn_charges,
                 molecule_file_options, include_type_labels, commandline_inputs):
        
        # Give access to inputs (update later on if command line over ride is given)
        self.commandline_inputs = commandline_inputs
        self.files = files
        self.parent_dir = parent_dir
        self.newfile = newfile
        self.atom_style = atom_style
        self.generate_map_file = generate_map_file
        self.write_rxn_mol2files = write_rxn_mol2files
        self.write_rxn_datafiles = write_rxn_datafiles
        self.write_moleculefiles = write_moleculefiles
        self.print_options = print_options
        self.map_near_edge_rxn_charges = map_near_edge_rxn_charges
        self.molecule_file_options = molecule_file_options
        self.include_type_labels = include_type_labels
        
        
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
        supported_tags = ['-files', '-dir', '-atomstyle', '-map', '-write-rxn-datafiles', '-write-rxn-mol2files',
                          '-write-moleculefiles', '-edge', '-molecule', '-type-labels', '-newfile']
        
        # set shortcut_tags mapping
        shortcut_tags = {'-f':'-files', '-d':'-dir', '-as':'-atomstyle', '-m':'-map' ,'-wmf':'-write-moleculefiles', '-wrd':'-write-rxn-datafiles',
                         '-wrm':'-write-rxn-mol2files', '-wmd':'-write-mol-datafiles', '-e':'-edge', '-mol':'-molecule', '-tl': '-type-labels',
                         '-nf':'-newfile'}
        
        # set default variables
        default_variables ={'-files':self.files, '-dir':self.parent_dir, '-atomstyle':self.atom_style, '-map':self.generate_map_file,
                            '-write-rxn-datafiles':self.write_rxn_datafiles, '-write-rxn-mol2files':self.write_rxn_mol2files, '-write-moleculefiles':self.write_moleculefiles,
                            '-edge': self.map_near_edge_rxn_charges, '-molecule': self.molecule_file_options, '-type-labels': self.include_type_labels,
                            '-newfile': self.newfile}
        
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
                print('WARNING override option   {:<24} not provided at the command line. Hard coded input being enforced: {}'.format(i, default_variables[i]))
        
                
        # Start changing python variables based on command line over rides
        print('\n\nCommand line run option override inputs confirmation (will confirm to user that the command line override does override the hard coded inputs):')
        # Function to convert T string to True and F string to False
        def T_F_string2boolean(tag, string):
            if string == 'F': return False
            elif string == 'T': return True
            else: print(f'\nERROR Provided tag: {tag} string is {string} which is not T or F. Use T or F string\n'); sys.exit()
            
        # Function to set -edge update variable accorindingly
        def edge_commandline_convert(tag, string):
            try:
                return int(string)
            except:
                return T_F_string2boolean(tag, string)
                
            
        # Fuction to get files dict from string
        def files_dict_from_srt(files_str):
            files_lst = files_str.split(','); files = {}; # { filetag : file }
            for entry in files_lst:
                entry = entry.split(':')
                
                # Warn and exit if entry is not two in length
                if len(entry) != 2:
                    print(f'ERROR -files or -f was used at command line with incorrect argument {"".join(entry)} in string'); sys.exit();
                
                # Find tag and file
                tag = entry[0]; file = entry[1];
                
                # Warn and exit if tag already exists in files
                if tag in files:
                    print(f'ERROR -files or -f was used at command line and multiple instances of tag: {tag} was provided (tags MUST BE UNIQUE)'); sys.exit();
                
                # add tag/file pair to files
                files[tag] = file
            return files
        
        
        ###############################################
        # set new -topo option and print confirmation #
        ###############################################  
        # set new -files option and print confirmation
        if tags['-files']:
            self.files = files_dict_from_srt(tags['-files'])
            print('Override confirmation for {:<22} Hard codeded input is being overridden with this input: {}'.format('-files', self.files))
          
        # set new -dir option and print confirmation
        if tags['-dir']:
            if tags['-dir'].isspace() or tags['-dir'] == '':
                print('-dir flag was used, with an empty directory string. The directory string must contain alpha-numeric characters')
            else:
                self.parent_dir = tags['-dir']
                print('Override confirmation for {:<22} Hard codeded input is being overridden with this input: {}'.format('-dir', self.parent_dir))
        
        # set new -atomstyle option and print confirmation
        if tags['-atomstyle']:
            self.atom_style = tags['-atomstyle']
            print('Override confirmation for {:<22} Hard codeded input is being overridden with this input: {}'.format('-atomstyle', self.atom_style))
            
        # set new -newfile option and print confirmation
        if tags['-newfile']:
            self.newfile = tags['-newfile']
            print('Override confirmation for {:<22} Hard codeded input is being overridden with this input: {}'.format('-newfile', self.newfile))
            
        # set new -map option and print confirmation
        if tags['-map']:
            self.generate_map_file = T_F_string2boolean('-map', (tags['-map']))
            print('Override confirmation for {:<22} Hard codeded input is being overridden with this input: {}'.format('-map', self.generate_map_file))
            
        # set new -edge option and print confirmation
        if tags['-edge']:
            self.map_near_edge_rxn_charges = edge_commandline_convert('-edge', (tags['-edge']))
            print('Override confirmation for {:<22} Hard codeded input is being overridden with this input: {}'.format('-edge', self.map_near_edge_rxn_charges))
            
        # set new -molecule option and print confirmation
        if tags['-molecule']:
            self.molecule_file_options = tags['-molecule'].split(',')
            print('Override confirmation for {:<22} Hard codeded input is being overridden with this input: {}'.format('-molecule', self.molecule_file_options))
            
        # set new -write-rxn-datafiles option and print confirmation
        if tags['-write-rxn-datafiles']:
            self.write_rxn_datafiles = T_F_string2boolean('-write-rxn-datafiles', (tags['-write-rxn-datafiles']))
            print('Override confirmation for {:<22} Hard codeded input is being overridden with this input: {}'.format('-write-rxn-datafiles', self.write_rxn_datafiles))
            
        # set new -write-rxn-datafiles option and print confirmation
        if tags['-write-rxn-mol2files']:
            self.write_rxn_mol2files = T_F_string2boolean('-write-rxn-mol2files', (tags['-write-rxn-mol2files']))
            print('Override confirmation for {:<22} Hard codeded input is being overridden with this input: {}'.format('-write-rxn-mol2files', self.write_rxn_mol2files))
            
        # set new -write-moleculefiles option and print confirmation
        if tags['-write-moleculefiles']:
            self.write_moleculefiles = T_F_string2boolean('-write-moleculefiles', (tags['-write-moleculefiles']))
            print('Override confirmation for {:<22} Hard codeded input is being overridden with this input: {}'.format('-write-moleculefiles', self.write_moleculefiles))
            
        # set new -type-labels option and print confirmation
        if tags['-type-labels']:
            self.include_type_labels = T_F_string2boolean('-type-labels', (tags['-type-labels']))
            print('Override confirmation for {:<22} Hard codeded input is being overridden with this input: {}'.format('-type-labels', self.include_type_labels))
        
        # print buffer
        print('\n\n')