# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
June 14th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

File contains command line options/man page and
command line over riding class to reset hard coded
inputs.

cell_builder run:
    python3 cell_builder.py -files detda_typed_IFF_merged.data:1,dgeba_typed_IFF_merged.data:2 -tl T -as full -tl T
    python3 cell_builder.py -files detda_typed_IFF_merged.data:1,dgeba_typed_IFF_merged.data:2 -tl T -as full -nf cell1 -dup 100 -ds 1.5 -rm T

"""

##############################
# Import Necessary Libraries #
##############################
import sys


def print_man_page(topofiles, force_field_joining, duplicate, distance_scale, newfile, atom_style, parent_directory, max_rotations,
                   reset_molids, unwrap_atoms_via_image_flags, include_type_labels, group_monomers_locally, seed, domain, maxtry, tolerance,
                   mixing_rule, boundary):
    
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
    print('\n\ncell_builder has been run with -opt or -man option to show the optional command line overrides available.')
    print('Command line option summary [-tag <tag-input>]:\n')
    print('python3 cell_builder.py [-files <files-string>] [-dir <new directory name>] [-atomstyle <atomstyle>] [-seed <int>]')
    print('                        [-newfile <new filename>] [-duplicate <int>] [-dist-scale <float>] [-type-labels <T|F>]')
    print('                        [-rx or -ry or -rz or -rall <float>] [-reset-molids <files|clusters|skip>] [-domain <string>]')
    print('                        [-unwrap <T|F>] [-offset <T|F>] [-ff-join <none|offset|merge>] [-tolerance <float or int>]')
    print('                        [-maxtry <int>] [-mixing <string>] [-boundary <string>] <-gui> <-opt>|<-man>')

    print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the')
    print('cell_builder.py file will be enforced. The command line inputs option gives more flexibility to the code')
    print('depending on IDE usage or command line usage. The user may specify as many or as few command line options')
    print('as desired, but if the command line option is used and not all options are provided by the user, the hard')
    print('coded inputs will be enforced, and the user will be warned.*')
    
    print('\nTo use command line override options there must be alteration between tags started with the - character and the')
    print('input value. Ordering of tag/input pairs do not matter. Example:')
    print('    python3 cell_builder.py -tag1 tag1-input  -tag2 tag2-input ...N-arbitrary tagsN/tagN-input pairs')
    
    # print lattice options
    print('\n\n')
    print('*******************')
    print('* Lattice Options *')
    print('*******************')

    # print -files option
    print('\n -files or -f <files string>   cell_builder variable: files    hard coded (converted from dict):')
    string_dict(topofiles)
    print('\n    Command line option to pass in files information that will then be converted into a python dictionary. This')
    print('    will allow the user to provide the files/qtys needed to run cell_builder at the command line. The <tag-input>')
    print('    MUST BE A SINGLE STRING with the following format (with comma delimiter and NO whitespaces):')
    print('        file1.data:qty1,file2.data:qty2,fileN.data:qtyN')
    print('    were qtyN sets the number of files to use. NOTE to use the tab completion you will have to provide whitespaces')
    print('    between filenames and then go back in the string and remove the whitespaces since <files string> must not')
    print('    contain any whitespaces. Example usege:')
    print('        python3 cell_builder.py  -files  molecule1.data:2,molecule2.data:1')
    
    # print -offset option
    print(f'\n -ff-join  or -ffj <none|offset|merge>   cell_builder variable: force_field_joining    hard coded: {force_field_joining}')
    print('    Command line option to defines how to handle the force field between multiple LAMMPS datafiles, where the following')
    print('    options exist:')
    print('      none   which assumes all LAMMPS coefficient types are the same in all the read in files and applies no offset to the')
    print('             files as they are being used to generate a large molecular system. This option should be used if your files where')
    print('             processed with LUNAR/bond_react_merge.py to ensure that the force field between the output system is consistent')
    print('             with the reaction templates.')
    print('      merge  which applies the merging processes present LUNAR/bond_react_merge.py to merge all coefficient types amongst all')
    print('             read in files. This option requires that all LAMMPS datafiles have the LUNAR/all2lmp.py style of comments.')
    print('      offset which applies an offset to each coefficient type in each file as it is read into cell_builder.py.')
    print('    Example usage:')
    print('        python3 cell_builder.py -ff-join merge')
    
    # print -dir option
    print(f'\n -dir or -d <new directory name>   cell_builder variable: parent_directory    hard coded: {parent_directory}')
    print('    Command line option to set new directory name to store all files bond_react_template_merge can write. The file')
    print('    management is based on the present working directory (pwd) and will store files one directory deep inside the')
    print('    pwd. The -dir option will specify the directory that is inside of pwd. This will allow for easier file ')
    print('    management and allows the user to merge multiple files for a certain system or research project and keep them')
    print('    all organized into a single directory.')
    print('    python variables*. Example usage:')
    print('        python3 cell_builder.py -dir EXAMPLE_TEST_SYSTEM_DIRECTORY')
    
    # print -atomstyle option
    print(f'\n -atomstyle or -as <atomstyle>   cell_builder variable: atom_style    hard coded: {atom_style}')
    print('    Command line option to set atom style in written LAMMPS .data file. Currently supported atom styles are full,')
    print('    charge, molecular. *NOTE: If another atomstyle is desired it can be added in the write_lmp.py script.* Example')
    print('    usage:')
    print('        python3 cell_builder.py -atomstyle full')
    
    # print -new-file option
    print(f'\n -newfile or -nf <new filename>   cell_builder variable: newfile   hard coded: {newfile}')
    print('    Command line option to set new filename. Example usage:')
    print('        python3 cell_builder.py -new-file CELL_DATAFILE')
            
    # print -duplicate option
    print(f'\n -duplicate or -dup <int> cell_builder variable: duplicate    hard coded: {duplicate}')
    print('    Command line option to set the system size by duplicating each file Ntimes. Example usage:')
    print('        python3 cell_builder.py -duplicate 100')
    
    # print -grp-mono option
    print(f'\n -grp-mono or -grp <T/F> cell_builder variable: group_monomers_locally    hard coded: {group_monomers_locally}')
    print('    Command line option to group the monomers locally to keep relative locations of hardener to resin close. Example usage:')
    print('        python3 cell_builder.py -grp-mono T')
                
    # print -dist-scale option
    print(f'\n -dist-scale or -ds <int> cell_builder variable: distance_scale    hard coded: {distance_scale}')
    print('    Command line option to set the size of each subcell that a molecule will be placed in. Each subcell')
    print('    dimension is set that it may hold the largest molecule from the read-in files in any orientation without')
    print('    passing the subcell dimensions. This is to avoid running molecules into other molecules when applying random')
    print('    rotations. distance_scale sets a multiplier value to increase or decrease (at your own risk) the dimension')
    print('    of each subcell and ultimately the spacing of all molecules. Example usage:')
    print('        python3 cell_builder.py -dist-scale 1.2')
    
    # print -seed option
    print(f'\n -seed or -s <int> cell_builder variable: seed    hard coded: {seed}')
    print('    Command line option to set the seed for generating random numbers to account for reproducibility concerns of')
    print('    generating random initial positions. If the seed value is zero the seed will default to using the current system')
    print('    time and if greater then zero, the seed will be applied to the random number generator. Example usage:')
    print('        python3 cell_builder.py -seed 12345')
        
    # print -reset-molids option
    print(f'\n -reset-molids or -rm <files|clusters|skip>   cell_builder variable: reset_molids    hard coded: {reset_molids}')
    print('    Command line option to reset molids based on different methods, where the following options are available:')
    print('        skip     will use the molIDs set in the files that are read in. If the file format supplied to cell_builder.py')
    print('                 does not have molIDs, every atom molID will default to one.')
    print()
    print('        files    will set molIDs based on the order the files are read into cell_builder.py, where every atom in the')
    print('                 first file will be assigned to molID one, every atom in the second file will be assigned to molID two')
    print('                 and so on.')
    print()
    print('        offset   will offset the molIDs in each file as each file is read in. If every atoms molID are the same in each')
    print('                 file the offset method and the files method will create the same result. If the file contains different')
    print('                 molIDs on different atoms and you wish to maintain the distinction, the offset is the method to use.')
    print()
    print('        clusters will perform a cluster analysis, where the “clusters” of atoms are determined via bonding connectivity,')
    print('                 where the criteria for atoms to be part of the same “cluster” is that the atoms must be linked by at')
    print('                 least one covalent bond to the cluster. The clusters are then sorted by the number of atoms, where the')
    print('                 largest number of atoms is identified as cluster one, then molIDs are incremented and are assigned to each')
    print('                 of the remaining clusters (i.e., the largest cluster of atoms will have molID 1, and then the smallest')
    print('                 cluster will have a molID of NCLUSTERS found).')
    print()
    print('    Depending on the different analysis and/or visualization the different reset_molid options may be useful, since they')
    print('    can be used to identity groups of atoms, which can be tedious depending on the type of modeling that is being used.')
    print('    Example usage:')
    print('        python3 cell_builder.py -reset-molids files')
    
    # print -rx or -ry or -rz or -rall
    print(f'\n -rx or -ry or -rz or -rall <float>   cell_builder variable: max_rotations    hard coded: {max_rotations}')
    print('    Command line option to adjust the maximum degree of rotation that any molecule or system will undergo during the')
    print('    random selection of phi, theta, and psi values (Degrees). The following -tags are available:')
    print('        -rx   <float> adjust maximum rotation in around x-direction')
    print('        -ry   <float> adjust maximum rotation in around y-direction')
    print('        -rz   <float> adjust maximum rotation in around z-direction')
    print('        -rall <float> adjust maximum rotation in around x, y, and z-directions with one tag')
    print('    Example usage:')
    print('        python3 cell_builder.py -rx 90 -ry 180 -rz 270')
    print('        python3 cell_builder.py -rall 360')
    
    # print -unwrap option
    print(f'\n -unwrap or -u <T|F>   cell_builder variable: unwrap_atoms_via_image_flags    hard coded: {unwrap_atoms_via_image_flags}')
    print('    Command line option to unwrap atoms via image flags when reading in the the files. This functionality assumes')
    print('    consistent image flags. The default should be to always unwrap atoms via image flags. There are a few scenerios that')
    print('    the atoms can stay "wrapped"; the datafiles have non-periodic molecules or special ReaxFF cases where bonds will be')
    print('    infered via inter-atomic distnaces of a periodic system. Example usage:')
    print('        python3 cell_builder.py -unwrap T')
    
    # print -type-labels option
    print(f'\n -type-labels or -tl <T|F>   cell_builder variable: include_type_labels    hard coded: {include_type_labels}')
    print('    Command line option to inlcude LAMMPS new type labels in written datafile to make bond/angle/dihedrals/improper')
    print('    types if applicable consistent between datafiles or molecule files (primary use for bond/react file generation)')
    print('    *NOTE: this option will also change how -write-bond-react molecule files are written to be consistent between')
    print('    the written .data file and the written molecule .moltemp file.* Example usage:')
    print('        python3 cell_builder.py -type-labels T')
    
    # print Random options
    print('\n\n')
    print('******************')
    print('* Random Options *')
    print('******************')
    
    # print -domain option
    print(f'\n -domain or -dn <string> cell_builder variable: domain    hard coded: {domain}')
    print('    Command line option to set the cubic lattice domain. The lattice points will always be generated about the 0, 0, 0')
    print('    position in x, y, and z. However, the number of lattice points in the x, y, and z directions can be set by the user')
    print('    or set to be cubic. The following options are available:')
    print('        cubic           which automatically determines the number of lattice points required based on the qty of files and the')
    print('                        duplicate variable.')
    print()
    print('        Ni x Nj x Nk    where Ni is the number of lattice points in the x-direction, Nj is the number of lattice points')
    print('                        in the y-direction, and Nk is the number of lattice points in the z-direction. Please note the')
    print('                        following about Ni x Nj x Nk:')
    print('                           If group_monomers_locally is False, Ni x Nj x Nk must be greater than duplicate*sum(qty of all')
    print('                                                               files), to allow for enough lattice points. If there is not')
    print('                                                               enough lattice points, the code will exit with an ERROR.')
    print('                           If group_monomers_locally is True,  Ni x Nj x Nk must be greater than duplicate value, to allow')
    print('                                                               for enough lattice points. If there is not enough lattice points,')
    print('                                                               the code will exit with an ERROR. Additionally, by default during')
    print('                                                               the initial grouping of monomers will occur with a cubic lattice.')
    print()
    print("        LxA x LyA x LzA where 'Lx' is the box size in the X-direction and 'A' differentiates it from Ni, 'Ly'")
    print("                        is the box size in the Y-direction and 'A' differentiates it from Nj, and 'Lz' is the")
    print("                        box size in the Z-direction and 'A' differentiates it from Nz. This method will randomly")
    print("                        place molecules (not on a lattice) and randomly rotate molecules. This method can be    ")
    print("                        computationally intensive. However it can create systems nearing densities of 0.55 g/cc ")
    print("                        in under a few minutes (depending on the random settings).                              ")
    print()
    print("        A x A x A       where the three 'A' means simulation box is being defined with files of QTY or ZEROs.   ")
    print('    Example usage:')
    print('        python3 cell_builder.py -domain 2x2x4')
    print('        python3 cell_builder.py -domain 10Ax20Ax30A')
    print('        python3 cell_builder.py -dn cubic')
    
    # print -maxtry option
    print(f'\n -maxtry or -mt <int>   cell_builder variable: maxtry   hard coded: {maxtry}')
    print('    Command line option to control the number of times to try to randomly insert a molecule into a system. After each')
    print('    insertion of a molecule, the next insertion becomes even more difficult as, the simulation cell is getting denser with')
    print('    each insertaion. Example usage:')
    print('        python3 cell_builder.py -maxtry 100')
    
    # print -tolerance option
    print(f'\n -tolerance or -tol <int or float>   cell_builder variable: tolerance   hard coded: {tolerance}')
    print("    Command line option to control how to check for overlaps in atoms during molecule insertion. If it is in int or float,")
    print("    it effects how to check for overlaps as follows:                                                         ")
    print("            tolerance = <float>, which means all atom diameters will be modeled as that <float> input        ")
    print("                                 provided.                                                                   ")
    print("            tolerance = <int>,   which sets the index of the sigma value in the LAMMPS Pair Coeff section    ")
    print("                                 of the LAMMPS datafile. For example Pair Coeffs are read from the LAMMPS    ")
    print("                                 data file as:                                                               ")
    print("                                     Pair Coeffs # lj/class2/coul/long                                       ")
    print("                                                                                                             ")
    print("                                     1  0.054  4.01 # [0.054, 4.01] -> index=1, sigam=4.01                   ")
    print("                                     2  0.054  3.90 # [0.054, 3.90] -> index=1, sigam=3.90                   ")
    print("                                     3  0.013  1.11 # [0.013, 1.11] -> index=1, sigam=1.11                   ")
    print("                                     :   :      :   :       :       :    :         :                         ")
    print('    Example usage:')
    print('        python3 cell_builder.py -tolerance 1   -mixing sixthpower')
    print('        python3 cell_builder.py -tolerance 2.0 -mixing tolerance')
    
    # print -mixing option
    print(f'\n -mixing or -mix <string>   cell_builder variable: mixing_rule   hard coded: {mixing_rule}')
    print("        Command line option to control how Pair Coeff LJ parameters are mixed if the tolerance variable is")
    print("        an <int>. The following strings are supported:                                           ")
    print("            'tolerance'   which means the tolerance variable is a <float> and to use the float variable to   ")
    print("                          check for atom overlaps. This option maybe needed for ReaxFF model generation, as  ")
    print("                          there are no Pair Coeffs in a ReaxFF LAMMPS datafile.                              ")
    print("            'geometric'   which means mix the i,j LJ parameters using geometric mixing rules (FFs like       ")
    print("                          DREIDING).                                                                         ")
    print("            'arithmetic'  which means mix the i,j LJ parameters using arithmetic mixing rules (FFs like      ")
    print("                          CHARMM).                                                                           ")
    print("            'sixthpower'  which means mix the i,j LJ parameters using sixthpower mixing rules (FFs like      ")
    print("                          PCFF).                                                                             ")
    print("            '-min'        NOTE: the '-min' ending can be appended to 'geometric' or 'arithmetic' or          ")
    print("                          'sixthpower' to create 'geometric-min' or 'arithmetic-min' or 'sixthpower-min',    ")
    print("                          which will multiple the mixed LJ-sigma values by 2^(1/6) to set the overlap        ")
    print("                          condition to place molecules with vdw energy at the LJ-minimum.                    ")  
    print("         The 'sixthpower' mixing rule is the most 'conservative' as it generates the largest mixed LJ sigma  ")
    print("         parameters and thus can ensure no overlapped atoms not matter what mixing rule ends up being        ")
    print("         applied in LAMMPS. Thus the 'sixthpower' mixing rule can be a good default. Examples:               ")
    print("             mixing_rule = 'tolerance'  # will user tolerance <float> to model all atom diameters the same.  ")
    print("             mixing_rule = 'sixthpower' # will combine LJ parameters using the sixthpower mixing rule to     ")
    print("                                        # model all atom diameters like in an MD simulation.                 ")
    print('    Example usage:')
    print('        python3 cell_builder.py -mixing sixthpower -tolerance 1')
    print('        python3 cell_builder.py -mixing tolerance  -tolerance 2.0')
    
    # print -boundary option
    string = '-'.join([i for i in boundary.split()])
    print(f'\n -boundary or -by <string> atom_typing variable: boundary    hard coded: {string}')
    print("        Command line option to control the boundary of the simulation cell, when inserting molecules. The")
    print("        boundary variable is set up like the LAMMPS 'boundary' command, where three flags are provided to")
    print("        set the x, y, or z boundary of the simulation cell. The flags are like LAMMPS flags: ")
    print("            p is periodic                                                                                    ")
    print("            f is non-periodic and fixed                                                                      ")
    print("        When the boundary is 'p p p' or full periodic, each image of each atom is checked, thus checking for ")
    print("        overlaps is more computationaly intensive. However allowing molecules to span the simulation cell    ")
    print("        'opens' more space to possible insert the molecule. This ultimately seems to make the code run time  ")
    print("        quicker when inserting molecules into a dense system as compared to a boundary of 'f f f' or a non   ")
    print("        periodic system. Examples:                                                                           ")
    print("            boundary = 'f-f-f' # non-periodic system                                                         ")
    print("            boundary = 'p-p-p' # fully-periodic system                                                       ")
    print("            boundary = 'p-f-f' # periodic in X-dir and non-periodic in Y- and Z-dir                          ")
    print("            boundary = 'f-f-p' # periodic in Z-dir and non-periodic in X- and Y-dir                          ")
    print('    To set the boundary flags at the commandline a "-" character must be supplied between each of the three flags')
    print('    (i.e. f-f-f or p-f-f or p-p-p or ...), such that there exists no whitespace in the "tag-input". Example usage:')
    print('        python3 atom_typing.py -boundary f-f-p')
    
    print('\n\n')
    print('****************')
    print('* Misc Options *')
    print('****************')
    
    # print -opt or -man option
    print('\n -opt or -man')
    print('    Will tell cell_builder to only print out avaiable command line options known as tagN and tagN-inputs. This is the')
    print('    only tag that doesnt require a tagN-input following the tag since the code will only look for if -opt is in')
    print('    command line input list.')
    
    # print -gui option
    print('\n -gui <no addition info required> auto_morse_bond_update variable: use_GUI')
    print('    Command line option flag to load GUI. Example usage:')
    print('        python3 auto_morse_bond_update.py -gui')
    
    # print requirements and dependencies
    print('\n\nRequirements and Dependencies:')
    print('    Requirements: python 3.7+ to guarantee ordering of dictionaries for generate_map_file algorithm')
    
    # print buffer and exit code
    print('\n\n')
    return

# Class to update inputs
class inputs:
    def __init__(self, topofiles, force_field_joining, duplicate, distance_scale, newfile, atom_style, parent_directory, max_rotations,
                 reset_molids, unwrap_atoms_via_image_flags, include_type_labels, group_monomers_locally, seed, domain, maxtry, tolerance,
                 mixing_rule, boundary, commandline_inputs):
        
        # Give access to inputs (update later on if command line over ride is given)
        self.commandline_inputs = commandline_inputs
        self.topofiles = topofiles
        self.parent_directory = parent_directory
        self.atom_style = atom_style
        self.include_type_labels = include_type_labels
        self.newfile = newfile
        self.duplicate = duplicate
        self.distance_scale = distance_scale
        self.reset_molids = reset_molids
        self.max_rotations = max_rotations
        self.unwrap_atoms_via_image_flags = unwrap_atoms_via_image_flags
        self.group_monomers_locally = group_monomers_locally
        self.seed = seed
        self.force_field_joining = force_field_joining
        self.domain = domain
        self.maxtry = maxtry
        self.tolerance = tolerance 
        self.mixing_rule = mixing_rule
        self.boundary = boundary
        
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
        supported_tags = ['-files', '-dir', '-atomstyle', '-type-labels', '-newfile', '-duplicate', '-dist-scale', '-reset-molids', '-rx', '-ry',
                          '-rz', '-rall', '-unwrap', '-grp-mono', '-seed', '-ff-join', '-domain', '-maxtry', '-tolerance', '-mixing', '-boundary']
        
        # set shortcut_tags mapping
        shortcut_tags = {'-f':'-files', '-d':'-dir', '-as':'-atomstyle', '-tl': '-type-labels', '-nf':'-newfile', '-dup':'-duplicate',
                         '-ds':'-dist-scale', '-rm':'-reset-molids', '-rx':'-rx', '-ry':'-ry', '-rz':'-rz', '-rall':'-rall', '-u':'-unwrap',
                         '-grp':'-grp-mono', '-s':'-seed', '-ffj':'-ff-join', '-dn':'-domain', '-mt':'-maxtry', '-tol':'-tolerance', '-mix':'-mixing',
                         '-by':'-boundary'}
        
        # set default variables
        default_variables ={'-files':self.topofiles, '-dir':self.parent_directory, '-atomstyle':self.atom_style, '-type-labels': self.include_type_labels,
                            '-newfile':self.newfile, '-duplicate':self.duplicate, '-dist-scale':self.distance_scale, '-reset-molids': self.reset_molids,
                            '-rall':self.max_rotations, '-rx':self.max_rotations['x'], '-ry':self.max_rotations['y'], '-rz':self.max_rotations['z'],
                            '-unwrap':self.unwrap_atoms_via_image_flags, '-grp-mono':self.group_monomers_locally, '-seed':self.seed,
                            '-ff-join':self.force_field_joining, '-domain':self.domain, '-maxtry':self.maxtry, '-tolerance':self.tolerance,
                            '-mixing':self.mixing_rule, '-boundary':self.boundary}
        
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
        rtags = [i for i in tags if i in ['-rx', '-ry', '-rz', '-rall'] and not tags[i]]
        for i in tags:
            if i in ['-rx', '-ry', '-rz', '-rall'] and rtags: continue
            if not tags[i]:
                print('WARNING override option   {:<24} not provided at the command line. Hard coded input being enforced: {}'.format(i, default_variables[i]))
        
                
        # Start changing python variables based on command line over rides
        print('\n\nCommand line run option override inputs confirmation (will confirm to user that the command line override does override the hard coded inputs):')
        # Function to convert T string to True and F string to False
        def T_F_string2boolean(tag, string):
            if string == 'F': return False
            elif string == 'T': return True
            else: print(f'\nERROR Provided tag: {tag} string is {string} which is not T or F. Use T or F string\n'); sys.exit()                
            
        # Fuction to get files dict from string
        def files_dict_from_srt(files_str):
            files_lst = files_str.split(','); files = {}; # { file : qty }
            for entry in files_lst:
                entry = entry.split(':')
                
                # Warn and exit if entry is not two in length
                if len(entry) != 2:
                    print(f'ERROR -files or -f was used at command line with incorrect argument {"".join(entry)} in string'); sys.exit();
                
                # Find tag and file
                file = entry[0]; qty = entry[1];
                
                try: qty = int(qty)
                except:
                    print(f'ERROR file: {file} qty: {qty} could not be converted to an int value')
                    sys.ext()
                
                # Warn and exit if tag already exists in files
                if file in files:
                    print(f'ERROR -files or -f was used at command line and multiple instances of file: {file} was provided (filenames MUST BE UNIQUE)'); sys.exit();
                
                # add tag/file pair to files
                files[file] = qty
            return files
        
        
        ###############################################
        # set new -topo option and print confirmation #
        ###############################################  
        # set new -files option and print confirmation
        if tags['-files']:
            self.topofiles = files_dict_from_srt(tags['-files'])
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-files', self.topofiles))
          
        # set new -dir option and print confirmation
        if tags['-dir']:
            if tags['-dir'].isspace() or tags['-dir'] == '':
                print('-dir flag was used, with an empty directory string. The directory string must contain alpha-numeric characters')
            else:
                self.parent_directory = tags['-dir']
                print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-dir', self.parent_directory))
        
        # set new -atomstyle option and print confirmation
        if tags['-atomstyle']:
            self.atom_style = tags['-atomstyle']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-atomstyle', self.atom_style))
            
        # set new -domain option and print confirmation
        if tags['-domain']:
            self.domain = tags['-domain']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-domain', self.domain))
            
        # set -new-file option and print confirmation
        if tags['-newfile']:
            self.newfile = tags['-newfile']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-newfile', self.newfile))
            
        # set new -duplicate option and print confirmation
        if tags['-duplicate']:
            try: self.duplicate = int(tags['-duplicate']) # try getting int
            except: 
                print('ERROR -duplicate input could not be converted to an int value')
                sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-duplicate', self.duplicate))
            
        # set new -dist-scale option and print confirmation
        if tags['-dist-scale']:
            try: self.distance_scale = float(tags['-dist-scale']) # try getting float
            except: 
                print('ERROR -dist-scale input could not be converted to a float value')
                sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-dist-scale', self.distance_scale))    
            
        # set new -seed option and print confirmation
        if tags['-seed']:
            try: self.seed = int(tags['-seed']) # try getting int
            except: 
                print('ERROR -seed input could not be converted to an integer value')
                sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-seed', self.seed))    
            
        # set new -reset-molids option and print confirmation
        if tags['-reset-molids']:
            self.reset_molids = tags['-reset-molids']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-reset-molids', self.reset_molids))
            
        # set new -ff-join option and print confirmation
        if tags['-ff-join']:
            self.force_field_joining = tags['-ff-join']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-ff-join', self.force_field_joining))
            
        # set new -grp-mono option and print confirmation
        if tags['-grp-mono']:
            self.group_monomers_locally = T_F_string2boolean('-grp-mono', (tags['-grp-mono']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-grp-mono', self.group_monomers_locally))
            
        # set new -rall option and print confirmation
        if tags['-rall']:
            try: rot = float(tags['-rall'])
            except: print(f'-rall: {tags["-rall"]} could not be converted to a float'); sys.exit()
            self.max_rotations['x'] = rot; self.max_rotations['y'] = rot; self.max_rotations['z'] = rot
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-rall', self.max_rotations))
            
        # set new -rx option and print confirmation
        if tags['-rx']:
            try: self.max_rotations['x'] = float(tags['-rx'])
            except: print(f'-rx: {tags["-rx"]} could not be converted to a float'); sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-rx', self.max_rotations))
            
        # set new -ry option and print confirmation
        if tags['-ry']:
            try: self.max_rotations['y'] = float(tags['-ry'])
            except: print(f'-ry: {tags["-ry"]} could not be converted to a float'); sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-ry', self.max_rotations))
            
        # set new -rz option and print confirmation
        if tags['-rz']:
            try: self.max_rotations['z'] = float(tags['-rz'])
            except: print(f'-rz: {tags["-rz"]} could not be converted to a float'); sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-rz', self.max_rotations))
            
        # set new -unwrap option and print confirmation
        if tags['-unwrap']:
            self.unwrap_atoms_via_image_flags = T_F_string2boolean('-type-labels', (tags['-unwrap']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-unwrap', self.unwrap_atoms_via_image_flags))
            
        # set new -type-labels option and print confirmation
        if tags['-type-labels']:
            self.include_type_labels = T_F_string2boolean('-type-labels', (tags['-type-labels']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-type-labels', self.include_type_labels))
            
        # set new -maxtry option and print confirmation
        if tags['-maxtry']:
            try: self.maxtry = int(tags['-maxtry']) # try getting int
            except: 
                print('ERROR -maxtry input could not be converted to an int value')
                sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-maxtry', self.maxtry))
            
        # set new -tolerance option and print confirmation
        if tags['-tolerance']:
            if '.' in tags['-tolerance']:
                self.tolerance = float(tags['-tolerance'])
            else:
                self.tolerance = int(tags['-tolerance'])
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-tolerance', self.tolerance))
            
        # set new -mixing option and print confirmation
        if tags['-mixing']:
            self.mixing_rule = tags['-mixing']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-mixing', self.mixing_rule))
            
        # set new -boundary option and print confirmation (ERROR checks will occur in the next step so only try float except set as input)
        if tags['-boundary']:
            self.boundary = ' '.join(tags['-boundary'].split('-'))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-boundary', self.boundary))
        
        # print buffer
        print('\n\n')