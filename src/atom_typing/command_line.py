#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
January 5th, 2023
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


def print_man_page(topofile, bondfile, parent_dir, newfile, ff_name, delete_atoms, mass_map, bondorder, maxbonded,
                   boundary, vdw_radius_scale, reset_charges, print_options, pdb_file, chargefile, include_comments_nta,
                   bonds_via_distance_override):

    
    # print general command line options
    print('\n\natom_typing has been run with -opt or -man option to show the optional command line overrides available. Command line')
    print('option summary [-tag <tag-input>]:')
    print('python3 atom_typing.py [-topo <topo-filename>] [-bond <bond-filename>] [-dir <new directory name>] [-newfile <string>]')
    print('                       [-ff <force field name>] [-reset-charges <T|F>] [-charge-file <Gasteiger-filename>]')
    print('                       [-nta-comments <T|F>] [-vdw-scale <float>] [-boundary <string>] [-bond-reset <T|F>]')
    print('                       [-del-method <mass or size>] [-del-crit<float>] <-gui> <-opt>|<-man>')
    print('                       [*NOTE: Not all options found in atom_typing.py are currently supported via command line overrides.*]')

    print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the atom_typing.py')
    print('file will be enforced. The command line inputs option gives more flexibility to the code depending on IDE usage or')
    print('command line usage. The user may specify as many or as few command line options as desired, but if the command line')
    print('option is used and not all options are provided by the user, the hard coded inputs will be enforced, and the user')
    print('will be warned.*')
    
    print('\nTo use command line override options there must be alteration between tags started with the - character and the')
    print('input value. Ordering of tag/input pairs do not matter. Example:')
    print('    python3 atom_typing.py -tag1 tag1-input  -tag2 tag2-input  -tag3 tag3-input  ...N-arbitrary tagsN/tagN-input pairs')
    
    # print options
    print('\n\nOptions:')
    
    # print -topo option
    print(f'\n -topo or -t <topo-filename>  atom_typing variable: topofile    hard coded: {topofile}')
    print('    Command line option to set topofile name of topo/coord file to be read into atom_typing for conversion. Currently')
    print('    supported topofile formats are:')
    print('         .mol =    .mol file in v2000 format (Info: atom positions, bonds | Set: center about (0,0,0), box set by')
    print('                   min/max atom extents, image flags zeroed)')
    print('         .sdf =    .sdf file in v2000 format (Info: atom positions, bonds | Set: center about (0,0,0), box set by')
    print('                   min/max atom extents, image flags zeroed)')
    print('         .mol2 =   .mol2 SYBYL file  (Info: atom positions, bonds | Set: center about (0,0,0), box set by')
    print('                   min/max atom extents, image flags zeroed)')
    print('         .pdb =    .pdb protein data bank file  (Info: atom positions, bonds | Set: center about (0,0,0), box set by')
    print('                   min/max atom extents, image flags zeroed)')
    print('         .smiles = .SMILES string placed in front of .smiles (Info: atom positions, bonds | Set: center about (0,0,0),')
    print('                   box set by min/max atom extents, image flags zeroed)')
    print('         .data  =  .data or dat file from LAMMPS (Info: atom positions, bonds, no center, box set by LAMMPS, image')
    print('                   flags set by LAMMPS)')
    print('    Example usage:')
    print('        python3 atom_typing.py -topo EXAMPLE_TOPO-FILE.data')
    
    # print -bond option
    print(f'\n -bond or -b <bond-filename>  atom_typing variable: bondfile    hard coded: {bondfile}')
    print('    Command line option to set bondfile name of ReaxFF bond order file to be read into atom_typing.py for')
    print('    generating bonds from a ReaxFF simulation. The code will read all time instances of bond orders and')
    print('    average them together. It will then apply a minimum bond order set by the bond types in bondorder dictionary')
    print('    to find bonds based on bond order cutoffs. After that the bonds will go through one more phase to check that')
    print('    the maximum number of bonded atoms to any element type does not excede what is listed in maxbonded dictionary.')
    print('    This option is for LAMMPS ReaxFF simulations ONLY. If you choose not to use this option for other atom typing')
    print('    set the variable to n.u. (meaning not used). If you use n.u. bonds will then be found via inter atomic distance')
    print('    searching set by a maximum distance of a combination of the vdw radii between bonding atoms. vdw_bond_scale')
    print('    will be a multiplier value to adjust this search distance. Example usage:')
    print('        python3 atom_typing.py -bond EXAMPLE_BOND-FILE.reaxc')
    
    # print -charge-file  option
    print(f'\n -charge-file  or -qf <Gasteiger-filename>  atom_typing variable: chargefile    hard coded: {chargefile}')
    print('    Command line option to set chargefile containing Gasteiger charge parameters that will be used with the Gasteiger')
    print('    charge method if reset_charges is True. Example usage:')
    print('        python3 atom_typing.py -charge-file   frc_files/Gasteiger_parameters.txt')

    # print -nta-comments  option
    print(f'\n -nta-comments  or -nc <T|F>  atom_typing variable: include_comments_nta    hard coded: {include_comments_nta}')
    print('    Command line option set the Boolean variable (T for True or F for False) to either write comments (T) to the new')
    print('    type assignment file or to not write comments in the new type assignment file (F). During atom typing the code will')
    print('    set comments for each atom (most of the time the comments are "Correctly found", but may differ depending on how the')
    print('    atom type was assigned). If you plan to manually edit the atom types or add the atomtype:NAME to the atom type the')
    print('    comments can quickly ruin the readability of the file. Please refer to the all2lmp.py chapter in the manual for the')
    print('    atomtype:NAME option. Example usage:')
    print('        python3 atom_typing.py -nta-comments T')
    
    # # print -dir option
    print(f'\n -dir or -d <new directory name>   atom_typing variable: parent_directory    hard coded: {parent_dir}')
    print('    Command line option to set new directory name to store all files atom_typing.py can write. The file')
    print('    management is based on the present working directory (pwd) and will store files one directory deep inside the')
    print('    pwd. The -dir option will specify the directory that is inside of pwd. This will allow for easier file ')
    print('    management and allows the user to merge multiple files for a certain system or research project and keep them')
    print('    all organized into a single directory.')
    print('    python variables*. Example usage:')
    print('        python3 atom_typing.py -dir EXAMPLE_TEST_SYSTEM_DIRECTORY')
    
    # print -newfile option
    print(f'\n -newfile or -nf <string>   atom_typing variable: newfile    hard coded: {newfile}')
    print('    Command line option to set the new output filename(s). The following options exist for using the newfile string')
    print('    for setting the output file basenames: ')
    print()
    print("      if newfile starts with ':' or ends with ':'")
    print("        The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix added to")
    print("        the file basename. The following are examples:")
    print("          Suffix (newfile = ':_merged'  and topofile = 'detda.data')")
    print("            basename = 'detda_merged', where the ':' character acts as a placeholder for the topofile basename.")
    print("          Prefix (newfile = 'merged-:'  and topofile = 'detda.data')")
    print("            basename = 'merged-detda', where the ':' character acts as a placeholder for the topofile basename.")
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
    print('        python3 atom_typing.py -newfile  :_typed')
    print('        python3 atom_typing.py -newfile ": typed"  # note the use of the " " bounding characters to pass a whitespace')
    
    # print -pdb option
    print(f'\n -pdb or -p <new file extension>   atom_typing variable: pdb_file    hard coded: {pdb_file}')
    print('    Command line option to write the optional .pdb file and what format to write it in. Currently supported')
    print('    pdb_file options and their meanings (Case matters):')
    print("        - 'skip'    will not write additional *_packmol.pdb and *_packmol.nta files")
    print("        - 'types'   will write additional *_packmol.pdb file but not *_packmol.nta, with the atom types set")
    print("                    in the atom name")
    print("        - 'typeIDs' will write additional *_packmol.pdb and *_packmol.nta files, with atomTypeIDs set in the")
    print("                    atom name column of the .pdb file and atom types set in the .nta file using the")
    print('                    "style type" method')
    print('    Example usage:')
    print('        python3 atom_typing.py -pdb types')
    
    # print -ff option
    print(f'\n -ff or -f <force field name>   atom_typing variable: ff_name    hard coded: {ff_name}')
    print('    Command line option to set force field name to assign correct atom types. Currently supported ff_names (Case matters):')
    print('       * Class2:')
    print('           - PCFF-IFF  most atom-types supported (FF-file:   frc_files/pcff_iff_v1_5_CNT_poly_solv.frc)')
    print('           - PCFF      most atom-types supported (FF-file:   frc_files/pcff.frc )')
    print('           - compass    all atom-types supported (FF-file:   frc_files/compass_published.frc)')
    print()
    print('       * Class1:')
    print('           - CVFF-IFF  most atom-types supported (FF-file:   frc_files/cvff_interface_v1_5.frc)')
    print('           - CVFF      most atom-types supported (FF-file:   frc_files/cvff.frc)')
    print('           - Clay-FF    all atom-types supported (FF-file:   frc_files/clayff.frc)')
    print('           - DREIDING  most atom-types supported (FF-file:   frc_files/all2lmp_dreiding.frc)')
    print('           - OPLS-AA    all atom-types supported (FF-file:   frc_files/oplsaa.frc {EXPEIRMENTAL})')
    print('    Example usage:')
    print('        python3 all2lmp.py -ff PCFF-IFF')
    
    # print -reset-charges option
    print(f'\n -reset-charges or -rq <T|F>   atom_typing variable: reset_charges    hard coded: {reset_charges}')
    print('    Command line option reset charges. This is performed by 1st finding the hybridization of each atom in the system')
    print('    via topological constraints or a mimimization of VSEPR angles onto each atom (when needed). Once the hybridization')
    print('    of each atom is known the code will assign Gasteiger parameters to each atom accordingly and then converge the')
    print('    per atom charge via the Gasteiger method. *NOTE: some hybridization might be set via comparison of angles, this')
    print('     means that for best results your system should be in the geometrical configuration associated with its "True"')
    print('    hybridized geometry (IE Sp1=linear, Sp2=trigonal planar, and Sp3=tetrahedral. Example usage:')
    print('        python3 atom_typing.py -reset-charges T')
    
    # print -vdw-scale option
    print(f'\n -vdw-scale or -vdw <float> atom_typing variable: vdw_radius_scale    hard coded: {vdw_radius_scale}')
    print('    Command line option set the vdw radii scale for finding bonds based on interatomic distances. Example usage:')
    print('        python3 atom_typing.py -vdw-scale 1.0')
    
    # print -boundary option
    string = '-'.join([i for i in boundary.split()])
    print(f'\n -boundary or -by <string> atom_typing variable: boundary    hard coded: {string}')
    print('    Command line option set the boundary for finding bonds based on interatomic distances. Three boundary flags')
    print('    must be provided to set the boundary conditions for finding bonds. Two flags are supported:')
    print('      - f is non-periodic and fixed')
    print('      -p is periodic')
    print('    To set the boundary flags at the commandline a "-" character must be supplied between each of the three flags')
    print('    (i.e. f-f-f or p-f-f or p-p-p or ...), such that there exists no whitespace in the "tag-input". Example usage:')
    print('        python3 atom_typing.py -boundary f-f-p')
    
    # print -bond-reset option
    print(f'\n -bond-reset or -br <T|F> atom_typing variable: bonds_via_distance_override    hard coded: {bonds_via_distance_override}')
    print('    Command line option set the over ride flag to reset bonds via interatomic distance searching. The bonds that will be found')
    print('    are dependant on the vdw_radius_scale, boundary, and maxbonded inputs. Example usage:')
    print('        python3 atom_typing.py -bond-reset T')
    
    # print -del-method option
    print(f'\n -del-method or -dm atom_typing variable: delete_atoms["method"]    hard coded: {delete_atoms["method"]}')
    print('    Command line option set the method of how to select which atoms will be deleted. The following methods are available:')
    print('       "mass"  will search for clusters via mass cut-off where anything less than the -del-crit will be')
    print('               removed before any analysis occurs such as setting charge or finding atom types.')
    print()
    print('       "size"  will search for clusters via number of atoms cut-off where anything less than the -del-crit')
    print('               will be removed before any analysis occurs such as setting charge or finding atom types.')
    print('    Example usage:')
    print('        python3 atom_typing.py -del-method mass  -del-crit 50.0')
    
    # print -del-crit option
    print(f'\n -del-crit or -dc atom_typing variable: delete_atoms["criteria"]    hard coded: {delete_atoms["criteria"]}')
    print('    Command line option set the criteria used to select which atoms will be deleted. The following methods are available:')
    print('       "mass"  will search for clusters via mass cut-off where anything less than the -del-crit will be')
    print('               removed before any analysis occurs such as setting charge or finding atom types.')
    print()
    print('       "size"  will search for clusters via number of atoms cut-off where anything less than the -del-crit')
    print('               will be removed before any analysis occurs such as setting charge or finding atom types.')
    print('    Example usage:')
    print('        python3 atom_typing.py -del-method mass  -del-crit 50.0')
    
    # print -opt or -man option
    print(f'\n -opt or -man   atom_typing variable: print_options    hard coded: {print_options}')
    print('    Will tell atom_typing to only print out avaiable command line options known as tagN and tagN-inputs. This is the')
    print('    only tag that doesnt require a tagN-input following the tag since the code will only look for if -opt is in')
    print('    command line input list. Example usage:')
    print('        python3 cluster_analysis.py -man')
    
    # print -gui option
    print('\n -gui <no addition info required> auto_morse_bond_update variable: use_GUI')
    print('    Command line option flag to load GUI. Example usage:')
    print('        python3 atom_typing.py -gui')
    
    # print requirements and dependencies
    print('\n\nRequirements and Dependencies:')
    print('    Requirements: python 3.7+ to guarantee ordering of dictionaries which stores most of the data')
    print('    Dependencies:')
    print('      * python rdkit module:')
    print('          - pip3 install rdkit (if pip manager is installed - adds functionality to read smiles strings)')
    
    # print buffer and exit code
    print('\n\n')
    return

# Class to update inputs
class inputs:
    def __init__(self, topofile, bondfile, parent_directory, newfile, ff_name, delete_atoms, mass_map, bondorder, maxbonded,
                 boundary, vdw_radius_scale, reset_charges, print_options, pdb_file, chargefile, include_comments_nta, 
                 bonds_via_distance_override, commandline_inputs):
        
        # Give access to inputs (update later on if command line over ride is given)
        self.commandline_inputs = commandline_inputs
        self.topofile = topofile
        self.bondfile = bondfile
        self.parent_directory = parent_directory
        self.newfile = newfile
        self.ff_name = ff_name
        self.delete_atoms = delete_atoms 
        self.mass_map = mass_map
        self.bondorder = bondorder
        self.maxbonded = maxbonded
        self.boundary = boundary
        self.vdw_radius_scale = vdw_radius_scale
        self.reset_charges = reset_charges
        self.print_options = print_options
        self.pdb_file = pdb_file
        self.chargefile = chargefile
        self.include_comments_nta = include_comments_nta
        self.bonds_via_distance_override = bonds_via_distance_override

        
        
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
        supported_tags = ['-topo', '-bond', '-dir', '-newfile', '-ff', '-reset-charges', '-vdw-scale', '-pdb', '-charge-file', '-nta-comments',
                          '-vdw-scale', '-boundary', '-bond-reset', '-del-method', '-del-crit']
        
        # set shortcut_tags mapping
        shortcut_tags = {'-t':'-topo', '-b':'-bond', '-d':'-dir', '-nf':'-newfile', '-f':'-ff', '-rq':'-reset-charges', '-vdw':'-vdw-scale', '-p':'-pdb',
                         '-qf':'-charge-file', '-nc':'-nta-comments', '-vdw':'-vdw-scale', '-br':'-bond-reset', '-by':'-boundary', '-dm': '-del-method', 
                         '-dc': '-del-crit'}
        
        # set default variables
        default_variables ={'-topo': self.topofile, '-bond': self.bondfile, '-dir': self.parent_directory, '-ext': self.newfile, '-ff': self.ff_name,
                            '-reset-charges': self.reset_charges, '-vdw-scale': self.vdw_radius_scale, '-pdb':self.pdb_file,
                            '-charge-file':self.chargefile, '-newfile':self.newfile, '-nta-comments':self.include_comments_nta,
                            '-boundary':self.boundary, '-bond-reset':self.bonds_via_distance_override, '-del-method': self.delete_atoms['method'],
                            '-del-crit':self.delete_atoms['criteria']}
        
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
        
        # set new -bond option and print confirmation
        if tags['-bond']:
            self.bondfile = tags['-bond']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-bond', self.bondfile))
            
        # set new -charge-file option and print confirmation
        if tags['-charge-file']:
            self.chargefile = tags['-charge-file']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-charge-file', self.chargefile))
            
        # set new -pdb option and print confirmation
        if tags['-pdb']:
            self.pdb_file = tags['-pdb']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-pdb', self.pdb_file))
                
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
            
        # set new -ff option and print confirmation
        if tags['-ff']:
            self.ff_name = tags['-ff']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-ff', self.ff_name))
        
        # set new -reset-charges option and print confirmation
        if tags['-reset-charges']:
            self.reset_charges = T_F_string2boolean('-reset-charges', (tags['-reset-charges']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-reset-charges', self.reset_charges))
            
        # set new -nta-comments option and print confirmation
        if tags['-nta-comments']:
            self.include_comments_nta = T_F_string2boolean('-nta-comments', (tags['-nta-comments']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-nta-comments', self.include_comments_nta))
            
        # set new -vdw-scale option and print confirmation (ERROR checks will occur in the next step so only try float except set as input)
        if tags['-vdw-scale']:
            try:
                self.vdw_radius_scale = float(tags['-vdw-scale'])
            except:
                self.vdw_radius_scale = tags['-vdw-scale']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-vdw-scale', self.vdw_radius_scale))
            
        # set new -boundary option and print confirmation (ERROR checks will occur in the next step so only try float except set as input)
        if tags['-boundary']:
            self.boundary = ' '.join(tags['-boundary'].split('-'))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-boundary', self.boundary))
            
        # set new -bond-reset option and print confirmation
        if tags['-bond-reset']:
            self.bonds_via_distance_override = T_F_string2boolean('-bond-reset', (tags['-bond-reset']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-bond-reset', self.bonds_via_distance_override))
            
        # set new -del-method option and print confirmation
        if tags['-del-method']:
            self.delete_atoms['method'] = tags['-del-method']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-del-method', self.delete_atoms['method']))
            
        # set new -del-crit option and print confirmation
        if tags['-del-crit']:
            try: 
                self.delete_atoms['criteria'] = float(tags['-del-crit'])
                print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-del-crit', self.delete_atoms['criteria']))
            except:
                print('Override FAILED for {:<18}, likely because tag input was not a float or an int: {}'.format('-del-crit', tags['-del-crit']))

        # print buffer
        print('\n\n')