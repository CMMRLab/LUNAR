#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
September 21st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

File contains command line options/man page and
command line over riding class to reset hard coded
inputs.

# Testing global inputs:
python3 sheet_builder.py -dir . -run sheet -bonds T -pbc T -r0 1.42
python3 sheet_builder.py -dir . -run sheet -bonds T -pbc T -r0 1.42 -t1 CP1 -t2 CP2 -t3 CP3 -t4 CP4

# full sheet test:
python3 sheet_builder.py -dir . -run sheet -bonds T -pbc T -r0 1.42 -sheet-name test1 -plane xy -sheet-edge zigzag -le 30 -lp 50 -stacking ABC -ss 10 -nl 6

# full symmetric tube test:
python3 sheet_builder.py -dir . -run symmetric-tube -bonds T -pbc T -r0 1.42 -stname sym-tube -sa z -tedge zigzag -sl 10 -sd 5 -ts 10 -nt 3

# full chiral tube test:
python3 sheet_builder.py -dir . -run chiral-tube -bonds T -pbc T -r0 1.42 -ctname chi-tube -n 10 -m 15 -ca z -cl 50
"""

##############################
# Import Necessary Libraries #
##############################
import sys


def print_man_page(sheet_basename, symmetric_tube_basename, chiral_tube_basename, run_mode, parent_directory, length_in_perpendicular, length_in_edgetype, sheet_edgetype, types,
                   bond_length, sheet_layer_spacing, sheet_nlayers, stacking, plane, tube_edgetype, tube_layer_spacing, symmetric_ntubes, symmetric_length, diameter, n, m,
                   chiral_length, symmetric_tube_axis, chiral_tube_axis, find_bonds, periodic_bonds, seed, functional_atoms, terminating_atoms, grafting_files, minimum_distance):

    
    # print general command line options
    print('\n\nsheet_builder.py has been run with -opt or -man option to show the optional command line overrides available. Command line')
    print('option summary [-tag <tag-input>]:')
    print('python3 sheet_builder.py [-dir <new directory name>] [-run-mode <string>] [-find-bonds <T|F>] [-pbc-bonds <T|F>] [-bond-length <float>]')
    print('                         [-type1, -type2, -type3, -type4 <string>] [-sheet-name <string>] [-plane <xy or xz or yz>] [-sheet-edge <armchair or zigzag>]')
    print('                         [-length-edge <float>] [-length-perp <float>] [-stacking <AA or AB or ABC>] [-sheet-spacing <float>] [-nlayers <int>]')
    print('                         [-sym-tube-name <string>] [-sym-axis <x or y or z>] [-tube-edge <armchair or zigzag>] [-sym-length <float>]')
    print('                         [-sym-diameter <float>] [-tube-spacing <float>] [-ntubes <int>] [-chi-tube-name <string>] [-n <int>] [-m <int>]')
    print('                         [-chi-axis <x or y or z>] [-chi-length <float>] [-fseed <int>] [-fatoms <string>] [-tatoms <string>] <-gui> <-opt>|<-man>')


    print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the sheet_builder.py')
    print('file will be enforced. The command line inputs option gives more flexibility to the code depending on IDE usage or')
    print('command line usage. The user may specify as many or as few command line options as desired, but if the command line')
    print('option is used and not all options are provided by the user, the hard coded inputs will be enforced, and the user')
    print('will be warned.*')
    
    print('\nTo use command line override options there must be alteration between tags started with the - character and the')
    print('input value. Ordering of tag/input pairs do not matter. Example:')
    print('    python3 sheet_builder.py -tag1 tag1-input  -tag2 tag2-input  -tag3 tag3-input  ...N-arbitrary tagsN/tagN-input pairs')
    
    # print global options
    print('\n\n')
    print('**********************************************')
    print('* Global Options that apply to all run modes *')
    print('**********************************************')
    
    # print -dir option
    print(f'\n -dir or -d <new directory name>   sheet_builder variable: parent_directory    hard coded: {parent_directory}')
    print('    Command line option to set new directory name to store all files atom_typing.py can write. The file')
    print('    management is based on the present working directory (pwd) and will store files one directory deep inside the')
    print('    pwd. The -dir option will specify the directory that is inside of pwd. This will allow for easier file ')
    print('    management and allows the user to merge multiple files for a certain system or research project and keep them')
    print('    all organized into a single directory.')
    print('    python variables*. Example usage:')
    print('        python3 sheet_builder.py -dir EXAMPLE_TEST_SYSTEM_DIRECTORY')
    
    # print -run-mode option
    print(f'\n -run-mode or -run <string>   sheet_builder variable: run_mode    hard coded: {run_mode}')
    print('    Command line option to set run mode to run the code in. Currently supported modes are:')
    print('        sheet          which will generate sp2 hybridized sheets. All settings for this mode can be found in')
    print('                       "mode: sheet" section of this file.')
    print('        symmetric-tube which will generate sp2 hybridized symmetric tubes in either armchair or zigzag')
    print('                       configuration. All settings for this mode can be found in "mode: symmetric-tube"')
    print('                       section of this file.')
    print('        chiral-tube   which will generate sp2 hybridized chiral tubes based.  All settings for this mode can')
    print('                      be found in "mode: chiral-tube" section of this file.')
    print('    Example usage:')
    print('        python3 sheet_builder.py -run-mode  sheet')
    
    # print -find-bonds option
    print(f'\n -find-bonds or -bonds <T|F>   sheet_builder variable: find_bonds    hard coded: {find_bonds}')
    print('    Command line option to control whether or not to find bonds based on interatomic distance searching.')
    print('    The following options are available:')
    print('        T for True, to find bonds ')
    print('        F for False, to not find bonds')
    print('    Example usage:')
    print('        python3 sheet_builder.py -find-bonds T')
    
    # print -pbc-bonds option
    print(f'\n -pbc-bonds or -pbc <T|F>   sheet_builder variable: periodic_bonds    hard coded: {periodic_bonds}')
    print('    Command line option to control whether bonds will pass through the periodic boundary or not.')
    print('    The following options are available:')
    print('        T for True, to find bonds that pass through the periodic boundary')
    print('        F for False, to find bonds that do not pass through the periodic boundary')
    print('    Example usage:')
    print('        python3 sheet_builder.py -pbc-bonds T')
    
    # print -bond-length option
    print(f'\n -bond-length or -r0 <float>   sheet_builder variable: bond_length    hard coded: {bond_length}')
    print('    Command line option to control the in-plane spacing of atoms when generating sheets or tubes based on the')
    print('    bond length between atoms. The bond length is supplied in units of Angstroms (usually around 1.42 A).')
    print('    Example usage:')
    print('        python3 sheet_builder.py -bond-length 1.42')
    
    # print -bond-length option
    print(f'\n -type1, -type2, -type3, -type4 or -t1, t2, t3, t4 <string>   sheet_builder variable: types    hard coded: {types}')
    print('    Command line option to control the "pattern" of atoms for generating sheets or tubes that contain multiple atom types')
    print('    or elements (for example in Boron-Nitride sheets or tubes). The "base unit" that is replicated to generate the sheets')
    print('    or tubes is describe in the figure below:')
    try:
        print()
        print("       type2___type3    replicated   2___3     2___3    types =  {1:'B', 2:'N', 3:'B', 4:'N'}   N___B     N___B    ")
        print("           /   \        --------->  1/   \ 4__1/   \ 4  ------------------------------------>  B/   \ N__B/   \ N  ")
        print("      type1     type4                \___/     \___/    |                                       \___/     \___/    ")
        print("                                     2   3     2   3    |                                       N   B     N   B    ")
        print("                                                        |                                                          ")
        print("                                                        |types = {1:'C', 2:'C', 3:'C', 4:'C'}   C___C     C___C     ")
        print("                                                        +------------------------------------> C/   \ C__C/   \ C   ")
        print("                                                                                                \___/     \___/     ")
        print("                                                                                                C   C     C   C     ")
        print()
        print("    Here are a few examples for different force:")
        print("        PCFF     -> types = {1:'cp',  2:'cp',  3:'cp',  4:'cp'}")
        print("        PCFF-IFF -> types = {1:'cg1', 2:'cg1', 3:'cg1', 4:'cg1'}")
        print("        DREIDING -> types = {1:'C_R', 2:'C_R', 3:'C_R', 4:'C_R'}")
        print("        ReaxFF   -> types = {1:'C',   2:'C',   3:'C',   4:'C'}")
        print("    NOTE: that when generating inputs to ReaxFF that you set the actually element symbol and not classical force field")
        print("    atom types. Also note that when generating inputs to ReaxFF, you should set find_bonds = False  which will allow ")
        print("    you to directly take the output LAMMPS datafile of this code and run in a LAMMPS ReaxFF simulation.")
        print()
        print("    Additionaly IFF's pi-electrons can be added to a system by setting each type as 'atomtype|pi-electron'")
        print("    where the '|' character seperates the atomtype from the pi-electron type. For example:")
        print()
        print("        types = {1:'cg1|cge', 2:'cg1|cge', 3:'cg1|cge', 4:'cg1|cge'} will generate:                             ")
        print("                                                            cge cge  cge cge                                     ")
        print("             2___3     2___3             cg1___cg1           |   |    |   |                                      ")
        print("            1/   \ 4__1/   \ 4   --->   cg1/   \ cg1   ---> cg1-cg1--cg1-cg1                                     ")
        print("             \___/     \___/               \___/             |   |    |   |                                      ")
        print("             2   3     2   3             cg1   cg1          cge cge  cge cge                                     ")
        print()
        print("        types = {1:'bbn|bbe', 2:'nbn|nbe', 3:'bbn|bbe', 4:'nbn|nbe'} will generate:                             ")
        print("                                                           bbe nbe  bbe nbe                                     ")
        print("             2___3     2___3             nbn___bbn          |   |    |   |                                      ")
        print("            1/   \ 4__1/   \ 4   --->   bbn/   \ nbn  ---> bbn-nbn--bbn-nbn                                     ")
        print("             \___/     \___/               \___/            |   |    |   |                                      ")
        print("             2   3     2   3             nbn   bbn         bbe nbe  bbe nbe                                     ")
    except:
        print('    Encountered Error loading "base unit" representation. Perhaps use sheet_builder.py GUI "Quick help" button')
        print('    to visualize the types to "base unit" representaion or open LUNAR/src/sheet_builder/command_line.py to see')
        print('    the representation.')
    print()
    print('    Final note, is that in the "sheet_builder.py" file, there is a charge dictionary, which is used to store the charge')
    print('    for each atom type that may be set in the in type1 - type4 command line tags. Thus the only way to adjust charge')
    print('    is to open the "sheet_builder.py" file and manually adjust the charge dictionary. Example usage:')
    print('        python3 sheet_builder.py -type1 cp -type2 cp -type3 cp -type4 cp')
    print('        python3 sheet_builder.py -type1 cg1|cge -type2 cg1|cge -type3 cg1|cge -type4 cg1|cge')
    
    # print adding functional and terminating atoms
    print('\n\n')
    print('*********************************************************************************************')
    print('* Options for adding different atoms to sheets and tubes (terminating or functionalization) *')
    print('*********************************************************************************************')
    
    # print -mdist option
    print(f'\n -mdist or -md <float>   sheet_builder variable: minimum_distance    hard coded: {minimum_distance}    ')
    print('   Command line option set the minimum distance value to, that is used when adding grafting files or      ')
    print('   functional groups. If the minimum_distance value is set to 0 (zero), no minimum distance constraint    ')
    print('   is applied. *NOTE: this peforms pairwise calculations that are periodic and non-periodic, which can    ')
    print('   be very computationally intensive depending on the size of the sheets. Currently there is no domain    ')
    print('   decomposition implentmented here.*                                                                     ')
    print('                                                                                                          ')
    print('   When the "minimum_distance" constraint is used for "grafting_files" files an additional keyword syntax ')
    print('   is available of "cylinder" or "cylinder:SF", where "cylinder" or "cylinder:SF" envokes an option to    ')
    print('   find out the smallest cylinder that would fit around the grafting fragment. The diameter of this       ')
    print('   cylinder is then used as the "minimum_distance" constraint. The "SF" in "cylinder:SF" is an optional   ')
    print('   keyword that sets a scale factor to multiply the diameter by to set the minimum distance constraint.   ')
    print('   For example say the fit cylinder had a length of 15.0 angstroms and a diameter of 10.0 angstroms. Then ')
    print('   the following is true:                                                                                 ')
    print('     cylinder     -> minimum_distance = 10.0                                                              ')
    print('     cylinder:2   -> minimum_distance = 20.0                                                              ')
    print('     cylinder:0.5 -> minimum_distance = 5.0                                                               ')
    print('                                                                                                          ')
    print('   If you attempt the use the "cylinder" or "cylinder:SF" option for "functional_groups" it will default  ')
    print('   the minimum distance to 0 (zero), as a single "line" of atoms will not have a cylinder to fit around   ')
    print('   them. Example usage:                                                                                   ')
    print('        python3 sheet_builder.py -mdist 6.0')
    
    # print -fseed option
    print(f'\n -fseed or -fs <int>   sheet_builder variable: seed    hard coded: {seed}')
    print('   Command line option to set a seed to the random number generate to define the random atoms the functional')
    print('   groups or grafting files will be added to. will be added to. If the seed value is set to 0 (zero), the ')
    print('   current system time from your computer is used to provide a seed to the random number generator. Example usage:')
    print('        python3 sheet_builder.py -fseed 12345')
    
    # print -fatoms option
    print(f'\n -fatoms or -fa <string>   sheet_builder variable: functional_atoms    hard coded: {functional_atoms}')
    print("   Command line option to set how functional atoms are added to the sheet or tube and what the functional")
    print("   group is. The string format is as follows:                                                               ")
    print('     BondingType<MaxPercent,Direction,molID>|Type1|Type2|TypeN, where "BondingType" is the atom type to add ')
    print('     the function group to, "MaxPercent" is a float or integer type to set the maximum percent of atoms to  ')
    print('     functionalize, "Direction" is the direction to stick the functional group, "molID" is the molecule     ')
    print('     identifier to attached the functional group to and the "|" character seperates types, and the "TypeN"  ')
    print('     sets the atom to add.                                                                                  ')
    print("                                                                                                            ")
    print('     The "Direction" character can be "+" or "-" with the following meanings:                               ')
    print('       Sheets:                                                                                              ')
    print('         + means point the functional groups in the positive direction                                      ')
    print('         - means point the functional groups in the positive direction                                      ')
    print('       Tubes:                                                                                               ')
    print('         + means point the functional groups point outwards relative to the circular cross-section (points  ')
    print('           outwards to the circumference)                                                                   ')
    print('         - means point the functional groups point inwards from the to the circular cross-section (points   ')
    print('           inwards to the center)                                                                           ')
    print("                                                                                                            ")
    print('     The "molID" is an integer value starting from 1 and going to nlayers for sheets and ntubes for tubes.  ')
    print('     The indexing of molIDs is as follows:                                                                  ')
    print('       Sheets: molID of 1 is the sheet that is the most negative in the normal direction (can be thought of ')
    print('               as the bottom sheet). Then each sequential sheet that is placed above have their molIDs      ')
    print('               incremented. For example a 3 layer stack molIDs would be:                                    ')
    print('                 1 -> bottom sheet                                                                          ')
    print('                 2 -> middle sheet                                                                          ')
    print('                 3 -> top sheet                                                                             ')
    print("                                                                                                            ")
    print('       Tubes:  molID of 1 is the center tube, then each sequential tube that is built outward have their    ')
    print('               molIDs incremented. For example a 3 concentric generation of tubes would have molIDs:        ')
    print('                 1 -> center tube                                                                           ')
    print('                 2 -> middle tube                                                                           ')
    print('                 3 -> outer most tube                                                                       ')
    print("                                                                                                            ")
    print('     The "*" character acts as a wildcard character to define which molIDs to add the functional groups to, ')
    print('     where setting the molID to "*", will randomly place the functional groups and any of the molIDs sheets ')
    print('     or tubes that have been built.                                                                         ')
    print('                                                                                                            ')
    print("   For example say types = {1:'C', 2:'C', 3:'C', 4:'C'} to generate a 3 layer set of carbon sheets          ")
    print('   (pointing the functional groups in the positive direction) or 3 concentric nanotubes (pointing the       ')
    print('   funcational groups in the outward direction) and the goal was to functionalize 5% of the carbon atoms    ')
    print("   with -OH functional group. Then the functional_atoms string would be 'C<5,+,3>|O|H', which would         ")
    print('   randomly add the -OH functional group to 5% of the C atoms to the molID of 3 (for sheets that is the top ')
    print('   sheet and for tubes that is the outer most tube)                                                         ')
    print("                                                                                                            ")
    print('   Additionaly the functional_atoms string can handle multiple BondingTypes by seperating them with the     ')
    print('   ";" character. So the generalized functional_atoms string becomes:                                       ')
    print('   TypeA<MaxPercentA,dirA,molidA>|TypeA1|TypeAN; TypeB<MaxPercentB,dirB,molidb>|TypeB1|TypeBN; ...          ')
    print('                                                                                                            ')
    print("   For example say types = {1:'B', 2:'N', 3:'B', 4:'N'} to generate a Boron-Nitride sheet or tube with      ")
    print('   alternating B/N atoms and the goal was to functionalize 10% of the Boron atoms with -OH functional       ')
    print('   group and to functionalize 20% of the Nitride atoms with -H functional group. Then the functional_atom   ')
    print("   string would be 'B<10,+,*>|O|H; N<20,+*>|H', which would randomly add the -OH functional group to 10%    ")
    print('   of the B atoms and add the -H functional group to 20% of the N atoms.                                    ')
    print("                                                                                                            ")
    print('   All examples above will place the atoms in a line along the orthagonal direction from the surface of     ')
    print('   tube, but say we wanted to added a functional group that resembles an epoxide ring (3 member ring with   ')
    print('   two carbons and 1 oxygen). Then we can add a "|" character to the end of the functional_atoms string.    ')
    print('   This method currently only works for adding a single atom functional group like oxygen to the sheets or  ')
    print('   tubes.                                                                                                   ')
    print('                                                                                                            ')
    print("     For example say types = {1:'C', 2:'C', 3:'C', 4:'C'} to generate a single carbon sheet or single       ")
    print('     nanotube and the goal was to functionalize 30% of the carbon atoms with the epoxide ring oxygen (that  ')
    print('     point in the negative direction of the sheet and inwards for the tube). Then the functional_atoms      ')
    print("     string would be 'C<30,-,*>|O|', where the last character is the '|' character. This will tell the code ")
    print('     to find a first neighbor from the random atom and center the oxygen atom between the first neighbor    ')
    print('     and itself. Finally, add two bonds to create the epoxide type ring. Note that each time the oxygen     ')
    print('     atom is added, it functionalizes two carbon atoms at a time. So say the sheet or tube had 100 carbon   ')
    print('     atoms and the functionalization MaxPercent was set to 30%, then only 15 oxygen atoms will be added     ')
    print('     (not 30).                                                                                              ')
    print('                                                                                                            ')
    print('   This option is also currently limited to only being able to add one kind of functional group to each     ')
    print('   type (1-4). So say you wanted to model a graphene sheet of all carbon atoms using the PCFF atom type     ')
    print('   "cp", but wanted to model functional groups of -O- and -OH. A work around is to use all2lmp.py atom      ')
    print('   naming scheme summarized as: "AtomType:Name", where the ":"-character provides a delimiter from the true ')
    print('   atom type and a name a user can supply. You can uniquely name each of the types (1-4) with a ":name".    ')
    print('                                                                                                            ')
    print("     For example you could set types {1:'cp:line', 2:'cp:line', 3:'cp:ring', 4:'cp:ring'} and then set the  ")
    print("     functionalization string as 'cp:line<30,+,*>|O|H; cp:ring<10,+*>|O|'. This will trick sheet_builder    ")
    print("     into recognizing the atoms as two different types 'cp:line' and 'cp:ring' which will allow you to add  ")
    print('     two different functional groups to the sheet, while maintaining an atom type all2lmp.py can recognize  ')
    print('     and automatically parameterize.		                                                                   ')
    print("                                                                                                            ")
    print("   If the functional_atoms entry/string is left blank, this option will not be envoked. Additionally,       ")
    print("   this option requires find_bonds to be True.                                                              ")
    print("                                                                                                            ")
    print('   Example usage:')
    print('    python3 sheet_builder.py -fatoms "C<30,+,*>|O|"')
    
    
    # print -gfiles option
    print(f'\n -gfiles or -gf <string>   sheet_builder variable: grafting_files    hard coded: {grafting_files}       ')
    print("   Command line option to set how grafting files are added to the sheet or tube and what the atomID(s) are ")
    print("   used from the grafting file as the atomID(s) to bond to the sheet or tube. The string format is as      ")
    print("   follows:                                                                                                ")
    print('    An entry to supply a string to set how grafting files are added to the sheet or tube and what the      ')
    print('    atomID(s) are used from the grafting file as the atomID(s) to bond to the sheet or tube. The string    ')
    print('    format is as follows:                                                                                  ')
    print('      BondingType<MaxPercent,Direction,molID><id1, id2>|FILENAME.EXT                                       ')
    print('    The "BondingType<MaxPercent,Direction,molID>" portion descirbing the "BondingType", "MaxPercent",      ')
    print('    "Direciton", and "molID" are the same as the functional_atoms and are described in that section of     ')
    print('    this page. Please refer to that section for a more in-depth understanding. FILENAME.EXT is a file that ')
    print('    defines the new atoms and bonds to graft onto the sheets or tubes. The following per-atom attributes   ')
    print('    are used from each file extension to set the atom type of the grafted atoms:                           ')
    print('      .pdb  the per-atom "atom_name" information from this file is used to set the atom type               ')
    print('      .mol  the per-atom "element" information from this file is used to set the atom type                 ')
    print('      .sdf  the per-atom "element" information from this file is used to set the atom type                 ')
    print('      .mol2 the per-atom "element" information from this file is used to set the atom type                 ')
    print('      .data the per-atom "type label" is the first attempt at being used, if that fails, the per-atom      ')
    print('            "comment from the Masses" section is the second attempt, and if that fails, the numeric LAMMPS ')
    print('            atomTypeID is used for the information from this file is used to set the atom typ .            ')
    print('                                                                                                           ')
    print('    The "<id1, id2>" portion allows users to set which atomID(s) are used to bond the grafting fragment to ')
    print('    the sheets by. You may supply one or  two atomIDs with the following meanings:                         ')
    print('      <id1>      will bond the grafting fragment via one covalent bond between the sheet or tube and the   ')
    print('                 grafting fragment                                                                         ')
    print('      <id1, id2> will bond the via two covalent bonds between the sheet or tube and the grafting fragment, ')
    print('                 making a ring between them (generating properly geometric rings like this is difficult -  ')
    print('                 therefore this method may not produce nice geometries and it is recommended to use a      ')
    print('                 "fix nve/limit 0.01" run to intialize a these system in LAMMPS).                          ')
    print("                                                                                                           ")
    print('   Example usage:')
    print('    python3 sheet_builder.py -fatoms "C<7,+,*><15>|EXAMPLES/sheet_builder/PBZ_graft.15.pdb"')
    
    # print -tatoms option
    print(f'\n -tatoms or -ta <string>   sheet_builder variable: terminating_atoms    hard coded: {terminating_atoms}  ')
    print("   Command line option to set how terminating atoms are added to the sheet or tube and what the termanting  ")
    print("   atoms are. This option requires that periodic_bonds is False, as this creates open valences on the       ")
    print('   "end" atoms of the sheet or tube. The string format is as follows:                                       ')
    print('     BondingType|Type1|Type2|TypeN, where "BondingType" is the atom type to add the terminating atoms to,   ')
    print('     the "|" character seperates types, and the "TypeN" sets the atom to add.                               ')
    print("                                                                                                            ")
    print("     For example say types = {1:'C', 2:'C', 3:'C', 4:'C'} to generate a carbon sheet or nanotube that is    ")
    print('     not periodically bonded and the goal was to terminate the open valences on the "edges" of the sheet or ')
    print("     tube with the -OH functional group. Then the termanting_atoms string would be 'C|O|H', which would     ")
    print('     terminate all "edge" atoms with the -OH group.                                                         ')
    print("                                                                                                            ")
    print("     Additionaly the terminating_atoms string can handle multiple BondingType's by seperating them with the ")
    print('     ";" character. So the generalized termanting_atoms string becomes:                                     ')
    print("       BondingTypeA|TypeA1|TypeAN; BondingTypeB|TypeB1|TypeBN; BondingTypeC|TypeC1|TypeCN; ...              ")
    print("                                                                                                            ")
    print("       For example say types = {1:'B', 2:'N', 3:'B', 4:'N'} to generate a Boron-Nitride sheet or tube with  ")
    print("       alternating B/N atoms that was not periodically bonded and the goal was to termanate the Boron atoms ")
    print("       with -H and to termanate the Nitride atoms with -OH, then the terminating_atoms string would be      ")
    print('       "B|H; N|O|H", which would termanate the "edge" Boron atoms with an -H and the "edge" Nitride atoms   ')
    print("       with a -OH group.                                                                                    ")
    print("                                                                                                            ")
    print("   If the terminating_atoms entry/string is left blank, this option will not be envoked. Additionally, this ")
    print("   option requires find_bonds to be True.                                                                   ")
    print('  Example usage:')
    print('    python3 sheet_builder.py -tatoms "C|O|H"')
    
    
    # print sheet options
    print('\n\n')
    print('*****************')
    print('* Sheet Options *')
    print('*****************')
    
    # print -sheet-name option
    print(f'\n -sheet-name or -sname <string>   sheet_builder variable: sheet_basename    hard coded: {sheet_basename}')
    print("    Command line option to set output basename of outputs for when running in 'sheet' mode. Example usage:")
    print('        python3 sheet_builder.py -run-mode sheet -sheet-name graphene')
    
    # print -plane option
    print(f'\n -plane or -p <xy or xz or yz>   sheet_builder variable: plane    hard coded: {plane}')
    print("    Command line option to set the plane to generate the sheet on. The available planes to generate sheets on are:")
    print("        'xy' which will generate a sheet on the XY plane")
    print("        'xz' which will generate a sheet on the XZ plane")
    print("        'yz' which will generate a sheet on the YZ plane")
    print("    The sheet size and shape will depend on the length_in_normal/length_in_edgetype and sheet_edgetype. Example usage:")
    print('        python3 sheet_builder.py -run-mode sheet -plane xy')
    
    # print -sheet-edge option
    print(f'\n -sheet-edge or -sedge <armchair or zigzag>   sheet_builder variable: sheet_edgetype    hard coded: {sheet_edgetype}')
    print("    Command line option to set the edge type of the sheet. The available edge types are:")
    print("        'armchair' which will generate a sheet in the armchair directions, with the length in the armchair")
    print("                   direction set by the length_in_edgetype variable.")
    print("        'zigzag'   which will generate a sheet in the zigzag directions, with the length in the zigzag")
    print("                   direction set by the length_in_edgetype variable.")
    print('    Example usage:')
    print('        python3 sheet_builder.py -run-mode sheet -sheet-edge armchair')
    
    # print -length-edge option
    print(f'\n -length-edge or -le <float>   sheet_builder variable: length_in_edgetype    hard coded: {length_in_edgetype}')
    print('    Command line option to control the size of the sheet. The length_in_edgetype will set the length in the')
    print('    edge type direction set by the sheet_edgetype variable. Example usage:')
    print('        python3 sheet_builder.py -run-mode sheet -length-edge 25')
    
    # print -length-perp option
    print(f'\n -length-perp or -lp <float>   sheet_builder variable: length_in_perpendicular    hard coded: {length_in_perpendicular}')
    print('    Command line option to control the size of the sheet. The length_in_perpendicular will set the length in the')
    print('    perpendicular direction from the edge type. Example usage:')
    print('        python3 sheet_builder.py -run-mode sheet -length-perp 15')
    
    # print -stacking option
    print(f'\n -stacking or -s <AA or AB or ABC>   sheet_builder variable: stacking    hard coded: {stacking}')
    print('    Command line option to set how sheets are stacked if sheet_nlayers is greater than one. The following stacking sequences are available:')
    print("        'AA'")
    print("        'AB'")
    print("        'ABC'")
    print("    where the schematic below shows AA, AB, and ABC stacking for six sheets:")
    print("        AA-stacking        AB-stacking        ABC-stacking                                                    ")
    print("        ___________        ___________        ___________                                                     ")
    print("        ___________           ___________         ___________                                                 ")
    print("        ___________        ___________                ___________                                             ")
    print("        ___________        ___________        ___________                                                     ")
    print("        ___________           ___________         ___________                                                 ")
    print("        ___________        ___________                ___________                                             ")
    print('    Example usage:')
    print('        python3 sheet_builder.py -run-mode sheet -stacking AB')
    
    # print -sheet-spacing option
    print(f'\n -sheet-spacing or -ss <float>   sheet_builder variable: sheet_layer_spacing    hard coded: {sheet_layer_spacing}')
    print('    Command line option to control the out-of-plane spacing of atoms when multiple layers are added. The layer spacing is supplied in units')
    print('    of Angstroms (usually around 3.354 A). Example usage:')
    print('        python3 sheet_builder.py -run-mode sheet -sheet-spacing 3.4')
    
    # print -sheet-spacing option
    print(f'\n -nlayers or -nl <int>   sheet_builder variable: sheet_nlayers    hard coded: {sheet_nlayers}')
    print('    Command line option to control the number of layers to generated. All models will be centered about (0, 0, 0) no matter the number of ')
    print('    layers generated. Example usage:')
    print('        python3 sheet_builder.py -run-mode sheet -nlayers 5')
    
    
    # print symmetric-tube options
    print('\n\n')
    print('**************************')
    print('* Symmetric-tube Options *')
    print('**************************')
    
    # print -sym-tube-name option
    print(f'\n -sym-tube-name or -stname <string>   sheet_builder variable: symmetric_tube_basename    hard coded: {symmetric_tube_basename}')
    print("    Command line option to set output basename of outputs for when running in 'symmetric-tube' mode. Example usage:")
    print('        python3 sheet_builder.py -run-mode symmetric-tube -sym-tube-name CNT')
    
    # print -sym-axis option
    print(f'\n -sym-axis or -sa <x or y or z>   sheet_builder variable: symmetric_tube_axis    hard coded: {symmetric_tube_axis}')
    print("    Command line option to set which axis the tube axis is aligned with. Example usage:")
    print('        python3 sheet_builder.py -run-mode symmetric-tube -sym-axis z')
    
    # print -tube-edge option
    print(f'\n -tube-edge or -tedge <armchair or zigzag>   sheet_builder variable: tube_edgetype    hard coded: {tube_edgetype}')
    print("    Command line option to set the edge type of the tube. Example usage:")
    print('        python3 sheet_builder.py -run-mode symmetric-tube -tube-edge zigzag')
    
    # print -sym-length option
    print(f'\n -sym-length or -sl <float>   sheet_builder variable: symmetric_length    hard coded: {symmetric_length}')
    print("    Command line option to control the length of the tube in the axial direction. Example usage:")
    print('        python3 sheet_builder.py -run-mode symmetric-tube -sym-length 15')
    
    # print -sym-length option
    print(f'\n -sym-diameter or -sd <float>   sheet_builder variable: diameter    hard coded: {diameter}')
    print("    Command line option to control the diameter of the tube. Example usage:")
    print('        python3 sheet_builder.py -run-mode symmetric-tube -sym-diameter 10')
    
    # print -tube-spacing option
    print(f'\n -tube-spacing or -ts <float>   sheet_builder variable: tube_layer_spacing    hard coded: {tube_layer_spacing}')
    print("    Command line option to control the out-of-plane spacing of atoms when multiple tubes are added. (i.e. when symmetric_ntubes is")
    print("    greater then 1). The layer spacing is supplied in units of Angstroms (usually around 3.354 A). Example usage:")
    print('        python3 sheet_builder.py -run-mode symmetric-tube -tube-spacing 3.4')
    
    # print -ntubes option
    print(f'\n -ntubes or -nt <int>   sheet_builder variable: symmetric_ntubes    hard coded: {symmetric_ntubes}')
    print("    Command line option to control the number of tubes to generate. All models will be centered about (0, 0, 0)")
    print("    no matter the number of tubes generated. When symmetric_ntubes is greater then 1, the diameter variable")
    print("    will control the inner most tube diameter and the N-outer layer tube diameters will be controlled by the")
    print("    symmetric_ntubes and tube_layer_spacing varaibles. Example usage:")
    print('        python3 sheet_builder.py -run-mode symmetric-tube -ntubes 3')
    
    
    # print chiral-tube options
    print('\n\n')
    print('***********************')
    print('* Chiral-tube Options *')
    print('***********************')
    
    # print -chi-tube-name option
    print(f'\n -chi-tube-name or -ctname <string>   sheet_builder variable: chiral_tube_basename    hard coded: {chiral_tube_basename}')
    print("    Command line option to set output basename of outputs for when running in 'chiral-tube' mode. Example usage:")
    print('        python3 sheet_builder.py -run-mode chiral-tube -chi-tube-name CNT')
    
    # print -sym-axis option
    print(f'\n -chi-axis or -ca <x or y or z>   sheet_builder variable: chiral_tube_axis    hard coded: {chiral_tube_axis}')
    print("    Command line option to set which axis the tube axis is aligned with. Example usage:")
    print('        python3 sheet_builder.py -run-mode chiral-tube -chi-axis z')
    
    # print -chi-length option
    print(f'\n -chi-length or -cl <float>   sheet_builder variable: chiral_length    hard coded: {chiral_length}')
    print('    Command line option to control the length of a chiral tube. NOTE the chiral indices m and n control the')
    print('    "base-tube" length where the chiral_length sets the desired length of a chiral tube, where the "base-tube"')
    print('    will be replicated enough times to get a tube close in length to chiral_length. Example usage:')
    print('        python3 sheet_builder.py -run-mode chiral-tube -chi-length 15')
    
    # print -n option
    print(f'\n -n or -n <int>   sheet_builder variable: n    hard coded: {n}')
    print("    Command line option to control the chiral index 'n'. Example usage:")
    print('        python3 sheet_builder.py -run-mode chiral-tube -n 5')
    
    # print -m option
    print(f'\n -m or -m <int>   sheet_builder variable: m    hard coded: {m}')
    print("    Command line option to control the chiral index 'm'. Example usage:")
    print('        python3 sheet_builder.py -run-mode chiral-tube -m 10')
    
    # print -m and -n info
    print('\n -n and -m additional info')
    print("    Symmetric tubes can be generated when:")
    print("        m = n, which generates an armchair tube")
    print("        m = 0, which generates a zigzag tube")
    print("    If it is desired to create armchair or zigzag tubes this code can be ran in mode = 'symmetric-tube', which")
    print("    provides easier control over the tube length and diameter. The 'symmetric-tube' mode also allows multiwall")
    print("    nanotubes to be generated with easy. Thus it is recommend only to use the 'chiral-tube' mode if you want")
    print("    to generate a chiral tube and use the 'symmetric-tube' to generate tubes in either armchair or zigzag")
    print("    configurations.")
    
    
    # print -opt or -man option
    print('\n\n')
    print('*************************')
    print('* Miscellaneous Options *')
    print('*************************')
    print('\n -opt or -man')
    print('    Will tell bond_react_merge_prep to only print out avaiable command line options known as tagN and tagN-inputs.')
    print('    This is the only tag that doesnt require a tagN-input following the tag since the code will only look for if')
    print('     -opt is in command line input list. Example usage:')
    print('        python3 sheet_builder.py -man')
        
    # print -gui option
    print('\n -gui <no addition info required> bond_react_merge_prep variable: use_GUI')
    print('    Command line option flag to load GUI. Example usage:')
    print('        python3 sheet_builder.py -gui')
    
    # print requirements and dependencies
    print('\n\nRequirements and Dependencies:')
    print('    Requirements: python 3.7+ to guarantee ordering of dictionaries which stores most of the data')
    print('    Dependencies: numpy if adding pi-electrons')
    
    # print buffer and exit code
    print('\n\n')
    return

# Class to update inputs
class inputs:
    def __init__(self, commandline_inputs, sheet_basename, symmetric_tube_basename, chiral_tube_basename, run_mode, parent_directory, length_in_perpendicular, length_in_edgetype,
                 sheet_edgetype, types, bond_length, sheet_layer_spacing, sheet_nlayers, stacking, plane, tube_edgetype, tube_layer_spacing, symmetric_ntubes, symmetric_length,
                 diameter, n, m, chiral_length, symmetric_tube_axis, chiral_tube_axis, find_bonds, periodic_bonds, seed, functional_atoms, terminating_atoms, grafting_files, minimum_distance):
        
        # Give access to inputs (update later on if command line over ride is given)
        self.commandline_inputs = commandline_inputs
        self.parent_directory = parent_directory
        self.run_mode = run_mode
        self.find_bonds = find_bonds
        self.periodic_bonds = periodic_bonds
        self.bond_length = bond_length
        self.types = types
        self.sheet_basename = sheet_basename
        self.plane = plane
        self.sheet_edgetype = sheet_edgetype
        self.length_in_perpendicular = length_in_perpendicular
        self.length_in_edgetype = length_in_edgetype
        self.stacking = stacking
        self.sheet_layer_spacing = sheet_layer_spacing
        self.sheet_nlayers = sheet_nlayers
        self.symmetric_tube_basename = symmetric_tube_basename
        self.symmetric_tube_axis = symmetric_tube_axis
        self.tube_edgetype = tube_edgetype
        self.symmetric_length = symmetric_length
        self.diameter = diameter
        self.tube_layer_spacing = tube_layer_spacing
        self.symmetric_ntubes = symmetric_ntubes
        self.chiral_tube_basename = chiral_tube_basename
        self.chiral_tube_axis = chiral_tube_axis
        self.chiral_length = chiral_length
        self.n = n
        self.m = m
        
        self.seed = seed     # '-fseed' or '-fs'
        self.minimum_distance = minimum_distance # '-mdist' or '-md'
        self.grafting_files = grafting_files   # '-gfiles' or '-gf'
        self.functional_atoms = functional_atoms   # '-fatoms' or '-fa'
        self.terminating_atoms = terminating_atoms # '-tatoms' or '-ta'
        
        
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
        supported_tags = ['-dir', '-run-mode', '-find-bonds', '-pbc-bonds', '-bond-length', '-type1', '-type2', '-type3', '-type4',
                          '-sheet-name', '-plane', '-sheet-edge', '-length-edge', '-length-perp', '-stacking', '-sheet-spacing', '-nlayers',
                          '-sym-tube-name', '-sym-axis', '-tube-edge', '-sym-length', '-sym-diameter', '-tube-spacing', '-ntubes',
                          '-chi-tube-name', '-chi-axis', '-chi-length', '-n', '-m', '-fseed', '-fatoms', '-tatoms','-gfiles', '-mdist']
        
        # set shortcut_tags mapping
        shortcut_tags = {'-d':'-dir', '-run':'-run-mode', '-bonds':'-find-bonds', '-pbc':'-pbc-bonds', '-r0':'-bond-length', '-t1':'-type1', '-t2':'-type2', '-t3':'-type3', '-t4':'-type4',
                         '-sname':'-sheet-name', '-p':'-plane', '-sedge':'-sheet-edge', '-le':'-length-edge', '-lp':'-length-perp', '-s':'-stacking', '-ss':'-sheet-spacing', '-nl':'-nlayers',
                         '-stname':'-sym-tube-name', '-sa':'-sym-axis', '-tedge':'-tube-edge', '-sl':'-sym-length','-sd':'-sym-diameter', '-ts':'-tube-spacing', '-nt':'-ntubes',
                         '-ctname':'-chi-tube-name', '-ca':'-chi-axis', '-cl':'-chi-length', '-n':'-n', '-m':'-m', '-fs':'-fseed', '-fa':'-fatoms', '-ta':'-tatoms',
                         '-gf':'-gfiles', '-md':'-mdist'}
        
        # set default variables
        default_variables ={'-dir':self.parent_directory, '-run-mode':self.run_mode, '-find-bonds':self.find_bonds, '-pbc-bonds':self.periodic_bonds,
                            '-bond-length':self.bond_length, '-type1':self.types[1], '-type2':self.types[2], '-type3':self.types[3], '-type4':self.types[4],
                            '-sheet-name':self.sheet_basename, '-plane':self.plane, '-sheet-edge':self.sheet_edgetype, '-length-edge':self.length_in_edgetype,
                            '-length-perp':self.length_in_perpendicular, '-stacking':self.stacking, '-sheet-spacing':self.sheet_layer_spacing, '-nlayers':self.sheet_nlayers,
                            '-sym-tube-name':self.symmetric_tube_basename, '-sym-axis':self.symmetric_tube_axis, '-tube-edge':self.tube_edgetype,
                            '-sym-length':self.symmetric_length, '-sym-diameter':self.diameter, '-tube-spacing':self.tube_layer_spacing, '-ntubes':self.symmetric_ntubes,
                            '-chi-tube-name':self.chiral_tube_basename, '-chi-axis':self.chiral_tube_axis, '-chi-length':self.chiral_length, '-n':self.n,'-m':self.m,
                            '-fseed':self.seed, '-fatoms':self.functional_atoms, '-tatoms':self.terminating_atoms,
                            '-gfiles':self.grafting_files, '-mdist':self.minimum_distance}
        
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
        #----------------#
        # Global options #
        #----------------#         
        # set new -dir option and print confirmation
        if tags['-dir']:
            if tags['-dir'].isspace() or tags['-dir'] == '':
                print('-dir flag was used, with an empty directory string. The directory string must contain alpha-numeric characters')
            else:
                self.parent_directory = tags['-dir']
                print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-dir', self.parent_directory))
                
        # set new -run-mode option and print confirmation
        if tags['-run-mode']:
            self.run_mode = tags['-run-mode']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-run-mode', self.run_mode))
                 
        # set new -find-bonds option and print confirmation
        if tags['-find-bonds']:
            self.find_bonds = T_F_string2boolean('-find-bonds', (tags['-find-bonds']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-find-bonds', self.find_bonds))
        
        # set new -pbc-bonds option and print confirmation
        if tags['-pbc-bonds']:
            self.periodic_bonds = T_F_string2boolean('-pbc-bonds', (tags['-pbc-bonds']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-pbc-bonds', self.periodic_bonds))
            
        # set new -bond-length option and print confirmation
        if tags['-bond-length']:
            self.bond_length = float(tags['-bond-length'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-bond-length', self.bond_length))
        
        # set new -type1 option and print confirmation
        if tags['-type1']:
            self.types[1] = tags['-type1']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-type1', self.types[1]))
            
        # set new -type2 option and print confirmation
        if tags['-type2']:
            self.types[2] = tags['-type2']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-type2', self.types[2]))
            
        # set new -type3 option and print confirmation
        if tags['-type3']:
            self.types[3] = tags['-type3']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-type3', self.types[3]))
            
        # set new -type4 option and print confirmation
        if tags['-type4']:
            self.types[4] = tags['-type4']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-type4', self.types[4]))
            
        #-------------------------------------------#
        # Termination and functionalization options #
        #-------------------------------------------#
        # set new -fseed option and print confirmation
        if tags['-fseed']:
            self.functional_seed = int(tags['-fseed'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-fseed', self.seed))
            
        # set new -mdist option and print confirmation
        if tags['-mdist']:
            try: self.minimum_distance = float(tags['-mdist'])
            except: self.minimum_distance = tags['-mdist']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-mdist', self.minimum_distance))
            
        # set new -fatoms option and print confirmation
        if tags['-fatoms']:
            self.functional_atoms = tags['-fatoms']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-fatoms', self.functional_atoms))
            
        # set new -gfiles option and print confirmation
        if tags['-gfiles']:
            self.grafting_files = tags['-gfiles']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-gfiles', self.grafting_files))
            
        # set new -tatoms option and print confirmation
        if tags['-tatoms']:
            self.terminating_atoms = tags['-tatoms']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-tatoms', self.terminating_atoms))
            
        #---------------#
        # Sheet options #
        #---------------#   
        # set new -sheet-name option and print confirmation
        if tags['-sheet-name']:
            self.sheet_basename = tags['-sheet-name']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-sheet-name', self.sheet_basename))
            
        # set new -plane option and print confirmation
        if tags['-plane']:
            self.plane = tags['-plane']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-plane', self.plane))
            
        # set new -sheet-edge option and print confirmation
        if tags['-sheet-edge']:
            self.sheet_edgetype = tags['-sheet-edge']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-sheet-edge', self.sheet_edgetype))
            
        # set new -length-edge option and print confirmation
        if tags['-length-edge']:
            self.length_in_edgetype = float(tags['-length-edge'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-length-edge', self.length_in_edgetype))
            
        # set new -length-perp option and print confirmation
        if tags['-length-perp']:
            self.length_in_perpendicular = float(tags['-length-perp'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-length-perp', self.length_in_perpendicular))
            
        # set new -stacking option and print confirmation
        if tags['-stacking']:
            self.stacking = tags['-stacking']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-stacking', self.stacking))
            
        # set new -sheet-spacing option and print confirmation
        if tags['-sheet-spacing']:
            self.sheet_layer_spacing = float(tags['-sheet-spacing'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-sheet-spacing', self.sheet_layer_spacing))
            
        # set new -nlayers option and print confirmation
        if tags['-nlayers']:
            self.sheet_nlayers = int(tags['-nlayers'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-nlayers', self.sheet_nlayers))
            
        #------------------------#
        # Symmetric-tube options #
        #------------------------#   
        # set new -sym-tube-name option and print confirmation
        if tags['-sym-tube-name']:
            self.symmetric_tube_basename = tags['-sym-tube-name']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-sym-tube-name', self.symmetric_tube_basename))
            
        # set new -sym-axis option and print confirmation
        if tags['-sym-axis']:
            self.symmetric_tube_axis = tags['-sym-axis']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-sym-axis', self.symmetric_tube_axis))
            
        # set new -tube-edge option and print confirmation
        if tags['-tube-edge']:
            self.tube_edgetype = tags['-tube-edge']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-tube-edge', self.tube_edgetype))
            
        # set new -sym-length option and print confirmation
        if tags['-sym-length']:
            self.symmetric_length = float(tags['-sym-length'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-sym-length', self.symmetric_length))
            
        # set new -sym-diameter option and print confirmation
        if tags['-sym-diameter']:
            self.diameter = float(tags['-sym-diameter'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-sym-diameter', self.diameter))
            
        # set new -tube-spacing option and print confirmation
        if tags['-tube-spacing']:
            self.tube_layer_spacing = float(tags['-tube-spacing'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-tube-spacing', self.tube_layer_spacing))
            
        # set new -ntubes option and print confirmation
        if tags['-ntubes']:
            self.symmetric_ntubes = int(tags['-ntubes'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-ntubes', self.symmetric_ntubes))
            
        #---------------------#
        # Chiral-tube options #
        #---------------------# 
        # set new -chi-tube-name option and print confirmation
        if tags['-chi-tube-name']:
            self.chiral_tube_basename = tags['-chi-tube-name']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-chi-tube-name', self.chiral_tube_basename))
            
        # set new -chi-axis option and print confirmation
        if tags['-chi-axis']:
            self.chiral_tube_axis = tags['-chi-axis']
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-chi-axis', self.chiral_tube_axis))
            
        # set new -chi-length option and print confirmation
        if tags['-chi-length']:
            self.chiral_length = float(tags['-chi-length'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-chi-length', self.chiral_length))

        # set new -n option and print confirmation
        if tags['-n']:
            self.n = int(tags['-n'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-n', self.n))
            
        # set new -m option and print confirmation
        if tags['-m']:
            self.m = int(tags['-m'])
            print('Override confirmation for {:<18} Hard coded input is being overridden with this input: {}'.format('-m', self.m))
        
        # print buffer
        print('\n\n')