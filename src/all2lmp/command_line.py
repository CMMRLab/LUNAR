#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.5
November 13th, 2024
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


def print_man_page(topofile, nta_file, frc_file, assumed, parent_dir, newfile, atom_style, ff_class, use_auto_equivalence, use_assumed_auto_fill, reset_molids,
                  reset_charges, write_txt_comments, write_bond_react, print_options, use_morse_bonds, include_type_labels, add2box, ignore_missing_parameters,
                  shift, rotate):

    
    # print general command line options
    print('\n\nall2lmp has been run with -opt or -man option to show the optional command line overrides available. Command line')
    print('option summary [-tag <tag-input>]:')
    print('python3 all2lmp.py [-topo <topo-filename>] [-nta <nta-filename>] [-frc <frc-filename>] [-asm <assumed-filename>]')
    print('                   [-dir <new directory name>] [-newfile <string>] [-atomstyle <atomstyle>] [-class <0|1|2|d|s1|s2|r>]')
    print('                   [-auto-equivs <T|F>] [-assumed-equivs <T|F>] [-reset-molids <T|F>] [-reset-charges <T|F>] <-gui>')
    print('                   [-write-comments <T|F>] [-write-bond-react <T|F>] [-morse-bond <T|F>] [-type-labels <T|F>] [-ignore <T|F>] ')
    print('                   [-add2box <float or int>] [-rx or -ry or -rz <float>] [-sx or -sy or -sz <float>] <-opt>|<-man>')

    print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the all2lmp.py')
    print('file will be enforced. The command line inputs option gives more flexibility to the code depending on IDE usage or')
    print('command line usage. The user may specify as many or as few command line options as desired, but if the command line')
    print('option is used and not all options are provided by the user, the hard coded inputs will be enforced, and the user')
    print('will be warned.*')
    
    print('\nTo use command line override options there must be alteration between tags started with the - character and the')
    print('input value. Ordering of tag/input pairs do not matter. Example:')
    print('    python3 all2lmp.py -tag1 tag1-input  -tag2 tag2-input  -tag3 tag3-input  ...N-arbitrary tagsN/tagN-input pairs')
    
    # print options
    print('\n\nOptions:')
    
    # print -topo option
    print(f'\n -topo or -t <topo-filename>  all2lmp variable: topofile    hard coded: {topofile}')
    print('    Command line option to set topofile name of topo/coord file to be read into all2lmp for conversion. Currently')
    print('    supported topofile formats are:')
    print('         .mol = .mol file in v2000 format (Info: atom positions, bonds | Set: center about (0,0,0), box set by')
    print('                 min/max atom extents, image flags zeroed)')
    print('         .sdf = .sdf file in v2000 format (Info: atom positions, bonds | Set: center about (0,0,0), box set by')
    print('                 min/max atom extents, image flags zeroed)')
    print('         .mol2 = .mol2 SYBYL file  (Info: atom positions, bonds | Set: center about (0,0,0), box set by')
    print('                 min/max atom extents, image flags zeroed)')
    print('         .pdb =  .pdb protein data bank file  (Info: atom positions, bonds | Set: center about (0,0,0), box set by')
    print('                 min/max atom extents, image flags zeroed)')
    print('         .data = .data or dat file from LAMMPS (Info: atom positions, bonds, no center, box set by LAMMPS, image')
    print('                 flags set by LAMMPS)')
    print('         .mdf = .mdf file from material studio (Info: atom positions, bonds, box set by file | Set: center about')
    print('                 (0,0,0),  atom wrapped, and image flags updated)')
    print('    See Materials Studio software below for msi2lmp conversion. The files will be determined by file extension, thus')
    print('    it is important files provided have the correct extension. When using VMD MDL MOL2 file you MUST save w/ the ')
    print('    the .mol2 file extension so all2lmp knows how to read/parse the file. Example usage:')
    print('        python3 all2lmp.py -topo EXAMPLE_TOPO-FILE.data')
    
    # print -nta option
    print(f'\n -nta or -n <nta-filename>   all2lmp variable: nta_file    hard coded: {nta_file}')
    print('    Command line option to set nta filename of file containing new type assignment (nta) file to be read into')
    print('    all2lmp for conversion. This file contains atom-ids and atom-types in columns, the atom-type in each row will')
    print('    be mapped onto the corresponding atom-id. The mapping between atom-ids and atom-types allows for the changing')
    print('    of atom-types or setting of atom-types. *NOTE: This code assumes correct atom-typing since the atom-types will')
    print('    set the atom, bond, angle, dihedral, improper coeff types and coeff parameters. See Materials Studio software')
    print('    below for msi2lmp conversion. Example usage:')
    print('        python3 all2lmp.py -nta EXAMPLE_NTA-FILE.nta')
    
    # Material Studio option
    print('\n -Biovias Materials Studio software support (msi)')
    print('    all2lmp v1.4 currently offers support to convert from Biovias Materials Studio software (msi). Currently')
    print('    another tool already exists for this purpose called msi2lmp that is available in the tool’s directory of all')
    print('    LAMMPS distributions, however the original msi2lmp has several limitations that all2lmp is trying to address.')
    print('    all2lmp uses the same inputs as the original msi2lmp code and that is a .car and .mdf file from material ')
    print('    studio. The .car and .mdf file can have material studio find atom-types via the forcite plugin and can also')
    print('    have material studio find charge via bond-incs (I have heard newer version of material studio has stopped')
    print('    supporting charge via bond-incs, but cannot confirm nor deny since I dont have direct access to material studio).')
    print('    This makes the use of material studio a nice route to build starting molecules however all2lmp comes with a set')
    print('    of atom-typing codes as well and all2lmp can reset charges via bond-incs as well. Current limitations for this')
    print('    tool are that the space group MUST BE P1! PBCs conditions are supported, and atoms will be wrapped (supports')
    print('    orthogonal and triclinic cells). If PBC flag is not provided box size is set to min/max of atom extents += 0.5')
    print('    angstroms for non-zero dims, if dim is zero the box size is set to +- 0.5 angstroms. Atoms will always be')
    print('    centered about 0, 0, 0, no user option currently exists to override this option. LAMMPS is better suited for')
    print('    this, and it makes wrapping atoms easier. To use all2lmp to convert from a material studio set of files to a')
    print('    LAMMPS datafile the following hard coded input variables for reading in material studio files will be used:')
    print('        topofile -> .mdf file (.mdf file contains bonding connectivity thus topofile variable seems fitting to use)')
    print('                     Info used from file: bonds')
    print('        nta_file -> .car file (.car file contains the atom types all2lmp will use as nta types thus nta_file variable')
    print('                     seems fitting to use.) Info used from file: atom-types, atom-positions, charge, and unit-cell')
    print('    Example command line override usage of all2lmp to convert from material studio .car and .mdf file using PCFF')
    print('    forcefield pareameters:')
    print('        python3 all2lmp.py -class 2 -topo benzene-class2b.mdf -nta benzene-class2b.car -frc PCFF.frc -ext PCFF')
    print('                                                             or with shortcut tags')
    print('        python3 all2lmp.py -c     2 -t    benzene-class2b.mdf -n   benzene-class2b.car -f   PCFF.frc -e   PCFF')
    
    # print -frc option
    print(f'\n -frc or -f <frc-filename>   all2lmp variable: frc_file    hard coded: {frc_file}')
    print('    Command line option to set frc filename of file containing force field parameters to be read into all2lmp for')
    print('    force field parameter assignment. Example usage:')
    print('        python3 all2lmp.py -frc EXAMPLE_FRC-FILE.frc')
    
    # print -add2box option
    print(f'\n -add2box or -a2b <float or in>   all2lmp variable: add2box    hard coded: {add2box}')
    print('    Command line option to set a floot or int that can be positive or negative that will be used to adjust each of the')
    print('    6-faces of the simulation cell box dimensions (angstroms). If the read-in topofile has none-zero image flags or')
    print('    image flags are dervived that are none-zero add2box will not be used since it is generally a bad idea to adjust')
    print('    the simulation cell dimensions if there are peirodic molecules spanning the current box. The box dimensions of each')
    print('    supported file format defaults are provided below:')
    print('        .mol or .sdf or .mol2 or .pdb or .car/.mdf (without box defined) will have image flags set to zero and simulation')
    print('        cell set at extent of atoms with a 0.5 angstrom buffer.')
    print()
    print('        .data or .car/.mdf (with box defined) will leave the simulation cell dimensions as is. A .data file will have the')
    print('        image flags left as is and a .car/.mdf (with box defined) will have the image flags derivied. ')
    print('    If add2box is set as zero no modification to the default simulation cell dimensions will occur. Example usage:')
    print('        python3 all2lmp.py -add2box 5')
    
    # print -ignore option
    print(f'\n -ignore or -i <T|F>   all2lmp variable: ignore_missing_parameters    hard coded: {ignore_missing_parameters}')
    print('    Command line option to ignore errors about missing parameters and provide a partially parameterized datafile if')
    print('    there are any missing parameters and if ignore_missing_parameters is True. If there are missing parameters and')
    print('    ignore_missing_parameters is False all parameters will be zeroed making them nonoperative in the written datafile.')
    print('    Please note that if ignore_missing_parameters is False and there are missing parameters all2lmp.py will prompt')
    print('    the user with a response during code execution to make this decision as well. However, if ignore_missing_parameters')
    print('    is True and there is missing parameters no prompt will be provided during all2lmp.py execution. Example usage:')
    print('        python3 all2lmp.py -ignore T')
    
    # print -asm option
    print(f'\n -asm or -a <assumed-filename>   all2lmp variable: assumed    hard coded: {assumed}')
    print('    Command line option to set assumed auto file coeffs filename. This file is unique to all2lmp and can be used')
    print('    to generalize bond/angle/dihedral/improper coeff types based on elements. The file contains section headers')
    print('    much like a True LAMMPS topofile, but under each header there are rows of elements that map onto generic atom')
    print('    types equivalences. If the code cannot find a fully parameterized coeff, this option allows for a level of')
    print('    generalization and assumption to automatically insert assumed types into the written topofile. The file')
    print('    provided gives examples and describes the layout of the file in the comments section. Comments in the file')
    print('    are led by the # character like LAMMPS and python and can help comment and organize to file or comment out')
    print('    headers to skip sections for testing. This option may be controversial to some and is meant for conversion from')
    print('    ReaxFF to other fix bond force fields where atom-type assignment is tricky. *NOTE: just providing the assumed')
    print('    auto fill filename is not enough to run the assumed auto fill option, you must also set use_assumed_auto_fill')
    print('    to True in the hard coded inputs section of the all2lmp.py script or using the [-assumed-equivs <T/F>] option')
    print('    at the command line.* Example usage:')
    print('        python3 all2lmp.py -asm EXAMPLE_ASSUMED-AUTO-FILL-FILE.coeffs')
    
    # print -dir option
    print(f'\n -dir or -d <new directory name>   all2lmp variable: parent_directory    hard coded: {parent_dir}')
    print('    Command line option to set new directory name to store all files bond_react_template_merge can write. The file')
    print('    management is based on the present working directory (pwd) and will store files one directory deep inside the')
    print('    pwd. The -dir option will specify the directory that is inside of pwd. This will allow for easier file ')
    print('    management and allows the user to merge multiple files for a certain system or research project and keep them')
    print('    all organized into a single directory.')
    print('    python variables*. Example usage:')
    print('        python3 all2lmp.py -dir EXAMPLE_TEST_SYSTEM_DIRECTORY')
    
    # print -newfile option
    print(f'\n -newfile or -nf <string>   all2lmp variable: newfile    hard coded: {newfile}')
    print('    Command line option to set the new output filename(s). The following options exist for using the newfile string')
    print('    for setting the output file basenames: ')
    print()
    print("      if newfile starts with ':' or ends with ':'")
    print("        The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix added to")
    print("        the file basename. The following are examples:")
    print("          Suffix (newfile = ':_IFF'  and topofile = 'detda.data')")
    print("            basename = 'detda_IFF', where the ':' character acts as a placeholder for the topofile basename.")
    print("          Prefix (newfile = 'IFF-:'  and topofile = 'detda.data')")
    print("            basename = 'IFF-detda', where the ':' character acts as a placeholder for the topofile basename.")
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
    print('        python3 all2lmp.py -newfile  :_IFF')
    print('        python3 all2lmp.py -newfile ": IFF"  # note the use of the " " bounding characters to pass a whitespace')
    
    # print -atomstyle option
    print(f'\n -atomstyle or -as <atomstyle>   all2lmp variable: atom_style    hard coded: {atom_style}')
    print('    Command line option to set atom style in written LAMMPS .data file. Currently supported atom styles are full,')
    print('    charge, molecular. *NOTE: If another atomstyle is desired it can be added in the write_lmp.py script.* Example')
    print('    usage:')
    print('        python3 all2lmp.py -atomstyle full')
    
    # print -class option
    print(f'\n -class or -c <0|1|2|d|s1|s2|i|ilmp>   all2lmp variable: ff_class    hard coded: {ff_class}')
    print('    Command line option to set the class of forcefield used. This will determine how the bonds/angles/diherdals/')
    print('    impropers are found, like wise how to coeffs are found and defined. There are 4 supported class types, 0, 1, 2')
    print('    and r, each of these MUST BE used with their correct .frc file. If the code detects the .frc file and class is')
    print('    inconsistent it will exit the code execution. Most cases have been tested, but it still is up to the user to')
    print('    make sure their inputs are consistent, and the outputs make sense. The class options are as follows with known')
    print('    FFs that use each class:')
    print('      0 = class0 (for OPLS-AA force field file in frc_files directory)')
    print('      1 = class1 (for cvff, clayff force field file in frc_files directory)')
    print('      2 = class2 (for PCFF-IFF, PCFF, compass force field file in frc_files directory)')
    print('      d = DREIDING (for all2lmp_dreiding.frc force field file in frc_files directory. Specfic to all2lmp ONLY!)')
    print()
    print('      s1 = skeleton datafile for class1; where skeleton means a complete LAMMPS datafile, but no coeffs in it (all')
    print('           N-body coeffs will have atomIDs in the Bonds, Angles, Dihedrals, and Impropers section match the ordering')
    print('           of the types in the coeff comment)')
    print()
    print('      s2 = skeleton datafile for class2; where skeleton means a complete LAMMPS datafile, but no coeffs in it (all')
    print('           N-body coeffs will have atomIDs in the Bonds, Angles, Dihedrals, and Impropers section match the ordering')
    print('           of the types in the coeff comment)')
    print()
    print('      i = for interatomic potentials like: ReaxFF, REBO, AIREBO, SNAP, ... (for .frc file specfic to all2lmp ONLY!')
    print('          contains elements and masses)')
    print("      ilmp = interatomic with the same meaning as 'i' ff_class, but the topofile is a LAMMPS datafile. When using 'ilmp'")
    print("               the atomTypeIDs set by the read-in LAMMPS datafile are maintained, where as when using 'i' and reading in")
    print("               a LAMMPS datafile, the atomTypeIDs are reset. This can be useful for when converting from a fix bond force")
    print("               field like PCFF to ReaxFF, where you want to keep the atomTypeID distinction based on the PCFF atom types.")
    print()
    print("  's1' and 's2' notes:")
    print("      Note, that when using either s1 or s2, the frc_file that is listed is not used and can be left as any frc file.")
    print("  'i' and 'ilmp' notes:")
    print("      Note that when using i, the nta_file that is listed is not used and can be left as any nta file. The element types and")
    print("      masses will be derived from the all2lmp_interatomic.frc file.")
    print('    Example usage:')
    print('        python3 all2lmp.py -class 2')

    # print -morse-bond option
    print(f'\n -morse-bond or -mb <T|F>  all2lmp variable: use_morse_bonds    hard coded: {use_morse_bonds}')
    print('    Command line option to use morse bond parameters over harmonic bond parameters (default is to use harmonic,')
    print('    thus invoking this command with T will enforce the morse bond parameters over the harmonic default). This')
    print('    option is only for class1 FFs only that have morse bond parameters in the .frc file! The use_auto_equivalence')
    print('    variable or override command -auto-equivs or -auto-e will also apply for class1 morse bond parameters. Example')
    print('    usage:')
    print('        python3 all2lmp.py -use_morse_bonds T')
    
    # print -auto-equivs option
    print(f'\n -auto-equivs or -auto-e <T|F>   all2lmp variable: use_auto_equivalence    hard coded: {use_auto_equivalence}')
    print('    Command line option to use auto-equivalent coeff types to supplement fully parameterized heuristic ones. T is for')
    print('    True and use option and F is for False. T will tell the code to use auto-equivalence types and F will skip usage.')
    print('    Example usage:')
    print('        python3 all2lmp.py -auto-equivs T')
    
    # print -assumed-equivs option
    print(f'\n -assumed-equivs or -assumed-e <T|F>   all2lmp variable: use_assumed_auto_fill    hard coded: {use_assumed_auto_fill}')
    print('    Command line option to use assumed auto fill equiv type mapping set by the -asm assumed auto file fill. T is')
    print('    for True and use option and F is for False. T will tell the code to use assumed auto fill option and F will')
    print('    skip usage. *NOTE: the usage and results of this option will be dependant on which elements to atom-type mapping')
    print('    you have in the -asm file read into the code. The usage of this option may be controversial to some and is meant')
    print('    for conversion from ReaxFF to other fix bond force fields where atom-type assignment is tricky. Example usage:')
    print('        python3 all2lmp.py -assumed-equivs T')
    
    # print -reset-molids option
    print(f'\n -reset-molids or -rm <T|F>   all2lmp variable: reset_molids    hard coded: {reset_molids}')
    print('    Command line option tp reset molids. This is done by finding all clustered atoms and then sorting them in size by')
    print('    number of atoms per cluster. The molids are then assigned based on largest cluster size to smallest cluster size.')
    print('    T is for True and use option and F is for False. T will tell the code to use find clusters (will add runtime to')
    print('    execution) and then assign molids to each atom based on which cluster it is contained in and F will skip usage.')
    print('    If this option is not used all molids will be set to 1. *NOTE: molids will only appear in the topofile if full')
    print('    atomstyle is used.* Example usage:')
    print('        python3 all2lmp.py -reset-molids T')
    
    # print -reset-charges option
    print(f'\n -reset-charges or -rq <T|F>   all2lmp variable: reset_charges    hard coded: {reset_charges}')
    print('    Command line option to reset charges. This is done by looping through the bonded atoms to a specfic atom and')
    print('    finding the bond-inc partial charge. The partial charge of each bonded atom is then summed together to set the')
    print('    total charge in that atom. T is for True and use option and F is for False. T will tell the code use bond-incs')
    print('    to set atom charge and F will skip usage. If usage is skipped the charge of the atom will come from the read')
    print('    topofile, for LAMMPS topofile this means it will be the charge of the original forcefield the file comes from')
    print('    and for .mol the code sets sets each atom charge as zeros. Example usage:')
    print('        python3 all2lmp.py -reset-charges T')
    
    # print -write-comments option
    print(f'\n -write-comments or -wc <T|F>   all2lmp variable: write_txt_comments    hard coded: {write_txt_comments}')
    print('    Command line option to write a .txt file of all comments all2lmp logs when finding and assigning forcefield')
    print('    parameter types. T is for True and use option and F is for False. T will tell the code to write the comments')
    print('    file which will be named as the final .data file, but with the .txt extension and F will skip usage. The final')
    print('    name will be set by the topofile name + "_" -ext option. Example usage:')
    print('        python3 all2lmp.py -write-comments T')
    
    # print -type-labels option
    print(f'\n -type-labels or -tl <T|F>   all2lmp variable: include_type_labels    hard coded: {include_type_labels}')
    print('    Command line option to inlcude LAMMPS new type labels in written datafile to make bond/angle/dihedrals/improper')
    print('    types if applicable consistent between datafiles or molecule files (primary use for bond/react file generation)')
    print('    *NOTE: this option will also change how -write-bond-react molecule files are written to be consistent between')
    print('    the written .data file and the written molecule .moltemp file.* Example usage:')
    print('        python3 all2lmp.py -type-labels T')
    
    # print -write-bond-react option
    print(f'\n -write-bond-react or -wbr <T|F>   all2lmp variable: write_bond_react    hard coded: {write_bond_react}')
    print('    Command line option to write LAMMPS molecule files to use for building bond react template files. T is for True')
    print('    and use option and F is for False. T will tell the code to write the molecule and coeffs file which will be')
    print('    named as the final .data file, but with the .moltemp and .ecoeffs extension and F will skip usage. The final')
    print('    name will be set by the topofile name + "_" -ext option. Example usage:')
    print('        python3 all2lmp.py -write-bond-react T')
    
    # print -rx or -ry or -rz
    print(f'\n -rx or -ry or -rz <float>   all2lmp variable: rotate    hard coded: {rotate}')
    print('    Command line option to control the optional rotation (in degrees) that can be applied to the molecular system about')
    print('    the x-axis, y-axis, and z-axis respectively. The values can be an int or a float value to set the rotation. To use')
    print('    this option the system must be none-periodic (i.e. all image flags are zero). Before the system is rotated, it is')
    print('    centered about (0,0,0), rotated, and then shifted back to its original location. After the system has been rotated,')
    print('    the simulation cell is redefined with the same amount of "x-padding", "y-padding", and "z-padding", between the atoms')
    print('    and the simulation cell before the rotation. The following -tags are available:')
    print('        -rx   <float> adjust the rotation in around x-axis')
    print('        -ry   <float> adjust the rotation in around y-axis')
    print('        -rz   <float> adjust the rotation in around z-axis')
    print('    Example usage:')
    print('        python3 all2lmp.py -rx 90 -ry 180 -rz 270')
    
    # print -sx or -sy or -sz
    print(f'\n -sx or -sy or -sz <float>   all2lmp variable: shift    hard coded: {shift}')
    print('    Command line option to control the optional shift to apply to a molecular systems. The values can be an int or a float value')
    print('    to set the shift in each direction. This shifts the atoms and the simulation cell. If -sx, -sy, and -sz  are set to ZERO, this')
    print('    option is not used. The following -tags are available:')
    print('        -sx   <float> adjust the shift about the x-axis')
    print('        -sy   <float> adjust the shift about the y-axis')
    print('        -sz   <float> adjust the shift about the z-axis')
    print('    Example usage:')
    print('        python3 all2lmp.py -sx 5 -sy 10 -sz 15')
    
    # print -opt or -man option
    print(f'\n -opt or -man   all2lmp variable: print_options    hard coded: {print_options}')
    print('    Will tell all2lmp to only print out avaiable command line options known as tagN and tagN-inputs. This is the')
    print('    only tag that doesnt require a tagN-input following the tag since the code will only look for if -opt is in')
    print('    command line input list. Example usage:')
    print('        python3 all2lmp.py -man')
        
    # print -gui option
    print('\n -gui <no addition info required> auto_morse_bond_update variable: use_GUI')
    print('    Command line option flag to load GUI. Example usage:')
    print('        python3 auto_morse_bond_update.py -gui')
    
    # print requirements and dependencies
    print('\n\nRequirements and Dependencies:')
    print('    Requirements: python 3.7+ to guarantee ordering of dictionaries which stores most of the data')
    print('    Dependencies: CORRECT ATOM-TYPING in either the .nta file or MSI files for conversion')
    
    # print buffer and exit code
    print('\n\n')
    return

# Class to update inputs
class inputs:
    def __init__(self, commandline_inputs, topofile, nta_file, frc_file, assumed, parent_dir, newfile, atom_style, ff_class, use_auto_equivalence,
                 use_assumed_auto_fill, reset_molids, reset_charges, write_txt_comments, write_bond_react, print_options,
                 use_morse_bonds, include_type_labels, add2box, ignore_missing_parameters, shift, rotate):
        
        # Give access to inputs (update later on if command line over ride is given)
        self.commandline_inputs = commandline_inputs
        self.topofile = topofile
        self.nta_file = nta_file
        self.frc_file = frc_file
        self.assumed = assumed
        self.parent_dir = parent_dir
        self.newfile = newfile
        self.atom_style = atom_style
        self.ff_class = ff_class
        self.use_auto_equivalence = use_auto_equivalence
        self.use_assumed_auto_fill = use_assumed_auto_fill
        self.reset_molids = reset_molids
        self.reset_charges = reset_charges
        self.write_txt_comments = write_txt_comments
        self.write_bond_react = write_bond_react
        self.print_options = print_options
        self.use_morse_bonds = use_morse_bonds
        self.include_type_labels = include_type_labels
        self.add2box = add2box
        self.ignore_missing_parameters = ignore_missing_parameters
        self.shift = shift
        self.rotate = rotate
        
        
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
        supported_tags = ['-topo', '-nta', '-frc', '-asm', '-dir', '-newfile', '-atomstyle', '-class', '-auto-equivs', '-assumed-equivs',
                          '-reset-molids', '-reset-charges', '-write-comments',  '-write-bond-react', '-morse-bond', '-type-labels',
                          '-add2box', '-ignore', '-rx', '-ry', '-rz', '-sx', '-sy', '-sz']
        
        # set shortcut_tags mapping
        shortcut_tags = {'-t':'-topo', '-n':'-nta', '-f':'-frc', '-a':'-asm', '-d':'-dir', '-nf':'-newfile', '-as':'-atomstyle', '-c':'-class',
                         '-auto-e':'-auto-equivs', '-assumed-e':'-assumed-equivs', '-rm':'-reset-molids', '-rq':'-reset-charges', '-i':'-ignore',
                         '-wc':'-write-comments', '-wbr':'-write-bond-react', '-mb': '-morse-bond', '-tl': '-type-labels', '-a2b':'-add2box',
                         '-rx':'-rx', '-ry':'-ry', '-rz':'-rz', '-sx':'-sx', '-sy':'-sy', '-sz':'-sz'}
        
        # set default variables
        default_variables ={'-topo': self.topofile, '-nta': self.nta_file, '-frc': self.frc_file, '-asm': self.assumed, '-dir':self.parent_dir, '-newfile': self.newfile,
                            '-atomstyle': self.atom_style, '-class': self.ff_class, '-auto-equivs': self.use_auto_equivalence, '-assumed-equivs': self.use_assumed_auto_fill, 
                            '-reset-molids': self.reset_molids, '-reset-charges': self.reset_charges, '-write-comments': self.write_txt_comments, '-ignore': self.ignore_missing_parameters,
                            '-write-bond-react': self.write_bond_react, '-morse-bond': self.use_morse_bonds, '-type-labels': self.include_type_labels,
                            '-add2box':self.add2box, '-rx':self.rotate['x'], '-ry':self.rotate['y'], '-rz':self.rotate['z'], '-sx':self.shift['x'], '-sy':self.shift['y'], '-sz':self.shift['z']}
        
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
        
        # set new -nta option and print confirmation
        if tags['-nta']:
            self.nta_file = tags['-nta']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-nta', self.nta_file))
            
        # set new -frc option and print confirmation
        if tags['-frc']:
            self.frc_file = tags['-frc']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-frc', self.frc_file))
            
        # set new -asm option and print confirmation
        if tags['-asm']:
            self.assumed = tags['-asm']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-asm', self.assumed))
            
        # set new -rx option and print confirmation
        if tags['-rx']:
            try: self.rotate['x'] = float(tags['-rx'])
            except: print(f'-rx: {tags["-rx"]} could not be converted to a float'); sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-rx', self.rotate))
            
        # set new -ry option and print confirmation
        if tags['-ry']:
            try: self.rotate['y'] = float(tags['-ry'])
            except: print(f'-ry: {tags["-ry"]} could not be converted to a float'); sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-ry', self.rotate))
            
        # set new -rz option and print confirmation
        if tags['-rz']:
            try: self.rotate['z'] = float(tags['-rz'])
            except: print(f'-rz: {tags["-rz"]} could not be converted to a float'); sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-rz', self.rotate))
            
        # set new -sx option and print confirmation
        if tags['-sx']:
            try: self.shift['x'] = float(tags['-sx'])
            except: print(f'-sx: {tags["-sx"]} could not be converted to a float'); sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-sx', self.shift))
            
        # set new -sy option and print confirmation
        if tags['-sy']:
            try: self.shift['y'] = float(tags['-sy'])
            except: print(f'-sy: {tags["-sy"]} could not be converted to a float'); sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-sy', self.shift))
            
        # set new -sz option and print confirmation
        if tags['-sz']:
            try: self.shift['z'] = float(tags['-sz'])
            except: print(f'-sz: {tags["-sz"]} could not be converted to a float'); sys.exit()
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-sz', self.shift))         
            
        # set new -add2box option and print confirmation
        if tags['-add2box']:
            self.add2box = float(tags['-add2box'])
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-add2box', self.add2box))
            
        # set new -ignore option and print confirmation
        if tags['-ignore']:
            self.ignore_missing_parameters = T_F_string2boolean('-ignore', (tags['-ignore']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-ignore', self.ignore_missing_parameters))
            
        # set new -dir option and print confirmation
        if tags['-dir']:
            if tags['-dir'].isspace() or tags['-dir'] == '':
                print('-dir flag was used, with an empty directory string. The directory string must contain alpha-numeric characters')
            else:
                self.parent_dir = tags['-dir']
                print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-dir', self.parent_dir))
            
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
        
        # set new -auto-equivs option and print confirmation
        if tags['-auto-equivs']:
            self.use_auto_equivalence = T_F_string2boolean('-auto-equivs', (tags['-auto-equivs']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-auto-equivs', self.use_auto_equivalence))
            
        # set new -assumed-equivs option and print confirmation
        if tags['-assumed-equivs']:
            self.use_assumed_auto_fill = T_F_string2boolean('-assumed-equivs', (tags['-assumed-equivs']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-assumed-equivs', self.use_assumed_auto_fill))
            
        # set new -reset-molids option and print confirmation
        if tags['-reset-molids']:
            self.reset_molids = T_F_string2boolean('-reset-molids', (tags['-reset-molids']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-reset-molids', self.reset_molids))
            
        # set new -reset-charges option and print confirmation
        if tags['-reset-charges']:
            self.reset_charges = T_F_string2boolean('-reset-charges', (tags['-reset-charges']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-reset-charges', self.reset_charges))
            
        # set new -type-labels option and print confirmation
        if tags['-type-labels']:
            self.include_type_labels = T_F_string2boolean('-type-labels', (tags['-type-labels']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-type-labels', self.include_type_labels))
            
        # set new -write-comments option and print confirmation
        if tags['-write-comments']:
            self.write_txt_comments = T_F_string2boolean('-write-comments', (tags['-write-comments']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-write-comments', self.write_txt_comments))
            
        # set new -write-bond-react option and print confirmation
        if tags['-write-bond-react']:
            self.write_bond_react = T_F_string2boolean('-write-bond-react', (tags['-write-bond-react']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-write-bond-react', self.write_bond_react))
            
        # set new -write-bond-react option and print confirmation
        if tags['-morse-bond']:
            self.use_morse_bonds = T_F_string2boolean('-morse-bond', (tags['-morse-bond']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-morse-bond', self.use_morse_bonds))
        
        # print buffer
        print('\n\n')