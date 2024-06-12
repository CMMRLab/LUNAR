# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.19
June 5th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
 
    ********************************************************
    * Requirements:                                        *
    *   python 3.7+                                        *
    *                                                      *
    * Dependencies:                                        *
    *   correct atom-typing                                *
    *                                                      *
    * Run methods:                                         *
    *   - IDE (manipulate variables and run from IDE)      *
    *   - GUI (manipulate variables and run from. Default  *
    *          settings set from this script)              *
    *   - command line (python3 all2lmp.py -man to get     *
    *                   command line override options and  *
    *                   examples)                          *
    *                                                      *
    ********************************************************
    
    https://acc2.ncbr.muni.cz/ from https://github.com/sb-ncbr/eem_parameters
    maybe a useful website for inserting charge into .mol2 files with a variety
    of different partial charge methods.
"""


##############
### Inputs ###
##############
###############################################################################################################################
# Python boolean variable to use GUI or not (True or False). All GUI settings are intialized from this script so adjusting    #
# this script will set the default settings that the GUI will load with. Please NOTE that when using the GUI to check the     #
# console or terminal print outs every time a system is run through the code because that is where the import information     #
# will be displayed. Examples:                                                                                                #
#    use_GUI = True  # Will launch GUI to manipulate variables and run from GUI                                               #
#    use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have to be from the IDE or command line #
#                                                                                                                             #
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account for the different     #
# ways individuals have their OS display settings setup. A GUI_zoom = 100, means that default settings are used, whereas a    #
# GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80 means decrease GUI by 20%. Examples:                        #
#   GUI_zoom = 100 # use default GUI size                                                                                     #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                                                 #
#                                                                                                                             #
# Update use_GUI and GUI_zoom as desired.                                                                                     #
###############################################################################################################################
use_GUI = True
GUI_zoom = 100


###############################################################################################################################
# Files to be read in. Four system files are assigned to the following variables:                                             #
#    topofile = the main file containing atom positions, element/type, and possibly bonds. Supported files:                   #
#               LAMMPS .data/.dat or .mol or .sdf or .mol2 or .pdb or MSI .mdf file formats (interpreted based on extension)  #
#                                                                                                                             #
#    nta_file = file containing atom types map. For MSI conversions this will come from the MSI .car file and performed       #
#               without any user intervention. For all other file conversions the nta_file will be a *.nta file with a        #
#               specfic format as such (read by headers):                                                                     #
#                  - "#" is a comment character and anything trailing after "#" will be ignored                               #
#                  - Currently supported headers with their formats and meanings:                                             #
#                      "style id" which will map atom types via atomIDs (compatibility .mol .mol2 .data .sdf)                 #
#                       1 cp   # will map cp atom type to topofile file atomID: 1                                             #
#                       2 hc   # will map hc atom type to topofile file atomID: 2                                             #
#                                                                                                                             #
#                      "style type" which will map atom types via atomTypeIDs (compatibility .data or manually adjust .pdb)   #
#                       1 cp   # will map cp atom type to topofile file atomTypeID: 1                                         #
#                       2 hc   # will map hc atom type to topofile file atomTypeID: 2                                         #
#                                                                                                                             #
#                      "equivs" which will map all types in column-1 to types set in column-2 (compatibility .mol .mol2 .data)#
#                       hn2  hn   # will map all atom types that are hn2 to atom type  hn                                     #
#                       hc   hpan # will map all atom types that are hc to atom type  hpan                                    #
#                                                                                                                             #
#                      "edge id" which adds atoms to bond-incs searches (compatibility .mol .mol2 .data .sdf)                 #
#                       1 cp hc  # will "extend" cp & hc onto edge atomID: 1 for setting charge via bond-incs                 #
#                       2 hc oc  # will "extend" hc & oc onto edge atomID: 2 for setting charge via bond-incs                 #
#                                                                                                                             #
#                      "charge id" which sets charge of atomIDs after bond-incs search (compatibility .mol .mol2 .data .sdf)  #
#                       1 1.111  # will set charge atomID: 1 to 1.1111 after bond-incs search                                 #
#                       2 2.222  # will set charge atomID: 2 to 2.2222 after bond-incs search                                 #
#                                                                                                                             #
#                      "charge type" which sets charge of atomTypeIDs after bond-incs (compatibility .data ONLY!)             #
#                       1 1.111  # will set charge atomTypeIDs: 1 to 1.1111 after bond-incs search                            #
#                       2 2.222  # will set charge atomTypeIDs: 2 to 2.2222 after bond-incs search                            #
#                                                                                                                             #
#                      "charge nta" which sets charge of atomTypeLabels after bond-incs (compatibility .mol .mol2 .data .sdf) #
#                      If using the "equivs" option the atom-type in this section corresponds to the 2nd-column atom-type in  #
#                      the "equivs" section (IE the atom-type that will be used to assign FF parameters).                     #
#                       hc 1.111  # will set charge atomTypeLabels: hc to 1.1111 after bond-incs search                       #
#                       cp 2.222  # will set charge atomTypeLabels: hc to 1.1111 after bond-incs search                       #
#                                                                                                                             #
#                      "neutralize system charge all" will make system charge nuetral by adding a fixed value to "all" atoms  #
#                       in the system                                                                                         #
#                                                                                                                             #
#                      "neutralize system charge zero" will make system charge nuetral by adding a fixed value to atoms that  #
#                       have a "zero" charge in the system                                                                    #
#                                                                                                                             #
#                      "neutralize system charge bond-inc" will make system charge nuetral by adding a fixed value to atoms   #
#                       that are defined as "bond-inc" atoms. "bond-inc" atoms will be "all" atoms in the system if not using #
#                       any of the optional "charge <opt>" headers. If using any of the "charge <opt>" headers, "bond-inc"    #
#                       atoms are defined as atoms that have not had charge set via the "charge <opt>" headers.               #
#                                                                                                                             #
#                      "neutralize system charge user-defined" will make system charge nuetral by adding a fixed value to     #
#                       atoms that are defined as "user-defined" atoms. "user-defined" atoms will be atoms in the system using#
#                       any of the "charge <opt>" headers and will ONLY have a fixed charge added to the atoms defined by the #
#                       "charge <opt>" headers.                                                                               #
#                                                                                                                             #
#               Either "style id" or "style type" is a mandatory header in the *.nta file with the correct mapping onto atom  #
#               info found in the read-in topofile. All other headers are optional and can be useful when attempting to add a #
#               charge to a variety of different situations. "edge id" is useful when creating bond/react templates since the #
#               summation of bond-incs of the missing atoms can be added via "edge id" header. "charge id", "charge type",    #
#               and "charge nta" can be useful when setting charges when that information is known. The "neutralize system    #
#               <opt>" is useful to ensure charge neutrality when combining the different charge methods available. If the    #
#               python reset_charges boolean variable is set to False, the system charge will remain unchanged by bond-incs,  #
#               but the "charge <opt>" and "neutralize system <opt>" will still function. Lending to more flexible charging   #
#               scenerios, that may or maynot be useful. If the python reset_charges boolean variable is set to True, the     #
#               system charge is guaranteed to be neutral since the "bond-inc" method will always produce a net zero charged  #
#               system. This means if charge is set via "bond-incs", you may skip usage of the "neutralize system <opt>" in   #
#               the nta_file.                                                                                                 #
#                                                                                                                             #
#    frc_file = the force field to search for energy coefficients after all bonds/angles/dihedrals/impropers are found and    #
#               typed. The file format is to be in an MSI *.frc file format. Currently, available *.frc files can be found in #
#               /frc_files directory and more can be added if/when needed.                                                    #
#                                                                                                                             #
#    assumed =  file containing general equivalence mapping of atom types to elements. This also has a specific format that   #
#               will not be addresed here and is mainly meant for ReaxFF to fix bond force field conversions. Most of the time#
#               you can ignore this as an option and leave it as some arbiturary python string. use_assumed_auto_fill must be #
#               True to try reading. If use_assumed_auto_fill = False the code will not try looking for the assumed file.     #
#                                                                                                                             #
# Example for a .mol file (found from $pwd/test directory) conversion to LAMMPS datafile w/ PCFF-IFF_v1.5:                    #
#     topofile = 'test/detda.mol'  # Contains atoms, bonds, and elements                                                      #
#     nta_file = 'test/detda.nta'  # Contains minimally the "style id" or "style type" header                                 #
#     frc_file = 'frc_files/pcff_iff_v1_5_CNT_poly_solv.frc' # PCFF-IFF_v1.5 *.frc parameter file                             #
#     assumed =  'arbitrary' # use_assumed_auto_fill = False so wont read (set as arbitrary string)                           #
#                                                                                                                             #
# Example for a .data file (found from $pwd/../ directory) conversion to LAMMPS datafile w/ cvff:                             #
#     topofile = '../detda.data' # Contains atoms, bonds, and atom types                                                      #
#     nta_file = '../detda.nta'  # Contains minimally the "style id" or "style type" header                                   #
#     frc_file = 'frc_files/cvff.frc' # cvff.frc parameter file                                                               #
#     assumed =  'general_assumed_equivs.coeffs' # use_assumed_auto_fill = True so read in user-defined assumed equivs file   #
#                                                                                                                             #
# Example for a MSI .mdf and .car file (found from $pwd directory) conversion to LAMMPS datafile w/ PCFF-IFF_v1.5:            #
#     topofile = 'furan.mdf'  # atom postions, atom types, charges, space group and simulation cell info                      #
#     nta_file = 'furan.car'  # bonding connectivty and atom types                                                            #
#     frc_file = 'frc_files/pcff_iff_v1_5_CNT_poly_solv.frc' # PCFF-IFF_v1.5 *.frc parameter file                             #
#     assumed =  'arbitrary' # use_assumed_auto_fill = False so wont read (set as arbitrary string)                           #
#                                                                                                                             #
# Example using the 'topofile' short cut for the nta_file for a MSI .mdf and .car file and a .mol and .nta file:              #
#     topofile = 'furan.mdf'  # atom postions, atom types, charges, space group and simulation cell info                      #
#     nta_file = 'topofile'   # nta_file will assume the basename of topofile with the .car ext (IE nta_file = furan.car)     #
#     topofile = 'furan.mol'  # Contains atoms, bonds, and elements                                                           #
#     nta_file = 'topofile'   # nta_file will assume the basename of topofile with the .car ext (IE nta_file = furan.nta)     #
#                                                                                                                             #
# If the nta_file variable is set to 'types_from_pdb.nta', all2lmp will look for atom type strings in the location of the     #
# atom name column of the .pdb file and then apply those atom types to the system. This method requires manual editing of the #
# .pdb file and is mainly meant for an integration method with packmol, since packmol will maintain the atom name column.     #
# Plese look in EXAMPLES/packmol_pdb_methods/atom_name_as_atom_type for an example water test of this method. Settings:       #
# 	topofile = 'EXAMPLES/packmol_pdb_methods/atom_name_as_atom_type/water_packmol_atom_name_as_atom_type.pdb'                 #
#	nta_file = 'types_from_pdb.nta'                                                                                           #
#                                                                                                                             #
# Update topofile, nta_file, frc_file, and assumed as needed.                                                                 #
###############################################################################################################################
topofile = 'EXAMPLES/EPON_862/atom_typing_Outputs/detda_typed.data'
nta_file = 'EXAMPLES/EPON_862/atom_typing_Outputs/detda_typed.nta'
frc_file = 'frc_files/pcff.frc'
assumed = 'frc_files/general_assumed_equivs.coeffs'


###############################################################################################################################
# add2box is a python float or int value that can be positive or negative that will be used to adjust each of the 6-faces of  #
# the simulation cell box dimensions (angstroms). If the read-in topofile has none-zero image flags or image flags are        #
# dervived that are none-zero add2box will not be used since it is generally a bad idea to adjust the simulation cell         #
# dimensions if there are peirodic molecules spanning the current box. The box dimensions of each supported file format       #
# defaults are provided below:                                                                                                #
#    .mol or .sdf or .mol2 or .pdb or .car/.mdf (without box defined) will have image flags set to zero and simulation cell   #
#    set at extent of atoms with a 0.5 angstrom buffer.                                                                       #
#                                                                                                                             #
#    .data or .car/.mdf (with box defined) will leave the simulation cell dimensions as is. A .data file will have the image  #
#     flags left as is and a .car/.mdf (with box defined) will have the image flags derivied.                                 #
#                                                                                                                             #
# If add2box is set as zero no modification to the default simulation cell dimensions will occur (update as needed).          #
###############################################################################################################################
add2box = 0.0


###############################################################################################################################
# Python variable (True or False) to ignore errors about missing parameters and provide a partially parameterized datafile if #
# there are any missing parameters and if ignore_missing_parameters is True. If there are missing parameters and              #
# ignore_missing_parameters is False all parameters will be zeroed making them nonoperative in the written datafile. Please   #
# note that if ignore is False and there are missing parameters all2lmp.py will prompt the user with a response during code   #
# execution to make this decision as well. However, if ignore_missing_parameters is True and there is missing parameters no   #
# prompt will be provided during all2lmp.py execution. Examples:                                                              #
#     ignore_missing_parameters = False  # if missing parameters exist the all2lmp will provide a prompt to decide how to     #
#                                          handle the situation                                                               #
#     ignore_missing_parameters = True   # if missing parameters exist all2lmp will provide a partially parameterized file    #
#                                          with zero's as the missing parameters.                                             #
#                                                                                                                             #
# Update ignore_missing_parameters as desired.                                                                                #
###############################################################################################################################
ignore_missing_parameters = True


###############################################################################################################################
# Python string variable type to set parent directory to store all new files this code will write. If the variable is left    #
# as an empty string, files will be written to the present working directory. Example to set parent_directory as present      #
# working directory:                                                                                                          #
#     parent_directory = '' or parent_directory = '.'                                                                         #
#                                                                                                                             #
# Example to set parent_directory as with path from topofile:                                                                 #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile'                                            #
#     files will be written to 'Furan_Resin/files/' since the 'topofile' string envokes this option.                          #
#                                                                                                                             #
# Example to set parent_directory as with path from topofile and build dirs from that location:                               #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/NEWDIR'                                     #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/NEWDIR' will build relative directories         #
#                                                                                                                             #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/../NEWDIR'                                  #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/../NEWDIR' will build relative directories      #
#                                                                                                                             #
# Example to set parent_directory to a location inside of present working dirctory:                                           #
#     parent_directory = 'inside'       # will save files in $pwd/inside                                                      #
#     parent_directory = 'inside/test'  # will save files in $pwd/inside/test                                                 #
#                                                                                                                             #
# Example to set parent_directory to a location outside of present working dirctory:                                          #
#     parent_directory = '../'      # will save files in $pwd/../                                                             #
#     parent_directory = '../test'  # will save files in $pwd/../test                                                         #
#                                                                                                                             #
# Update parent_directory as desired.                                                                                         #
###############################################################################################################################
#parent_directory = 'topofile/../all2lmp_Outputs'
parent_directory = 'topofile'


###############################################################################################################################
# Is a Python string variable in which to set the new output filename(s). The following options exist for using the newfile   #
# string for setting the output file basenames:                                                                               #
#                                                                                                                             #
#   if newfile starts with ':' or ends with ':'                                                                               #
#     The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix added to the file    #
#     basename. The following are examples:                                                                                   #
#       Suffix (newfile = ':_IFF'  and topofile = 'detda.data')                                                               #
#         basename = 'detda_IFF', where the ':' character acts as a placeholder for the topofile basename.                    #
#       Prefix (newfile = 'IFF-:'  and topofile = 'detda.data')                                                               #
#         basename = 'IFF-detda', where the ':' character acts as a placeholder for the topofile basename.                    #
#     Recommended usage: common and safe as this method is safe and output filename(s) carry similar names to input           #
#     filename(s).                                                                                                            #
#                                                                                                                             #
#   if newfile == 'ANYTEXT'                                                                                                   #
#     The output filename(s) will be set as 'ANYTEXT'. For example:                                                           #
#       newfile = 'detda_renamed' and topofile = 'detda.data')                                                                #
#         basename = 'detda_renamed'                                                                                          #
#     Recommended usage: occasional and safe as output filename(s) may no longer carry similar names as output filename(s),   #
#     but is safe as output filename(s) will not overwrite input filename(s)                                                  #
#                                                                                                                             #
#   if newfile == ''                                                                                                          #
#     The output filename(s) will be the same as the input filename(s). This can be a dangerous option as you may             #
#     inadvertently overwrite a file and then must assume the file to contain certain information, but it contains other      #
#     information.                                                                                                            #
#     Recommended usage: rare and dangerous as this could lead to using incorrect files if not paying attention very          #
#     carefully.                                                                                                              #
#                                                                                                                             #
# Update files newfile as needed.                                                                                             #
###############################################################################################################################
newfile = ':_IFF'


###############################################################################################################################
# Python string variable to set atom style format of the written data file. Currently supported atom styles:                  #
# charge, molecular, full, angle, bond, atomic, dipole, dpd, or line. Example:                                                #
#     atom_style = 'full' # will set atom style in written data file as full                                                  #
#                                                                                                                             #
# Update atom_style as desired.                                                                                               #
###############################################################################################################################
atom_style = 'full'


###############################################################################################################################
# Python string or int variable to set what type of forcefield (FF) to apply to the system. *NOTE the FF type set here MUST   #
# be consistent with the information inside the read in FF set by the frc_file variable.* The following options are available:#
#    0   = class0      (int data type - FF files: opls-AA {EXPEIRMENTAL})                                                     #
#    1   = class1      (int data type - FF files: cvff-IFF, cvff, clayff)                                                     #
#    2   = class2      (int data type - FF files: PCFF-IFF, PCFF, compass)                                                    #
#   'd'  = DREIDING    (str data type - FF file: all2lmp_dreiding.frc)                                                        #
#   'i'  = interatomic (str data type - FF files: all2lmp all2lmp interatomic specfic .frc file, for potentials like ReaxFF,  #
#                      REBO, AIREBO, SNAP, ... {contains elements and mass info, simulation cell info, insertion of extra     #
#                      elements as atom types to make reaxFF files consistant with one another if different elements are in   #
#                      each file})                                                                                            #
#   's1' = skeleton datafile for class1; where skeleton means a complete LAMMPS datafile, but no coeffs in it (all N-body     #
#          coeffs will have atomIDs in the Bonds, Angles, Dihedrals, and Impropers section match the ordering of the types in #
#          the coeff comment)                                                                                                 #
#   's2' = skeleton datafile for class2; where skeleton means a complete LAMMPS datafile, but no coeffs in it (all N-body     #
#          coeffs will have atomIDs in the Bonds, Angles, Dihedrals, and Impropers section match the ordering of the types in #
#          the coeff comment)                                                                                                 #
#                                                                                                                             #
# Examples:                                                                                                                   #
#    ff_class = 1    # sets FF application for class1 FF's                                                                    #
#    ff_class = 'i'  # sets FF application for interatomic potentials like: ReaxFF, REBO, AIREBO, SNAP, ...                   #
#    ff_class = 's2' # generates a class2 LAMMPS datafile with no energy coeffs, but in class2 format                         #
#                                                                                                                             #
# Update ff_class as required.                                                                                                #
###############################################################################################################################
ff_class = 2


###############################################################################################################################
# Python dictionary with the keys set as 'x', 'y', and 'z'  which controls the optional shift to apply to a molecular systems.#
# The values of the dictionary are int or float values to set the shift in each direction. This shifts the atoms and the      #
# simulation cell. If 'x', 'y', and 'z'  are set to ZERO, this option is not used.                                            #
#                                                                                                                             #
# Update shift as required.                                                                                                   #
###############################################################################################################################
shift = {'x': 0.0, 'y': 0.0, 'z': 0.0}


###############################################################################################################################
# Python dictionary with the keys set as 'x', 'y', and 'z'  which controls an optional rotation (in degrees) that can be      #
# applied to the molecular system about the x-axis, y-axis, and z-axis respectively. The values of the dictionary are int or  #
# float values to set the rotation. To use this option the system must be none-periodic (i.e. all image flags are zero).      #
# Before the system is rotated, it is centered about (0,0,0), rotated, and then shifted back to its original location. After  #
# the system has been rotated, the simulation cell is redefined with the same amount of "x-padding", "y-padding", and         #
# "z-padding", between the atoms and the box before the rotation. If 'x', 'y', and 'z'  are set to ZERO, this option is not   #
# used.                                                                                                                       #
#                                                                                                                             #
# Update rotate as required.                                                                                                  #
###############################################################################################################################
rotate = {'x': 0.0, 'y': 0.0, 'z': 0.0}


###############################################################################################################################
# Python boolean variable to use auto-equivalences to supplement fully parameterized interactions with heuristic ones (True   #
# or False). *NOTE usage will also result in searches for bond-incs, crossterm r0's, and crossterm theta0's to be supplemented#
# as well.*                                                                                                                   #
#                                                                                                                             #
# Examples:                                                                                                                   #
#    use_auto_equivalence = True  # Will try finding fully parameterized interations and supplement w/auto-equivs if needed.  #
#    use_auto_equivalence = False # Will try finding fully parameterized interations and skip auto-equiv.                     #
#                                                                                                                             #
# Update use_auto_equivalence as desired (Default should be True).                                                            #
###############################################################################################################################
use_auto_equivalence = True


###############################################################################################################################
# Python boolean variable to use morse bonds over harmonic bonds for Class1 or DRIEDING FF's ONLY! use_auto_equivalence still #
# applies and will search auto-equivalent forms for the morse bond, should they exist in the FF *.frc file (True or False).   #
#                                                                                                                             #
# Examples:                                                                                                                   #
#    use_morse_bonds = True   # Will set bond-coeffs as morse bond definition for Class1 or DRIEDING FF's                     #
#    use_morse_bonds = False  # Will set bond-coeffs as standard harmonic definition for Class1 or DRIEDING FF's              #
#                                                                                                                             #
# Update use_morse_bonds as desired.                                                                                          #
###############################################################################################################################
use_morse_bonds = False


###############################################################################################################################
# Python boolean variable to reset molids for molecules in system (True or False). Examples:                                  #
#    reset_molids = True  # Will find all molecules in system and assign new molids to them                                   #
#    reset_molids = False # Will leave molids as defaults that have been read-in (If they DO NOT EXIST like in .mol or .mol2  #
#                           files all atoms will be assigned to molid: 1)                                                     #
#                                                                                                                             #
# Update reset_molids as desired (Default should be True).                                                                    #
###############################################################################################################################
reset_molids = True


###############################################################################################################################
# Python boolean variable to reset charges via the bond-increments (bond-incs) method. If the read-in frc file does not have  #
# bond-incs for the bonding pair a zero will be summed and the user will be warned (True or False). Examples:                 #
#    reset_charges = True  # Will attempt finding bond-incs for all bonded atoms and sum to find per/atom charge              #
#    reset_charges = False # Will keep charge of read-in file, if charges don't exist they will be set to zeros               #
#                                                                                                                             #
# Update reset_charges as desired (Default should be True).                                                                   #
###############################################################################################################################
reset_charges = True


###############################################################################################################################
# Python boolean variable to use assumed auto fill coeffs to fill in coeffs that could not be found, which uses coeff types   #
# specified in the assumed file (True or False). This typically should be reserved for reaxFF to fix bond forcefiel(s)        #
# applications since there are some uncertainties in atom-types for the conversion. Examples:                                 #
#    use_assumed_auto_fill = True # Will read the file from the assumed variable and apply the general element mapping to     #
#                                   coeffs that could not be found.                                                           #
#    use_assumed_auto_fill = False # Will skip usage of this option.                                                          #
#                                                                                                                             #
# Update use_assumed_auto_fill as desired (Default should be False).                                                          #
###############################################################################################################################
use_assumed_auto_fill = False


###############################################################################################################################
# Python boolean variable to write a verbose *.txt file with a long list of comments for parameters. Typically used for debug #
# info (True or False). Examples:                                                                                             #
#    write_txt_comments = True  # Will write verbose comment file (useful to keep track of parameters).                       #
#    write_txt_comments = False # Will skipping writing of verbose comment file.                                              #
#                                                                                                                             #
# Update write_txt_comments as desired (Default should be True).                                                              #
###############################################################################################################################
write_txt_comments = True


###############################################################################################################################
# Python boolean variable to assign LAMMPS's new Type Labels section(s) to the written LAMMPS data file (True or False).      #
# Examples:                                                                                                                   #
#    include_type_labels = False # Will write the "standard" format of the LAMMPS data file                                   #
#    include_type_labels = True  # Will write the Type Labels section(s) in the LAMMPS data file                              #
#                                                                                                                             #
# Update include_type_labels as desired (Default should be False).                                                            #
###############################################################################################################################
include_type_labels = False


###############################################################################################################################
# Python boolean variable to write files for fix bond/react templates (True or False). The include_type_labels option is      #
# integrated with this option as well. *NOTE: this will allow bond/react template creation via the old manual method that     #
# CMMRL group use to do. There is a secondary code called bond_react_merge.py which will read all data files of a system and  #
# automates the mapping of the coefficient Types and also can auto-generate the reaction map file. Thus only use this option  #
# if you want to perform manual bond/react template creation. Otherwise, set it to False and use bond_react_template_merge*.  #
# Examples:                                                                                                                   #
#    write_bond_react = True  # Will write a *.ecoeffs (energy coeffs file) and *.lmpmol (LAMMPS molecule file)               #
#    write_bond_react = False # Will skip writing of a *.ecoeffs and *.moltemp                                                #
#                                                                                                                             #
# Update write_bond_react as desired (Default should be False).                                                               #
###############################################################################################################################
write_bond_react = False


###############################################################################################################################
# Python boolean variable to print options avaible for the code at the command line interface and exit (True or False). True  #
# will print options and exit code execution and False will allow to run normally. The same printouts can be accessed with    #
# the command line override -opt or -man command line inputs. Examples                                                        #
#     print_options = True   # Will print out command line manual and exit                                                    #
#     print_options = False  # Will allow code to run to completion                                                           #
#                                                                                                                             #
# Update print_options as desired.                                                                                            #
###############################################################################################################################
print_options = False




############################################
### Import main from src.all2lmp and run ###
############################################
if __name__ == "__main__":  
    # Import main from src.all2lmp.main
    from src.all2lmp.main import main
    import sys
    
    ###########################################################################################################################
    # Set IDE/command line mode option:                                                                                       #
    # - if run in IDE commandline_inputs lst will have zero entries. Which means use hard coded inputs from inputs section.   #
    #                                                                                                                         #
    # - if run in command line, but commandline_inputs still has zero entries, this means use hard coded inputs from inputs.  #
    #                                                                                                                         #
    # - else run in command line and commandline_inputs has entries determine what entries will be overrided from hard coded  #
    #   inputs from inputs section. This overide will occur in in the inputs section called Command Line Override. *NOTE: if  #
    #   user does not specify certain command line options, the default settings from inputs section will be enforced. A      #
    #   message will tell user if inputs from inputs section have been enforced. The user can find all command line options   #
    #   by running:                                                                                                           #
    #       python3 all2lmp.py -opt                                                                                           #
    #              or                                                                                                         #
    #       python3 all2lmp.py -man                                                                                           #
    #   in the terminal and the code will print out all command line option and terminate before any  further analysis is done#
    ###########################################################################################################################
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: use_GUI = False
    if use_GUI or '-gui' in commandline_inputs:
        from src.all2lmp.GUI import all2lmp_GUI
        print('\n\n\nall2lmp is currently running in GUI mode, where all GUI inputs are intialized from all2lmp.\n\n\n')
        all2lmp_GUI(topofile, nta_file, frc_file, assumed, parent_directory, newfile, atom_style, ff_class, use_auto_equivalence, use_assumed_auto_fill,
                    reset_molids, reset_charges, write_txt_comments, write_bond_react, print_options, use_morse_bonds, include_type_labels, add2box, 
                    ignore_missing_parameters, shift, rotate, GUI_zoom, [])   
    else:        
        # If not commandline_inputs the code is running in IDE/command line mode with out command line inputs
        # this means that all inputs to the code are handled in the Inputs section below this block of code.
        if not commandline_inputs:
            print('\n\n\nall2lmp is currently running in IDE mode or command line mode with no command line inputs.')
            print('This means all inputs used for the code comes from the all2lmp.py script in the inputs section.\n\n\n')
            
        # Run main all2lmp classes/functions (and get parameters class)
        main(topofile, nta_file, frc_file, assumed, parent_directory, newfile, atom_style, ff_class, use_auto_equivalence, use_assumed_auto_fill, reset_molids,
             reset_charges, write_txt_comments, write_bond_react, print_options, use_morse_bonds, include_type_labels, add2box, ignore_missing_parameters,
             shift, rotate, commandline_inputs)
    
