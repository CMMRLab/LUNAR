# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.9
January 5th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


    **********************************************************
    * Requirements:                                          *
    *   python 3.7+                                          *
    *                                                        *
    * Dependencies:                                          *
    *   python rdkit module:                                 *
    *    - pip3 install rdkit (if pip manager is installed)  *
    *    - not required (only adds functionality to read     *
    *      smiles strings in this LAT code only)             *
    *                                                        *
    * Run methods:                                           *
    *   - IDE (manipulate variables and run from IDE)        *
    *   - GUI (manipulate variables and run from. Default    *
    *          settings set from this script.)               *
    *   - command line (python3 atom_typing.py -man to get   *
    *                   command line override options and    *
    *                   examples)                            *
    *                                                        *
    **********************************************************
    
    https://acc2.ncbr.muni.cz/ from https://github.com/sb-ncbr/eem_parameters
    maybe a useful website for inserting charge into .mol2 files with a variety
    of different partial charge methods.
"""


##############
### Inputs ###
##############
#############################################################################################################################
# Python boolean variable to use GUI or not (True or False). All GUI settings are intialized from this script so adjusting  #
# this script will set the default settings that the GUI will load with. Please NOTE that when using the GUI to check the   #
# console or terminal print outs every time a system is run through the code because that is where the import information   #
# will be displayed. Examples:                                                                                              #
#   use_GUI = True  # Will launch GUI to manipulate variables and run from GUI                                              #
#   use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have to be from the IDE or command line#
#                                                                                                                           #
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account for the different   #
# ways individuals have their OS display settings setup. A GUI_zoom = 100, means that default settings are used, whereas a  #
# GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80 means decrease GUI by 20%. Examples:                      #
#   GUI_zoom = 100 # use default GUI size                                                                                   #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                                               #
#                                                                                                                           #
# Update use_GUI and GUI_zoom as desired.                                                                                   #
#############################################################################################################################
use_GUI = True
GUI_zoom = 100


#############################################################################################################################
# Files to be read in. Two system files are assigned to the following variables:                                            #
#   topofile = the main file containing atom positions, element/type, and possibly bonds. Supported files: LAMMPS .data/.dat#
#              or .mol or .sdf or .mol2 or .pdb file formats. If rdkit is install topofile can be a smiles string with a    #
#              .smiles extension (interpreted/read based on extension)                                                      #
#   bondfile = optional file that contains ReaxFF bonding information to create bonds for a ReaxFF system to convert to a   #
#              fix bond force field. If the topofile contains bonding info in it set bondfile = 'n.u' (short for 'not used')#
#              and the code WILL NOT attempt to read bondfile if topofile has bonding info in it, but bondfile variable     #
#              still needs to exist, hence 'n.u.' string.                                                                   #
#                                                                                                                           #
# Example for .mol/.data files with bonds in them:                                                                          #
#     topofile = 'detda.mol'  # Contains atoms, bonds, and elements                                                         #
#     bondfile = 'n.u.'       # Set bondfile variable to not used                                                           #
#                                                                                                                           #
#     topofile = 'detda.data' # Contains atoms, bonds, and atomtypes                                                        #
#     bondfile = 'n.u.'       # Set bondfile variable to not used                                                           #
#                                                                                                                           #
# Example for .data file for ReaxFF with a bond order file                                                                  #
#     topofile = 'SMP_10_mix_test.data'  # Contains atoms and atomtypes                                                     #
#     bondfile = 'SMP_10_mix_test.reaxc' # Contains bondorder info to create bonds                                          #
#                                                                                                                           #
# Example for .data file for ReaxFF with a bond order file w/ bondorder filename set with the 'topofile.ext' shortcut       #
#     topofile = 'SMP_10_mix_test.data'  # Contains atoms and atomtypes                                                     #
#     bondfile = 'topofile.reaxc' # Will make bondfile = 'SMP_10_mix_test.reaxc'                                            #
#                                                                                                                           #
# Example for if rdkit is installed use smiles extension add detda smiles code.                                             #
#     topofile = 'CCC1=CC(=C(C(=C1N)C)N)CC.smiles' # .smiles extension will be removed                                      #
#     bondfile = 'n.u.'                            # Set bondfile variable to not used                                      #
#                                                                                                                           #
# Paths can be added to file to read. Such as:                                                                              #
#     topofile = 'test/detda.mol'  # will read from $pwd/test directory                                                     #
#     topofile = '../detda.mol'    # will read from $pwd/../ directory                                                      #
#                                                                                                                           #
# Update topofile and bondfile as needed for the systems that you want to perform atom typing of.                           #
#############################################################################################################################
topofile = 'EXAMPLES/EPON_862/atom_typing_Inputs/detda.mol'
bondfile = 'n.u.'


#############################################################################################################################
# Python string variable type to set parent directory to store all new files this code will write. If the variable is left  #
# as an empty string, files will be written to the present working directory. Example to set parent_directory as present    #
# working directory:                                                                                                        #
#     parent_directory = '' or parent_directory = '.'                                                                       #
#                                                                                                                           #
# Example to set parent_directory as with path from topofile:                                                               #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile'                                          #
#     files will be written to 'Furan_Resin/files/' since the 'topofile' string envokes this option.                        #
#                                                                                                                           #
# Example to set parent_directory as with path from topofile and build dirs from that location:                             #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/NEWDIR'                                   #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/NEWDIR' will build relative directories       #
#                                                                                                                           #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/../NEWDIR'                                #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/../NEWDIR' will build relative directories    #
#                                                                                                                           #
# Example to set parent_directory to a location inside of present working dirctory:                                         #
#     parent_directory = 'inside'       # will save files in $pwd/inside                                                    #
#     parent_directory = 'inside/test'  # will save files in $pwd/inside/test                                               #
#                                                                                                                           #
# Example to set parent_directory to a location outside of present working dirctory:                                        #
#     parent_directory = '../'      # will save files in $pwd/../                                                           #
#     parent_directory = '../test'  # will save files in $pwd/../test                                                       #
#                                                                                                                           #
# Update parent_directory as desired.                                                                                       #
#############################################################################################################################
#parent_directory = 'topofile/../atom_typing_Outputs'
parent_directory = 'topofile'


#############################################################################################################################
# Python string variable to set the file that contains the Gasteiger charge parameters that will be used with the Gasteiger #
# charge method if reset_charges is True. Examples                                                                          #
#    chargefile = 'frc_files/Gasteiger_parameters.txt'  # Will read Gasteiger parameters from the selected file             #
#                                                                                                                           #
# Update parent_directory as desired.                                                                                       #
#############################################################################################################################
chargefile = 'frc_files\Gasteiger_parameters.txt'


#############################################################################################################################
# Is a Python string variable in which to set the new output filename(s). The following options exist for using the newfile #
# string for setting the output file basenames:                                                                             #
#                                                                                                                           #
#   if newfile starts with ':' or ends with ':'                                                                             #
#     The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix added to the file  #
#     basename. The following are examples:                                                                                 #
#       Suffix (newfile = ':_typed'  and topofile = 'detda.data')                                                           #
#         basename = 'detda_typed', where the ':' character acts as a placeholder for the topofile basename.                #
#       Prefix (newfile = 'typed-:'  and topofile = 'detda.data')                                                           #
#         basename = 'typed-detda', where the ':' character acts as a placeholder for the topofile basename.                #
#     Recommended usage: common and safe as this method is safe and output filename(s) carry similar names to input         #
#     filename(s).                                                                                                          #
#                                                                                                                           #
#   if newfile == 'ANYTEXT'                                                                                                 #
#     The output filename(s) will be set as 'ANYTEXT'. For example:                                                         #
#       newfile = 'detda_renamed' and topofile = 'detda.data')                                                              #
#         basename = 'detda_renamed'                                                                                        #
#     Recommended usage: occasional and safe as output filename(s) may no longer carry similar names as output filename(s), #
#     but is safe as output filename(s) will not overwrite input filename(s)                                                #
#                                                                                                                           #
#   if newfile == ''                                                                                                        #
#     The output filename(s) will be the same as the input filename(s). This can be a dangerous option as you may           #
#     inadvertently overwrite a file and then must assume the file to contain certain information, but it contains other    #
#     information.                                                                                                          #
#     Recommended usage: rare and dangerous as this could lead to using incorrect files if not paying attention very        #
#     carefully.                                                                                                            #
#                                                                                                                           #
# Update files newfile as needed.                                                                                           #
#############################################################################################################################
newfile = ':_typed'


#############################################################################################################################
# Set ff_name to determine atom-types for that specific ff. Currently supported ff_names (Case matters):                    #
# * Class2:                                                                                                                 #
#    - PCFF-IFF  most atom-types supported (FF-file:   frc_files/pcff_iff_v1_5_CNT_poly_solv.frc)                           #
#    - PCFF      most atom-types supported (FF-file:   frc_files/pcff.frc )                                                 #
#    - compass    all atom-types supported (FF-file:   frc_files/compass_published.frc)                                     #
#                                                                                                                           #
# * Class1:                                                                                                                 #
#    - CVFF-IFF  most atom-types supported (FF-file:   frc_files/cvff_interface_v1_5.frc)                                   #
#    - CVFF      most atom-types supported (FF-file:   frc_files/cvff.frc)                                                  #
#    - Clay-FF    all atom-types supported (FF-file:   frc_files/clayff.frc)                                                #
#    - DREIDING  most atom-types supported (FF-file:   frc_files/all2lmp_dreiding.frc)                                      #
#    - OPLS-AA    all atom-types supported (FF-file:   frc_files/oplsaa.frc {EXPEIRMENTAL})                                 #
#                                                                                                                           #
# General typing:                                                                                                           #
#    - general:0  will set atom type as elementRINGnb (IE C in 6 member ring TYPE = C63)                                    #
#    - general:1  will set atom type as elementRINGnb-(1st-neighs) and (1st-neighs) = (count:elementRINGnb,...)             #
#    - general:2  will set atom type as elementRINGnb-(1st-neighs)-(2nd-neighs) and (ith-neighs) = (count:elementRINGnb,...)#
#    - general:3  will set atom type as elementRINGnb-(1st)-(2nd)-(3rd) and (ith-neighs) = (count:elementRINGnb,...)        #
#    - general:4  will set atom type as elementRINGnb-(1st)-(2nd)-(3rd)-(4th) and (ith-neighs) = (count:elementRINGnb,...)  #
#    NOTES: When using the general atom typing module the code will set an "equivs" header with all the unqiue atom type    #
#    strings in column-1. The user will have to manually add column-2 to set the equivalent atom type that they would like  #
#    to use. all2lmp.py will then automatically convert all atom types in the "type id" section to the mapped atom types    #
#    set in the "equivs" section. This will allow the LAT codes to be more extendable for building models when atom typing  #
#    modules have not yet been built for a specific force field.                                                            #
#                                                                                                                           #
# Examples:                                                                                                                 #
#    ff_name = 'PCFF-IFF'  # Will load PCFF-IFF atom typing module to set PCFF-IFF atom types                               #
#    ff_name = 'PCFF'      # Will load PCFF atom typing module to set PCFF atom types                                       #
#    ff_name = 'DREIDING'  # Will load DREIDING atom typing module to set DREIDING atom types                               #
#    ff_name = 'general:1' # Will load general atom typing module to and will use atom and 1st neighs to set atom types     #
#                                                                                                                           #
# Update ff_name as desired.                                                                                                #
#############################################################################################################################
ff_name = 'PCFF-IFF'


#############################################################################################################################
# Python boolean variable to reset charges with Gasteiger method. The code will first perform a hybridization analysis to   #
# set correct Gasteiger parameters (True or False).                                                                         #
#                                                                                                                           #
# This operation uses the algorthim from the following paper:                                                               #
#    - Iterative Equalization of Oribital Electronegatiity A Rapid Access to Atomic Charges                                 #
#                                                                                                                           #
# Then in all2lmp set reset_charges as False to keep Gasteiger charges from this phase of the model building process. This  #
# process will be performed after delete atoms so if you delete any atoms they will not play are role in the charge         #
# calculation. Examples:                                                                                                    #
#    reset_charges = True   # Will find hybridization state, set Gasteiger parameters and converge them                     #
#    reset_charges = False  # Will keep charges as is from read in file or if not present set the to zeros                  #
#                                                                                                                           #
# Update reset_charges as desired.                                                                                          #
#############################################################################################################################
reset_charges = False


#############################################################################################################################
# Python dictionary to control cluster size that will be kept. Anything cluster less than this size will not  show up in    #
# the written *.data file nor the written *.nta file (atomIDs will be sorted and reset to make them contiguous). If you     #
# want to use this code to keep all atoms, set this variable to zero. This will allow the code to find all clusters and all #
# clusters. Dictionary keys/values and their meanings:                                                                      #
# delete_atoms = {'method': METHOD-OPTION ('mass' or 'size' string)                                                         #
#                 'criteria': CRITERIA-OPTION (float or int value) }                                                        #
#                                                                                                                           #
# where METHOD-OPTION sets the search method for determining cluster quanitiy to keep. Supported methods:                   #
#    'mass' will search for clusters via mass cut-off where anything less than the CRITERIA-OPTION will be                  #
#           removed before any analysis occurs such as setting charge or finding atom types.                                #
#    'size' will search for clusters via number of atoms cut-off where anything less than the CRITERIA-OPTION               #
#           will be removed before any analysis occurs such as setting charge or finding atom types.                        #
#                                                                                                                           #
# The reason for this option is if you are trying to convert a ReaxFF simulation to a fix bond force field, it is possible  #
# that the ReaxFF simulation has small by-products that you may want to remove and this is the easiest stage to remove the  #
# possible by-products. For topofile extension(s) like a .mol or .mol2 or .smiles the 'criteria' likely should be set to    #
# zero to keep all atoms. Examples:                                                                                         #
#    delete_atoms = {'method': 'mass', 'criteria': 100 } # Will remove clusters ligher then 100 amu                         #
#    delete_atoms = {'method': 'mass', 'criteria': 0 }   # Will keep all clusters                                           #
#    delete_atoms = {'method': 'size', 'criteria': 0 }   # Will keep all clusters                                           #
#    delete_atoms = {'method': 'size', 'criteria': 5 }   # Will remove smaller clusters than 5 atoms                        #
#                                                                                                                           #
# Update delete_atoms as desired.                                                                                           #
#############################################################################################################################
delete_atoms = {'method': 'mass', # 'size' or 'mass'
                'criteria': 0.0} # if 'size' natoms if 'mass' mass in amu


#############################################################################################################################
# Set pdb_file to determine if any additional files will be written by atom_typing.py. These options exist to help with     #
# interfacing with packmol, but may also be found to be useful outside of interfacing with packmol. Currently supported     #
# pdb_file options and their meanings (Case matters):                                                                       #
#    - 'skip'    will not write additional *_packmol.pdb and *_packmol.nta files                                            #
#    - 'types'   will write additional *_packmol.pdb file but not *_packmol.nta, with the atom types set in the atom name   #
#    - 'typeIDs' will write additional *_packmol.pdb and *_packmol.nta files, with atomTypeIDs set in the atom name column  #
#                of the .pdb file and atom types set in the .nta file using the "style type" method                         #
#                                                                                                                           #
# Update pdb_file as desired.                                                                                               #
#############################################################################################################################
pdb_file = 'skip'


#############################################################################################################################
# Python boolean variable to print options avaible for the code at the command line interface and exit (True or False).     #
# True will print options and exit code execution and False will allow to run normally. The same printouts can be accessed  #
# with the command line override -opt or -man command line inputs. Examples:                                                #
#     print_options = True   # Will print out command line manual and exit                                                  #
#     print_options = False  # Will allow code to run to completion                                                         #
#                                                                                                                           #
# Update print_options as desired.                                                                                          #
#############################################################################################################################
print_options = False


#############################################################################################################################
# The mass_map dictionary is now a "global" dictionary stored in src/masses.py. The purpose of this was to simplify adding  #
# new elements, where the new elements can now be applied to every code that uses the mass_map. If you get an "ERROR Not    #
# all masses in ... are in the mass_map dictionary.", you will now have to open src/masses.py and update the mass_map       #
# dictionary found in that file.                                                                                            #
#############################################################################################################################
import src.masses as masses
mass_map = masses.mass_map 


#############################################################################################################################
# Specific inputs for determining bonds if bonds are not provided in the read-in file. atom_typing.py currently supports:   #
#     1.) Bonds via ReaxFF bond order file:                                                                                 #
#           Which will read in a ReaxFF bond order file in the format of the original format of the 1st stand-alone ReaxFF  #
#           code developed by Adri Van Duin's group via the LAMMPS command:                                                 #
#               fix             1 all reax/c/bonds 100 ${myid}.reaxc                                                        #
#                                     or                                                                                    #
#               fix             1 all reaxff/bonds 100 ${myid}.reaxff                                                       #
#           This option can read bond order files with multiple time steps of bond order info and average the bond orders   #
#           together to generate a time-averaged bond order. This means it is up to the user to log as many timesteps of    #
#           data as desired and it is also up to the user to understand the effects of temperature and dynamics present     #
#           during bond order sampling. Once the bond orders are averaged together a complex bond determination algorithm   #
#           will determine viable bonds by applying minimum bond order cut-offs and max bonded cut-offs both of which can   #
#           be adjusted in the following variables respectively; bondorder (python dictionary with key as tuple pairs for   #
#           each bond type and value as float minimum bond order cut-off) and maxbonded (python dictionary with key as      #
#           element symbols and value as maxbonded integer cut-off).                                                        #
#                                                                                                                           #
#    2.) Bonds via interatomic distance searching:                                                                          #
#          Which will find bonds by looping through all atomIDs and search for interatomic distances that are within the    #
#          vdw waals radius cutoff (cut-offs used can be found in /src/bonds_via_distance.py/vdw_radius dictionary). The    #
#          "scale" of the vdw radius cut-off can be adjusted by the vdw_radius_scale variable (python float value to set    #
#          the multiplication value of the vdw radius). The bonds will be searched via the boundary conditions supplied by  #
#          the boundary variable (python string with 3-white space separated boundary flags). Boundary flags are set up     #
#          like LAMMPS where 'p' is periodic and 'f' is not periodic. If boundary = 'p p p' all directions are assumed to   #
#          be periodic and bonds will be searched across all faces of the simulation cell. Similarly, any combinations of   #
#          'p' or 'f' provided will set the search space for bonds via interatomic distances. Lastly, a max bonded cutoff   #
#          set by the maxbonded dictionary will be applied and will reject the longest bonds if the number of bonded atoms  #
#          per element type goes over the max bonded value set in the python dictionary. NOTE: The usage of this option     #
#          will be automatically enforced if the read-in topofile does not contain any bonds. This means if you wish to     #
#          use this for a ReaxFF atom typing then set bondfile = 'n.u.' and just search the interatomic distances to find   #
#          the most reasonable bonds. Future codes will be written to read-in files without bonds such as .cif, .xyz,       #
#          POSCAR, and so on which will utilize this functionality to offer support for a greater variety of input files as #
#          well. The 'p' flag in any direction is currently limited to orthogonal boxes ONLY. The code will exit if the     #
#          topofile has no bonds and 'p' is in any location in the boundary string and the box is not orthogonal.           #
#                                                                                                                           #
# Below you will find the python variables discussed above and can manipulate them as needed. The ordering is random and    #
# the bondfile variable is above next to the topofile variable. Adjust each setting as needed and keep in mind which method #
# from above is being used to generate the bonds. The comment will start with:                                              #
#     ReaxFF: which will signify usage for only using ReaxFF bond order files                                               #
#     Distance: which will signify usage for only finding bonds via interatomic distances                                   #
#     ReaxFF/Distance: which will signify usabge for either ReaxFF or Distance based usage                                  #
#                                                                                                                           #
# *NOTE: if performing a ReaxFF atom typing, it is still recommended to use a bond order file due to the run time in python #
# to search interatomtic distances as well as a bond order file can be written to over small periods of time (like 1ps) to  #
# perform time averaging to capture thermal fluctuations. The bonds by distance option is still a valid method, but for     #
# systems with many atoms (above 10,000) it is very slow in python and no time averaging of atomic positions can be found   #
# in a static snapshot of the LAMMPS datafile.                                                                              #
#############################################################################################################################

# ReaxFF: Possbile bonding configurations bond order cut-offs. This bond order cut-off will be used for the time-averaged 
# bond order for each bond type. This gives the ability to adjust the minimum bond order cut-off for each bond type. If the
# bond does not exist in this bondorder dictionary, the 'unknown' key bondorder cut-off will be enforce. Add onto the 
# dictionary as needed, the code will adjust based on the info provided in the bondorder and maxbonded dictionary. The code
# will re-order the dictionary bonding pair element keys before use so only one ordering of a bond pair needs to be given 
# (IE forward and reverse ordering of C-H bond only needs either be listed as ('C','H') or ('H','C') and the code will look
# for both forward and reverse orderings). Update as needed.

#          like-paired     pairing of all other configurations
bondorder = {('C','C'):    0.3, ('C','H'):    0.3, ('C','O'):    0.3, ('C','N'):    0.3, ('C','S'):    0.3,   # C-other  bonds
             ('C','F'):    0.3, ('C', 'Si'):  0.3, ('C', 'D'):   0.3, ('C', 'Cl'):  0.3, ('C', 'Ca'):  0.3,   # C-other  bonds
             ('C', 'Br'):  0.3, ('C', 'P'):   0.3,  ('C', 'Al'): 0.3, ('C', 'Mg'):  0.3, ('C', 'Li'):  0.3,   # C-other  bonds
             ('C', 'Fe'):  0.3, ('C', 'Na'):  0.3,                                                            # C-other  bonds
             
             ('H','H'):    0.3, ('H','O'):    0.3, ('H','N'):    0.3, ('H','S'):    0.3, ('H','F'):    0.3,   # H-other  bonds
             ('H', 'Si'):  0.3, ('H', 'D'):   0.3, ('H', 'Cl'):  0.3, ('H', 'Ca'):  0.3, ('H', 'Br'):  0.3,   # H-other  bonds
             ('H', 'P'):   0.3, ('H', 'Al'):  0.3, ('H', 'Mg'):  0.3, ('H', 'Li'):  0.3, ('H', 'Fe'):  0.3,   # H-other  bonds
             ('H', 'Na'):  0.3,                                                                               # H-other  bonds
             
             ('O','O'):    0.3, ('O','N'):    0.3, ('O','S'):    0.3, ('O','F'):    0.3, ('O', 'Si'):  0.3,   # O-other  bonds
             ('O', 'D'):   0.3, ('O', 'Cl'):  0.3, ('O', 'Ca'):  0.3, ('O', 'Br'):  0.3, ('O', 'P'):   0.3,   # O-other  bonds
             ('O', 'Al'):  0.3, ('O', 'Mg'):  0.3, ('O', 'Li'):  0.3, ('O', 'Fe'):  0.3, ('O', 'Na'):  0.3,   # O-other  bonds
             
             ('N','N'):    0.3, ('N','S'):    0.3, ('N','F'):    0.3, ('N', 'Si'):  0.3, ('N', 'D'):   0.3,   # N-other  bonds
             ('N', 'Cl'):  0.3, ('N', 'Ca'):  0.3, ('N', 'Br'):  0.3, ('N', 'P'):   0.3, ('N', 'Al'):  0.3,   # N-other  bonds
             ('N', 'Mg'):  0.3, ('N', 'Li'):  0.3, ('N', 'Fe'):  0.3, ('N', 'Na'):  0.3,                      # N-other  bonds
             
             ('S','S'):    0.3, ('S','F'):    0.3, ('S', 'Si'):  0.3, ('S', 'D'):   0.3, ('S', 'Cl'):  0.3,   # S-other  bonds
             ('S', 'Ca'):  0.3, ('S', 'Br'):  0.3, ('S', 'P'):   0.3, ('S', 'Al'):  0.3, ('S', 'Mg'):  0.3,   # S-other  bonds
             ('S', 'Li'):  0.3, ('S', 'Fe'):  0.3, ('S', 'Na'):  0.3,                                         # S-other  bonds
             
             ('F','F'):    0.3, ('F', 'Si'):  0.3, ('F', 'D'):   0.3, ('F', 'Cl'):  0.3, ('F', 'Ca'):  0.3,   # F-other  bonds
             ('F', 'Br'):  0.3, ('F', 'P'):   0.3, ('F', 'Al'):  0.3, ('F', 'Mg'):  0.3, ('F', 'Li'):  0.3,   # F-other  bonds
             ('F', 'Fe'):  0.3, ('F', 'Na'):  0.3,                                                            # F-other  bonds
             
             ('Si', 'Si'): 0.3, ('Si', 'D'):  0.3,  ('Si', 'Cl'):0.3, ('Si', 'Ca'): 0.3, ('Si', 'Br'): 0.3,   # Si-other bonds
             ('Si', 'P'):  0.3, ('Si', 'Al'): 0.3,  ('Si', 'Mg'):0.3, ('Si', 'Li'): 0.3, ('Si', 'Fe'): 0.3,   # Si-other bonds
             ('Si', 'Na'): 0.3,                                                                               # Si-other bonds
             
             ('Al', 'Al'): 0.3, ('Al', 'Mg'): 0.3, ('Al', 'Li'): 0.3, ('Al', 'Fe'): 0.3, ('Al', 'Na'): 0.3,   # Al-other bonds
             
             ('Mg', 'Mg'): 0.3, ('Mg', 'Li'): 0.3, ('Mg', 'Fe'): 0.3, ('Mg', 'Na'): 0.3,                      # Mg-other bonds
             
             ('Li', 'Li'): 0.3, ('Li', 'Fe'): 0.3, ('Li', 'Na'): 0.3,                                         # Li-other bonds
             
             ('Fe', 'Fe'): 0.3, ('Fe', 'Na'): 0.3,                                                            # Fe-other bonds
             
             ('Na', 'Na'): 0.3,                                                                               # Na-other bonds
             
             ('D', 'D'):   0.3, ('D', 'Cl'):  0.3, ('D', 'Ca'):  0.3, ('D', 'Br'):  0.3,  ('D', 'P'):  0.3,   # D-other  bonds
             'unknown':    0.3}                                                                               # unknown  bonds

# ReaxFF/Distance: maxbonded cut-off based on element type to reduce the number of bonded atoms to each element type in
# ReaxFF systems or systems where bonds must be inferred via interatomic distances. Update as needed.
maxbonded = {'C':4, 'H':1, 'O':2, 'N':4, 'S':2, 'F':1, 'Si':4, 'Xe':0, 'Ne':0, 'Kr':0, 'He':0, 'D':1, 'Cl':1,
             'Ca':1, 'Br':1, 'Ar':0, 'P':5, 'Al': 6, 'Mg': 6, 'Li': 6, 'Fe': 6, 'Na': 1, 'K':0, 'Cs':0, 'Ba':0,
             'Sr':0, 'Pb':0}


# Distance: vdw_radius_scale will set the float value multiplier to adjust the vdw radius for interatomic distance searching.
# A value of 1.0 sets the vdw as it is, a value of 1.1 adds 10% more distance onto the vdw radius of each element, ...
vdw_radius_scale = 1.1


# Distance: Boundary to set the type of boundary to search for bonds. Setup like LAMMPS with the following meaning:
#   p is periodic
#   f is non-periodic and fixed
# This will be set by a python string with white space between each x, y, or z box-face flags. (IE 'p p p' sets all box sides
# to be periodic; 'f f f' sets all box sides to be non-periodic; 'p f f' sets the two x-box sides to be periodic and leaves 
# the y/z-box sides as non-periodic ...). *NOTE: there MUST BE white space between each direction flag and there MUST BE 
# 3-direction flags of 'f' or 'p' in any combination. Also searching for all periodic faces (27 total images) is much more 
# computationally expensive. The interatomic distance search code is time-optimized as much as Josh has been able to make it,
# but if suffers from pythons slow nested looping issues limitations ...
boundary = 'p p p'

# Override option to find bonds via interatomic distances. bonds_via_distance is usually only used if no bonds currenly exist
# in the readin topofile however this override will allow you to search for bonds via distance even if there are already bonds 
# defined in the topofile. This option is useful if you want to generate a CNT or graphene sheet in VMD and want the structure to
# be periodically bonded. You may then adjust boundary to set which directions to search for periodic bonds of (True or False).
# When reading in a .mol or .mol2 file the atoms are centered about 0, 0, 0 and the simulation cell dimensions are set at the 
# extent of atoms in each direction + 0.5 angstroms on all edges (+1 in X-dir, +1 in Y-dir, and +1 in Z-dir) so atoms near one
# edge of the of the simulation cell will only be ~1 angstroms from the other edge of the simulation cell. Thus the bonds found
# to be periodic should be of high qualitiy. 
bonds_via_distance_override = False # Usualy should be False except for Rare cases described above




################################################
### Import main from src.atom_typing and run ###
################################################
if __name__ == "__main__":  
    # Import main from src.atom_typing.main
    from src.atom_typing.main import main
    import sys

    ####################################################################################################
    # Set IDE/command line mode option:                                                                #
    # - if run in IDE commandline_inputs lst will have zero entries. Which means use hard coded inputs #
    #   from inputs section.                                                                           #
    #                                                                                                  #
    # - if run in command line, but commandline_inputs still has zero entries, this means use hard     #
    #   coded inputs from inputs section.                                                              #
    #                                                                                                  #
    # - else run in command line and commandline_inputs has entries determine what entries will be     #
    #   overrided from hard coded inputs from inputs section. This overide will occur in in the inputs #
    #   section called Command Line Override. *NOTE: if user does not specify certain command line     #
    #   options, the default settings from inputs section will be enforced. A message will tell user   #
    #   if inputs from inputs section have been enforced. The user can find all command line options   #
    #   by running:                                                                                    #
    #       python3 atom_typing.py -opt                                                                #
    #              or                                                                                  #
    #       python3 atom_typing.py -man                                                                #
    #   in the terminal and the code will print out all command line option and terminate before any   #
    #   further analysis is done                                                                       #
    ####################################################################################################
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: use_GUI = False
    if use_GUI or '-gui' in commandline_inputs:
        from src.atom_typing.GUI import atom_typing_GUI
        print('\n\n\natom_typing is currently running in GUI mode, where all GUI inputs are intialized from atom_typing.\n\n\n')
        atom_typing_GUI(topofile, bondfile, parent_directory, newfile, ff_name, delete_atoms, mass_map, bondorder, maxbonded, boundary,
                        vdw_radius_scale, reset_charges, print_options, [], bonds_via_distance_override, pdb_file, chargefile, GUI_zoom)
    else:
        # If not commandline_inputs the code is running in IDE/command line mode with out command line inputs
        # this means that all inputs to the code are handled in the Inputs section below this block of code.
        if not commandline_inputs:
            print('\n\n\natom_typing is currently running in IDE mode or command line mode with no command line inputs.')
            print('This means all inputs used for the code comes from the atom_typing.py script in the inputs section.\n\n\n')
        
        # Run main atom_typing
        main(topofile, bondfile, parent_directory, newfile, ff_name, delete_atoms, mass_map, bondorder, maxbonded, boundary,
             vdw_radius_scale, reset_charges, print_options, commandline_inputs, bonds_via_distance_override, pdb_file, chargefile)