# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.9
April 2nd, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


    **************************************************************
    * Requirements:                                              *
    *   python 3.7+                                              *
    *                                                            *
    * Dependencies:                                              *
    *   LAMMPS datafile bond coeffs must be in class2 or         *
    *   harmonic format                                          *
    *                                                            *
    *   python matplotlib module:                                *
    *    - pip3 install matplotlib (if pip manager is installed) *
    *                                                            *
    * Run methods:                                               *
    *   - IDE (manipulate variables and run from IDE)            *
    *   - GUI (manipulate variables and run from. Default        *
    *          settings set from this script)                    *
    *   - command line (python3 auto_morse_bond_update.py -man   *
    *                   to get command line override options and *
    *                   examples)                                *
    *                                                            *
    *                                                            *
    **************************************************************
"""

##############
### Inputs ###
##############
###########################################################################################################
# Python boolean variable to use GUI or not (True or False). All GUI settings are intialized from this    #
# script so adjusting this script will set the default settings that the GUI will load with. Please NOTE  #
# that when using the GUI to check the console or terminal print outs every time a system is run through  #
# the code because that is where the import information will be displayed. Examples:                      #
#   use_GUI = True  # Will launch GUI to manipulate variables and run from GUI                            #
#   use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have to be from the  #
#                     IDE or command line                                                                 #
#                                                                                                         #
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account   #
# for the different ways individuals have their OS display settings setup. A GUI_zoom = 100, means that   #
# default settings are used, whereas a GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80    #
# means decrease GUI by 20%. Examples:                                                                    #
#   GUI_zoom = 100 # use default GUI size                                                                 #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                             #
#                                                                                                         #
# Update use_GUI and GUI_zoom as desired.                                                                 #
###########################################################################################################
use_GUI = True
GUI_zoom = 100


###########################################################################################################
# LAMMPS datafile with bonds in it (REQUIRED!). The LAMMPS datafile can have new LAMMPS "type labels" in  #
# it. The only supported atom styles are full, charge, or molecular. More atom styles can be coded if     #
# needed. NOTE: This code will only works for LAMMPS datafiles in class2 or harmonic bond coeff format    #
# due to ordering of the parameters in the bond coeffs.                                                   #
###########################################################################################################
topofile = 'EXAMPLES/auto_morse_bond/detda_typed_IFF.data'


###########################################################################################################
# Python string variable to set the file that contains the bond typing rules and the corresponding Morse  #
# bond dissociation energy coefficients. The morsefile currently has most of the bond typing rules set for#
# all elemental/hybridization bonding configurations of the elements set in the mass_map dictionary. The  #
# morsefile may be added onto or adjusted as desired to users liking, however, the file is meant to be as #
# comprehensive as possible and should rarely need to be modified.                                        #
###########################################################################################################
morsefile = 'frc_files/Morse_parameters.txt'


###########################################################################################################
# Option to zero effected crossterms. Will help stabilize IFF-R completely (True or False). This option   #
# techincally changes the force field from that defined as IFF-R, since IFF-R does not do anything about  #
# crossterms that use bond lengths to compute them. Therefore usage of this option requires a recognition #
# of this and the associated force field name will be IFF-RX0 standing for Reactive Interface Force Field #
# w/crossterms (X) zerod (0). This option only is for class2 FFs, if your file is in class1/harmonic      #
# format this option can be True or False (since class1 has no xterms, this option does nothing).         #
###########################################################################################################
zero_effected_xterms = False


###########################################################################################################
# Is a Python string variable in which to set the new output filename(s). The following options exist for #
# using the newfile string for setting the output file basenames:                                         #
#                                                                                                         #
#   if newfile starts with ':' or ends with ':'                                                           #
#     The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix  #
#     added to the file basename. The following are examples:                                             #
#       Suffix (newfile = ':_morse_bond'  and topofile = 'detda.data')                                    #
#         basename = 'detda_morse_bond', where the ':' character acts as a placeholder for the topofile   #
#                     basename.                                                                           #
#       Prefix (newfile = 'morse_bond-:'  and topofile = 'detda.data')                                    #
#         basename = 'morse_bond-detda', where the ':' character acts as a placeholder for the topofile   #
#                     basename.                                                                           #
#     Recommended usage: common and safe as this method is safe and output filename(s) carry similar      #
#     names to input filename(s).                                                                         #
#                                                                                                         #
#   if newfile == 'ANYTEXT'                                                                               #
#     The output filename(s) will be set as 'ANYTEXT'. For example:                                       #
#       newfile = 'detda_renamed' and topofile = 'detda.data')                                            #
#         basename = 'detda_renamed'                                                                      #
#     Recommended usage: occasional and safe as output filename(s) may no longer carry similar names as   #
#     output filename(s), but is safe as output filename(s) will not overwrite input filename(s)          #
#                                                                                                         #
#   if newfile == ''                                                                                      #
#     The output filename(s) will be the same as the input filename(s). This can be a dangerous option as #
#     you may inadvertently overwrite a file and then must assume the file to contain certain information,#
#     but it contains other information.                                                                  #
#     Recommended usage: rare and dangerous as this could lead to using incorrect files if not paying     #
#     attention very carefully.                                                                           #
#                                                                                                         #
# Update files newfile as needed.                                                                         #
###########################################################################################################
newfile = ':_morse_bond'


###########################################################################################################
# Parameters that affect Morse bond fit onto quartic bonds. alpha_specs sets the range of alpha parameters#
# to test for minimizing the morse bond onto the harmonic bond. alpha_scale sets a multipleir value to the#
# fit alpha parameter found via minimization of eulcidean distances. alpha_scale = 1.0 means use exact    #
# minimized value, 1.1 means increase the alpha paramter by 10%, 1.2 means increase the alpha paramter by #
# 20%. ... radius_specs sets the bond length spectrum to plot over (NOTE should always start at zero and  #
# minimally go to 6 angstroms in 0.01 increments). This code now also supports the auto-generation of a   #
# psuedo-LAMMPS input script for generating fix bond/break fixes. The Rmax in fix bond/break can be       #
# adjusted with bondbreak_scale, which sets a a multiplier value to the r0 of each bond coeff to find the #
# Rmax value for fix bond/break. bondbreak_scale = 2.0 means set Rmax 200% of orginal bond length.        #
########################################################################################################### 
radius_specs = {'start': 0.0, 'end': 8.0, 'increment': 0.01} # usually: {'start':0.0, 'end': 8.0, 'increment': 0.01}
alpha_specs  = {'start': 1.0, 'end': 3.5, 'increment': 0.1} # usually: {'start':1.0, 'end': 3.5, 'increment': 0.1}
alpha_scale = 1.0 # usually: 1.0 (increase to create narrower well and decrease to widen well)
bondbreak_scale = 2.0 # usually: 1.7-2.0
include_rcut = False # option to have the bond_coeffs set as ID d0 r0 alpha rcut     # where rcut will be used within LAMMPS to shift the morse potential


###########################################################################################################
# Python string variable type to set parent directory to store all new files this code will write. If the #
# variable is left as an empty string, files will be written to the present working directory. Example to #
# set parent_directory as present working directory:                                                      #
#     parent_directory = '' or parent_directory = '.'                                                     #
#                                                                                                         #
# Example to set parent_directory as with path from topofile:                                             #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile'                        #
#     files will be written to 'Furan_Resin/files/' since the 'topofile' string envokes this option.      #
#                                                                                                         #
# Example to set parent_directory as with path from topofile and build dirs from that location:           #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/NEWDIR'                 #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/NEWDIR' will build relative #
#     directories                                                                                         #
#                                                                                                         #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/../NEWDIR'              #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/../NEWDIR' will build       #
#     relative directories                                                                                #
#                                                                                                         #
# Example to set parent_directory to a location inside of present working dirctory:                       #
#     parent_directory = 'inside'       # will save files in $pwd/inside                                  #
#     parent_directory = 'inside/test'  # will save files in $pwd/inside/test                             #
#                                                                                                         #
# Example to set parent_directory to a location outside of present working dirctory:                      #
#     parent_directory = '../'      # will save files in $pwd/../                                         #
#     parent_directory = '../test'  # will save files in $pwd/../test                                     #
#                                                                                                         #
# Update parent_directory as desired.                                                                     #
###########################################################################################################
parent_directory = 'topofile'


###########################################################################################################
# Python string variable to set atom style format of the written data file. Currently supported atom      #
# styles: charge, molecular, full, angle, bond, atomic, dipole, dpd, or line. Example:                    #
#     atom_style = 'full' # will set atom style in written data file as full                              #
#                                                                                                         #
# Update atom_style as desired.                                                                           #
###########################################################################################################
atom_style = 'full'


###########################################################################################################
# Python int variable to set what type of forcefield the system is to initialize zeroes and update once   #
# parameters are found.The following options are available:                                               #
#    1   = class1   (int data type - FF files: cvff, clayff, DREIDING, opls-AA, CHARMM)                   #
#    2   = class2   (int data type - FF files: PCFF-IFF, PCFF, compass, CFF91)                            #
#                                                                                                         #
# Examples:                                                                                               #
#     ff_class = 1 # Will find bond coeffs based on LAMMPS bond style harmonic                            #
#     ff_class = 2 # Will find bond coeffs based on LAMMPS bond style class2                              #
#                                                                                                         #
# Update ff_class as required.                                                                            #
###########################################################################################################
ff_class = 2


###########################################################################################################
# Minimum bond length of bond coeff to update from quadratic/quartic to morse bond potential (Angstrom).  #
###########################################################################################################
min_bond_length = 1.2


###########################################################################################################
# Bond coeffs TypeIDs to skip over and leave as harmonic (if empty no coeffs will be skipped over).       #
###########################################################################################################
coeffs2skip = []


###########################################################################################################
# The mass_map dictionary is now a "global" dictionary stored in src/masses.py. The purpose of this was   #
# to simplify adding new elements, where the new elements can now be applied to every code that uses the  #
# mass_map. If you get an "ERROR Not all masses in ... are in the mass_map dictionary.", you will now     #
# have to open src/masses.py and update the mass_map dictionary found in that file.                       #
###########################################################################################################
import src.masses as masses
mass_map = masses.mass_map 


###########################################################################################################
# Commands for creating Files that this code can write (True or False responses).                         #
###########################################################################################################
files2write = {'write_datafile' : True,   # File containing extracted clusters w/ bonds (name *NEWFILE.data)
               'write_pdffile'  : True,   # File containing plotted morse vs harmonic fits (name *NEWFILE.pdf)
               'write_bondbreak': True,   # File containing fix bond/break fixes with Rmax set (name *NEWFILE.script)
               }


###########################################################################################################
# Python boolean variable to assign LAMMPS's new Type Labels section(s) to the written LAMMPS data file   #
# (True or False). NOTE This option ONLY works for datafiles with all2lmp style comments or datafiles     #
# that already have type labels defined in them or a combination of the two. Examples:                    #
#    include_type_labels = False # Will write the "standard" format of the LAMMPS data file               #
#    include_type_labels = True  # Will write the Type Labels section(s) in the LAMMPS data file          #
#                                                                                                         #
# Update include_type_labels as desired (Default should be False).                                        #
###########################################################################################################
include_type_labels = False




#####################################################
### Import main from src.bond_react_merge and run ###
#####################################################
if __name__ == "__main__":  
    # Import main from src.auto_morse_bond.main
    from src.auto_morse_bond.main import main
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
    #       python3 auto_morse_bond_update.py -opt                                                     #
    #              or                                                                                  #
    #       python3 auto_morse_bond_update.py -man                                                     #
    #   in the terminal and the code will print out all command line option and terminate before any   #
    #   further analysis is done                                                                       #
    ####################################################################################################
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: use_GUI = False
    if use_GUI or '-gui' in commandline_inputs:
        from src.auto_morse_bond.GUI import auto_morse_bond_GUI
        #from src.auto_morse_bond.GUI_simple import auto_morse_bond_GUI
        print('\n\n\nauto_morse_bond_update is currently running in GUI mode, where all GUI inputs are intialized from auto_morse_bond_update.\n\n\n')
        auto_morse_bond_GUI(topofile, morsefile, parent_directory, newfile, mass_map, min_bond_length, coeffs2skip,
                            radius_specs, alpha_specs, alpha_scale, files2write, atom_style, zero_effected_xterms,
                            bondbreak_scale, ff_class, include_type_labels, include_rcut, GUI_zoom)
    else:
        # Run main auto_morse_bond classes/functions (and get modified m class)
        main(topofile, morsefile, parent_directory, newfile, mass_map, min_bond_length, coeffs2skip,
             radius_specs, alpha_specs, alpha_scale, files2write, atom_style, zero_effected_xterms,
             bondbreak_scale, ff_class, include_type_labels, include_rcut, commandline_inputs=commandline_inputs)