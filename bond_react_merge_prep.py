# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.4
April 2nd, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    ***********************************************************
    * Requirements:                                           *
    *   python 3.7+                                           *
    *                                                         *
    * Dependencies:                                           *
    *   correct atom-typing                                   *
    *                                                         *
    * Run methods:                                            *
    *   - IDE (manipulate variables and run from IDE)         *
    *   - GUI (manipulate variables and run from. Default     *
    *          settings set from this script)                 *
    *   - command line (python3 bond_react_merge_prep.py -man *
    *                   to get command line override options  *
    *                   and examples)                         *
    *                                                         *
    ***********************************************************


bond_react_merge_prep.py is a standalone code to format any
LAMMPS datafile to be compatable with bond_react_merge.py,
which assumes a comment style that all2lmp.py writes. This
code will read a LAMMPS .data file and a .cta (comment type
assigment file) and then finds the required comments and 
writes a new LAMMPS datafile that is fully compatable
with bond_react_merge.py. This code could be useful to
build bond/react templates of an already simulated system or
if the orginal datafiles being used to build a bond/react
simulaiton did not come from all2lmp.py.
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
#   use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have be from the IDE #
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
# Files to be read in (Change these variables to match your information). The topofile variable has to be #
# a LAMMPS .data/.dat and the cta_file has a specific format not discussed here (see the examples folder).#
###########################################################################################################
topofile = 'EXAMPLES/bond_react_merge_prep/smp20230330_IFF_rlxd0_testv3_8proc_1ts.data'
cta_file = 'EXAMPLES/bond_react_merge_prep/smp20230330_IFF_rlxd0_testv3_8proc_1ts.cta'


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
# Is a Python string variable in which to set the new output filename(s). The following options exist for #
# using the newfile string for setting the output file basenames:                                         #
#                                                                                                         #
#   if newfile starts with ':' or ends with ':'                                                           #
#     The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix  #
#     added to the file basename. The following are examples:                                             #
#       Suffix (newfile = ':_cta'  and topofile = 'detda.data')                                           #
#         basename = 'detda_cta', where the ':' character acts as a placeholder for the topofile          #
#                     basename.                                                                           #
#       Prefix (newfile = 'cta-:'  and topofile = 'detda.data')                                           #
#         basename = 'cta-detda', where the ':' character acts as a placeholder for the topofile          #
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
newfile = ':_cta'


###########################################################################################################
# Python string variable to set atom style format of the written data file. Currently supported atom      #
# styles: charge, molecular, full, angle, bond, atomic, dipole, dpd, or line. Example:                    #
#     atom_style = 'full' # will set atom style in written data file as full                              #
#                                                                                                         #
# Update atom_style as desired.                                                                           #
###########################################################################################################
atom_style = 'full'


###########################################################################################################
# There are cases were LAMMPS datafiles may be built where not all coeffIDs are used by atomIDs in the    #
# current LAMMPS datafile. bond_react_merge_prep.py will be default assign 'N/A' for atom CoeffIDs,       #
# 'N/A  N/A' for bond coeff IDs, etc ... for any coeffIDs in which no atomIDs currently used the CoeffIDs.#
# bond_react_merge.py is compatible with the 'N/A' comments where if they exist in the file read by       #
# bond_react_merge.py will just not assign any atomID lists to these coeffs, however for cleanliness      #
# purposes it may be desired to remove these coeffIDs. This option will remove the unused CoeffIDs.       #
#                                                                                                         #
# Update rm_unused_coeffs as desired.                                                                     #
###########################################################################################################
rm_unused_coeffs = True




##########################################################
### Import main from src.bond_react_merge_prep and run ###
##########################################################
if __name__ == "__main__":  
    # Import main from src.bond_react_merge_prep.main
    from src.bond_react_merge_prep.main import main
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
    #       python3 bond_react_merge_prep.py -opt                                                                             #
    #              or                                                                                                         #
    #       python3 bond_react_merge_prep.py -man                                                                             #
    #   in the terminal and the code will print out all command line option and terminate before any  further analysis is done#
    ###########################################################################################################################
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: use_GUI = False
    if use_GUI or '-gui' in commandline_inputs:
        print('\n\n\nbond_react_merge_prep is currently running in GUI mode, where all GUI inputs are intialized from bond_react_merge_prep.\n\n\n')
        from src.bond_react_merge_prep.GUI import bond_react_merge_prep_GUI 
        bond_react_merge_prep_GUI(topofile, cta_file, newfile, atom_style, parent_directory, rm_unused_coeffs, GUI_zoom)
    else:
        # Run main bond_react_merge_prep classes/functions (and get parameters class)
        m = main(topofile, cta_file, newfile, atom_style, parent_directory, rm_unused_coeffs, commandline_inputs=commandline_inputs)