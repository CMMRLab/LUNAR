# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
April 1st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


    *********************************************************
    * Requirements:                                         *
    *   python 3.7+                                         *
    *                                                       *
    * Dependencies:                                         *
    *   N/A                                                 *
    *                                                       *
    * Run methods:                                          *
    *   - IDE (manipulate variables and run from IDE)       *
    *   - GUI (manipulate variables and run from. Default   *
    *          settings set from this script)               *
    *   - command line (python3 atom_removal.py -man to get *
    *                   command line override options and   *
    *                   examples)                           *
    *                                                       *
    *                                                       *
    *                                                       *
    *********************************************************
"""

##############
### Inputs ###
##############
#################################################################################################################
# Python boolean variable to use GUI or not (True or False). All GUI settings are intialized from this script so#
# adjusting this script will set the default settings that the GUI will load with. Please NOTE that when using  #
# the GUI to check the console or terminal print outs every time a system is run through the code because that  #
# is where the import information will be displayed. Examples:                                                  #
#   use_GUI = True  # Will launch GUI to manipulate variables and run from GUI                                  #
#   use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have to be from the IDE    #
#                                                                                                               #
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account for the #
# different ways individuals have their OS display settings setup. A GUI_zoom = 100, means that default settings#
# are used, whereas a GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80 means decrease GUI by 20%.#
# Examples:                                                                                                     #
#   GUI_zoom = 100 # use default GUI size                                                                       #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                                   #
#                                                                                                               #
# Update use_GUI and GUI_zoom as desired.                                                                       #
#################################################################################################################
use_GUI = True
GUI_zoom = 100


#################################################################################################################
# Input variables for atom removal:                                                                             #
#    topofile = LAMMPS .data with atom styles as full, charge, or molecular.                                    #
#    atoms2remove = [1, 2, 3] list of atomIDs, TypeIDs, cluster-size or cluster-mass cutoff to remove from the  #
#                   system. NOTE for cluster-size or cluster-mass only a single value maybe supplied in list.   #
#   method                                                                                                      #
#     'atomIDs' will identify atoms based on atomIDs                                                            #
#     'TypeIDs' will identify atoms based on atom TypeIDs                                                       #
#     'cluster-mass' will perform cluster analysis and identify atoms based on a cutoff value set in            #
#                    atoms2remove, where all cluster mass less than or equal to the cutoff value will be removed#                                                                      #
#     'cluster-size' will perform cluster analysis and identify atoms based on a cutoff value set in            #
#                    atoms2remove, where all cluster sizes (number of atoms) less than or equal to the cutoff   #
#                    value will be removed                                                                      #
#################################################################################################################
topofile = 'EXAMPLES/atom_removal_test/pre1_typed_IFF_GT_merged.data'
method = 'atomIDs'
atoms2remove = [1, 2, 32, 33]


#################################################################################################################
# Is a Python string variable in which to set the new output filename(s). The following options exist for using #
# the newfile string for setting the output file basenames:                                                     #
#                                                                                                               #
#   if newfile starts with ':' or ends with ':'                                                                 #
#     The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix added  #
#     to the file basename. The following are examples:                                                         #
#       Suffix (newfile = ':_typed'  and topofile = 'detda.data')                                               #
#         basename = 'detda_typed', where the ':' character acts as a placeholder for the topofile basename.    #
#       Prefix (newfile = 'typed-:'  and topofile = 'detda.data')                                               #
#         basename = 'typed-detda', where the ':' character acts as a placeholder for the topofile basename.    #
#     Recommended usage: common and safe as this method is safe and output filename(s) carry similar names to   #
#     input filename(s).                                                                                        #
#                                                                                                               #
#   if newfile == 'ANYTEXT'                                                                                     #
#     The output filename(s) will be set as 'ANYTEXT'. For example:                                             #
#       newfile = 'detda_renamed' and topofile = 'detda.data')                                                  #
#         basename = 'detda_renamed'                                                                            #
#     Recommended usage: occasional and safe as output filename(s) may no longer carry similar names as output  #
#     filename(s), but is safe as output filename(s) will not overwrite input filename(s)                       #
#                                                                                                               #
#   if newfile == ''                                                                                            #
#     The output filename(s) will be the same as the input filename(s). This can be a dangerous option as you   #
#     may inadvertently overwrite a file and then must assume the file to contain certain information, but it   #
#     contains other information.                                                                               #
#     Recommended usage: rare and dangerous as this could lead to using incorrect files if not paying attention #
#     very carefully.                                                                                           #
#                                                                                                               #
# Update files newfile as needed.                                                                               #
#################################################################################################################
newfile = ':_rm_atoms'


#################################################################################################################
# Python string variable to set atom style format of the written data file. Currently supported atom styles:    #
# charge, molecular, full, angle, bond, atomic, dipole, dpd, or line. Example:                                  #
#     atom_style = 'full' # will set atom style in written data file as full                                    #
#                                                                                                               #
# Update atom_style as desired.                                                                                 #
#################################################################################################################
atom_style = 'full'


#################################################################################################################
# Python string variable type to set parent directory to store all new files this code will write. If the       #
# variable is left as an empty string, files will be written to the present working directory. Example to set   #
# parent_directory as present working directory:                                                                #
#     parent_directory = '' or parent_directory = '.'                                                           #
#                                                                                                               #
# Example to set parent_directory as with path from topofile:                                                   #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile'                              #
#     files will be written to 'Furan_Resin/files/' since the 'topofile' string envokes this option.            #
#                                                                                                               #
# Example to set parent_directory as with path from topofile and build dirs from that location:                 #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/NEWDIR'                       #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/NEWDIR' will build relative       #
#     directories                                                                                               #
#                                                                                                               #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/../NEWDIR'                    #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/../NEWDIR' will build relative    #
#     directories                                                                                               #
#                                                                                                               #
# Example to set parent_directory to a location inside of present working dirctory:                             #
#     parent_directory = 'inside'       # will save files in $pwd/inside                                        #
#     parent_directory = 'inside/test'  # will save files in $pwd/inside/test                                   #
#                                                                                                               #
# Example to set parent_directory to a location outside of present working dirctory:                            #
#     parent_directory = '../'      # will save files in $pwd/../                                               #
#     parent_directory = '../test'  # will save files in $pwd/../test                                           #
#                                                                                                               #
# Update parent_directory as desired.                                                                           #
#################################################################################################################
parent_directory = 'topofile'


#################################################################################################################
# Python boolean variable to assign LAMMPS's new Type Labels section(s) to the written LAMMPS data file (True or#
# False). NOTE This option ONLY works for datafiles with all2lmp style comments or datafiles that already have  #
# type labels defined in them or a combination of the two. Examples:                                            #
#    include_type_labels = False # Will write the "standard" format of the LAMMPS data file                     #
#    include_type_labels = True  # Will write the Type Labels section(s) in the LAMMPS data file                #
#                                                                                                               #
# Update include_type_labels as desired (Default should be False).                                              #
#################################################################################################################
include_type_labels = True




#################################################
### Import main from src.atom_removal and run ###
#################################################
if __name__ == "__main__":  
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
    #       python3 atom_removal.py -opt                                                               #
    #              or                                                                                  #
    #       python3 atom_removal.py -man                                                               #
    #   in the terminal and the code will print out all command line option and terminate before any   #
    #   further analysis is done                                                                       #
    ####################################################################################################
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: use_GUI = False
    if use_GUI or '-gui' in commandline_inputs:
        print('\n\n\natom_removal is currently running in GUI mode, where all GUI inputs are intialized from atom_removal.\n\n\n')
        from src.atom_removal.GUI import atom_removal_GUI 
        atom_removal_GUI(topofile, parent_directory, newfile, atom_style, atoms2remove, include_type_labels, GUI_zoom, method=method)
    else:
        # Import main from src.atom_removal.main and run
        from src.atom_removal.main import main
        main(topofile, parent_directory, newfile, atom_style, atoms2remove, include_type_labels, method=method, commandline_inputs=commandline_inputs)