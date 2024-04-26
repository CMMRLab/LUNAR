# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.5
April 1st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

This code is meant for IFF/IFF-R ONLY. It will add in the virtual
pi-electrons for IFF/IFF-R graphite type systems, will reset charges
and will also convert some coeff types to the cg1/cge types.

    *********************************************************
    * Requirements:                                         *
    *   python 3.7+                                         *
    *                                                       *
    * Run methods:                                          *
    *   - IDE (manipulate variables and run from IDE)       *
    *   - GUI (manipulate variables and run from. Default   *
    *          settings set from this script)               *
    *   - command line (python3 convert2graphite.py -man to *
    *                   get command line override options   *
    *                   and examples)                       *
    *                                                       *
    *********************************************************
"""



##############
### Inputs ###
##############
##############################################################################################################
# Python boolean variable to use GUI or not (True or False). All GUI settings are intialized from this script#
# so adjusting this script will set the default settings that the GUI will load with. Please NOTE that when  #
# using the GUI to check the console or terminal print outs every time a system is run through the code      #
# because that is where the import information will be displayed. Examples:                                  #
#   use_GUI = True  # Will launch GUI to manipulate variables and run from GUI                               #
#   use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have to be from the     #
#                     IDE or command line                                                                    #
#                                                                                                            #
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account for  #
# the different ways individuals have their OS display settings setup. A GUI_zoom = 100, means that default  #
# settings are used, whereas a GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80 means decrease#
# GUI by 20%. Examples:                                                                                      #
#   GUI_zoom = 100 # use default GUI size                                                                    #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                                #
#                                                                                                            #
# Update use_GUI and GUI_zoom as desired.                                                                    #
##############################################################################################################
use_GUI = True
GUI_zoom = 100


##############################################################################################################
# LAMMPS topofile with bonds in it (REQUIRED!). The LAMMPS datafile can have new LAMMPS "type labels" in     #
# it. The only supported atom styles are full, charge, or molecular. More atom styles can be coded if needed.#
##############################################################################################################
#topofile ='EXAMPLES/add_pi_electrons/Fullerene-C60_typed_PCFF-IFF.data'
topofile = 'EXAMPLES/add_pi_electrons/detda_typed_IFF.data'


##############################################################################################################
# Python list variable to specify atomTypeIDs to either convert to "cg1" charge or "cg1" consistent          #
# parameters or to add "cge" virtual pi electron to. Note if the read-in topofile has comments of the atom   #
# types in the Masses section you may also specify a string type of the atom type, which will be converted   #
# to the corresponding atomTypeID interanlly. You may also mix and match int TypeIDs and str atom type       #
# comments in types2convert. Examples:                                                                       #
#    types2convert = [4, 6]  # Will set atomTypeIDs 4 and 6 to operate on (Usually TypeIDs for "cp" or "c5") #
#    types2convert = [4]     # Will set atomTypeID 4 to operate on (Usually TypeIDs for "cp" or "c5")        #
#    types2convert = [4,'cp']# Will set atomTypeID 4 and 'cp' atom type 'cp' to operate on                   #
#    types2convert = ['cp']  # Will set atom type 'cp' to operate on                                         #
#    types2convert = ['cp','c5'] # Will set atom type 'cp' and 'c5' to operate on                            #
#                                                                                                            #
# Update types2convert as desired.                                                                           #
##############################################################################################################
types2convert = ['cp']


##############################################################################################################
# Is a Python string variable in which to set the new output filename(s). The following options exist for    #
# using the newfile string for setting the output file basenames:                                            #
#                                                                                                            #
#   if newfile starts with ':' or ends with ':'                                                              #
#     The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix     #
#     added to the file basename. The following are examples:                                                #
#       Suffix (newfile = ':_pi_electrons'  and topofile = 'detda.data')                                     #
#         basename = 'detda_pi_electrons', where the ':' character acts as a placeholder for the topofile    #
#                     basename.                                                                              #
#       Prefix (newfile = 'pi_electrons-:'  and topofile = 'detda.data')                                     #
#         basename = 'pi_electrons-detda', where the ':' character acts as a placeholder for the topofile    #
#                     basename.                                                                              #
#     Recommended usage: common and safe as this method is safe and output filename(s) carry similar names   #
#     to input filename(s).                                                                                  #
#                                                                                                            #
#   if newfile == 'ANYTEXT'                                                                                  #
#     The output filename(s) will be set as 'ANYTEXT'. For example:                                          #
#       newfile = 'cnt_renamed' and topofile = 'cnt.data')                                                   #
#         basename = 'cnt_renamed'                                                                           #
#     Recommended usage: occasional and safe as output filename(s) may no longer carry similar names as      #
#     output filename(s), but is safe as output filename(s) will not overwrite input filename(s)             #
#                                                                                                            #
#   if newfile == ''                                                                                         #
#     The output filename(s) will be the same as the input filename(s). This can be a dangerous option as    #
#     you may inadvertently overwrite a file and then must assume the file to contain certain information,   #
#     but it contains other information.                                                                     #
#     Recommended usage: rare and dangerous as this could lead to using incorrect files if not paying        #
#     attention very carefully.                                                                              #
#                                                                                                            #
# Update newfile as needed.                                                                                  #
##############################################################################################################
newfile = ':_pi_electrons'


##############################################################################################################
# Python string variable to set how to handle charge on "compounds" (materials that are not either pure      #
# graphite or CNT or fullerene, where there are other atom types bonded to the aromatic carbon atoms). The   #
# IFF charge model using the virtual pi-electrons was originally formulated only for pure graphitic systems  #
# and thus when trying to use the pi-electron charge with system that have different atom types other than   #
# aromatic carbons, the total system charge will become none-neutral. One method around this is to enforce   #
# the system charge to be neutral via the net_zero_charge option, however this method can have adverse       #
# effects such as completely removing charge from all none-aromatic atoms or arbitrarily scaling all none-   #
# aromatic atoms by the same charge value to achieve charge neutrality. This option provides a greater       #
# level of control to enforce charge neutrality. If your system is pure graphite or CNT or fullerene this    #
# option can be set to any supported options as none of the options effect this type of system. The          #
# following options exist:                                                                                   #
#   'none'                   which does not apply any constraint to how neighbor charges are handled. If     #
#                            there are any first neighbors bonded to the aromatic carbon atoms, this method  #
#                            may result in a none-charge neutral system.                                     #
#                                                                                                            #
#   'check-neighbors'        which checks that all neighbors are aromatic and only places the pi-electron if #
#                            all neighbors are aromatic.                                                     #
#                                                                                                            #
#   'accumulate-carbon'      which will accumulate any residual charge into the carbon atom, to ensure the   #
#                            "local grouping of atoms" stays charge neutral and thus the entire system will  #
#                            remain charge neutral.                                                          #
#                                                                                                            #
#   'accumulate-pi-electron' which will accumulate any residual charge into the pi-electron atoms, to ensure #
#                            the "local grouping of atoms" stays charge neutral and thus the entire system   #
#                            will remain charge neutral.                                                     #
#                                                                                                            #
#   'accumulate-neighbor'    which will accumulate any residual charge into the first neighboring atom, to   #
#                            ensure the "local grouping of atoms" stays charge neutral and thus the entire   #
#                            system will remain charge neutral.                                              #
#                                                                                                            #
# Update neighbor_charge_constraint as needed.                                                               #
##############################################################################################################
neighbor_charge_constraint = 'none'


##############################################################################################################
# Python boolean variable to convert all bonds, angles, dihedrals, impropers, and crossterms that have       #
# TypeIDs in types2convert to "cg1" consistent parameters (True or False). Examples:                         #
#    convert2cg1 = False  # Will keep all bonds, angles, dihedrals, impropers, and crossterms as is.         #
#    convert2cg1 = True   # Will find all bonds, angles, dihedrals, impropers, and crossterms that only have #
#                           "cg1" atom types (set by types2convert list) in them and update to "cg1"         #
#                           consistent parameters.                                                           #
#                                                                                                            #
# Update convert2cg1 as desired.                                                                             #
##############################################################################################################
convert2cg1 = False


##############################################################################################################
# Python boolean variable to add pi electrons to atom TypeIDs in types2convert (True or False). Usage  will  #
# result in the addition of virtual pi electron atom, new bonds, and new angles. The atom types, bond types, #
# and angle types will also be updated to be consistent with PCFF-IFF v1.5. Lastly, charges will also be     #
# updated and existing atoms and set to the new pi electrons. Examples:                                      #
#    add_pi_electrons = True  # Will add pi electrons and updated all coeffs and charges                     #
#    add_pi_electrons = False # Will not add pi electrons                                                    #
#                                                                                                            #
# Update add_pi_electrons as desired.                                                                        #
##############################################################################################################
add_pi_electrons = True


##############################################################################################################
# Python boolean variable to reset charge of atom types in types2convert to update charge of "cp" atom type  #
# to a "cg1" atom type (True or False). Examples:                                                            #
#    reset_charges = False  # Will keep charges as is from read in file.                                     #
#    reset_charges = True   # Will reset the charge of all atom types listed in types2convert to be          #
#                             consistent with "cg1" charge. *NOTE convert2cg1 will not automatically         #
#                             reset charge and add_pi_eletrons will automatically reset charge, so use       #
#                             accordingly.*                                                                  #
#                                                                                                            #
# Update reset_charges as desired.                                                                           #
##############################################################################################################
reset_charges = False


##############################################################################################################
# Python boolean variable to make the system charge neutral after adding in pi-electrons via add_pi_electrons#
# boolean variable or after resetting cg1 charge via reset_charges boolean variable. If the system is a      #
# simple graphite sheet(s) or CNT then the system charge will already be charge neutral and this option can  #
# be left as False. However, if you have functional groups attached to graphite sheet(s) or CNTs and the     #
# system was charged via another method or the system is amorphous/glassy carbon from ReaxFF there will be a #
# residual non-zero net charge of the system and this option shohuld be set as True. Examples:               #
#    net_zero_charge = False # Will leave system charge as is depending on add_pi_electrons or reset_charges #
#    net_zero_charge = True  # Will make system system charge net zero avoiding charge manipulation of cg1   #
#                              and cge atom types by adding a fixed charge value to all other atom types in  #
#                              the system                                                                    #
#                                                                                                            #
# Update net_zero_charge as desired.                                                                         #
##############################################################################################################
net_zero_charge =  False


##############################################################################################################
# Python string variable to set atom style format of the written data file. Currently supported atom styles: #
# charge, molecular, or full. Example:                                                                       #
#     atom_style = 'full' # will set atom style in written data file as full                                 #
#                                                                                                            #
# Update atom_style as desired.                                                                              #
##############################################################################################################
atom_style = 'full'


##############################################################################################################
# Python string variable type to set parent directory to store all new files this code will write. If the    #
# variable is left as an empty string, files will be written to the present working directory. Example to    #
# set parent_directory as present working directory:                                                         #
#     parent_directory = '' or parent_directory = '.'                                                        #
#                                                                                                            #
# Example to set parent_directory as with path from topofile:                                                #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile'                           #
#     files will be written to 'Furan_Resin/files/' since the 'topofile' string envokes this option.         #
#                                                                                                            #
# Example to set parent_directory as with path from topofile and build dirs from that location:              #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/NEWDIR'                    #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/NEWDIR' will build relative    #
#     directories                                                                                            #
#                                                                                                            #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/../NEWDIR'                 #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/../NEWDIR' will build          #
#     relative directories                                                                                   #
#                                                                                                            #
# Example to set parent_directory to a location inside of present working dirctory:                          #
#     parent_directory = 'inside'       # will save files in $pwd/inside                                     #
#     parent_directory = 'inside/test'  # will save files in $pwd/inside/test                                #
#                                                                                                            #
# Example to set parent_directory to a location outside of present working dirctory:                         #
#     parent_directory = '../'      # will save files in $pwd/../                                            #
#     parent_directory = '../test'  # will save files in $pwd/../test                                        #
#                                                                                                            #
# Update parent_directory as desired.                                                                        #
##############################################################################################################
parent_directory = 'topofile'


##############################################################################################################
# Python boolean variable to assign LAMMPS's new Type Labels section(s) to the written LAMMPS data file      #
# (True or False). NOTE This option ONLY works for datafiles with all2lmp style comments or datafiles that   #
# already have type labels defined in them or a combination of the two. Examples:                            #
#    include_type_labels = False # Will write the "standard" format of the LAMMPS data file                  #
#    include_type_labels = True  # Will write the Type Labels section(s) in the LAMMPS data file             #
#                                                                                                            #
# Update include_type_labels as desired (Default should be False).                                           #
##############################################################################################################
include_type_labels = False




#####################################################
### Import main from src.add_pi_electrons and run ###
#####################################################
if __name__ == "__main__":  
    # Import main from src.add_pi_electrons.main
    from src.add_pi_electrons.main import main
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
    #       python3 add_pi_electrons.py -opt                                                                                  #
    #              or                                                                                                         #
    #       python3 add_pi_electrons.py -man                                                                                  #
    #   in the terminal and the code will print out all command line option and terminate before any  further analysis is done#
    ###########################################################################################################################
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: use_GUI = False
    if use_GUI or '-gui' in commandline_inputs:
        from src.add_pi_electrons.GUI import add_pi_electrons_GUI 
        print('\n\n\nadd_pi_electrons is currently running in GUI mode, where all GUI inputs are intialized from add_pi_electrons.\n\n\n')
        add_pi_electrons_GUI(topofile, types2convert, atom_style, reset_charges, net_zero_charge, convert2cg1,
                             add_pi_electrons, parent_directory, newfile, include_type_labels, neighbor_charge_constraint, GUI_zoom)
    else:
        # Run main convert2graphite
        main(topofile, types2convert, atom_style, reset_charges, net_zero_charge, convert2cg1, add_pi_electrons,
             parent_directory, newfile, include_type_labels, neighbor_charge_constraint, commandline_inputs=commandline_inputs)







