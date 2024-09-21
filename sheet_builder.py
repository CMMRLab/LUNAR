# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
September 21st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

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
##################################################
### Global Inputs: apply to all modes sections ###
##################################################
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
# Python string variable type to set parent directory to store all new files this code will write. If the    #
# variable is left as an empty string, files will be written to the present working directory. Example to    #
# set parent_directory as present working directory:                                                         #
#     parent_directory = '' or parent_directory = '.'                                                        #
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
parent_directory = 'EXAMPLES/sheet_builder'


##############################################################################################################
# Python string variable type to set the mode to run the code in. Currently supported modes are:             #
#     'sheet'          which will generate sp2 hybridized sheets. All settings for this mode can be found in #
#                      "mode: sheet" section of this file.                                                   #
#                                                                                                            #
#     'symmetric-tube' which will generate sp2 hybridized symmetric tubes in either armchair or zigzag       #
#                      configuration. All settings for this mode can be found in "mode: symmetric-tube"      #
#                      section of this file.                                                                 #
#                                                                                                            #
#     'chiral-tube'   which will generate sp2 hybridized chiral tubes based.  All settings for this mode can #
#                     be found in "mode: chiral-tube" section of this file.                                  #
#                                                                                                            #
# Update mode as desired.                                                                                    #
##############################################################################################################
run_mode = 'sheet'


##############################################################################################################
# Python Boolean variable to control whether or not to find bonds based on interatomic distance searching.   #
# If find_bonds is True bonds will be found via interatomic distance searching and bonds will be written in  #
# the output LAMMPS datafile. If find_bonds is False no bonds will be found and no bonds will be written in  #
# the output LAMMPS datafile. When finding bonds the bonds can be set to be periodic or non-periodic, see    #
# the periodic_bonds flag. NOTE find_bonds must be True to generate an outputted .nta file.                  #
#                                                                                                            #
# Update find_bonds as desired.                                                                              #
##############################################################################################################
find_bonds = True


##############################################################################################################
# Python Boolean variable to control whether bonds will pass through the periodic boundary or not. If the    #
# Boolean is True, when bonds are determined via interatomic distances bonds will be found that pass through #
# a periodic boundary. If the Boolean is False no periodic bonds will be found.                              #
#                                                                                                            #
# Update periodic_bonds as desired.                                                                          #
##############################################################################################################
periodic_bonds = False


##############################################################################################################
# Python float variable to control the in-plane spacing of atoms when generating sheets or tubes based on    #
# the bond length between atoms. The bond length is supplied in units of Angstroms (usually around 1.42 A).  #
#                                                                                                            #
# Update bond_length as desired.                                                                             #
##############################################################################################################
bond_length = 1.42


##############################################################################################################
# Python dictionary with four keys 1,2,3,4 and four values that are strings. The strings will assign an atom #
# type or element in a specfic location of the ring. This will allow users to control the "pattern" of atoms #
# for generating sheets or tubes that contain multiple atom types or elements (for example in Boron-Nitride  #
# sheets or tubes). The "base unit" that is replicated to generate the sheets or tubes is describe in the    #
# figure below:                                                                                              #
#                                                                                                            #
#  type2___type3    replicated   2___3    2___3    types =  {1:'B', 2:'N', 3:'B', 4:'N'}   N___B    N___B    #
#      /   \        --------->  1/   \4__1/   \4   ------------------------------------>  B/   \N__B/   \N   #
# type1     type4                \___/    \___/    |                                       \___/    \___/    #
#                                2   3    2   3    |                                       N   B    N   B    #
#                                                  |                                                         #
#                                                  |types =  {1:'B', 2:'N', 3:'N', 4:'B'}   N___N    N___N   #
#                                                  +------------------------------------>  B/   \B__B/   \B  #
#                                                  |                                        \___/    \___/   #
#                                                  |                                        N   N    N   N   #
#                                                  |                                                         #
#                                                  |types =  {1:'N', 2:'B', 3:'B', 4:'N'}   B___B    B___B   #
#                                                  +------------------------------------>  N/   \N__N/   \N  #
#                                                  |                                        \___/    \___/   #
#                                                  |                                        B   B    B   B   #
#                                                  |                                                         #
#                                                  |types =  {1:'C', 2:'C', 3:'C', 4:'C'}   C___C    C___C   #
#                                                  +------------------------------------>  C/   \C__C/   \C  #
#                                                                                           \___/    \___/   #
#                                                                                           C   C    C   C   #
#                                                                                                            #
# The values of strings will set the atom type or element to repeat in the "pattern". An .nta file will be   #
# written such that the outputs of this code can be read into LUNAR/all2lmp.py to assign a specifc force     #
# field to the atomic positions generated. Thus the strings in the values should be set to the force field   #
# specfic atom types that you plan to apply using LUNAR/all2lmp.py. Here are a few examples for different    #
# force fields for a carbon system:                                                                          #
#    PCFF     -> types = {1:'cp',  2:'cp',  3:'cp',  4:'cp'}                                                 #
#    PCFF-IFF -> types = {1:'cg1', 2:'cg1', 3:'cg1', 4:'cg1'}                                                #
#    DREIDING -> types = {1:'C_R', 2:'C_R', 3:'C_R', 4:'C_R'}                                                #
#    ReaxFF   -> types = {1:'C',   2:'C',   3:'C',   4:'C'}                                                  #
# NOTE: that when generating inputs to ReaxFF that you set the actually element symbol and not classical     #
# force field atom types. Also note that when generating inputs to ReaxFF, you should set find_bonds = False #
# which will allow you to directly take the output LAMMPS datafile of this code and run in a LAMMPS ReaxFF   #
# simulation.                                                                                                #
#                                                                                                            #
# Additionaly IFF's pi-electrons can be added to a system by setting each type as 'atomtype|pi-electron'     #
# where the '|' character seperates the atomtype from the pi-electron type. For example:                     #
#                                                                                                            #
#     types = {1:'cg1|cge', 2:'cg1|cge', 3:'cg1|cge', 4:'cg1|cge'} will generate:                            #
#                                                     cge cge  cge cge                                       #
#         2___3    2___3            cg1___cg1          |   |    |   |                                        #
#        1/   \4__1/   \4   --->   cg1/   \cg1   ---> cg1-cg1--cg1-cg1                                       #
#         \___/    \___/              \___/            |   |    |   |                                        #
#         2   3    2   3            cg1   cg1         cge cge  cge cge                                       #
#                                                                                                            #
#     types = {1:'bbn|bbe', 2:'nbn|nbe', 3:'bbn|bbe', 4:'nbn|nbe'} will generate:                            #
#                                                     bbe nbe  bbe nbe                                       #
#         2___3    2___3            nbn___bbn          |   |    |   |                                        #
#        1/   \4__1/   \4   --->   bbn/   \nbn   ---> bbn-nbn--bbn-nbn                                       #
#         \___/    \___/              \___/            |   |    |   |                                        #
#         2   3    2   3            nbn   bbn         bbe nbe  bbe nbe                                       #
#                                                                                                            #
#                                                                                                            #
# Update types as desired.                                                                                   #
##############################################################################################################
#types = {1: 'C',       2:'C',        3:'C',        4:'C'}
#types = {1: 'B',       2: 'N',       3: 'B',       4: 'N'}
#types = {1: 'B',       2: 'N',       3: 'N',       4: 'B'}
#types = {1: 'cp',      2: 'cp',      3: 'cp',      4: 'cp'}
#types = {1: 'bbn|bbe', 2: 'nbn|nbe', 3: 'bbn|bbe', 4: 'nbn|nbe'}
types = {1: 'C', 2: 'C', 3: 'C', 4: 'C'}



##############################################################################################################
# Python dictionary to set charges based on atom type strings. The dictionary is set up as:                  #
#     charges = {'atomtype': charge }  # where 'atomtype' is a string defining the atom type and charge is a #
#                                      # float value defining the charge of the atom type.                   #
#                                                                                                            #
# Note that if an atom type is provided in the types dictionary, but a charge does not exist for that atom   #
# type in the charge dictionary, the charge will be set to zero.                                             #
#                                                                                                            #
# Update charges as desired.                                                                                 #
##############################################################################################################
charges = {'cg1':  0.250, 'cge': -0.125,  # IFF cg1/cge charge (graphitic-carbon/pi-electron)
           'bbn':  0.650, 'bbe': -0.100,  # IFF bbn/bbe charge (boron/pi-electron)
           'nbn': -0.250, 'nbe': -0.100,  # IFF nbn/nbe charge (nitride/pi-electron)
           'cp' :  0.000                  # PCFF aromatic carbon
           }

##############################################################################################################
# Python dictionary to set masses and element symbols. The dictionary is set up as:                          #
#     masses = {'atomtype': [mass, 'element'] }  # where 'atomtype' is a string defining the atom type and   #
#                                                # mass is a float and element is a string.                  #
# If a system attempts to be built and the mass is not in masses, the mass and element pair will attempt to  #
# be found from LUNAR/src/masses.py/mass_map dictionary and if it does not exist there the mass will be set  #
# to 0.000 and the elemet will be set to ''.                                                                 #
#                                                                                                            #
# Update masses as desired.                                                                                  #
##############################################################################################################
masses = {'cg1': [10.011150, 'C'],
          'cge': [1.0000000, 'H'],
          'bbn': [8.8110000, 'B'],
          'bbe': [1.0000000, 'H'],
          'nbn': [12.006700, 'N'],
          'nbe': [1.0000000, 'H'],
          'cp':  [12.011150, 'C'],
          'C':   [12.011000, 'C'],
          'B':   [10.810000, 'B'],
          'N':   [14.007000, 'N'],
          'O':   [15.999400, 'O'],
          'H':   [1.0080000, 'H']
          }


###################################################################################################
### Settings for adding different atom types (terminating of functional) to the sheets or tubes ###
###################################################################################################
##############################################################################################################
# functional_atoms                                                                                           #
#   An entry to supply a string to set how functional atoms are added to the sheet or tube and what the      #
#   functional group is. The string format is as follows:                                                    #
#     BondingType<MaxPercent>|Type1|Type2|TypeN, where "BondingType" is the atom type to add the function    #
#     group to, "MaxPercent" is a float or integer type to set the maximum percent of atoms to functionalize #
#     the "|" character seperates types, and the "TypeN" sets the atom to add.                               #
#                                                                                                            #
#     For example say types = {1:'C', 2:'C', 3:'C', 4:'C'} to generate a carbon sheet or nanotube and the    #
#     goal was to functionalize 5% of the carbon atoms with -OH functional group. Then the functional_atoms  #
#     string would be 'C<5>|O|H', which would randomly add the -OH functional group to 5% of the C atoms.    #
#                                                                                                            #
#     Additionaly the functional_atoms string can handle multiple BondingType's by seperating them with the  #
#     ";" character. So the generalized functional_atoms string becomes:                                     #
#       BondingTypeA<MaxPercentA>|TypeA1|TypeAN; BondingTypeB<MaxPercentB>|TypeB1|TypeBN; ...                #
#                                                                                                            #
#       For example say types = {1:'B', 2:'N', 3:'B', 4:'N'} to generate a Boron-Nitride sheet or tube with  #
#       alternating B/N atoms and the goal was to functionalize 10% of the Boron atoms with -OH functional   #
#       group and to functionalize 20% of the Nitride atoms with -H functional group. Then the functional_   #
#       atoms string would be 'B<10>|O|H; N<20>|H', which would randomly add the -OH functional group to 10% #
#       of the B atoms and add the -H functional group to 20% of the N atoms.                                #
#                                                                                                            #
#     All examples above will place the atoms in a line along the orthagonal direction from the surface of   #
#     tube, but say we wanted to added a functional group that resembles an epoxide ring (3 member ring with #
#     two carbons and 1 oxygen). Then we can add a "|" character to the end of the functional_atoms string.  #
#     This method currently only works for adding a single atom functional group like oxygen to the sheets   #
#     or tubes.                                                                                              #
#                                                                                                            #
#       For example say types = {1:'C', 2:'C', 3:'C', 4:'C'} to generate a carbon sheet or nanotube and the  #
#       goal was to functionalize 30% of the carbon atoms with the epoxide ring oxygen. Then the             #
#       functional_atoms string would be 'C<30>|O|', where the last character is the "|" character. This     #
#       will tell the code to find a first neighbor from the random atom and center the oxygen atom between  #
#       the first neighbor and itself. Finally, add two bonds to create the epoxide type ring. Note that     #
#       each time the oxygen atom is added, it functionalizes two carbon atoms at a time. So say the sheet   #
#       or tube had 100 carbon atoms and the functionalization MaxPercent was set to 30%, then only 15       #
#       oxygen atoms will be added (not 30).                                                                 #
#                                                                                                            #
#   If the functional_atoms entry/string is left blank, this option will not be envoked. Additionally,       #
#   this option requires find_bonds to be True.                                                              #
#                                                                                                            #
# functional_seed                                                                                            #
#   An entry to supply a seed to the random number generate to define the random atoms the functional groups #
#   will be added to. If the seed value is set to ZERO, the current system time from your computer is used   #
#   to provide a seed to the random number generator.                                                        #
#                                                                                                            #
# Update functional_atoms and functional_seed as desired.                                                    #
##############################################################################################################
#functional_atoms = 'cg1<10>|O|; C<50>|O|'
functional_atoms = ''
functional_seed = 12345


##############################################################################################################
# terminating_atoms                                                                                          #
#   An entry to supply a string to set how terminating atoms are added to the sheet or tube and what the     #
#   termanting atoms are. This option requires that periodic_bonds is False, as this creates open valences   #
#   on the "end" atoms of the sheet or tube. The string format is as follows:                                #
#     BondingType|Type1|Type2|TypeN, where "BondingType" is the atom type to add the terminating atoms to,   #
#     the "|" character seperates types, and the "TypeN" sets the atom to add.                               #
#                                                                                                            #
#     For example say types = {1:'C', 2:'C', 3:'C', 4:'C'} to generate a carbon sheet or nanotube that is    #
#     not periodically bonded and the goal was to terminate the open valences on the "edges" of the sheet or #
#     tube with the -OH functional group. Then the termanting_atoms string would be 'C|O|H', which would     #
#     terminate all "edge" atoms with the -OH group.                                                         #
#                                                                                                            #
#     Additionaly the terminating_atoms string can handle multiple BondingType's by seperating them with the #
#     ";" character. So the generalized termanting_atoms string becomes:                                     #
#       BondingTypeA|TypeA1|TypeAN; BondingTypeB|TypeB1|TypeBN; BondingTypeC|TypeC1|TypeCN; ...              #
#                                                                                                            #
#       For example say types = {1:'B', 2:'N', 3:'B', 4:'N'} to generate a Boron-Nitride sheet or tube with  #
#       alternating B/N atoms that was not periodically bonded and the goal was to termanate the Boron atoms #
#       with -H and to termanate the Nitride atoms with -OH, then the terminating_atoms string would be      #
#       'B|H; N|O|H', which would termanate the "edge" Boron atoms with an -H and the "edge" Nitride atoms   #
#       with a -OH group.                                                                                    #
#                                                                                                            #
#   If the terminating_atoms entry/string is left blank, this option will not be envoked. Additionally, this #
#   option requires find_bonds to be True.                                                                   #
#                                                                                                            #
# Update terminating_atoms as desired.                                                                       #
##############################################################################################################
#terminating_atoms = 'cg1|H; C|H; cg1|H'
terminating_atoms = ''



#####################
### mode: 'sheet' ###
#####################
##############################################################################################################
# Python string variable type to set output basename of outputs for when running in 'sheet' mode.            #
#                                                                                                            #
# Update sheet_basename as desired.                                                                          #
##############################################################################################################
sheet_basename = 'sheet'


##############################################################################################################
# Python string variable type to set the plane to generate the sheet on. The available planes to generate    #
# sheets on are:                                                                                             #
#     'xy' which will generate a sheet on the XY plane                                                       #
#     'xz' which will generate a sheet on the XZ plane                                                       #
#     'yz' which will generate a sheet on the YZ plane                                                       #
# The sheet size and shape will depend on the length_in_normal/length_in_edgetype and sheet_edgetype.        #
#                                                                                                            #
# Update plane as desired.                                                                                   #
##############################################################################################################
plane = 'xy'


##############################################################################################################
# Python string variable type to set the edge type of the sheet. The available edge types are:               #
#     'armchair' which will generate a sheet in the armchair directions, with the length in the armchair     #
#                direction set by the length_in_edgetype variable.                                           #
#     'zigzag'   which will generate a sheet in the zigzag directions, with the length in the zigzag         #
#                direction set by the length_in_edgetype variable.                                           #
#                                                                                                            #
# Update sheet_edgetype as desired.                                                                          #
##############################################################################################################
sheet_edgetype = 'armchair'


##############################################################################################################
# Python float variables to control the size of the sheet. The length_in_edgetype will set the length in the #
# edge type direction set by the sheet_edgetype variable and the length_in_perpendicular sets the            #
# perpendicular length from the edge type. Here is the meaning of length_in_edgetype and                     #
# length_in_perpendicular for different plane and sheet_edgetype settings:                                   #
#     plane = 'xy' and sheet_edgetype = 'zigzag'                                                             #
#       length_in_perpendicular   = length in X-direction                                                    #
#       length_in_edgetype = length in Y-direction                                                           #
#                                                                                                            #
#     plane = 'xz' and sheet_edgetype = 'zigzag'                                                             #
#       length_in_perpendicular   = length in X-direction                                                    #
#       length_in_edgetype = length in Z-direction                                                           #
#                                                                                                            #
#     plane = 'yz' and sheet_edgetype = 'zigzag'                                                             #
#       length_in_perpendicular   = length in Y-direction                                                    #
#       length_in_edgetype = length in Z-direction                                                           #
#                                                                                                            #
#     plane = 'xy' and sheet_edgetype = 'armchair'                                                           #
#       length_in_perpendicular   = length in Y-direction                                                    #
#       length_in_edgetype = length in X-direction                                                           #
#                                                                                                            #
#     plane = 'xz' and sheet_edgetype = 'armchair'                                                           #
#       length_in_perpendicular   = length in X-direction                                                    #
#       length_in_edgetype = length in Z-direction                                                           #
#                                                                                                            #
#     plane = 'yz' and sheet_edgetype = 'armchair'                                                           #
#       length_in_perpendicular   = length in Z-direction                                                    #
#       length_in_edgetype = length in Y-direction                                                           #
#                                                                                                            #
# The units of length_in_perpendicular and length_in_edgetype is in angstroms.                               #
#                                                                                                            #
# Update length_in_perpendicular and length_in_edgetype as desired.                                          #
##############################################################################################################
length_in_perpendicular = 20.0
length_in_edgetype = 10.0


##############################################################################################################
# Python string variable type to set how sheets are stacked if sheet_nlayers is greater than one. The        #
# following stacking sequences are available:                                                                #
#  'AA'                                                                                                      #
#  'AB'                                                                                                      #
#  'ABC'                                                                                                     #
# where the schematic below shows AA, AB, and ABC stacking for six sheets:                                   #
#                                                                                                            #
#      AA-stacking        AB-stacking        ABC-stacking                                                    #
#      ___________        ___________        ___________                                                     #
#      ___________           ___________         ___________                                                 #
#      ___________        ___________                ___________                                             #
#      ___________           ___________     ___________                                                     #
#      ___________        ___________            ___________                                                 #
#      ___________           ___________             ___________                                             #
#                                                                                                            #
# The "shift" in AB- and ABC-stacking is the bond_length. Also note for AB- and ABC-stacking, the            #
# simulation cell will be set by the "top layer" and the rest of the layers will be wrapped to fit within a  #
# periodic simulation cell.                                                                                  #
#                                                                                                            #
# Update stacking as desired.                                                                                #
##############################################################################################################
stacking = 'AB'


##############################################################################################################
# Python float variable to control the out-of-plane spacing of atoms when multiple layers are added. The     #
# layer spacing is supplied in units of Angstroms (usually around 3.354 A).                                  #
#                                                                                                            #
# Update sheet_layer_spacing as desired.                                                                     #
##############################################################################################################
sheet_layer_spacing = 3.354


##############################################################################################################
# Python int variable to control the number of layers to generated. All models will be centered about (0, 0, #
# 0) no matter the number of layers generated.                                                               #
#                                                                                                            #
# Update sheet_nlayers as desired.                                                                           #
##############################################################################################################
sheet_nlayers = 1




##############################
### mode: 'symmetric-tube' ###
##############################
##############################################################################################################
# Python string variable type to set output basename of outputs for when running in 'symmetric-tube' mode.   #
#                                                                                                            #
# Update symmetric_tube_basename as desired.                                                                 #
##############################################################################################################
symmetric_tube_basename = 'tube-symmetric'


##############################################################################################################
# Python string variable type to set which axis the tube axis is aligned with. The available axis are:       #
#     'x' which will build the tube in the X-direction                                                       #
#     'y' which will build the tube in the Y-direction                                                       #
#     'Z' which will build the tube in the Z-direction                                                       #
#                                                                                                            #
# Update axis as desired.                                                                                    #
##############################################################################################################
symmetric_tube_axis = 'z'


##############################################################################################################
# Python string variable type to set the edge type of the tube. The available edge types are:                #
#     'armchair' which will generate a tube with the "circumference" being in the armchair direction         #
#     'zigzag'   which will generate a tube with the "circumference" being in the zigzag direction           #
#                                                                                                            #
# Update tube_edgetype as desired.                                                                           #
##############################################################################################################
tube_edgetype = 'armchair'


##############################################################################################################
# Python float variables to control the size of the tube. The symmetric_length variable controls the length  #
# of the tube in the axial direction and diameter controls the diameter of the tube.                         #
#                                                                                                            #
# Update symmetric_length and diameter as desired.                                                           #
##############################################################################################################
symmetric_length = 25.0
diameter = 12.0


##############################################################################################################
# Python float variable to control the out-of-plane spacing of atoms when multiple tubes are added (i.e.     #
# when symmetric_ntubes is greater then 1). The layer spacing is supplied in units of Angstroms (usually     #
# around 3.354 A)                                                                                            #
#                                                                                                            #
# Update tube_layer_spacing as desired.                                                                      #
##############################################################################################################
tube_layer_spacing = 3.354


##############################################################################################################
# Python int variable to control the number of tubes to generate. All models will be centered about (0, 0, 0)#
# no matter the number of tubes generated. When symmetric_ntubes is greater then 1, the diameter variable    #
# will control the inner most tube diameter and the N-outer layer tube diameters will be controlled by the   #
# symmetric_ntubes and tube_layer_spacing varaibles.                                                         #
#                                                                                                            #
# Update ntubes as desired.                                                                                  #
##############################################################################################################
symmetric_ntubes = 1




###########################
### mode: 'chiral-tube' ###
###########################
##############################################################################################################
# Python string variable type to set output basename of outputs for when running in 'chiral-tube' mode.      #
#                                                                                                            #
# Update chiral_tube_basename as desired.                                                                    #
##############################################################################################################
chiral_tube_basename = 'tube-chiral'


##############################################################################################################
# Python string variable type to which axis the tube axis is aligned with. The available axis are:           #
#     'x' which will build the tube in the X-direction                                                       #
#     'y' which will build the tube in the Y-direction                                                       #
#     'Z' which will build the tube in the Z-direction                                                       #
#                                                                                                            #
# Update axis as desired.                                                                                    #
##############################################################################################################
chiral_tube_axis = 'z'


##############################################################################################################
# Python float variables to control the length of a chiral tube. NOTE the chiral indices m and n control the #
# "base-tube" length where the chiral_length sets the desired length of a chiral tube, where the "base-tube" #
# will be replicated enough times to get a tube close in length to chiral_length.                            #
#                                                                                                            #
# Update chiral_length as desired.                                                                           #
##############################################################################################################
chiral_length = 20.0


##############################################################################################################
# Python int variables to control the chiral indices of m and n. All physical properties of a chiral single  #
# wall tube depend on the chiral indices. Symmetric tubes can be generated when:                             #
#    m = n, which generates an armchair tube                                                                 #
#    m = 0, which generates a zigzag tube                                                                    #
# If it is desired to create armchair or zigzag tubes this code can be ran in mode = 'symmetric-tube', which #
# provides easier control over the tube length and diameter. The 'symmetric-tube' mode also allows multiwall #
# nanotubes to be generated with easy. Thus it is recommend only to use the 'chiral-tube' mode if you want   #
# to generate a chiral tube and use the 'symmetric-tube' to generate tubes in either armchair or zigzag      #
# configurations.                                                                                            #
#                                                                                                            #
# Update m and n as desired.                                                                                 #
##############################################################################################################
n = 5
m = 10




##################################################
### Import main from src.sheet_builder and run ###
##################################################
if __name__ == "__main__":  
    # Import main from src.sheet_builder.main
    from src.sheet_builder.main import main
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
    #       python3 sheet_builder.py -opt                                                                                     #
    #              or                                                                                                         #
    #       python3 sheet_builder.py -man                                                                                     #
    #   in the terminal and the code will print out all command line option and terminate before any  further analysis is done#
    ###########################################################################################################################
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: use_GUI = False
    if use_GUI or '-gui' in commandline_inputs:
        from src.sheet_builder.GUI import sheet_builder_GUI 
        print('\n\n\nsheet_builder is currently running in GUI mode, where all GUI inputs are intialized from sheet_builder.\n\n\n')
        sheet_builder_GUI(sheet_basename, symmetric_tube_basename, chiral_tube_basename, run_mode, parent_directory, length_in_perpendicular, length_in_edgetype,
                          sheet_edgetype, types, bond_length, sheet_layer_spacing, sheet_nlayers, stacking, plane, tube_edgetype, tube_layer_spacing,
                          symmetric_ntubes, symmetric_length, diameter, n, m, chiral_length, symmetric_tube_axis, chiral_tube_axis, find_bonds, periodic_bonds,
                          charges, masses, functional_seed, functional_atoms, terminating_atoms, GUI_zoom)
    else:
        main(sheet_basename, symmetric_tube_basename, chiral_tube_basename, run_mode, parent_directory, length_in_perpendicular, length_in_edgetype, sheet_edgetype, types,
             bond_length, sheet_layer_spacing, sheet_nlayers, stacking, plane, tube_edgetype, tube_layer_spacing, symmetric_ntubes, symmetric_length, diameter, n, m,
             chiral_length, symmetric_tube_axis, chiral_tube_axis, find_bonds, periodic_bonds, charges, masses, functional_seed, functional_atoms, terminating_atoms, 
             commandline_inputs=commandline_inputs)