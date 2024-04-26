# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.7
April 11th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    ************************************************************
    * Requirements:                                            *
    *   python 3.7+                                            *
    *                                                          *
    *   If unwrapping atoms via image flags, the image flags   *
    *   MUST BE CONSISTENT.                                    *
    *                                                          *
    *  ALL READ IN FILES MUST HAVE ALL ENERGY COEFFS MERGED,   *
    *  since no resetting of coeffIDs are performed. All energy*
    *  parameters will be mapped from the last file in the     *
    *  files dictionary and all coeffIDs will be left as is.   *
    *  bond_react_merge.py can be used to merge all coeffs,    *
    *  even if you are not planning on performing a bond/react *
    *  simulation by just using the dataN file-tag to generate *
    *  the inputs to this code.                                *
    *                                                          *
    * Run methods:                                             *
    *   - IDE (manipulate variables and run from IDE)          *
    *   - GUI (manipulate variables and run from. Default      *
    *          settings set from this script)                  *
    *   - command line (python3 cell_builder.py -man to get    *
    *                   command line override options and      *
    *                   examples)                              *
    *                                                          *
    ************************************************************
"""
###############################################################################################################
# Python boolean variable to use GUI or not (True or False). All GUI settings are intialized from this script #
# so adjusting this script will set the default settings that the GUI will load with. Please NOTE that when   #
# using the GUI to check the console or terminal print outs every time a system is run through the code       #
# because that is where the import information will be displayed. Examples:                                   #
#   use_GUI = True  # Will launch GUI to manipulate variables and run from GUI                                #
#   use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have to be from the      #
#                     IDE or command line                                                                     #
#                                                                                                             #
# Python variable to set the maxfiles the GUI will load with (maxfiles must be an integer value). Please NOTE #
# that the GUI supports file overloading when running, so maxfiles only sets the intial number of files and   #
# more will be dynamically added if more files are provide while the GUI is loaded. Examples:                 #
#   maxfiles = 4 # GUI will load with 4 files slots and the GUI will automatically add more when "overloaded" #
#   maxfiles = 8 # GUI will load with 8 files slots and the GUI will automatically add more when "overloaded" #
#                                                                                                             #
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account for   #
# the different ways individuals have their OS display settings setup. A GUI_zoom = 100, means that default   #
# settings are used, whereas a GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80 means decrease #
# GUI by 20%. Examples:                                                                                       #
#   GUI_zoom = 100 # use default GUI size                                                                     #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                                 #
#                                                                                                             #
# Update use_GUI, GUI_zoom, and maxfiles as desired.                                                          #
###############################################################################################################
use_GUI = True
GUI_zoom = 100
maxfiles = 4


###############################################################################################################
# Python dictionary variable to hold all files to generate the system. The keys of the dictionary are the     #
# filenames (with or without file paths) and the values are integer numbers representing the quantiy of each  #
# file. The values that represent the quantity values of each file can set the stoichiometric ratio of all    #
# files (NOTE you only have to specify the lowest quanity value and then use the duplicate varaible to scale  #
# up the number of atoms). Dictionary format:                                                                 #
#     files = {filename (with path if desired) : qty }                                                        #
#                                                                                                             #
# Examples:                                                                                                   #
#    files = {'detda_typed_IFF_merged.data' : 1, # read-in 1 file of detda molecule to mix with 2 of dgeba    #
#             'dgeba_typed_IFF_merged.data' : 2} # molecules, thus setting the stoichiometric ratio.          #
#                                                                                                             #
# Update files as desired.                                                                                    #
###############################################################################################################
files = {'EXAMPLES/cell_builder/system1_detda_typed_IFF_merged_qty=1.data' : 1,
         'EXAMPLES/cell_builder/system1_dgebf_typed_IFF_merged_qty=2.data' : 2,
         'EXAMPLES/cell_builder/system1_detda_dgebf_grouping_example.script': 0}


###############################################################################################################
# Python string to define how to handle the force field between multiple LAMMPS datafiles, where the following#
# options exist:                                                                                              #
#   'none'  which assumes all LAMMPS coefficient types are the same in all the read in files. and applies no  #
#           offset to the files as they are being used to generate a large molecular system. This option      #
#           should be used if your files where processed with LUNAR/bond_react_merge.py to ensure that the    #
#           force field between the output system is consistent with the reaction templates.                  #
#  'merge'  which applies the merging processes present LUNAR/bond_react_merge.py to merge all coefficient    #
#           types amongst all read in files. This option requires that all LAMMPS datafiles have the LUNAR    #
#           all2lmp.py style of comments.                                                                     #
#  'offset' which applies an offset to each coefficient type in each file as it is read into cell_builder.py. #
###############################################################################################################
force_field_joining = 'none'


###############################################################################################################
# Python int variable to set the system size by duplicating each file Ntimes. The file quantities set by the  #
# files dictionary will set the stoichiometric ratio of the system and the duplicate variable will be the     #
# multiplier increase the system size while maintain the correct stoichiometric ratio. The minimum interger   #
# value is 1 and the maximum is as large as a system as you wish to randomize and build. Examples:            #
#     distance_scale = 10   #  Will increase number of atoms of read in files 10  times.                      #
#     distance_scale = 100  #  Will increase number of atoms of read in files 100 times.                      #
#                                                                                                             #
# Update duplicate as desired.                                                                                #
###############################################################################################################
duplicate = 10


###############################################################################################################
# Python float variable to set the size of each subcell that a molecule will be placed in. Each subcell       #
# dimension is set that it may hold the largest molecule from the read-in files in any orientation without    #
# passing the subcell dimensions. This is to avoid running molecules into other molecules when applying random#
# rotations. distance_scale sets a multiplier value to increase or decrease (at your own risk) the dimension  #
# of each subcell and ultimately the spacing of all molecules. Examples:                                      #
#     distance_scale = 1.0 # Each subcell is the minimum size to hold the largest molecule in th read-in files#
#     distance_scale = 1.2 # Each subcell is 20% larger then required and will space molecules apart more.    #
#     distance_scale = 0.8 # Each subcell is 20% smaller then required and will reduce the space between      #
#                            molecules (may create overlapped atoms). If max_rotations are set to ZEROES all  #
#                            molecules will keep the same orientation and thus can be "packed" together with  #
#                            a reduced risk of overlapping atoms and likely is the only case where you would  #
#                            set this less then 1.0.                                                          #
#                                                                                                             #
# *NOTE: sometimes setting distance_scale to 1.0 may still cause molecules to have very high intial energies  #
# that may cause LAMMPS to crash with a "Out of range atoms" ERROR. If this happens, increase the distance    #
# scale variable to ensure that the intial energy of the system are not large during the first few timesteps. #
#                                                                                                             #
# Update distance_scale as desired (Default should be 1.2).                                                   #
###############################################################################################################
distance_scale = 1.2


###############################################################################################################
# Is a Python Boolean (True or False) to group the monomers in the files dictionary locally and then randomly #
# place and rotate the local grouping of monomers. This is useful if you want to keep close proximity of a    #
# harder and resin molecules to allow for quicker crosslinking simulations or if you would like to add in the #
# artifact that every harder and resin start in a relatively close position  to help achieve 100% conversion. #
#                                                                                                             #
# Update group_monomers_locally as desired.                                                                   #
###############################################################################################################
group_monomers_locally = False


###############################################################################################################
# Is a Python string variable in which to set the cubic lattice domain. The lattice points will always be     #
# generated about the 0, 0, 0 position in x, y, and z. However, the number of lattice points in the x, y, and #
# z directions can be set by the user or set to be cubic. The following options are available:                #
#   'cubic'        which automatically determines the number of lattice points required based on the qty of   #
#                  files and the duplicate variable.                                                          #
#	'Ni x Nj x Nk' where 'Ni' is the number of lattice points in the x-direction, 'Nj' is the number of       #
#                  lattice points in the y-direction, and 'Nk' is the number of lattice points in the         #
#                  z-direction. Please note the following about 'Ni x Nj x Nk':                               #
#                      If group_monomers_locally is False, 'Ni x Nj x Nk' must be greater than                #
#                                                           duplicate*sum(qty of all files), to allow for     #
#                                                           enough lattice points. If there is not enough     #
#                                                           lattice points, the code will exit with an ERROR. #
#                      If group_monomers_locally is True, 'Ni x Nj x Nk' must be greater than duplicate value,#
#                                                          to allow for enough lattice points. If there is not#
#                                                          enough lattice points, the code will exit with an  #
#                                                          ERROR. Additionally, by default during the initial #
#                                                          grouping of monomers will occur with a 'cubic'     #
#                                                          lattice.                                           #
#                                                                                                             #
# Examples:                                                                                                   #
#     domain = 'cubic' # will generate enough lattice points in a cube to place molecules on                  #
#     domain = '2x2x8' # will generate 2x2x8 number of lattice points, 2 in x-dir, 2 in y-dir, and 8 in z-dir #
#     domain = '6x8x2' # will generate 6x8x2 number of lattice points, 6 in x-dir, 8 in y-dir, and 2 in z-dir #
#                                                                                                             #
# Update domain as desired.                                                                                   #
###############################################################################################################
domain = 'cubic'


###############################################################################################################
# Python dictionary to store and adjust the maximum degree of rotation that any molecule or system will       #
# undergo during the random selection of phi, theta, and psi values (Degrees). Format:                        #
#     max_rotations = {'x': float or int, 'y': float or int, 'z': float or int} # Keys are strings of x, y, z #
#                                                                                                             #
# Examples:                                                                                                   #
#     max_rotations = {'x': 90,  'y': 180, 'z': 270} # max rotation about x=90, max rotation about y=180, ... #
#     max_rotations = {'x': 90,  'y': 45,  'z': 100} # max rotation about x=90, max rotation about y=45, ...  #
#     max_rotations = {'x':  0,  'y': 0,   'z': 0  } # keeps all molecules in the orginal orientations.       #
#                                                                                                             #
# Update max_rotations as desired (Default should be {'x': 360, 'y': 360, 'z': 360}).                         #
###############################################################################################################
max_rotations = {'x': 360.0, 'y': 360.0, 'z': 360.0}


###############################################################################################################
# Python positive integer value or ZERO to set a seed for generating random numbers to account for            #
# reproducibility concerns of generating random intial positions. If the seed value is zero the seed will     #
# default to using the current system time. Examples:                                                         #
#     seed = 0     # not setting a seed and will default to current system time                               #
#     seed = 12345 # setting an exact seed to allow for reproducibile postions and orientations               #
#                                                                                                             #
# Update seed as desired.                                                                                     #
###############################################################################################################
seed = 12345


###############################################################################################################
# Python boolean variable to unwrap atoms via image flags when reading in the the files. This functionality   #
# assumes consistent image flags. The default should be to always unwrap atoms via image flags. There are a   #
# few scenerios that the atoms can stay "wrapped"; the datafiles have non-periodic molecules or special ReaxFF#
# cases where bonds will be infered via inter-atomic distnaces of a periodic system. Examples:                #
#    unwrap_atoms_via_image_flags = True  # Will unwrap periodic atoms when reading in file.                  #
#    unwrap_atoms_via_image_flags = False # Will not unwrap periodic atoms when reading in file.              #
#                                                                                                             #
# Update unwrap_atoms_via_image_flags as desired (Default should be True).                                    #
###############################################################################################################
unwrap_atoms_via_image_flags = True


###############################################################################################################
# Python string variable type to set new file name(s).                                                        #
# Example for newfile = 'cell_building'                                                                       #
#     filename(s) =   cell_building.data, cell_building.mol2, and cell_building.txt will be generated         #
# Update newfile as desired.                                                                                  #
###############################################################################################################
newfile = 'system1_cell_replicate'


###############################################################################################################
# Python string variable to set atom style format of the written data file. Currently supported atom styles:  #
# charge, molecular, full, angle, bond, atomic, dipole, dpd, or line. Example:                                #
#     atom_style = 'full' # will set atom style in written data file as full                                  #
#                                                                                                             #
# Update atom_style as desired.                                                                               #
###############################################################################################################
atom_style = 'full'


###############################################################################################################
# Python string variable type to set parent directory to store all new files this code will write. If the     #
# variable is left as an empty string, files will be written to the present working directory. Examples:      #
# Example to set parent_directory as present working directory:                                               #
#     parent_directory = '' or parent_directory = '.'                                                         #
#                                                                                                             #
# Example to set parent_directory as with path from 1st topofile in the files dictionary:                     #
#     files = {Furan_Resin/files/furan.data': 2} and  parent_directory = 'topofile'                           #
#     files will be written to 'Furan_Resin/files/' since the 'topofile' string envokes this option.          #
#                                                                                                             #
# Example to set parent_directory as with path from 1st topofile in the files dictionary and build dirs from  #
# that location:                                                                                              #
#     files = {Furan_Resin/files/furan.data' : 2} and  parent_directory = 'topofile/NEWDIR'                   #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/NEWDIR' will build relative     #
#     directories                                                                                             #
#                                                                                                             #
#     files = {Furan_Resin/files/furan.data' : 1} and  parent_directory = 'topofile/../NEWDIR'                #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/../NEWDIR' will build relative  #
#     directories                                                                                             #
#                                                                                                             #
# Example to set parent_directory to a location inside of present working dirctory:                           #
#     parent_directory = 'inside'       # will save files in $pwd/inside                                      #
#     parent_directory = 'inside/test'  # will save files in $pwd/inside/test                                 #
#                                                                                                             #
# Example to set parent_directory to a location outside of present working dirctory:                          #
#     parent_directory = '../'      # will save files in $pwd/../                                             #
#     parent_directory = '../test'  # will save files in $pwd/../test                                         #
#                                                                                                             #
# Update parent_directory as desired.                                                                         #
###############################################################################################################
parent_directory = 'topofile'


###############################################################################################################
# Python boolean variable to assign LAMMPS's new Type Labels section(s) to the written LAMMPS data file (True #
# or False). NOTE This option ONLY works for datafiles with all2lmp style comments or datafiles that already  #
# have type labels defined in them or a combination of the two. Examples:                                     #
#    include_type_labels = False # Will write the "standard" format of the LAMMPS data file                   #
#    include_type_labels = True  # Will write the Type Labels section(s) in the LAMMPS data file              #
#                                                                                                             #
# Update include_type_labels as desired (Default should be False).                                            #
###############################################################################################################
include_type_labels = False


###############################################################################################################
# Is a Python string variable with the following options and their meanings:                                  #
#	'skip'     will use the molIDs set in the files that are read in. If the file format supplied to          #
#              cell_builder.py does not have molIDs, every atom molID will default to one.                    #
#                                                                                                             #
#	'files'    will set molIDs based on the order the files are read into cell_builder.py, where every atom   #
#              in the first file will be assigned to molID one, every atom in the second file will be assigned#
#              to molID two and so on.                                                                        #
#                                                                                                             #
#   'offset'   will offset the molIDs in each file as each file is read in. If every atoms molID are the same #
#              in each file the 'offset' method and the 'files' method will create the same result. If the    #
#              file contains different molIDs on different atoms and you wish to maintain the distinction,    #
#              the offset is the method to use.                                                               #
#                                                                                                             #
#	'clusters' will perform a cluster analysis, where the "clusters" of atoms are determined via bonding      #
#              connectivity, where the criteria for atoms to be part of the same "cluster" is that the atoms  #
#              must be linked by at least one covalent bond to the cluster. The clusters are then sorted by   #
#              the number of atoms, where the largest number of atoms is identified as cluster one, then      #
#              molIDs are incremented and are assigned to each of the remaining clusters (i.e., the largest   #
#              cluster of atoms will have molID 1, and then the smallest cluster will have a molID of         #
#              NCLUSTERS found).                                                                              #
#                                                                                                             #
# Depending on the different analysis and/or visualization the different reset_molid options may be useful,   #
# since they can be used to identity groups of atoms, which can be tedious depending on the type of modeling  #
# that is being used. Examples:                                                                               #
#    reset_molids = 'skip'     # every duplicated atom will maintain the molID in the read in file            #
#    reset_molids = 'files'    # every atom in each file will be assigned to a unique molID                   #
#    reset_molids = 'offset'   # each time a file is read the molIDs will be offset incrementally             #
#    reset_molids = 'clusters' # molIDs will be consistent with the clusters of atoms                         #
#                                                                                                             #
# Update reset_molids as desired (Default should be 'clusters').                                              #
###############################################################################################################
reset_molids = 'clusters'




#################################################
### Import main from src.cell_builder and run ###
#################################################
if __name__ == "__main__":  
    # Import main from src.cell_builder.main
    from src.cell_builder.main import main
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
    #       python3 cell_builder.py -opt                                                               #
    #              or                                                                                  #
    #       python3 cell_builder.py -man                                                               #
    #   in the terminal and the code will print out all command line option and terminate before any   #
    #   further analysis is done                                                                       #
    ####################################################################################################
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: use_GUI = False
    if use_GUI or '-gui' in commandline_inputs:
        from src.cell_builder.GUI import cell_builder_GUI
        print('\n\n\ncell_builder is currently running in GUI mode, where all GUI inputs are intialized from cell_builder.\n\n\n')
        cell_builder_GUI(files, force_field_joining, duplicate, distance_scale, newfile, atom_style, parent_directory, max_rotations, reset_molids,
                         unwrap_atoms_via_image_flags, include_type_labels, group_monomers_locally, seed, domain, GUI_zoom, nfiles=maxfiles)
    else:   
        # Run main cell_builder
        main(files, force_field_joining, duplicate, distance_scale, newfile, atom_style, parent_directory, max_rotations,
             reset_molids, unwrap_atoms_via_image_flags, include_type_labels, group_monomers_locally, seed, domain,
             commandline_inputs=commandline_inputs)
