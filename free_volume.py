# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.12
April 14, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    **********************************************************
    * Requirements:                                          *
    *   python 3.7+                                          *
    *                                                        *
    * Dependencies:                                          *
    *   python tqdm module:                                  *
    *    - pip3 install tqdm (if pip manager is installed)   *
    *   python numpy module (if run_mode = 'numpy'):         *
    *    - pip3 install numpy (if pip manager is installed)  *
    *   python numba module: (if run_mode = 'numba'):        *
    *    - pip install numba (if pip manager is installed)   *
    *                                                        *
    * Run methods:                                           *
    *   - IDE (manipulate variables and run from IDE)        *
    *   - GUI (manipulate variables and run from. Default    *
    *          settings set from this script)                *
    *   - command line (python3 free_volume.py -man to get   *
    *                   command line override options and    *
    *                   examples)                            *
    *                                                        *
    *                                                        *
    **********************************************************
    
"""
##################################################################################################################
# Python boolean variable to use GUI or not (True or False). All GUI settings are intialized from this script so #
# adjusting this script will set the default settings that the GUI will load with. Please NOTE that when using   #
# the GUI to check the console or terminal print outs every time a system is run through the code because that is#
# where the import information will be displayed. Examples:                                                      #
#   use_GUI = True  # Will launch GUI to manipulate variables and run from GUI                                   #
#   use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have to be from the IDE or  #
#                     command line                                                                               #
#                                                                                                                #
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account for the  #
# different ways individuals have their OS display settings setup. A GUI_zoom = 100, means that default settings #
# are used, whereas a GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80 means decrease GUI by 20%. #
# Examples:                                                                                                      #
#   GUI_zoom = 100 # use default GUI size                                                                        #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                                    #
#                                                                                                                #
# Update use_GUI and GUI_zoom as desired.                                                                        #
##################################################################################################################
use_GUI = True
GUI_zoom = 100


##################################################################################################################
# LAMMPS datafile with or without bonds in it. The LAMMPS datafile can have new LAMMPS "type labels" in it. The  #
# only supported atom styles are full, charge, or molecular. More atom styles can be coded if needed. The LAMMPS #
# datafile CAN ONLY have orthogonal simulation cell.                                                             #
#                                                                                                                #
# Update topofile as desired.                                                                                    #
##################################################################################################################
topofile = 'EXAMPLES/free_volume/propanal_natoms_1000_typed_IFF.data'


##################################################################################################################
# Set the max voxel size in angstroms to discretize the 3D Simulation cell by. Be aware that small voxels impose #
# a very large computational cost.                                                                               #
#                                                                                                                #
# max_voxel_size as desired.                                                                                     #
##################################################################################################################
max_voxel_size = 1.0


##################################################################################################################
# Python float or int value that will set the probe size in insert into the voxelated grid. If the probe diameter#
# is zero the atom volume calculation is exactly the vdw volume of the atoms. Whereas if the probe diameter is   #
# greater than zero the calculation starts to approach the Positron annihilation lifetime spectroscopy (PALS)    #
# method and the probe diameter supplied to free_volume.py can be thought as the probe size of Positron          #
# annihilation lifetime spectroscopy (PALS) method.                                                              #
#                                                                                                                #
# Additionally, if probe_diameter  is set to ‘min-voxel’ the probe diameter will be updated to the minimum       #
# dimension of a voxel found in your system and dependent on the max_voxel_size variable. Examples:              #
#     probe_diameter = 0.0 # use atom volume as the exact vdw volume                                             #
#     probe_diameter = 1.1 # use a larger probe to compare results to PALS expeirmental results                  #
#     probe_diameter = 'min-voxel' # uses the minimum voxel dimension set for your system                        #
#                                                                                                                #
# Update probe_diameter as desired.                                                                              #
##################################################################################################################
probe_diameter = 1.1


##################################################################################################################
# Python string variable type to set parent directory to store all new files this code will write. If the        #
# variable is left as an empty string, files will be written to the present working directory. Example to set    #
# parent_directory as present working directory:                                                                 #
#     parent_directory = '' or parent_directory = '.'                                                            #
#                                                                                                                #
# Example to set parent_directory as with path from topofile:                                                    #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile'                               #
#     files will be written to 'Furan_Resin/files/' since the 'topofile' string envokes this option.             #
#                                                                                                                #
# Example to set parent_directory as with path from topofile and build dirs from that location:                  #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/NEWDIR'                        #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/NEWDIR' will build relative        #
#     directories                                                                                                #
#                                                                                                                #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/../NEWDIR'                     #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/../NEWDIR' will build relative     #
#     directories                                                                                                #
#                                                                                                                #
# Example to set parent_directory to a location inside of present working dirctory:                              #
#     parent_directory = 'inside'       # will save files in $pwd/inside                                         #
#     parent_directory = 'inside/test'  # will save files in $pwd/inside/test                                    #
#                                                                                                                #
# Example to set parent_directory to a location outside of present working dirctory:                             #
#     parent_directory = '../'      # will save files in $pwd/../                                                #
#     parent_directory = '../test'  # will save files in $pwd/../test                                            #
#                                                                                                                #
# Update parent_directory as desired.                                                                            #
##################################################################################################################
parent_directory = 'topofile'


##################################################################################################################
# Boundary to set the type of boundary to search for bonds. Setup like LAMMPS with the following meaning:        #
#   p is periodic                                                                                                #
#   f is non-periodic and fixed                                                                                  #
# This will be set by a python string with white space between each x, y, or z box-face flags. (IE 'p p p' sets  #
# all box sides to be periodic; 'f f f' sets all box sides to be non-periodic; 'p f f' sets the two x-box sides  #
# to be periodic and leaves the y/z-box sides as non-periodic ...). *NOTE: there MUST BE white space between     #
# each direction flag and there MUST BE 3-direction flags of 'f' or 'p' in any combination. Also searching for   #
# all periodic faces (27 total images) is much more computationally expensive. The interatomic distance search   #
# code is time-optimized as much as Josh has been able to make it, but if suffers from pythons slow nested       #
# looping limitations ...                                                                                        #
##################################################################################################################
boundary = 'p p p'


##################################################################################################################
# Run mode to set how the code computes atom to voxel distances. All run modes required the tqdm module, whereas #
# some of the more advanced run modes requires numpy and/or numba modules. The following run mode optiohs are    #
# available (order from slowest to quickest):                                                                    #
#     'stl'       runs on the standard python libraby (slowest) and requires the "tqdm" module.                  #
#     'numpy'     uses numpy vectorization to speedup the code (medium speed) and requires the "tqdm" and "numpy"#
#                 modules.                                                                                       #
#     'numba'     uses numpy looped version compiled to machine code via the numba.njit methodolgy (quickest) and#
#                 requires the "tqdm", "numpy", and "numba" modules.                                             #
#     'numba-p'   uses numpy looped version compiled to machine code and runs in parallel via the numba.njit and #
#                 parallel=True argument and requires the "tqdm", "numpy", and "numba" modules. numba will       #
#                 control the number of processors that is informed based on your hardware setup.                #
#     'stl-dd'    same meaning as 'stl', but uses a domain decomposition and linked cell list algorithm instead  #
#                 of the default algorithm.                                                                      #
#     'numba-dd'  same meaning as 'numba', but uses a domain decomposition and linked cell list algorithm        #
#                 instead of the default algorithm.                                                              #
#     'numba-ddp' same meaning as 'numba-p', but uses a domain decomposition and linked cell list algorithm      #
#                 instead of the default algorithm.                                                              #
#     'CUDA'      uses numpy looped version compiled to machine code via the numba.cuda.jit methodolgy (quickest)#
#                 and requires the "tqdm", "numpy", and "numba" modules. Will excecute the slow parts of the code#
#                 on the GPU (NOTE must have an NVIDIA GPU for this method to work on your machine).             #
#     'CUDA-dd'   uses numpy looped version compiled to machine code via the numba.cuda.jit methodolgy (quickest)#
#                 and requires the "tqdm", "numpy", and "numba" modules. Will excecute the slow parts of the code#
#                 on the GPU (NOTE must have an NVIDIA GPU for this method to work on your machine).             #
#     'dd'        variants uses a domain decomposition, where the subdomains sizes are determined via the        #
#                 following equation:                                                                            #
#                   subdomain = = round(2.1*max(system_vdw_radii), 2) + 2*max_voxel_size + probe_diameter        #
#                   restrictions that simulation cell must at least be 3*subdomain in each direction to use any  #
#                   of the ‘dd’ variants.                                                                        #
#                                                                                                                #
# Examples:                                                                                                      #
#    run_mode = 'stl'     # runs on the standard python libray + usage of the tqdm progress bars                 #
#    run_mode = 'numpy'   # runs on the standard python libray + numpy vectorization + tqdm progress bars        #
#    run_mode = 'numba'   # runs on the standard python libray + numpy arrays + numba njit + tqdm progress bars  #
#    run_mode = 'numba-p' # runs on the standard python libray + numpy arrays + numba njit + tqdm progress bars  #
#                                                                                                                #
# Update run_mode as desired (If numba is installed, default should be 'numba-ddp', since it is the quickest).   #
##################################################################################################################
run_mode = 'CUDA-dd'


##################################################################################################################
# Commands for setting the number of threads per block when parallelizing on the GPU. Depending on the number of #
# atoms in your system and the simulation cell size that is being discretized based on max_voxel_size variable,  #
# the number of atoms and number of voxels may differ. Therefore to optimize the calculation there is support    #
# for selecting the number of thread per block that gets assigned to each parallel run on the GPU. The blocks    #
# per grid is then computed as such:                                                                             #
#   blocks_per_grid_atoms  = math.ceil( natoms  / CUDA_threads_per_block_atoms )                                 #
#   blocks_per_grid_voxels = math.ceil( nvoxels / CUDA_threads_per_block_voxels )                                #
# ultimately setting the grid size that CUDA will use to parallize the code. The blocks per grid should be based #
# on doubling multiples of 8 and maximize at 1024 (ie. 8, 16, 32, 64, 128, 256, 512, 1024), however some GPUs    #
# may max out lower then 1024 threads per block.                                                                 #
#                                                                                                                #
# Examples:                                                                                                      #
#   CUDA_threads_per_block_atoms = 128    # Will use 128 threads per blocks for atoms parallilzation             #
#   CUDA_threads_per_block_voxels = 256   # Will use 256 threads per blocks for voxels parallilzation            #
#                                                                                                                #
# On the developers GPU, it was found that a good starting point for setting the blocks per grid is roughly every#
# 150 atoms add 1 to the threads per block for atoms and then select the closests doubling multiple of 8's and   #
# for smaller max_voxel_sizes (<0.5), double the threads per blocks for atoms to get the threads per blocks for  #
# voxels. Example:                                                                                               #
#   System: 1200 atoms                                                                                           #
#   Voxel size: < 0.5                                                                                            #
#   CUDA_threads_per_block_atoms  = 1200/150 = 8                                                                 #
#   CUDA_threads_per_block_voxels = 2*8 = 16                                                                     #
#                                                                                                                #
# It was also found some of the voxel based task are still quicker on the CPU when comparing to different        #
# CUDA_threads_per_block_voxels, so a shortcut was devised to select running the code on th CPU over the GPU if  #
# the run_mode is 'CUDA' or 'CUDA-dd', where is CUDA_threads_per_block_voxels is set to 0 (zero), the code will  #
# perform the voxel based tasks (such as voxel generation and free volume voxel connectivity on the CPU instead).#
#                                                                                                                #
# Update CUDA_threads_per_block_atoms and CUDA_threads_per_block_voxels as desired.                              #
##################################################################################################################
CUDA_threads_per_block_atoms = 8
CUDA_threads_per_block_voxels = 16


##################################################################################################################
# Commands for creating Files that this code can write (True or False responses). NOTE depending on system volume#
# some of the *.data files that this code writes can get very large (IE a 1000 atom system with in gaseous state #
# using a max_voxel_size = 0.1 has been seen to generate 2GB and larger files).                                  #
##################################################################################################################
files2write = {'write_atoms_free': False,  # File containing the voxels assigned to the atoms and the free volume (name *_atoms_free.data)
               'write_bonds_free': False,  # File containing atoms/bonds and the free volume voxels (name *_bonds_free.data)
               'write_atoms_only': False,  # File containing the voxels assigned to the atoms with the LAMMPS atomTypeIDs set as element types (name *_atoms_only.data)
               'write_free_only' : False,  # File containing the voxels assigned to thethe free volume (name *_free_only.data)
               'write_all_voxels': False,  # File containing voxels that were generated before assigning any to the free volume or atom volume (name *_voxels_only.data)
			   'write_spat_dis-x': False,  # File containing the spatial free volume distribution in X-dir (name *__spatial_distribution_direction_x.csv)
               'write_spat_dis-y': False,  # File containing the spatial free volume distribution in Y-dir (name *__spatial_distribution_direction_y.csv)
               'write_spat_dis-z': False,  # File containing the spatial free volume distribution in Z-dir (name *__spatial_distribution_direction_z.csv)
               }

##################################################################################################################
# Python boolean variable to set the usage of computing the free volume voxel connectivity and then computing    #
# free volume distributions. True or False, where True turns on the analysis and False shuts it off. This option #
# also works with the boundary option and will compute the connectivity of the voxels based on the boundary.     #
##################################################################################################################
compute_free_volume_distributions = False # Seems to work okay, but leave as False for the most part


##################################################################################################################
# The mass_map dictionary is now a "global" dictionary stored in src/masses.py. The purpose of this was to       #
# simplify adding new elements, where the new elements can now be applied to every code that uses the mass_map.  #
# If you get an "ERROR Not all masses in ... are in the mass_map dictionary.", you will now have to open         #
# src/masses.py and update the mass_map dictionary found in that file.                                           #
##################################################################################################################
import src.masses as masses
mass_map = masses.mass_map   


##################################################################################################################
# Is a Python string value to tell free_volume.py which vdw radii to use for setting the occupied vdw volume.    #
# The following methods are supported:                                                                           #
#    ‘dict’ where the used vdw radii will come from the vdw_radius dictionary.                                   #
#    ‘class1’ where the used vdw radii will come from the read-in topofile and the topofile uses the 12-6 LJ     #
#             potential and parameter ordering in file is [epsilon sigma].                                       #
#	 ‘class2’ where the used vdw radii will come from the read-in topofile and the topofile uses the 9-6 LJ      #
#             potential and parameter ordering in file is [epsilon sigma].                                       #
# Examples:                                                                                                      #
# vdw_method = 'dict'   # will get vdw radii from vdw_radius dict                                                #
# vdw_method = 'class1' # will get vdw radii from LAMMPS datafile assume Pair Coeffs section exists use 12-6 LJ  #
# vdw_method = 'class2' # will get vdw radii from LAMMPS datafile assume Pair Coeffs section exists use 9-6 LJ   #
#                                                                                                                #
# Update vdw_method as desired.                                                                                  #
##################################################################################################################
vdw_method = 'dict'


##################################################################################################################
# vdw radius to search for which voxels are filled by some part of an atom.                                      #
#                                                                                                                #
# Ref1 = Batsanov, Stepan S. "Van der Waals radii of elements." Inorganic materials 37.9 (2001): 871-885.        #
# Ref2 = https://en.wikipedia.org/wiki/Van_der_Waals_radius                                                      #
# The goal is to set via Ref1, but if not found in Ref1 used wikipedia Ref2                                      #
##################################################################################################################
vdw_radius = {'C':  1.70, # Ref1
              'H':  1.20, # Ref1
              'O':  1.55, # Ref1
              'N':  1.60, # Ref1
              'S':  1.80, # Ref1
              'F':  1.50, # Ref1
              'Si': 2.10, # Ref1
              'Xe': 2.16, # Ref2
              'Ne': 1.54, # Ref2
              'Kr': 2.02, # Ref2
              'He': 1.40, # Ref2
              'D':  2.40, # Ref1 (2*H vdw radius)
              'Cl': 1.80, # Ref1
              'Ca': 2.40, # Ref1
              'Br': 1.90, # Ref1
              'Ar': 1.88, # Ref2
              'P':  1.90, # Ref1
              'Al': 2.10, # Ref1
              'Mg': 2.20, # Ref1
              'Li': 2.20, # Ref1
              'Fe': 2.05, # Ref1
              'Na': 2.40, # Ref1
              'K':  2.80, # Ref1
              'Cs': 3.00, # Ref1
              'Ba': 2.70, # Ref1
              'Sr': 2.55, # Ref1
              'Pb': 2.30, # Ref1
              'Mo': 2.45,
              }



###################################
### Import needed files and run ###
###################################
if __name__ == "__main__":  
    import src.free_volume.main as main
    import sys
    
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: use_GUI = False
    if use_GUI or '-gui' in commandline_inputs:
        print('\n\n\nfree_volume is currently running in GUI mode, where all GUI inputs are intialized from free_volume.\n\n\n')
        from src.free_volume.GUI import free_volume_GUI
        free_volume_GUI(topofile, max_voxel_size, mass_map, vdw_radius, boundary, parent_directory,
                        compute_free_volume_distributions, files2write, run_mode, probe_diameter, 
                        vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels, GUI_zoom)
    else:
        main.main(topofile, max_voxel_size, mass_map, vdw_radius, boundary, parent_directory,
                  compute_free_volume_distributions, files2write, run_mode, probe_diameter, 
                  vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels, commandline_inputs)
