# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.6
December 6th, 2023
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


#################################################################################################################
# Python string variable type to set directory where LAMMPS datafiles are stored. The code will change to this  #
# directory to read ALL files in this directory with the *.data extension. This means ONLY have LAMMPS *.data   #
# files in this directory you want to read into this code to automate the molecule logging process. Examples:   #
#     files_directory = 'EXMAPLES/auto_cluster_analysis' # Will read all *.data files in this directory         #
#                                                                                                               #
# Update files_directory as desired.                                                                            #
#################################################################################################################
files_directory = 'Manuscript/EPON_862/Replicate_packmol_liq'

##################################################################################################################
# Set the max voxel size in angstroms to discretize the 3D Simulation cell by. NOTE small voxels impose a very   #
# large computational cost. For testing of the code it is recommend to have a small system size (less then 1000  #
# atoms) and to set max_voxel_size = 1 angstrom.                                                                 #
#                                                                                                                #
# max_voxel_size as desired.                                                                                     #
##################################################################################################################
max_voxel_size = 0.25


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
parent_directory = 'topofile/CUDA_voxel_0.25_probe_1.1'


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
#     'dd'        variants uses a domain decomposition, where the subdomains sizes are determined via the        #
#                 following equation:                                                                            #
#                   subdomain = = round(2.1*max(system_vdw_radii), 2) + 2*max_voxel_size + probe_diameter        #
#                   restrictions that simulation cell must at least be 3*subdomain in each direction to use any  #
#                   of the ‘dd’ variants.                                                                        #
#                                                                                                                #
# Relative execution time table for the four different run modes using the following settings:                   #
#    topofile = 'EXAMPLES/free_volume/propanal_natoms_1000_typed_IFF.data'; probe_diameter = 0.0;                #
#    boundary = 'p p p'; max_voxel_size = 0.1; compute_free_volume_distributions = False;                        #
#          +-----------------------------------------------------------------------------------------+           #
#          | run mode  |   stl    |  numpy   |  numba   | numba-p  |  stl-dd  | numba-dd | numba-ddp |           #
#          |-----------+----------+----------+----------+----------+----------+----------+-----------|           #
#          |   stl     |    1.0   |    0.09  |    0.02  |  0.01    |  0.02    |  0.02    |  0.02     |           #
#          |  numpy    |   10.7   |    1.0   |    0.22  |   0.1    |  2.18    |  0.28    |  0.25     |           #
#          |  numba    |   49.4   |    4.63  |    1.0   |   0.5    |  10.1    |  1.31    |  1.17     |           #
#          | numba-p   |   95.4   |    8.9   |    1.9   |   1.0    |  19.4    |  2.50    |  2.23     |           #
#          | stl-dd    |   4.89   |    0.5   |    0.1   |   0.05   |  1.0     |  0.13    |  0.12     |           #
#          | numba-dd  |   37.8   |    3.6   |    0.8   |   0.4    |  7.7     |  1.0     |  0.89     |           #
#          | numba-ddp |   42.3   |    3.9   |    0.9   |   0.5    |  8.6     |  1.12    |  1.0      |           #
#          +-----------------------------------------------------------------------------------------+           #
# Note that numba-p is around 100x quicker then the stl run mode and sbout 1.9x quicker then the numba run mode  #
# for this test case, but numba-p maybe slower then numba for some systems depending on the overhead cost of     #
# parallelizing the code. Additionaly the first time you run numba and numba-p on your machine some of the run   #
# time is spent compiling the code, where sequential runs the code will be compliled. Therefore when using the   #
# numba or numpa-p run methods for the first time on your computer it is recommended to use a small system and   #
# a large voxel size to compile to code on your machine and then you can use larger systems and smaller voxel    #
# sizes.                                                                                                         #
#                                                                                                                #
# Examples:                                                                                                      #
#    run_mode = 'stl'     # runs on the standard python libray + usage of the tqdm progress bars                 #
#    run_mode = 'numpy'   # runs on the standard python libray + numpy vectorization + tqdm progress bars        #
#    run_mode = 'numba'   # runs on the standard python libray + numpy arrays + numba njit + tqdm progress bars  #
#    run_mode = 'numba-p' # runs on the standard python libray + numpy arrays + numba njit + tqdm progress bars  #
#                                                                                                                #
# Update run_mode as desired (If numba is installed, default should be 'numba', since it is the quickest).       #
##################################################################################################################
run_mode = 'numba-ddp'
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
CUDA_threads_per_block_atoms = 128
CUDA_threads_per_block_voxels = 256


##################################################################################################################
# Commands for creating Files that this code can write (True or False responses). NOTE depending on system volume#
# some of the *.data files that this code writes can get very large (IE a 1000 atom system with in gaseous state #
# using a max_voxel_size = 0.1 has been seen to generate 2GB and larger files).                                  #
##################################################################################################################
files2write = {'write_atoms_free':  False,  # File containing the voxels assigned to the atoms and the free volume (name *_atoms_free.data)
                'write_bonds_free': False,  # File containing atoms/bonds and the free volume voxels (name *_bonds_free.data)
                'write_atoms_only': False,  # File containing the voxels assigned to the atoms with the LAMMPS atomTypeIDs set as element types (name *_atoms_only.data)
                'write_free_only' : False,  # File containing the voxels assigned to thethe free volume (name *_free_only.data)
                'write_all_voxels': False,  # File containing voxels that were generated before assigning any to the free volume or atom volume (name *_voxels_only.data)
			    'write_spat_dis-x': True,  # File containing the spatial free volume distribution in X-dir (name *__spatial_distribution_direction_x.csv)
                'write_spat_dis-y': True,  # File containing the spatial free volume distribution in Y-dir (name *__spatial_distribution_direction_y.csv)
                'write_spat_dis-z': True,  # File containing the spatial free volume distribution in Z-dir (name *__spatial_distribution_direction_z.csv)
                }

##################################################################################################################
# Python boolean variable to set the usage of computing the free volume voxel connectivity and then computing    #
# free volume distributions. True or False, where True turns on the analysis and False shuts it off. This option #
# also works with the boundary option and will compute the connectivity of the voxels based on the boundary.     #
##################################################################################################################
compute_free_volume_distributions = False # Seems to work okay, but leave as False for the most part
compute_free_volume_distributions = True


##################################################################################################################
# Python dictionary to set element symbol from the mass of the atom if reading in a LAMMPS .data file. Each key  #
# in the dictionary signifies an element and must have a value of a list with at least one mass in the list. No  #
# examples are given since this is self-explanatory, however if a read file does cannot find the information     #
# needed from this dictionary and error will be issued and the code will exit.                                   #
#                                                                                                                #
#        element: [list of masses to identify element type from datafile                                         #
#        symbol :  minimally each list needs at least 1 mass in them]                                            #
##################################################################################################################
mass_map = {'C':   [12.0000, 12.01115, 12.01100, 12.0112, 10.01115],
            'H':   [1.0080, 1.00797, 1.00782, 1.0000],
            'O':   [15.9990, 15.99940, 15.99491, 15.994],
            'N':   [14.0000, 14.00670],
            'S':   [32.06400, 32.06],
            'F':   [18.9984, 18.998400],
            'Si':  [28.086, 28.085000, 28.085500, 28.06],
            'Xe':  [131.30000],
            'Ne':  [20.18300],
            'Kr':  [83.80000],
            'He':  [4.00300],
            'D':   [2.01400],
            'Cl':  [35.45300, 35.450000],
            'Ca':  [40.08000],
            'Br':  [79.90900, 79.904000],
            'Ar':  [39.94400],
            'P':   [30.97380, 30.973700],
            'Al':  [26.98154],
            'Mg':  [24.30500],
            'Li':  [6.941000],
            'Fe':  [55.84700],
            'Na':  [22.99000],
            'K':   [39.10],
            'Cs':  [132.9100],
            'Ba':  [137.3300],
            'Sr':  [87.6200],
            'Pb':  [207.2000],
            'Mo':  [95.95, 95.94]
             }  


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
vdw_method = 'class2'


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
    import src.io_functions as io_functions
    import src.free_volume.main as main
    import time
    import os

    ##########################
    # Setting up directories # 
    ##########################
    # Find present working directory and Find/create paths to store code results
    pwd = os.getcwd(); 
    path = os.path.join(pwd, files_directory);
    
    # Check if path exists. IF not create.
    if not os.path.isdir(path): 
        os.makedirs(path, exist_ok=True)
        
    # Change the current working directory to path to write all new files to outputs directory
    os.chdir(path)
    print('Moved to {} directory to find LAMMPS datafiles to iterate through'.format(path))
    
    #################################################################
    # Get all topofiles from files_directory and run free_volume.py #
    #################################################################
    start_time = time.time()
    topofiles = sorted([file for file in os.listdir(path) if file.endswith('.data')])
    for n, topofile in enumerate(topofiles, 1):
        log=io_functions.LUNAR_logger()
        main.main(topofile, max_voxel_size, mass_map, vdw_radius, boundary, parent_directory, compute_free_volume_distributions,
                  files2write, run_mode, probe_diameter, vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels, [], log=log)
        print(f'{n} of {len(topofiles)} completed')
        cool_time = 5
        print(f'Starting {cool_time} second cooling period before moving to next file.')
        time.sleep(cool_time)
        print(f'Finished {cool_time} second cool break, moving to next file.')
        
    execution_time = (time.time() - start_time)
    print('Execution time in seconds: ' + str(execution_time))
