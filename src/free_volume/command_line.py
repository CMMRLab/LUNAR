#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
Sept 29th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

File contains command line options/man page and
command line over riding class to reset hard coded
inputs.
"""

##############################
# Import Necessary Libraries #
##############################
import sys


def print_man_page(topofile, max_voxel_size, boundary, parent_directory, compute_free_volume_distributions, run_mode, probe_diameter,
                   vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels):

    # print general command line options
    print('\n\nfree_volume.py has been run with -opt or -man option to show the optional command line overrides available. Command line option summary [-tag <tag-input>]:')
    print('python3 free_volume.py [-topo <topo-filename>] [-dir <new directory name>] [-voxel <float>] [-pbc <p-p-p or f-f-f or p-f-f or ...>]')
    print("                       [-probe-diameter <float or 'min-voxel'>] [-run-mode <stl or stl-dd or numpy or numba or numba-p or numba-dd")
    print("                       or numba-ddp or CUDA or CUDA-dd>] [-cuda-atoms <int>] [-cuda-voxels <int>] <-opt>|<-man> [*NOTE: Not all")
    print("                       options found in free_volume.py are currently supported via command line overrides.*]")
    
    
    print('\n*NOTE: If the command line input options are not used the hard coded inputs in the inputs section of the atom_typing.py')
    print('file will be enforced. The command line inputs option gives more flexibility to the code depending on IDE usage or')
    print('command line usage. The user may specify as many or as few command line options as desired, but if the command line')
    print('option is used and not all options are provided by the user, the hard coded inputs will be enforced, and the user')
    print('will be warned.*')
    
    print('\nTo use command line override options there must be alteration between tags started with the - character and the')
    print('input value. Ordering of tag/input pairs do not matter. Example:')
    print('    python3 free_volume.py -tag1 tag1-input  -tag2 tag2-input  -tag3 tag3-input  ...N-arbitrary tagsN/tagN-input pairs')
    
    # print options
    print('\n\nOptions:')
    
    # print -topo option
    print(f'\n -topo or -t <topo-filename>  free_volume variable: topofile    hard coded: {topofile}')
    print('    Command line option to set topofile name of topo/coord file to be read into atom_typing for conversion. Currently')
    print('    supported topofile formats are:')
    print('         .data = .data or dat file from LAMMPS (Info: atom positions, bonds, no center, box set by LAMMPS, image')
    print('                 flags set by LAMMPS)')
    print('    Example usage:')
    print('        python3 free_volume.py -topo EXAMPLE_TOPO-FILE.data')
    
    # # print -dir option
    print(f'\n -dir or -d <new directory name>   free_volume variable: parent_directory    hard coded: {parent_directory}')
    print('    Command line option to set new directory name to store all files atom_typing.py can write. The file')
    print('    management is based on the present working directory (pwd) and will store files one directory deep inside the')
    print('    pwd. The -dir option will specify the directory that is inside of pwd. This will allow for easier file ')
    print('    management and allows the user to merge multiple files for a certain system or research project and keep them')
    print('    all organized into a single directory.')
    print('    python variables*. Example usage:')
    print('        python3 free_volume.py -dir EXAMPLE_TEST_SYSTEM_DIRECTORY')
    
    # print -pbc option
    print(f'\n -pbc or -p <p-p-p or f-f-f or p-f-f or ...>   free_volume variable: boundary    hard coded: {boundary}')
    print('    Command line option to set boundary conditions. NOTE when setting at the command line there must not be any')
    print('    white space, thus the x, y, and z face boundaries must be strung together with the "-" character (IE p-p-p or')
    print('    f-f-f or p-f-f or ....). Example usage:')
    print('        python3 free_volume.py -pbc  p-p-p')
    
    # print -run-mode option
    print(f'\n -run-mode or -run <stl or numpy or numba>   free_volume variable: run_mode    hard coded: {run_mode}')
    print('    Command line option to set Run mode to set how the code computes atom to voxel distances. All run modes')
    print('    required the tqdm module, whereas some of the more advanced run modes requires numpy and/or numba')
    print('    modules. The following run mode optiohs are available (order from slowest to quickest):')
    print('        stl        runs on the standard python libraby (slowest)')
    print('        numpy      uses numpy vectorization to speedup the code (medium speed)')
    print('        numba      uses numpy looped version compiled to machine code via the numba.njit methodolgy (quickest)')
    print('        numba-p    uses numpy looped version compiled to machine code via the numba.njit methodolgy and runs in parallel (quickest)')
    print('        numba-dd   uses numpy looped version compiled to machine code via the numba.njit methodolgy applying a domain decompositions (quickest)')
    print('        numba-ddp  uses numpy looped version compiled to machine code via the numba.njit methodolgy  applying a domain decompositions and runs in parallel (quickest)')
    print('        CUDA       uses numpy looped version compiled to machine code for NVIDIA GPU via the numba.cuda.jit methodolgy (quickest)')
    print('        CUDA-dd    uses numpy looped version compiled to machine code for NVIDIA GPU via the numba.cuda.jit methodolgy  applying a domain decompositions (quickest)')
    print('    Example usage:')
    print('        python3 free_volume.py -run-mode  numba')
    
    # print -cuda-atoms option
    print(f'\n -cuda-atoms or -ca <int>   free_volume variable: CUDA_threads_per_block_atoms    hard coded: {CUDA_threads_per_block_atoms}')
    print('     Command line option to set the number of threads per block for atom parallilzation operations on the GPU. The values should be')
    print('     be doubling multiples of 8 (i.e. 8, 16, 32, 64, 128, 256, 512, 1024). Example usage:')
    print('        python3 free_volume.py -cuda-atoms 128')
    
    # print -cuda-voxels option
    print(f'\n -cuda-voxels or -cv <int>   free_volume variable: CUDA_threads_per_block_voxels    hard coded: {CUDA_threads_per_block_voxels}')
    print('     Command line option to set the number of threads per block for voxel parallilzation operations on the GPU. The values should be')
    print('     be doubling multiples of 8 (i.e. 8, 16, 32, 64, 128, 256, 512, 1024). Additionally is set to zero the voxel operations such as')
    print('     generating the the voxels or the voxel connecitivty analysis will be performed on the CPU instead of the GPU. Example usage:')
    print('        python3 free_volume.py -cuda-voxels 128')
    
    # print -vdw-method option
    print(f'\n -vdw-method or -vdw <dict or class1 or class2>   free_volume variable: run_mode    hard coded: {vdw_method}')
    print('    Command line option to tell free_volume.py which vdw radii to use for setting the occupied vdw volume.')
    print('    The following methods are supported:')
    print('        ‘dict’ where the used vdw radii will come from the vdw_radius dictionary.')
    print('        ‘class1’ where the used vdw radii will come from the read-in topofile and the topofile uses the 12-6 LJ')
    print('                 potential and parameter ordering in file is [epsilon sigma].')
    print('        ‘class2’ where the used vdw radii will come from the read-in topofile and the topofile uses the 9-6 LJ')
    print('                 potential and parameter ordering in file is [epsilon sigma].')
    print('    Example usage:')
    print('        python3 free_volume.py -vdw-method dict')

    # print -voxel option
    print(f'\n -voxel or -v <float>   free_volume variable: max_voxel_size    hard coded: {max_voxel_size}')
    print('    Command line option to set the max voxel size. Example usage:')
    print('        python3 free_volume.py -voxel  0.5')
    
    # print -probe-diameter option
    print(f"\n -probe-diameter or -pd <float or 'min-voxel'>   free_volume variable: probe_diameter    hard coded: {probe_diameter}")
    print('    Command line option to set the -probe-diameter. Example usage:')
    print('        python3 free_volume.py -probe-diameter  0.5')
    
    # print -free-vol-dist option
    print(f'\n -free-vol-dist or -fve <T|F>   free_volume variable: compute_free_volume_distributions    hard coded: {compute_free_volume_distributions}')
    print('    Command line option use compute the free volume voxel connectivty and distributions. T is for')
    print('    True and use option and F is for False. T will tell the code to use compute free volume voxel')
    print('    connectivity and F will skip usage. Example usage:')
    print('        python3 free_volume.py -free-vol-dist T')
    
    # print requirements and dependencies
    print('\n\nRequirements and Dependencies:')
    print('    Requirements: python 3.7+ to guarantee ordering of dictionaries which stores most of the data')
    
    # print buffer and exit code
    print('\n\n')
    return

# Class to update inputs
class inputs:
    def __init__(self, topofile, max_voxel_size, boundary, parent_directory, compute_free_volume_distributions, run_mode, probe_diameter,
                 vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels, commandline_inputs,):
        
        # Give access to inputs (update later on if command line over ride is given)
        self.commandline_inputs = commandline_inputs
        self.topofile = topofile
        self.parent_directory = parent_directory
        self.boundary = boundary
        self.max_voxel_size = max_voxel_size
        self.compute_free_volume_distributions = compute_free_volume_distributions
        self.run_mode = run_mode
        self.probe_diameter = probe_diameter
        self.vdw_method = vdw_method
        self.CUDA_threads_per_block_atoms = CUDA_threads_per_block_atoms
        self.CUDA_threads_per_block_voxels =CUDA_threads_per_block_voxels
        
        
        # Check that the given command line inputs are even for alternating tags/tag-inputs
        if (len(commandline_inputs) % 2) != 0:
            print('\nERROR command line option over ride used but odd number of command line arguments provided (missing tag or tag-input value).\n')
            sys.exit()
        
        # Check that tags are alternating between tag and tag-input
        nonalternating_tags = [commandline_inputs[i] for i in range(0, len(commandline_inputs), 2) if commandline_inputs[i][0] != '-']
        if nonalternating_tags:
            print(f'\nERROR tags are not alernating between tag and tag-input. Incorrect tag(s) list: {str(nonalternating_tags)} in {str(commandline_inputs)}\n')
            sys.exit()
        
        
        # Check that tag is supported and log if tag from the command line
        # set supported tags
        supported_tags = ['-topo', '-dir', '-pbc', '-voxel', '-free-vol-dist', '-run-mode', '-probe-diameter', '-vdw-method', '-cuda-atoms', '-cuda-voxels']
        
        # set shortcut_tags mapping
        shortcut_tags = {'-t':'-topo', '-d':'-dir', '-p':'-pbc', '-v':'-voxel', '-fve':'-free-vol-dist', '-run':'-run-mode', '-pd':'-probe-diameter',
                         '-vdw':'-vdw-method', '-ca':'-cuda-atoms', '-cv': '-cuda-voxels'}
        
        # set default variables
        default_variables ={'-topo': self.topofile, '-dir': self.parent_directory, '-pbc': self.boundary, '-voxel': self.max_voxel_size,
                            '-free-vol-dist':self.compute_free_volume_distributions, '-run-mode': self.run_mode, '-probe-diameter':self.probe_diameter,
                            '-vdw-method':self.vdw_method, '-cuda-atoms':self.CUDA_threads_per_block_atoms, '-cuda-voxels': self.CUDA_threads_per_block_voxels}
        
        # set tag/tag-input pair as empty string and update
        tags = {i:'' for i in supported_tags}
        
        # Update tags and raise exception if user provided unsupported tag
        for i in range(0, len(commandline_inputs), 2):
            user_tag = commandline_inputs[i]
            tag_input = commandline_inputs[i+1]
            
            # Find if shortcut was used and update to full tag if in shortcuts
            if user_tag in shortcut_tags: 
                user_tag = shortcut_tags[user_tag]
            
            # 1st check make sure it is a tag
            if user_tag not in supported_tags:
                print(f'\nERROR requesting unsupported command line tag   {str(user_tag)}\n')
                sys.exit()
                
            # If 1st check clears add input to tags
            tags[user_tag] = tag_input
            
        # Loop through tag_checks and warn user hard coded variables will be enforced
        print('\n\nCommand line run option override checks (will warn if command line run option is used but not all options are provided at the command line):')
        for i in tags:
            if not tags[i]:
                print('WARNING override option   {:<18} not provided at the command line. Hard coded input being enforced: {}'.format(i, default_variables[i]))
        
                
        # Start changing python variables based on command line over rides
        print('\n\nCommand line run option override inputs confirmation (will confirm to user that the command line override does override the hard coded inputs):')
        def T_F_string2boolean(tag, string):
            if string == 'F': return False
            elif string == 'T': return True
            else: print(f'\nERROR Provided tag: {tag} string is {string} which is not T or F. Use T or F string\n'); sys.exit()
        
        
        ###############################################
        # set new -topo option and print confirmation #
        ###############################################
        if tags['-topo']:
            self.topofile = tags['-topo']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-topo', self.topofile))
            
        # set new -cuda-atoms option and print confirmation (ERROR checks will occur in the next step so only try int except set as input)
        if tags['-cuda-atoms']:
            try: self.CUDA_threads_per_block_atoms = int(tags['-cuda-atoms'])
            except: self.CUDA_threads_per_block_atoms = tags['-cuda-atoms']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-cuda-atoms', self.CUDA_threads_per_block_atoms))
            
        # set new -cuda-voxels option and print confirmation (ERROR checks will occur in the next step so only try int except set as input)
        if tags['-cuda-voxels']:
            try: self.CUDA_threads_per_block_voxels = int(tags['-cuda-voxels'])
            except: self.CUDA_threads_per_block_voxels = tags['-cuda-voxels']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-cuda-voxels', self.CUDA_threads_per_block_voxels))
                
        # set new -dir option and print confirmation
        if tags['-dir']:
            if tags['-dir'].isspace() or tags['-dir'] == '':
                print('-dir flag was used, with an empty directory string. The directory string must contain alpha-numeric characters')
            else:
                self.parent_directory = tags['-dir']
                print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-dir', self.parent_directory))
            
        # set new -pbc option and print confirmation
        if tags['-pbc']:
            self.boundary= ' '.join(tags['-pbc'].split('-'))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-pbc', self.boundary))
            
        # set new -run-mode option and print confirmation
        if tags['-run-mode']:
            self.run_mode = tags['-run-mode']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-run-mode', self.run_mode))
            
        # set new -vdw-method option and print confirmation
        if tags['-vdw-method']:
            self.vdw_method = tags['-vdw-method']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-vdw-method', self.vdw_method))
            
        # set new -voxel option and print confirmation (ERROR checks will occur in the next step so only try float except set as input)
        if tags['-voxel']:
            try: self.max_voxel_size = float(tags['-voxel'])
            except: self.max_voxel_size = tags['-voxel']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-voxel', self.max_voxel_size))
            
        # set new -probe-diameter option and print confirmation (ERROR checks will occur in the next step so only try float except set as input)
        if tags['-probe-diameter']:
            try: self.probe_diameter = float(tags['-probe-diameter'])
            except: self.probe_diameter = tags['-probe-diameter']
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-probe-diameter', self.probe_diameter))
            
        # set new -free-vol-dist option and print confirmation
        if tags['-free-vol-dist']:
            self.compute_free_volume_distributions = T_F_string2boolean('-free-vol-dist', (tags['-free-vol-dist']))
            print('Override confirmation for {:<18} Hard codeded input is being overridden with this input: {}'.format('-free-vol-dist', self.compute_free_volume_distributions))

        # print buffer
        print('\n\n')