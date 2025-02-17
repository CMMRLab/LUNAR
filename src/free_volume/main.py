# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.11
April 3rd, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.free_volume.command_line as command_line
import src.free_volume.out2console as out2console
import src.free_volume.write_lmp as write_lmp
import src.io_functions as io_functions
import src.read_lmp as read_lmp
import glob
import time
import sys
import os



############################################
# Main function to analyze the free volume #
############################################
def main(topofile, max_voxel_size, mass_map, vdw_radius, boundary, parent_directory, compute_free_volume_distributions,
         files2write, run_mode, probe_diameter, vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels,
         commandline_inputs, log=None):
    
    # Set up Tristan's "array" analysis using recursion
    if not os.path.isfile(str(topofile)):
        if log is None: log = io_functions.LUNAR_logger()
        log.configure(level='production', print2console=True, write2log=True)
        files = glob.glob(topofile); array_time = time.time()
        if files:
            for n, file in enumerate(files, 1):
                log.clear_all()
                log.out('\n\nUsing array input option:')
                log.out(' - topofile     : {}'.format(topofile))
                log.out(' - matched file : {}'.format(file))
                log.out(' - elapsed time : {:.2f} (seconds)'.format(time.time() - array_time))
                log.out(' - progress     : {} of {} ({:.2f}%)'.format(n, len(files), 100*(n/len(files))))
                if file.endswith('voxels_only.data') or file.endswith('atoms_only.data') or file.endswith('free_only.data'):
                    log.warn(f' - WARNING matched file {file} has free_volume extension. Skipping file to avoid processing a free_volume output.')
                    continue
                elif file.endswith('atoms_free.data') or file.endswith('bonds_free.data'):
                    log.warn(f' - WARNING matched file {file} has free_volume extension. Skipping file to avoid processing a free_volume output.')
                    continue
                try: # we dont want crashes to exit this loop
                    main(file, max_voxel_size, mass_map, vdw_radius, boundary, parent_directory, compute_free_volume_distributions,
                         files2write, run_mode, probe_diameter, vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels,
                         commandline_inputs, log)
                except: pass
            print('\a') # Alert
        else: log.error(f'ERROR topofile: {topofile} unwrapped zero files or does not exist')
        return 
    else:
        # To remove recursion, un-indedent this all the way to the return and delete the lines above
        start_time = time.time()
        
        # Configure log (default is level='production', switch to 'debug' if debuging)
        if log is None:
            log = io_functions.LUNAR_logger()
        log.configure(level='production')
        #log.configure(level='debug')
        
        # set version and print starting information to screen
        version = 'v1.11 / 3 April 2024'
        log.out(f'\n\nRunning free_volume {version}')
        log.out(f'Using Python version {sys.version}')
        
        #-----------------------#
        # Command Line Override #
        #-----------------------#
        # if -opt or -man option is in commandline_inputs print options and stop code execution
        if '-opt' in commandline_inputs or  '-man' in commandline_inputs:
            # call man page and exit if '-opt' or '-man' is provided at the command line
            command_line.print_man_page(topofile, max_voxel_size, boundary, parent_directory, compute_free_volume_distributions, run_mode, probe_diameter,
                                        vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels)
            sys.exit()
        
        #---------------------------------------------------------------------------------#
        # if -opt option is NOT in commandline_inputs start assiging over ride parameters #
        #---------------------------------------------------------------------------------#
        if '-opt' not in commandline_inputs and commandline_inputs:    
            # call inputs for commandline over rides
            over_rides = command_line.inputs(topofile, max_voxel_size, boundary, parent_directory, compute_free_volume_distributions, run_mode, probe_diameter,
                                             vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels, commandline_inputs)
            
            # Set new inputs from over_rides class
            topofile = over_rides.topofile
            parent_directory = over_rides.parent_directory
            max_voxel_size = over_rides.max_voxel_size
            boundary = over_rides.boundary
            compute_free_volume_distributions = over_rides.compute_free_volume_distributions
            run_mode = over_rides.run_mode
            probe_diameter = over_rides.probe_diameter
            vdw_method = over_rides.vdw_method
            CUDA_threads_per_block_atoms = over_rides.CUDA_threads_per_block_atoms
            CUDA_threads_per_block_voxels = over_rides.CUDA_threads_per_block_voxels
        
        #--------------------------------------------------#
        # Define run mode to find voxels to atom distances #
        #--------------------------------------------------#
        if run_mode == 'stl':
            import src.free_volume.voxels_standard as voxels
            log.out('Running in stl mode (standard library mode)')
        elif run_mode == 'stl-dd':
            import src.free_volume.voxels_standard_domain as voxels
            log.out('Running in stl mode (standard library mode w/domain decomposition optimization)')
        elif run_mode == 'numpy':
            import src.free_volume.voxels_vectorized as voxels
            log.out('Running in numpy mode (numpy vectorized mode)')
        elif run_mode in ['numba', 'numba-p']:
            import src.free_volume.voxels_compiled as voxels
            log.out('Running in numba mode (compiled numpy arrays mode)')
        elif run_mode in ['numba-dd', 'numba-ddp']:
            import src.free_volume.voxels_compiled_domain as voxels
            log.out('Running in numba mode (compiled numpy arrays mode w/domain decomposition optimization)')  
        elif 'CUDA-dd' in run_mode:
            import src.free_volume.voxels_CUDA_domain as voxels
            log.out('Running in numba CUDA-dd mode (some parts will still run on CPU)')
        elif 'CUDA' in run_mode:
            import src.free_volume.voxels_CUDA as voxels
            log.out('Running in numba CUDA mode (some parts will still run on CPU)')
        else: log.error(f'ERROR input run_mode = {run_mode} not supported')
    
        #------------------#
        # Find free volume #
        #------------------#
        # Read lammps datafile as m
        if os.path.isfile(topofile):
            m = read_lmp.Molecule_File(topofile, method='forward', sections=['Atoms', 'Bonds'])
            log.out(f'Read in {m.filename} LAMMPS datafile')
        else: log.error(f'ERROR lammps datafile: {topofile} does not exist')
        log.out(f'System boundary {boundary}')
        
        # Add elements to m
        m.elements = set()
        for i in m.atoms:
            # Find element from mass_map and add to atom instance
            atom = m.atoms[i]; mass = m.masses[atom.type].coeffs[0];
            try:
                atom.element = [i for i in mass_map if mass in mass_map[i]][0]
                m.elements.add(atom.element)
            except: log.error(f'Not all masses in {topofile} are in the mass_map dictionary. Failed mass = {mass}')
        for i in m.masses:
            mass = m.masses[i].coeffs[0]
            try: 
                element = [i for i in mass_map if mass in mass_map[i]][0]
                m.masses[i].element = element 
                try: 
                    pair = m.pair_coeffs[i];
                    pair.element = element
                except: pass
            except: log.out(f'Not all masses in {topofile} are in the mass_map dictionary. Failed mass = {mass}')
    
                
    
        # Generate voxels
        m.CUDA_threads_per_block_atoms = CUDA_threads_per_block_atoms
        m.CUDA_threads_per_block_voxels = CUDA_threads_per_block_voxels
        v = voxels.generate(m, max_voxel_size, boundary, vdw_radius, probe_diameter, vdw_method, run_mode, log)
        probe_diameter = v.probe_diameter
        
        # Find free volume
        grid = voxels.analyze(m, v, vdw_radius, boundary, max_voxel_size, compute_free_volume_distributions, mass_map, run_mode, files2write, probe_diameter, vdw_method, log)
        
        
        #-------------------------------------#
        # Setting up directories and where to #
        # write final files and results to    #
        #-------------------------------------#
        # Find present working directory
        pwd = os.getcwd()
        
        # Find/create paths to store code results
        path = os.path.join(pwd, parent_directory)
        
        # If parent_directory == 'topofile' use topofile path as path
        if 'topofile' in parent_directory:
            log.out('\nUsing path from topofile to set parent_directory ...')
            path = io_functions.get_dir_from_topofile(topofile, parent_directory)
        
        # Check if path exists. IF not create.
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
            
        # Change the current working directory to path
        # to write all new files to outputs directory
        os.chdir(path)
        
        
        #----------------------------------#
        # Find basename and write outfiles #
        #----------------------------------#
        basename = os.path.basename(topofile)
        if basename.endswith('data.gz'):
            basename = basename.replace('.gz', '') # remove .gz suffix
        basename = basename[:basename.rfind('.')]
        #string_boundary = '-'.join(boundary.split())
        #basename = '{}_boundary_{}_voxel_size_{}_probe_diameter_{}_run_mode_{}_vdw_method_{}'.format(basename, string_boundary, max_voxel_size, probe_diameter, run_mode, vdw_method)
        
        # Write LAMMPS pseudo files
        if files2write['write_all_voxels']: write_lmp.voxels_only(m, v, basename, version, log)
        if files2write['write_atoms_only']: write_lmp.atoms_only(m, grid, basename, version, log)
        if files2write['write_atoms_free']: write_lmp.atoms_free(m, grid, compute_free_volume_distributions, basename, version, log)
        if files2write['write_bonds_free']: write_lmp.bonds_free(m, grid, compute_free_volume_distributions, basename, version, log)
        if files2write['write_free_only']: write_lmp.free_only(m, grid, compute_free_volume_distributions, basename, version, log)
    
        # Script run time
        execution_time = (time.time() - start_time)
        
        # log to console
        out2console.out(m, v, grid, execution_time, boundary, max_voxel_size, compute_free_volume_distributions, basename, version, run_mode, files2write, probe_diameter, vdw_method, log)
        
        # Print file locations
        log.out(f'\n\nAll outputs can be found in {path} directory')
        
        # Print completion of code
        log.out('\n\nNormal program termination\n\n')
        
        # Show number of warnings and errors
        log.out_warnings_and_errors()
        
        # write log
        log.write_logged(basename+'.log.lunar')
        
        # Change back to the intial directory to keep directory free for deletion
        os.chdir(pwd)
        return