# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.15
July 14, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.free_volume.dump_splitting as dump_splitting
import src.free_volume.array_logging as array_logging
import src.free_volume.command_line as command_line
import src.free_volume.out2console as out2console
import src.free_volume.write_lmp as write_lmp
import src.glob_wildcards as glob_wildcards
import src.io_functions as io_functions
import src.read_dump as read_dump
import src.read_lmp as read_lmp
import traceback
import glob
import time
import sys
import os



############################################
# Main function to analyze the free volume #
############################################
def main(topofile, dumpfile, dump_settings, max_voxel_size, mass_map, vdw_radius, boundary, parent_directory, compute_free_volume_distributions,
         files2write, run_mode, probe_diameter, vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels,
         commandline_inputs, log=None):
    
    #-----------------------#
    # Command Line Override #
    #-----------------------#
    # if -opt or -man option is in commandline_inputs print options and stop code execution
    if '-opt' in commandline_inputs or  '-man' in commandline_inputs:
        # call man page and exit if '-opt' or '-man' is provided at the command line
        command_line.print_man_page(topofile, dumpfile, dump_settings, max_voxel_size, boundary, parent_directory, compute_free_volume_distributions, run_mode, probe_diameter,
                                    vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels)
        sys.exit()
    
    #---------------------------------------------------------------------------------#
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters #
    #---------------------------------------------------------------------------------#
    if '-opt' not in commandline_inputs and commandline_inputs:    
        # call inputs for commandline over rides
        over_rides = command_line.inputs(topofile, dumpfile, dump_settings, max_voxel_size, boundary, parent_directory, compute_free_volume_distributions, run_mode, probe_diameter,
                                         vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels, commandline_inputs)
        
        # Set new inputs from over_rides class
        topofile = over_rides.topofile
        dumpfile = dump_settings
        dump_settings = dump_settings
        parent_directory = over_rides.parent_directory
        max_voxel_size = over_rides.max_voxel_size
        boundary = over_rides.boundary
        compute_free_volume_distributions = over_rides.compute_free_volume_distributions
        run_mode = over_rides.run_mode
        probe_diameter = over_rides.probe_diameter
        vdw_method = over_rides.vdw_method
        CUDA_threads_per_block_atoms = over_rides.CUDA_threads_per_block_atoms
        CUDA_threads_per_block_voxels = over_rides.CUDA_threads_per_block_voxels
        
    #---------------------------#
    # Set up dump file analysis #
    #---------------------------#
    if os.path.isfile(str(topofile)) and os.path.isfile(str(dumpfile)):
        if log is None: log = io_functions.LUNAR_logger()
        log.configure(level='production', print2console=True, write2log=True)
        
        # Read topofile and dump file
        mol = read_lmp.Molecule_File(topofile, method='forward', sections=['Atoms', 'Bonds'])
        log.out(f'Read in {mol.filename} LAMMPS datafile')
        
        dump = read_dump.read_dumpfile(dumpfile)
        log.out(f'Read in {dump.filename} LAMMPS dump')
        
        # Get the steps in dump and corresponding section numbers
        step_lst = sorted( dump.steps.keys() )
        sections = [i+1 for i in range(len(step_lst))]
        
        # User settings from dump_settings. The start and end values may be
        # numeric bounds or the strings "start" and "end", which select the
        # first and last available timestep or section, respectively.
        #
        # Examples:
        #   dump_settings = 'start=start; end=end; nevery=1; style=step'
        #   dump_settings = 'start=0; end=100000; nevery=10; style=step'
        #   dump_settings = 'start=start; end=end; nevery=1; style=section'
        #   dump_settings = 'start=1; end=10; nevery=2; style=section'
        dump_dict = dump_splitting.get_misc_setting(dump_settings)
        
        style  = dump_dict.get('style',  'step')
        nevery = dump_dict.get('nevery', 1)
        start  = dump_dict.get('start',  'start')
        end    = dump_dict.get('end',    'end')
       
        # Resolve generic start/end strings using either the available
        # timestep values or the sequential dump section numbers.
        if style == 'step':
            if start == 'start': 
                start = min(step_lst)
            if end == 'end':
                end = max(step_lst)
        elif style == 'section':
            if start == 'start':
                start = min(sections)
            if end == 'end': 
                end = max(sections)
        else:
            log.error('ERROR - dump splitting: style="{}" is not a dump splitting supported style. Use "step" or "section"'.format(style))
        
        if not isinstance(nevery, int):
            log.error('ERROR - dump splitting: nevery="{}" is not an int.'.format(nevery))
        if not isinstance(start, (int,float)):
            log.error('ERROR - dump splitting: start="{}" is not an int or a float.'.format(start))
        if not isinstance(end, (int,float)):
            log.error('ERROR - dump splitting: end="{}" is not an int or a float.'.format(end))
            
        # Pick the steps to use for analysis
        steps_in_bounds = []
        for section, step in zip(sections, step_lst):
            if style == 'step':
                gauge_value = step
        
            elif style == 'section':
                gauge_value = section
        
            if start <= gauge_value <= end:
                steps_in_bounds.append(step)
        
        # Apply nevery after bounds filtering
        steps = steps_in_bounds[::nevery]
        
        # Incrementally write sections from the dump file to a datafile
        array_file = 'free_volume_dump_logging'
        processed_files, failed_files = [], []
        array_data = {} # {'filename': output_dict}
        loop_time = time.time()
        pwd_loop  = os.getcwd()
        path, m, grid = None, None, None
        loop_dumpfile, loop_dump_settings = '', ''
        for n, step in enumerate(steps, 1):
            datafile = '{}.{}.data'.format(dumpfile, step)
            os.chdir(pwd_loop) # in-case of a crash 
            log.clear_all()
            log.out('\n\nUsing dump splitting option:')
            log.out(' - topofile     : {}'.format(topofile))
            log.out(' - dumpfile     : {}'.format(dumpfile))
            log.out(' - step in dump : {}'.format(step))
            log.out(' - elapsed time : {:.2f} (seconds)'.format(time.time() - loop_time))
            log.out(' - progress     : {} of {} ({:.2f}%)'.format(n, len(steps), 100*(n/len(steps))))
            
            # Try converting to a datafile
            try:
                box, atoms = dump.parse_atoms_and_box(step)
                dump_splitting.topo_dump_to_data(mol, box, atoms, datafile)
            except:
                stack_trace_string = traceback.format_exc()
                log.out(f'\nAn error occurred: {stack_trace_string}')
                log.out(f'Skipping step={step} from dumpfile="{dumpfile}", because the step failed to be converted to a datafile.')
                failed_files.append( datafile )
                continue
            
            try: # we dont want crashes to exit this loop
                m, grid = main(datafile, loop_dumpfile, loop_dump_settings, max_voxel_size, mass_map, vdw_radius, boundary, parent_directory, compute_free_volume_distributions,
                               files2write, run_mode, probe_diameter, vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels,
                               commandline_inputs, log)
                
                # Generate dict of outputs
                output_dict = {}
                output_dict['step'] = step
                output_dict['lx'] = float(m.lx)
                output_dict['ly'] = float(m.ly)
                output_dict['lz'] = float(m.lz)
                output_dict['density'] = float(m.density)
                output_dict['cell_volume'] = float(grid.simulation_volume)
                output_dict['atom_volume'] = float(grid.atom_volume)
                output_dict['free_volume'] = float(grid.free_volume)
                output_dict['percent_free_volume'] = float(grid.pfree_volume)
                output_dict['elapsed_time_seconds'] = time.time() - loop_time
                array_data[datafile] = output_dict
                processed_files.append( datafile )
                path = os.path.normpath(m.path)
            except: 
                stack_trace_string = traceback.format_exc()
                log.out(f'\nAn error occurred: {stack_trace_string}')
                failed_files.append( datafile )
                
        # Log values for array processing
        if array_data and path not in [None, '']:                 
            # Change to location to write file, then change back to pwd
            log.out('\n\nWriting files in: {}'.format(path))
            log.out(' - array_file={}'.format(f'{array_file}.csv'))
            os.chdir(pwd_loop) # in-case of a crash  
            os.chdir(path)   
            
            # Write array data
            array_logging.output(array_file, array_data)
            
            # Finalize array files
            os.chdir(pwd_loop)
        
        # Log any failed files
        if failed_files:
            log.out('\n\nFiles that failed to be analyzed for unknown reasons:')
            for file in failed_files:
                log.out(' - {}'.format(file))
                
        # Finalize array processing
        execution_time = (time.time() - loop_time)
        log.out('\n\nDump analysis time in seconds: ' + str(execution_time))
        print('\a') # Alert    
        return m, grid

    #---------------------------------------------------#
    # Set up Tristan's "array" analysis using recursion #
    #---------------------------------------------------#
    elif not os.path.isfile(str(topofile)):
        if log is None: log = io_functions.LUNAR_logger()
        log.configure(level='production', print2console=True, write2log=True)
        #log.configure(level='production', print2console=False, write2log=True)
        files = glob.glob(topofile); array_time = time.time()
        if files:
            # Try sorting files
            try:
                from natsort import natsorted
                files = natsorted(files)
                log.out('Read in xrdfiles were sorted using natsort, thus they should be iterated through')
                log.out('in a sequence that represents evolution or progression and logged in that order.')
            except:
                files = sorted(files)
                log.warn('WARNING natsort was not installed and read in xrdfiles were sorted using pythons')
                log.out('sorted function. The file iteration sequence may be unorder from a evolution or')
                log.out('progression stand point. Depending on the names given to the read in files. You may')
                log.out('intall natsort with:  pip3 install natsort   (if pip manager is installed)')
            
            # Setup loop data structures
            array_file = 'free_volume_array_logging'
            processed_files, failed_files = [], []
            array_data = {} # {'filename': output_dict}
            loop_time = time.time()
            pwd_loop  = os.getcwd()
            path, m, grid = None, None, None
            loop_dumpfile, loop_dump_settings = '', ''
            for n, file in enumerate(files, 1):
                os.chdir(pwd_loop) # in-case of a crash 
                log.clear_all()
                log.out('\n\nUsing array input option:')
                log.out(' - topofile     : {}'.format(topofile))
                log.out(' - matched file : {}'.format(file))
                log.out(' - elapsed time : {:.2f} (seconds)'.format(time.time() - array_time))
                log.out(' - progress     : {} of {} ({:.2f}%)'.format(n, len(files), 100*(n/len(files))))
                # if io_functions.check_outfile_existance(file, ':', parent_directory, filetype='topofile'):
                #     log.warn(f' - WARNING matched file {file} already has been processed and was skipped')
                #     continue
                # if file.endswith('voxels_only.data') or file.endswith('atoms_only.data') or file.endswith('free_only.data'):
                #     log.warn(f' - WARNING matched file {file} has free_volume extension. Skipping file to avoid processing a free_volume output.')
                #     continue
                # elif file.endswith('atoms_free.data') or file.endswith('bonds_free.data'):
                #     log.warn(f' - WARNING matched file {file} has free_volume extension. Skipping file to avoid processing a free_volume output.')
                #     continue
                try: # we dont want crashes to exit this loop
                    m, grid = main(file, loop_dumpfile, loop_dump_settings, max_voxel_size, mass_map, vdw_radius, boundary, parent_directory, compute_free_volume_distributions,
                                   files2write, run_mode, probe_diameter, vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels,
                                   commandline_inputs, log)
                    
                    # Get wildcards
                    wildcards = glob_wildcards.get_glob_wildcards(str(topofile), file)
                    
                    # Generate dict of outputs
                    output_dict = {}
                    output_dict['wildcards'] = wildcards
                    output_dict['lx'] = float(m.lx)
                    output_dict['ly'] = float(m.ly)
                    output_dict['lz'] = float(m.lz)
                    output_dict['density'] = float(m.density)
                    output_dict['cell_volume'] = float(grid.simulation_volume)
                    output_dict['atom_volume'] = float(grid.atom_volume)
                    output_dict['free_volume'] = float(grid.free_volume)
                    output_dict['percent_free_volume'] = float(grid.pfree_volume)
                    output_dict['elapsed_time_seconds'] = time.time() - loop_time
                    array_data[file] = output_dict
                    processed_files.append( file )
                    path = os.path.normpath(m.path)
                except: 
                    stack_trace_string = traceback.format_exc()
                    log.out(f'\nAn error occurred: {stack_trace_string}')
                    failed_files.append( file )
                    
            # Log values for array processing
            if array_data and path not in [None, '']:                 
                # Change to location to write file, then change back to pwd
                log.out('\n\nWriting files in: {}'.format(path))
                log.out(' - array_file={}'.format(f'{array_file}.csv'))
                os.chdir(pwd_loop) # in-case of a crash  
                os.chdir(path)   
                
                # Write array data
                array_logging.output(array_file, array_data)
                
                # Finalize array files
                os.chdir(pwd_loop)
            
            # Log any failed files
            if failed_files:
                log.out('\n\nFiles that failed to be analyzed for unknown reasons:')
                for file in failed_files:
                    log.out(' - {}'.format(file))
                    
            # Finalize array processing
            execution_time = (time.time() - loop_time)
            log.out('\n\nArray time in seconds: ' + str(execution_time))
            print('\a') # Alert
        
        else: log.error(f'ERROR topofile: {topofile} unwrapped zero files or does not exist')
        return m, grid
    else:
        # To remove recursion, un-indedent this all the way to the return and delete the lines above
        start_time = time.time()
        
        # Configure log (default is level='production', switch to 'debug' if debuging)
        if log is None:
            log = io_functions.LUNAR_logger()
        log.configure(level='production')
        #log.configure(level='debug')
        
        # set version and print starting information to screen
        version = 'v1.14 / 14 July 2026'
        log.out(f'\n\nRunning free_volume: {version}')
        log.out(f'Using Python version: {sys.version}')
        log.out(f'Using Python executable: {sys.executable}')
        
        #------------------#
        # Find free volume #
        #------------------#
        # Read lammps datafile as m
        if os.path.isfile(topofile):
            m = read_lmp.Molecule_File(topofile, method='forward', sections=['Atoms', 'Bonds'])
            log.out(f'Read in {m.filename} LAMMPS datafile')
        else: log.error(f'ERROR lammps datafile: {topofile} does not exist')
        log.out(f'System boundary {boundary}')
        
        #--------------------------------------------------#
        # Define run mode to find voxels to atom distances #
        #--------------------------------------------------#
        if run_mode == 'stl':
            import src.free_volume.voxels_standard as voxels
            log.out('Running in stl mode (standard library mode)')
            if m.xy != 0 or m.xz != 0 or m.yz != 0:
                log.error('ERROR simulaiton cell is NOT orthogonal and "stl" does not support triclinic. Switch to "numba-dd" or "numba-ddp" or "CUDA".')
        
        elif run_mode == 'stl-dd':
            import src.free_volume.voxels_standard_domain as voxels
            log.out('Running in stl mode (standard library mode w/domain decomposition optimization)')
            if m.xy != 0 or m.xz != 0 or m.yz != 0:
                log.error('ERROR simulaiton cell is NOT orthogonal and "stl" does not support triclinic. Switch to "numba-dd" or "numba-ddp" or "CUDA".')
        
        elif run_mode == 'numpy':
            import src.free_volume.voxels_vectorized as voxels
            log.out('Running in numpy mode (numpy vectorized mode)')
            if m.xy != 0 or m.xz != 0 or m.yz != 0:
                log.error('ERROR simulaiton cell is NOT orthogonal and "stl" does not support triclinic. Switch to "numba-dd" or "numba-ddp" or "CUDA".')
        
        elif run_mode in ['numba-dd', 'numba-ddp']:
            import src.free_volume.voxels_compiled_domain as voxels
            log.out('Running in numba mode (compiled numpy arrays mode w/domain decomposition optimization)')  
        
        elif 'CUDA-dd' in run_mode:
            import src.free_volume.voxels_CUDA_domain as voxels
            log.out('Running in numba CUDA-dd mode (some parts will still run on CPU)')
            
            # CUDA Checks
            from numba import cuda
            log.out(f' * cuda.is_available() = {cuda.is_available()}')
            try:
                dev = cuda.get_current_device()
                log.out(f' * CUDA device: {dev.name}')
            except Exception as e:
                log.out(f'Failed to get CUDA device: {e}')
                log.error('ERROR unable to initialize CUDA context. Use run_mode "numba-dd" or "numba-ddp"')
        
        else: log.error(f'ERROR input run_mode = {run_mode} not supported')

        
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
        m.path = path
        
        
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
        m = out2console.out(m, v, grid, execution_time, boundary, max_voxel_size, compute_free_volume_distributions, basename, version, run_mode, files2write, probe_diameter, vdw_method, log)
        
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
        return m, grid