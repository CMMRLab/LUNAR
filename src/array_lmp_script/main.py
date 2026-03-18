# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
March 18, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.array_lmp_script.auto_script_template as auto_script_template
import src.array_lmp_script.mode_from_script as mode_from_script
import src.array_lmp_script.lmp_data_metrics as lmp_data_metrics
import src.io_functions as io_functions
from itertools import product
import numpy as np
import shutil
import time
import glob
import copy
import sys
import os




#############################################################
# Function to import file from path useage:                 #
#  module = import_file('path/to/file/default_density.py')  #
#  print(module.mode) # print mode dict                     #
#############################################################
def import_file(path):
    from importlib import util
    root = os.path.basename(path)
    root = root[:root.rfind('.')]
    spec = util.spec_from_file_location(root, path)
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module  


# Nutter eval() so people cant do evil things
def nuttered_eval(string, globals_dict=None, locals_dict=None):
    # Setup the globals namespace to limit scope of what eval() can do
    allowed_builtins = ['min','max','sum','abs','len','map','range','reversed']
    copied_builtins = globals()['__builtins__'].copy()
    if globals_dict is None: globals_dict = {}
    if locals_dict is None: locals_dict = {}
    globals_dict['__builtins__'] = {key:copied_builtins[key] for key in allowed_builtins}
    globals_dict['np'] = globals()['np']
    return eval(str(string), globals_dict, locals_dict)


# Function to update an array dict using other keys values from that same dict. For example:
#     batch_array = {'$<ncores>': 2,
#                    '$<inv_cores>': '1/$<ncores>'}
#
# The '1/$<ncores>' value is calling the '$<ncores>' key and thus the return dictionary would be:
#     batch_array = {'$<ncores>': 2,
#                    '$<inv_cores>': 0.5}
def eval_values_from_keys(array_dict):
    globals_dict = {i:j for i, j in array_dict.items()}
    for key, value in array_dict.items():
        try:
            for key1, value1 in globals_dict.items():
                if key1 in str(value):
                    value = str(value).replace(key1, str(value1))
            value = nuttered_eval(value, globals_dict)
            array_dict[key] = value
        except: pass
    return array_dict



# Define main function
def main(mode, get_mode_from_script=True, log=None):
    start_time = time.time()
    
    
    # Configure log (default is level='production', switch to 'debug' if debuging)
    if log is None:
        log = io_functions.LUNAR_logger()
    log.configure(level='production')
    #log.configure(level='debug', print2console=False, write2log=False)
    

    # Set version and print starting information to screen
    version = 'v1.2 / 18 March 2026'
    log.out(f'\n\nRunning array_lmp_script: {version}')
    log.out(f'Using Python version: {sys.version}')
    log.out(f'Using Python executable: {sys.executable}')
    
    
    # Read current lmp_script content
    lmp_script = mode['lmp_script']
    root_filename = os.path.basename(lmp_script)
    root_directory = os.path.dirname(lmp_script)
    lmp_encoding = mode_from_script.get_encoding(lmp_script)
    lmp_script_content = ''
    with open(lmp_script, 'r', encoding=lmp_encoding, errors='ignore') as f:
        lmp_script_content = f.read()
        
    # Find first line and last line that starts width '#MODE:' to 
    # stop find and replace operations at these lines.
    lines = lmp_script_content.split('\n')
    mode_lines, script_lines = dict(), dict()
    for n, line in enumerate(lines):
        if line.startswith('#MODE:'):
            mode_lines[n] = line
        else: script_lines[n] = line
    
    # Try parsing a mode from the lmp_script
    if get_mode_from_script:
        tmp_mode = mode_from_script.parse(lmp_script, log)
        if tmp_mode: mode = tmp_mode
        # print( list(mode.keys()) )
        # print( mode['copy_files'] )
        
        
    # Check to see if there is a batch script.  If the file does not 
    # exist, generate a minimal string to represent the batch file.
    batch_script = mode['batch_script']
    if mode['lmp_script_find'] in batch_script:
        batch_script = batch_script.replace(mode['lmp_script_find'], lmp_script)
    if os.path.exists(batch_script):
        batch_encoding = mode_from_script.get_encoding(batch_script)
        with open(batch_script, 'r', encoding=batch_encoding, errors='ignore') as f:
            batch_script_content = f.read()
    else:
        batch_script_content = ''
        
    # Get the auto script path and name
    auto_script = mode['auto_script']
    if mode['lmp_script_find'] in auto_script:
        root_name, root_ext = os.path.splitext(root_filename) 
        auto_script = auto_script.replace(mode['lmp_script_find'], root_name)
    if not auto_script.endswith('.sh'): auto_script += '.sh'
    auto_script_path = os.path.normpath(auto_script)
        
    
    # Change the current working directory to path
    # to write all new files to outputs directory
    pwd = os.getcwd()
    path = os.path.normpath(os.path.join(pwd, root_directory))
    os.chdir(path)
    
    
    # Update values that call other keys in the batch array
    lmp_array = eval_values_from_keys(mode['lmp_array']) 
    batch_array = eval_values_from_keys(mode['batch_array']) 
    
    # 1st pass: Generate basic lmp_values
    lmp_keys = list(lmp_array.keys())
    lmp_values = []
    for i in lmp_array.values():
        try: lst = list(nuttered_eval(i))
        except: lst = [str(i)]
        lmp_values.append( lst )
    
    # Start iterating through the products
    log.out('\n\n\nIterating through products: ')
    lmp_product = list(product(*lmp_values))
    build_directory = mode['build_directory']
    auto_submit = [] # [(submit-dir, submit-script)]
    copy_files = [[i for i in mode['copy_files']] for _ in range(len(lmp_product))] # Need a deep copy to avoid mutation issues
    for i, cartesian in enumerate(lmp_product):
        build_script = copy.deepcopy(mode['build_script'])
        build_dir    = copy.deepcopy(build_directory)
        new_content  = copy.deepcopy(lmp_script_content)
        new_lines    = copy.deepcopy(script_lines)
        files        = copy.deepcopy(copy_files[i])
        for key, value in zip(lmp_keys, cartesian):
            build_dir = build_dir.replace(key, str(value))
            build_script = build_script.replace(key, str(value))
            new_content = new_content.replace(key, str(value))
            for n, file in enumerate(files):
                files[n] = file.replace(key, str(value))
            for n, line in new_lines.items():
                new_lines[n] = line.replace(key, str(value))
                
        # enforce no whitespaces
        build_dir = build_dir.replace(' ', '_') 
        build_script = build_script.replace(' ', '_') 
                
            
        # Build the new directory
        log.out('')
        log.out('Starting {} of {}'.format(i+1, len(lmp_product)))
        log.out(' - Build Directory: {}'.format(build_dir))
        submit_dir = build_dir
        build_dir = os.path.normpath(os.path.join(path, build_dir))
        if not os.path.isdir(build_dir): 
            os.makedirs(build_dir, exist_ok=True)
            
            
        # Move files from root directory to newly built directory
        files_find_and_replace = {} # {'copy_files[N]':filename}
        if files:
            log.out(' - Moving files from root directory: {}'.format(root_directory))
            for n, file in enumerate(files):
                if not os.path.isfile(str(file)):
                    copy_paths = glob.glob(file)
                else:
                    copy_paths = [file]
                if copy_paths:
                    for cfile in copy_paths:
                        if os.path.isfile(str(file)):
                            find_pattern = 'copy_files[{}]'.format(n)
                            files_find_and_replace[find_pattern] = cfile
                        full_path = os.path.abspath(cfile)
                        if not os.path.isfile(full_path): continue
                        if os.path.basename(cfile) == os.path.basename(build_script): continue
                        if os.path.basename(cfile) == os.path.basename(batch_script): continue
                        try:
                            shutil.copy(full_path, build_dir)
                            log.out('   * {}'.format(os.path.basename(cfile)))
                        except: 
                            log.out('   * FAILED to copy: {}'.format(os.path.basename(cfile)))
                            
        # Find any 'copy_files[N]' strings in new_content and replace accordingly
        for key, value in files_find_and_replace.items():
            new_content = new_content.replace(key, str(value))
        
        # Write the new LAMMPS script
        if build_script in ['', ' ']:
            build_script = root_filename
        if mode['lmp_script_find'] in build_script:
            root_name, root_ext = os.path.splitext(root_filename) 
            build_script = build_script.replace(mode['lmp_script_find'], root_name) #+ root_ext
        log.out(' - Writing: {}'.format(os.path.basename(build_script)))
        new_script_path = os.path.normpath(os.path.join(build_dir, build_script))
        with open(new_script_path, 'w', encoding=lmp_encoding, newline=None) as f: 
            try: # Try to regenerate new_content from mode/script lines to avoid mutating a mode
                lines = list(mode_lines.keys()) + list(new_lines.keys())
                lines = sorted(lines)
                rebuilt_new_content = ''
                for n in lines:
                    if n in mode_lines:
                        rebuilt_new_content += '{}\n'.format(mode_lines[n])
                    if n in new_lines:
                        rebuilt_new_content += '{}\n'.format(new_lines[n])
                new_content = rebuilt_new_content
            except: pass
            f.write(new_content)
        
        # Write batch script if user provided a template
        if batch_script_content:
            batch_script = '{}.sh'.format(build_script)
            log.out(' - Writing: {}'.format(os.path.basename(batch_script)))
            batch_script_path = os.path.normpath(os.path.join(build_dir, batch_script))
            auto_submit.append( (submit_dir, os.path.basename(batch_script)) )
            
            new_batch_script_content = batch_script_content
            for key, value in batch_array.items():
                value = str(value)
                if mode['build_script_find'] in value:
                    value = value.replace(mode['build_script_find'], os.path.basename(build_script))
                new_batch_script_content = new_batch_script_content.replace(key, str(value))
                
            with open(batch_script_path, 'w', encoding=batch_encoding, newline='\n') as f: 
                f.write(new_batch_script_content)
            

    # Write the auto submit script 
    if batch_script_content and auto_script:            
        log.out('\n\nWriting auto_script="{}"'.format(auto_script))
        auto_script_template.write(auto_script_path, auto_submit, version, submit=mode['auto_submit'])
    elif auto_script:
        log.out('\n\nCould not write auto_script="{}", as there was not a valid batch script'.format(auto_script))
    else:
        log.out('\n\nIf both a batch template script and auto_script name are supplied an auto_submit script could be generated.')
        
        
    # Write outputs to copy and paste into LAMMPS
    # script to define mode from the LAMMPS script
    lines = mode_from_script.mode2lines(mode)
    log.out('\n\nTo recreate these outputs from a mode defined in a LAMMPS script')
    log.out('copy and paste the following comments into that LAMMPS script:\n#')
    for line in lines:
        if line == '\n': line = ''
        line = '#MODE: {}'.format(line)
        log.out(line)
        
        
    # Print completion of code
    log.out('\n\nNormal program termination\n\n')
    
    # Script run time
    execution_time = (time.time() - start_time)
    log.out('Execution time in seconds: ' + str(execution_time))
    
    # Show number of warnings and errors
    log.out_warnings_and_errors()
    
    # write log
    log.write_logged(root_filename+'.log.lunar')
    
    # Change back to the intial directory to keep directory free for deletion
    os.chdir(pwd)
    return
