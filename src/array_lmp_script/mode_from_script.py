# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
October 21, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.io_functions as io_functions
import os


def get_encoding(path):
    # Detect encoding
    try:
        import chardet
        with open(path, 'rb') as f:
            result = chardet.detect(f.read())
            encoding = result['encoding']
    except: encoding = 'utf-8'
    return encoding





# Function to parse out the mode dictionary from a LAMMPS
# script comment section that starts with '#MODE:'
def parse(lmp_script, log=None):
    
    # Parse the lines that start with '#MODE:',
    # to get the mode from the lmp_script file
    encoding = get_encoding(lmp_script)
    lines = []
    with open(lmp_script, 'r', encoding=encoding, errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#MODE:'):
                code = line.replace('#MODE:', '')
                code = code.strip()
                lines.append(code)
                

    # Try getting the mode from the pasred mode lines
    basename = os.path.basename(lmp_script)
    mode, success = {}, False
    if lines:
        if log is not None: log.out(f'\n\nGetting mode from: {basename} ')
        try:
            # Safely use exec() with no built-ins to get
            # the mode dictionary from the script
            globals_dict = {"__builtins__": {}}  # no built-ins
            locals_dict = {} 
            exec('\n'.join(lines), globals_dict, locals_dict)
            mode = locals_dict['mode']
            success = True
        except: 
            if log is not None:
                log.out(f' - Failed to parse mode from: {basename}')
    
    # Check that the mode has all needed keys
    if success:
        required = ['auto_script', 'build_script', 'build_directory',
                    'auto_submit', 'lmp_array', 'batch_array']
        if set(required) <= mode.keys():
            mode['lmp_script'] = lmp_script
            if log is not None: log.out(f' - Successfully parsed mode from: {basename}')
        else: mode = {}
    return mode


# Function to format a mode dictionary to a clean printable string
def mode2lines(mode):
    lines = []
    lmp_array = mode['lmp_array']
    if len(lmp_array) > 0:
        string = 'lmp_array = {'
        for n, find in enumerate(lmp_array):
            if n == 0:
                string += "'{}' : '{}'".format(str(find), str(lmp_array[find]))
            else:
                string = "'{}' : '{}'".format(str(find), str(lmp_array[find]))
                
            if n+1 < len(lmp_array): string += ','
            else: string += '}'
            
            if n == 0: 
                out = string
            else:
                out = "{:>13}{}".format('', string)
            lines.append(out)
    else: lines.append('lmp_array = {}')
    lines.append('\n')
    
    batch_array = mode['batch_array']
    if len(batch_array) > 0:
        string = 'batch_array = {'
        for n, find in enumerate(batch_array):
            if n == 0:
                string += "'{}' : '{}'".format(str(find), str(batch_array[find]))
            else:
                string = "'{}' : '{}'".format(str(find), str(batch_array[find]))
                
            if n+1 < len(batch_array): string += ','
            else: string += '}'
            
            if n == 0: 
                out = string
            else:
                out = "{:>15}{}".format('', string)
            lines.append(out)
    else: lines.append('batch_array = {}')
    lines.append('\n')
    
    copy_files = mode['copy_files']
    if len(copy_files) > 0:
        string = 'copy_files = ['
        for n, file in enumerate(copy_files):
            if n == 0:
                string += "'{}'".format(str(file))
            else:
                string = "'{}'".format(str(file))
                
            if n+1 < len(copy_files): string += ','
            else: string += ']'
            
            if n == 0: 
                out = string
            else:
                out = "{:>14}{}".format('', string)
            lines.append(out)
    else: lines.append('copy_files = []')
    lines.append('\n')
    
    lines.append("{}{}: '{}',".format('mode = {', "'lmp_script'", io_functions.path_to_string(mode['lmp_script']) ))
    lines.append("{:>8}{}: '{}',".format('', "'batch_script'", io_functions.path_to_string(mode['batch_script']) ))
    lines.append("{:>8}{}: '{}',".format('', "'auto_script'", mode['auto_script'] ))
    lines.append("{:>8}{}: '{}',".format('', "'build_script'", mode['build_script'] ))
    lines.append("{:>8}{}: '{}',".format('', "'build_directory'", mode['build_directory'] ))
    lines.append("{:>8}{}: '{}',".format('', "'lmp_script_find'", mode['lmp_script_find'] ))
    lines.append("{:>8}{}: '{}',".format('', "'build_script_find'", mode['build_script_find'] ))
    lines.append("{:>8}{}: '{}',".format('', "'auto_submit'", mode['auto_submit'] ))
    lines.append("{:>8}{}: {},".format('', "'lmp_array'", 'lmp_array' ))
    lines.append("{:>8}{}: {},".format('', "'batch_array'", 'batch_array'))
    lines.append("{:>8}{}: {}{}".format('', "'copy_files'", 'copy_files', '}' ))
    
    # for line in lines:
    #     print(line)
    return lines