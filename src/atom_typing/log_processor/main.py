# -*- coding: utf-8 -*-
"""
Revision 1.0
February 20, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.atom_typing.log_processor.read_atom_typing as read_atom_typing
import src.atom_typing.log_processor.keywords as keywords
import src.glob_wildcards as glob_wildcards
import src.io_functions as io_functions
import matplotlib.pyplot as plt
import glob
import time
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


##################################################################
# Function meant for expanding a table find all logger keywords  #
# in that table. This is used when a keyword in logger is set to #
# 'expand-table:TABLENAME', where 'TABLENAME' could be:          #
#    'rings'           -> 'expand-table:rings'                   #
#    'clusters'        -> 'expand-table:clusters'                #
#    'hybridizations'  -> 'expand-table:hybridizations'          #
#    'ringed_clusters' -> 'expand-table:ringed_clusters'         #
#    :                 -> :                                      #
##################################################################
def expand_table(table, dict_name):
    all_table = []
    for valueID0 in table:
        for valueID1 in table[valueID0]:
            # ring's dictionary may be nested  "3-deep"
            if isinstance(table[valueID0][valueID1], dict):
                for valueID2 in table[valueID0][valueID1]:
                    if isinstance(valueID0, str):
                        index0 = "'{}'".format(valueID0)
                    else: index0 = "{}".format(valueID0)
                    if isinstance(valueID1, str):
                        index1 = "'{}'".format(valueID1)
                    else: index1 = "{}".format(valueID1)
                    if isinstance(valueID2, str):
                        index2 = "'{}'".format(valueID2)
                    else: index2 = "{}".format(valueID2)
                    all_table.append("{}[{}][{}][{}]".format(dict_name, index0, index1, index2))
            else: # every other dictionary will be nested "2-deep"
                if isinstance(valueID0, str):
                    index0 = "'{}'".format(valueID0)
                else: index0 = "{}".format(valueID0)
                if isinstance(valueID1, str):
                    index1 = "'{}'".format(valueID1)
                else: index1 = "{}".format(valueID1)
                all_table.append("{}[{}][{}]".format(dict_name, index0, index1))
    return all_table

################################################################
# Function to parse logger string for keywords. This function: #
#  - checks and removes any duplicates                         #
#  - expands tables via the 'expand-table:TABLENAME' string    #
################################################################
def parse_logger(string, logfiles, log):
    # String newline characters from start and end of logger and the replace
    # all newline characters with a ',' incase users forget to place a comma 
    # at the end of a line
    string = string.strip()
    string = string.replace('\n', ',') 
    
    # Start parsing
    header = []
    for word in string.split(','):
        if word == '': continue
        word = word.strip()
        # use expand table method
        if word.startswith('expand-table:'):
            table = word.split(':')[-1].strip()
            if not table: continue
            for file in logfiles:
                table_dict = getattr(logfiles[file], table)
                try:
                    all_table = expand_table(table_dict, table)
                    for keyword in all_table:
                        if keyword not in header:
                            header.append(keyword)
                except: log.warn(f'WARNING {word} failed for some reasons')
                    
        # simply append word if not already logged
        else:
            if word not in header and word != '':
                header.append(word)
            else: log.warn(f'WARNING {word} was duplicated in logger. Using first instance of {word}')
    return header


#########################################################################
### Main function to perform all atom-typing logfile processing tasks ###
#########################################################################
def main(mode, savefig=True, dpi=300, log_clear=True, log=None): 
    start_time = time.time()
    
    #------------------------------------------------------------------------------#
    # Configure log (default is level='production', switch to 'debug' if debuging) #
    #------------------------------------------------------------------------------#
    if log is None: log = io_functions.LUNAR_logger()
    log.configure(level='production', print2console=True, write2log=True)
    if log_clear: log.clear_all()
    
    #--------------------------------------#
    # Read in all log files into logs list #
    #--------------------------------------#
    files = glob.glob(mode['logfile']); logfiles = []
    if files:
        logfiles = {} # {filename: logfile object}
        for file in files:
            root = file[:file.rfind('.')]
            if root == mode['newfile']:
                log.warn(f'WARNING matched newfile {mode["newfile"]} name matched with a basename log file and was skipped')
                continue
            
            # Read logfile and perform some wildcard magic and get the most numeric value from the wildcards
            logfile = read_atom_typing.logfile(file)
            wildcards = glob_wildcards.get_glob_wildcards(mode['logfile'], file)
            numeric = glob_wildcards.get_numeric_wildcard(wildcards)
            logfile.filename = os.path.basename(file)
            logfile.wildcards = wildcards
            logfile.numeric = numeric 
            logfiles[file] = logfile
        
    else: log.warn(f'WARNING logfile: {mode["logfile"]} unwrapped zero files')
    
    #------------------------#
    # Sort files accordingly #
    #------------------------#
    # Sort files in ascending order
    sorted_files = list(logfiles.keys())
    direction = mode['sorting_direction']
    if mode['sorting_method'] == 'sort':
        log.out(f'Sorting {len(sorted_files)} loaded files in {direction} direction, with Python "sorted()" function ...')
        sorted_files = sorted(sorted_files)
    elif mode['sorting_method'] == 'natsort':
        log.out(f'Sorting {len(sorted_files)} loaded files in {direction} direction, with natsorts "natsorted()" function ...')
        from natsort import natsorted
        sorted_files = natsorted(sorted_files)
    elif mode['sorting_method'] == 'numericsort':
        log.out(f'Sorting {len(sorted_files)} loaded files in {direction} direction, with using the numeric option ...')
        sorted_files = sorted(sorted_files)
        numeric_sort = {filename:logfiles[filename].numeric for filename in sorted_files} # {'filename' : numeric-value}
        numeric_sort = dict(sorted(numeric_sort.items(), key=lambda x: x[1])) # [0=keys;1=values]
        sorted_files = list(numeric_sort.keys())
    elif 'wildcards' in mode['sorting_method'] and '[' in  mode['sorting_method'] and ']' in  mode['sorting_method']:
        sorted_files = sorted(sorted_files); wildcard_sort = {} # {'filename' : wildcard-value}
        for value, filename in enumerate(sorted_files):
            index = glob_wildcards.get_value_between_characters(mode['sorting_method'])
            if index and isinstance(index, int):
                try: value = logfiles[filename].wildcards[index]
                except: pass
            wildcard_sort[filename] = value
        wildcard_sort = dict(sorted(wildcard_sort.items(), key=lambda x: x[1])) # [0=keys;1=values]
        sorted_files = list(wildcard_sort.keys())
        log.out(f'Sorting {len(sorted_files)} loaded files in {direction} direction, with using the wildcards[INDEX] option ...')
    else:
        log.error(f'ERROR sorting method {mode["sorting_method"]} is not supported')
        
    # Reverse the list if user wants descending order
    if direction == 'descending':
        sorted_files.reverse()
        
    #------------------------------------------------------------------#
    # Setup the globals namespace to limit scope of what eval() can do #
    #------------------------------------------------------------------#
    allowed_builtins = ['min','max','sum','abs','len','map','range','reversed']
    copied_builtins = globals()['__builtins__'].copy()
    globals_dict = {'__builtins__': {key:copied_builtins[key] for key in allowed_builtins}}
        
    #-------------------------------------------#    
    # Parse logger string and get values to log #
    #-------------------------------------------#
    log.out('Checking logger for consistency ...')
    header = parse_logger(mode['logger'], logfiles, log)
    
    # Go through and evaluate headers
    values = []; alias = 'logfile' # need to set alias for eval()
    for file in sorted_files:
        logfile = logfiles[file]; logfile_values = []
        for string in header:
            string = string.rstrip()
            string = string.strip()            
            value = None
            if string in keywords.shortcuts:
                string = keywords.shortcuts[string]
            if string.startswith('CharYield:'):
                try:
                    mass_initial = float(string.split(':')[-1])
                    value = round(100*(logfile.mass/mass_initial), 4)
                except: pass
            else:
                string2eval = '{}.{}'.format(alias, string)
                globals_dict[alias] = logfile
                try: value = eval(string2eval, globals_dict)
                except: pass
            logfile_values.append(value)
        values.append(logfile_values)
        
    #-----------------------------------------------------------------------#
    # Setting up directories and where to  write final files and results to #
    #-----------------------------------------------------------------------#
    # Find present working directory
    pwd = os.getcwd()
    
    # Find/create paths to store code results
    path = os.path.join(pwd, mode['parent_directory'])
    
    # If parent_directory == 'topofile' use topofile path as path
    if 'logfile' in mode['parent_directory']:
        log.out('Using path from 1st file in glob file list to set parent_directory ...')
        first_file = list(logfiles.keys())[0]
        
        # get_dir_from_topofile() uses 'topofile' instead of 'logfile'. So need to replace
        parent_directory = mode['parent_directory']
        parent_directory = parent_directory.replace('logfile', 'topofile') 
        path = io_functions.get_dir_from_topofile(first_file, parent_directory)
        
    # Check if path exits. IF not create
    if not os.path.isdir(path): 
        os.makedirs(path, exist_ok=True)
        
    # Change the current working directory to path so all files get written to that directory
    os.chdir(path)
        
    #------------------------#
    # Write logged .csv file #
    #------------------------#
    try:
        with open(mode['newfile']+'.csv', 'w') as f:
            f.write('{}\n'.format(', '.join(header)))
            for logfile_values in values:
                row = []
                for value in logfile_values:
                    if isinstance(value, (float, int)):
                        value = str(value)
                    elif isinstance(value, (list, tuple)):
                        value = ' '.join([str(i) for i in value])
                    else: value = str(value)
                    row.append(value)
                f.write('{}\n'.format(', '.join(row)))
    except: log.warn(f'WARNING could not save {mode["newfile"]+".csv"}. Likely opened by another program')
            
    #-------------------------------#
    # Invert data and generate plot #
    #-------------------------------#
    # Invert data
    index_map = {n:i for n, i in enumerate(header)}
    data = {i:[] for i in header}
    for logfile_values in values:
        for n, value in enumerate(logfile_values):
            data[index_map[n]].append(value)
    
    # Set up x and y plotting
    dependent_variables = mode['dependent_variables']
    independent_variable = mode['independent_variable']
    if dependent_variables and independent_variable:
        log.out('\n\nSetting up plotting option ...')
        if 'logger' in dependent_variables:
            yvalues = header
        else: 
            yvalues = parse_logger(dependent_variables, logfiles, log)
        
        # Start plotting
        fig, ax = plt.subplots(figsize=(6, 4))
        for string in yvalues:
            xvalue = independent_variable.strip()
            yvalue = string.strip()
            if xvalue == yvalue: continue
            try:
                x = data[xvalue]; y = data[yvalue]
                if len(x) != len(y):
                    continue
                elif not all(isinstance(i, (int, float)) for i in x) or not all(isinstance(i, (int, float)) for i in y):
                    continue
                else:
                    ax.plot(x, y, '-', lw=2, label=string)
            except:
                log.warn(f'WARNING X={independent_variable} or Y={string} could not be plotted. This may happen for a few reasons:')
                log.out('  - Perhaps something was mispelled')
                log.out('  - Perhaps it wasnt in the logger')
                log.out('  - Perhaps values are non-numeric')
                log.out('  - Perhaps length of X- and Y- arrays are different')
        try:
            ax.set_xlabel(independent_variable)
            ax.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
        except: pass
        plt.show()
        if savefig:
            try: fig.savefig(mode['newfile']+'.jpeg', dpi=dpi)
            except: log.warn(f'WARNING could not save {mode["newfile"]+".jpeg"}. Likely opened by another program')


    #--------------------------------------#
    # Wrap up the execution of this script #
    #--------------------------------------#    
    
    # Print file locations
    log.out(f'\n\nAll outputs can be found in {path} directory')
    
    # Print completion of code
    log.out('\n\nNormal program termination\n\n')
    
    # Script run time
    execution_time = (time.time() - start_time)
    log.out('Execution time in seconds: ' + str(execution_time))
    
    # Show number of warnings and errors
    log.out_warnings_and_errors()
    
    # Write log file
    try: log.write_logged(mode['newfile']+'.log.lunar')
    except: log.warning(f'WARNING could not save {mode["newfile"]+".log.lunar"}. Likely opened by another program')
            
    # Change back to the intial directory to keep directory free for deletion
    os.chdir(pwd)
    return header, values, data