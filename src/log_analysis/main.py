# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.5
April 28, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.log_analysis.vectorized_string_compute as vectorized_string_compute
import src.log_analysis.whittaker_smoother as whittaker_smoother
import src.log_analysis.signal_processing as signal_processing
import src.log_analysis.rfr_mechanical as rfr_mechanical
import src.log_analysis.rfr_thermal as rfr_thermal
import src.log_analysis.misc_funcs as misc_funcs
import src.log_analysis.read_log as read_log
import src.io_functions as io_functions
import src.glob_wildcards as glob_wildcards
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from cycler import cycler
import numpy as np
import traceback
import copy
import glob
import math
import time
import os


#colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#color_index = 0

# Color wheel defined by matplotlib:
#   colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
# However, we can construct our own color wheel to prioritize the colors we want first and we can 
# have way more colors defined than what matplotlib defines. Below is Josh's preferred color wheel
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple', 'tab:red', 'tab:gray','tab:olive', 'tab:cyan', 'tab:pink', 'teal',
          'crimson', 'lime', 'tomato',  'blue', 'orange', 'green', 'purple', 'red', 'gray', 'olive', 'cyan', 'pink', 'tab:brown']

# Get the colors from the 'tab20c' colormap
#colors = plt.cm.tab20c.colors

# Set the default color cycle using rcParams
plt.rcParams['axes.prop_cycle'] = cycler(color=colors)


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


#----------------------------#
# Array processing functions #
#----------------------------#
# Function that parses the Misc column in the stack, but useful for other parsing
def get_misc_setting(misc, use_eval=True):
    # Setup the globals namespace to limit scope of what eval() can do
    allowed_builtins = ['min','max','sum','abs','len','map','range','reversed']
    copied_builtins = globals()['__builtins__'].copy()
    globals_dict = {}
    globals_dict['__builtins__'] = {key:copied_builtins[key] for key in allowed_builtins}
    
    # Parse misc string
    setting = {} # {keyword:float or int or Boolean}
    tmp1 = misc.split(';')
    for tmp2 in tmp1:
        tmp3 = tmp2.split('=')
        if len(tmp3) >= 2:
            i = tmp3[0].strip()
            if use_eval:
                try: j = eval(tmp3[1], globals_dict)
                except: j = str(tmp3[1])
            else: j = str(tmp3[1])
            setting[i] = j
    return setting


#----------------------------#
# Array processing functions #
#----------------------------#
# Function to perform averaging over multiple wildcard files and then
# update the mode and files structs to analyze the single file.
def wildcard_averaging(mode, files, log):
    # Parse:
    #   sections = "1,2,3-6; wavg=True"
    sections = mode['sections']
    sections_only, additional_settings = sections.split(';', 1)
    sections_only       = sections_only.strip()
    additional_settings = additional_settings.strip()
    settings = get_misc_setting(additional_settings)

    # If user wants wildcard averaging go through and collect all series to average together
    wildcard_averaging = settings.get('wavg', False)
    if not wildcard_averaging: # by-pass if user did ask for this
        mode['sections'] = sections_only
        return mode, files
    
    # Start averaging together files
    import numpy as np
    xdata_column = mode['xdata']
    ydata_column = mode['ydata']
    xcompute     = mode['xcompute']
    ycompute     = mode['ycompute']
    keywords     = mode['keywords']
    sections     = sections_only
    xdata_series, ydata_series, processed_files = [], [], []
    for logfile in files:
        if not os.path.isfile(logfile): continue
    
        logged = read_log.file(logfile, keywords=keywords, pflag=log.print2console)
        data = logged.get_data(sections, remove_duplicates=True, pflag=log.print2console) # {column-name:[lst of data]}
        
        if xcompute:
            xdata, message = vectorized_string_compute.thermo_data(data, xcompute)
            log.out('  xcompute={}'.format(xcompute))
            log.out('  {}'.format(message))
        else: xdata = data[xdata_column]
        
        if ycompute:
            ydata, message = vectorized_string_compute.thermo_data(data, ycompute)
            log.out('  ycompute={}'.format(ycompute))
            log.out('  {}'.format(message))
        else: ydata = data[ydata_column]
            
        xdata = np.array(xdata, dtype=float)
        ydata = np.array(ydata, dtype=float)
        
        # Ensure increasing x for interpolation
        order = np.argsort(xdata)
        xdata = xdata[order]
        ydata = ydata[order]
        xdata_series.append( xdata )
        ydata_series.append( ydata )
        processed_files.append( logfile )
        
    if not xdata_series:
        log.warn('WARNING wildcard averaging found no valid files.')
        return mode, files

    # Common overlapping x-domain
    xmin = max(x.min() for x in xdata_series)
    xmax = min(x.max() for x in xdata_series)
    if xmin >= xmax:
        log.warn('WARNING wildcard averaging failed because X-ranges do not overlap.')
        return mode, files
    
    # Use the shortest series length to avoid inventing extra resolution
    npts = min(len(x) for x in xdata_series)
    xavg = np.linspace(xmin, xmax, npts)
    
    # Interpolate all y-series onto common x-grid
    y_interp = []
    for xdata, ydata in zip(xdata_series, ydata_series):
        y_interp.append(np.interp(xavg, xdata, ydata))
    
    yarr = np.vstack(y_interp)
    yavg = np.mean(yarr, axis=0)
    if yarr.shape[0] > 1:
        ystd = np.std(yarr, axis=0, ddof=1)
    else:
        ystd = np.zeros_like(yavg)
    
    # Write averaged LAMMPS-style log file in path of last processed file
    last_file = processed_files[-1]
    outdir = os.path.dirname(os.path.abspath(last_file))
    
    array_file = mode.get('array_file', '')
    if array_file:
        basename = array_file
    else:
        basename = 'wildcard_average'
    
    avg_logfile = os.path.join(outdir, '{}.log.lammps'.format(basename))
    if xcompute: xname = 'xavg'
    else: xname = xdata_column
    if ycompute: yname = 'yavg'
    else: yname = ydata_column
    with open(avg_logfile, 'w') as f:
        f.write('LAMMPS log file generated by LUNAR/log_analysis.py wildcard averaging\n')
        f.write('Averaged {} files:\n'.format(len(processed_files)))
        for file in files:
            if file in processed_files:
                f.write(' - {:>8}: {}\n'.format('Averaged', file))
            else:
                f.write(' - {:>8}: {}\n'.format('Failed', file))
        f.write('\n\n')
        
        f.write('{:>10} {:>24} {:>24} {:>24}\n'.format('Step', xname, yname, yname + '_std'))
        for i, (x, y, s) in enumerate(zip(xavg, yavg, ystd)):
            f.write('{:>10d} {:>24.16f} {:>24.16f} {:>24.16f}\n'.format(i, x, y, s))
    
        f.write('Loop time of 0 on 1 procs for {} steps\n'.format(len(xavg)))
    
    log.out(' - wildcard averaged logfile={}'.format(avg_logfile))
    
    # Update mode and files for the new averaged file
    mode['sections'] = 'all'
    mode['xdata']    = xname
    mode['ydata']    = yname
    mode['logfile']  = avg_logfile
    mode['xcompute'] = '' # reset any computes since this is already handled
    mode['ycompute'] = '' # reset any computes since this is already handled
    files = [avg_logfile]
    return mode, files


#----------------------------#
# Array processing functions #
#----------------------------#
# Partial name matching for wildcard switching
def partial_name(compute):
    # Strip function_name and parentheses
    stripped = ''
    left_parentheses = False
    right_parentheses_removed = compute.strip()
    if right_parentheses_removed.endswith(')'):
        right_parentheses_removed = right_parentheses_removed[:-1]
    for i in right_parentheses_removed:
        if i == '(':
            left_parentheses = True
            continue
        if left_parentheses:
            stripped += i
    
    # Generate map
    partial_to_variable = {} # {'partial_name':['variable', {'keyword':'key_var'}], ... }
    split = stripped.split(',')
    for pair in split:
        if not pair: continue
        
        # We need to support semi-colon stringing operations
        # to tie back to a specific analysis run
        semi_colon = pair.split(';', 1)
        before = semi_colon[0]
        after  = semi_colon[1:]
        if after:
            after = ' '.join(after)
            after_semi_colon = after.split(';')
            analysis_dict = {}
            for string in after_semi_colon:
                equate = string.split('=')
                if len(equate) == 2:
                    keyword = equate[0].strip()
                    key_var = equate[1].strip()
                    analysis_dict[keyword] = key_var
        else: analysis_dict = {}
        
        # Determine partial_to_variable
        paired = before.split('=')
        if len(paired) == 2:
            partial = paired[0].strip()
            variable = paired[1].strip()
            partial_to_variable[partial] = [variable, analysis_dict]
    return partial_to_variable


#----------------------------#
# Array processing functions #
#----------------------------#
# Wildcard matching for wildcard switching
def wild_match(compute):
    # Strip function_name and parentheses
    stripped = ''
    left_parentheses = False
    right_parentheses_removed = compute.strip()
    if right_parentheses_removed.endswith(')'):
        right_parentheses_removed = right_parentheses_removed[:-1]
    for i in right_parentheses_removed:
        if i == '(':
            left_parentheses = True
            continue
        if left_parentheses:
            stripped += i
    
    # Generate map
    wildcard_to_variable = {} # {('wildcard[i]', 'match'):['variable', {'keyword':'key_var'}]}
    split = stripped.split(',')
    for pair in split:
        if not pair: continue
    
        # We need to support semi-colon stringing operations
        # to tie back to a specific analysis run
        semi_colon = pair.split(';', 1)
        before = semi_colon[0]
        after  = semi_colon[1:]
        if after:
            after = ' '.join(after)
            after_semi_colon = after.split(';')
            analysis_dict = {}
            for string in after_semi_colon:
                equate = string.split('=')
                if len(equate) == 2:
                    keyword = equate[0].strip()
                    key_var = equate[1].strip()
                    analysis_dict[keyword] = key_var
        else: analysis_dict = {}
        
        # Determine partial_to_variable
        paired = before.split('=')
        if len(paired) == 3:
            wildcard = paired[0].strip()
            match    = paired[1].strip()
            variable = paired[2].strip()
            wildcard_to_variable[(wildcard, match)] = [variable, analysis_dict]
    return wildcard_to_variable


#----------------------------#
# Array processing functions #
#----------------------------#
# Write out the csv file for array processing
def array_csv(outputs, abouts, wildcards, basename, log, max_len=10):
    # sub_analysis to skip
    sub2skip = ['xdata-raw', 'ydata-raw', 'xdata-clean', 'ydata-clean', 'xdata', 'ydata',
                'lo_bounds', 'hi_bounds', 'lo_line', 'hi_line', 'lo_reg_x', 'hi_reg_x', 'lo_reg_y', 'hi_reg_y']
    
    # Find all unique columns
    all_wildcards = set() # {Nwildcards1, Nwildcards2, ...}
    all_outputs   = []    # ['analyis|sub-analysis',  'analyis|sub-analysis[index]', ...]
    all_abouts    = {}    # {output-name:about-text}
    for file in outputs:
        all_wildcards.add(len(wildcards[file]))
        for analysis in outputs[file]:
            for sub_analysis in outputs[file][analysis]:
                if sub_analysis in sub2skip: continue
                column = '{}|{}'.format(analysis, sub_analysis)    
                data = outputs[file][analysis][sub_analysis]
                try: about = abouts[file][analysis][sub_analysis]
                except: about = ''
                if isinstance(data, (list, tuple)):
                    if len(data) > max_len: continue
                    for n, value in enumerate(data):
                        if value is None: continue
                        list_column = '{}[{}]'.format(column, n)
                        if list_column not in all_outputs:
                            all_outputs.append(list_column)
                else:
                    if data is None: continue
                    if column not in all_outputs:
                        all_outputs.append(column)
                
                # Log here in-case of a continue due to data being
                # X-data or Y-data (or similar such as smoothed data)
                all_abouts[column] = about

    # Start making unique data structure to hold all possible combinations
    # of wildards and computes - initialize all with None's in case value
    # is missing.
    all_outputs = sorted(all_outputs, key=lambda s: s[0])
    columns = ['wildcards[{}]'.format(i) for i in range(max(all_wildcards))]
    columns.extend( list(all_outputs) )
    csv_data = {} # {'filename':{'wildcards[0]':'wildcard', 'Modulus|b1-raw':value, ... Nanalysis*Nsub-analysis},  ... Nfiles}
    for file in outputs:
        # Initialize all columns
        csv_data[file] = {i:str(None) for i in columns}
        
        # Upate wildcards
        for n, wildard in enumerate(wildcards[file]):
            column = 'wildcards[{}]'.format(n)
            csv_data[file][column] = wildard
        
        # Update analysis|sub-analysis data
        for analysis in outputs[file]:
            for sub_analysis in outputs[file][analysis]:
                if sub_analysis in sub2skip: continue
                column = '{}|{}'.format(analysis, sub_analysis)    
                data = outputs[file][analysis][sub_analysis]
                if isinstance(data, (list, tuple)):
                    if len(data) > max_len: continue
                    for n, value in enumerate(data):
                        if value is None: continue
                        list_column = '{}[{}]'.format(column, n)
                        value = str(value)
                        value = value.replace(',', '-')
                        csv_data[file][list_column] = value
                else:
                    if data is None: continue
                    value = str(data)
                    value = value.replace(',', '-')
                    csv_data[file][column] = value

        
    # Write output info
    try:
        output_name = '{}_output.csv'.format(basename)
        with open(output_name, 'w') as f:
            # Join with comma's and write titles
            titles = ', '.join(['filename'] + columns)
            f.write('{}\n'.format(titles))
            
            # Write rows
            for file in csv_data:
                data = csv_data[file]
                outputs = [str(file)]
                for column in columns:
                    outputs.append(str(data[column]))
                f.write('{}\n'.format(', '.join(outputs)))
        log.out(f'Successfully wrote: {output_name}')
    except: 
        stack_trace_string = traceback.format_exc()
        log.out(f'ERROR could not write {output_name}. File is likely open in another program.')
        log.out('Traceback:')
        log.out(stack_trace_string)

        
    # Write about info
    try:
        about_name = '{}_about.txt'.format(basename)
        with open(about_name, 'w') as f:
            f.write('This file contains details about the array_file header naming. Please note that if the about\n')
            f.write('details is a list, the actual column name in the corresponding *_output.csv file will be \n')
            f.write('indexed like a python list (starting at index 0). For example:\n')
            f.write('\n')
            f.write('  Calling "Regression Fringe Response Modulus" and setting the Name as RFR-mechanical\n')
            f.write('  the yield point computed via derivatives will have a columns:\n')
            f.write('    RFR-mechanical|yield_point_derivative[0] -> X-location of yield point\n')
            f.write('    RFR-mechanical|yield_point_derivative[1] -> Y-location of yield point\n')
            f.write('\n')
            f.write('  Since the about statement defines this value as a list:\n')
            f.write('    List of floats found from yield point determination, using derivative methods.\n')
            f.write('    Order [X-yp, Y-yp] or [None, None] if not using yp. NOTE: SHIFTED based on\n')
            f.write('    shift method\n')
            
            if all_abouts:
                longest_key = max(all_abouts, key=len)
                max_len = len(longest_key)
                top_bottom = (max_len + 4)*'*'
                text_len = 2*max_len
                for key in all_abouts:
                    about = all_abouts[key]                
                    header = '* {:^{width}} *'.format(key, width=max_len)
                    wrapped = '\n'.join(about[i:i+text_len].lstrip() for i in range(0, len(about), text_len))
                    f.write('\n\n')
                    f.write('{}\n'.format(top_bottom))
                    f.write('{}\n'.format(header))
                    f.write('{}\n'.format(top_bottom))
                    f.write('{}\n'.format(wrapped))
        log.out(f'Successfully wrote: {about_name}')
    except: 
        stack_trace_string = traceback.format_exc()
        log.out(f'ERROR could not write {about_name}. File is likely open in another program.')
        log.out('Traceback:')
        log.out(stack_trace_string)
    return

#---------------------------------#
# Array processing final function #
#---------------------------------#
def run(mode, plot=True, savefig=True, dpi=300, log=None, log_clear=False, ding=True):
    if log is None:
        log = io_functions.LUNAR_logger()
        log.configure(level='production')

    # Set up Tristan's "array" analysis using recursion
    if not os.path.isfile(str(mode['logfile'])):
        files = glob.glob(mode['logfile'])
        array_time = time.time(); analyzed = None; 
        figures   = {} # {'figure-description':[fig_instance, ..., Nprocessed_files]}
        wildcards = {} # {filename : [wildcard1, wildcard2, ...]}
        outputs   = {} # {filename : dict-of-outputs}
        abouts    = {} # {filename : dict-of-outputs}
        dirnames  = [] # [dirname1, dirname2, ...]
        processed_files, failed_files = [], []
        
        # Check if user wants all matched files to be averaged together based on the sections format of
        #   sections = "1,2,3-6; wavg=True; other_possible_settings=1"
        # where additional settings appear after the first semi-colon character (;) 
        if ';' in mode['sections']:
            mode, files = wildcard_averaging(mode, files, log)
        if files:
            # Try sorting files
            try:
                from natsort import natsorted
                files = natsorted(files)
                log.out('Read in logfiles were sorted using natsort, thus they should be iterated through')
                log.out('in a sequence that represents evolution or progression and logged in that order.')
            except:
                files = sorted(files)
                log.warn('WARNING natsort was not installed and read in logfiles were sorted using pythons')
                log.out('sorted function. The file iteration sequence may be unorder from a evolution or')
                log.out('progression stand point. Depending on the names given to the read in files. You may')
                log.out('intall natsort with:  pip3 install natsort   (if pip manager is installed)')
            
            outer_mode = copy.deepcopy(mode)
            for n, file in enumerate(files, 1):
                plt.close('all')
                log.clear_all()
                log.out('\n\nUsing array input option:')
                log.out(' - logfile       : {}'.format(mode['logfile']))
                log.out(' - matched file  : {}'.format(file))
                log.out(' - elapsed time  : {:.2f} (seconds)'.format(time.time() - array_time))
                log.out(' - progress      : {} of {} ({:.2f}%)'.format(n, len(files), 100*(n/len(files))))
                if file.endswith('.jpeg') or file.endswith('.log.lunar') or file.endswith('.eps'):
                    log.warn(f' - WARNING matched file {file} has an extension that does not make sense to process. Skipping file')
                    continue
                
                # We will be mutating the mode, so make a deep copy to avoid
                # weird mutation issues by making incremental modifications
                inner_mode = copy.deepcopy(outer_mode)
                inner_mode['logfile'] = file
                wilcards_lst = glob_wildcards.get_glob_wildcards(mode['logfile'], file)
                
                # Check if user wants to use an X-compute and Y-compute to set 
                # X-data and Y-data columns based on info in the logfile name
                xcompute = inner_mode['xcompute'].strip()
                ycompute = inner_mode['ycompute'].strip()
                if xcompute.startswith('partial_name(') and ycompute.startswith('partial_name('):
                    try:
                        xpartial_to_variable = partial_name(xcompute)
                        ypartial_to_variable = partial_name(ycompute)
                        same_keys = set(xpartial_to_variable.keys()) == set(ypartial_to_variable.keys())
                        if same_keys:
                            parital_matched_key = None
                            for key in xpartial_to_variable:
                                if key in file:
                                    parital_matched_key = key
                                    break
                            
                            # Set up variables
                            if parital_matched_key is not None:
                                x_variable, x_analysis_dict = xpartial_to_variable[parital_matched_key]
                                y_variable, y_analysis_dict = ypartial_to_variable[parital_matched_key]
                                log.out(' - partial_name(): {} -> {}'.format(xcompute, x_variable))
                                log.out(' - partial_name(): {} -> {}'.format(ycompute, y_variable))
                                
                                # Update the xdata, ydata, xcompute, and ycompute keys
                                compute_chars = ['${', '}']
                                if all(c in x_variable for c in compute_chars):
                                    inner_mode['xcompute'] = x_variable
                                else:
                                    inner_mode['xcompute'] = ''
                                    inner_mode['xdata'] = x_variable
                                    
                                if all(c in y_variable for c in compute_chars):
                                    inner_mode['ycompute'] = y_variable
                                else:
                                    inner_mode['ycompute'] = ''
                                    inner_mode['ydata'] = y_variable
                                
                                # Update the analysis columns
                                analysis_dict = {**x_analysis_dict, **y_analysis_dict}
                                for n, (method, xlo, xhi, misc, name) in enumerate(inner_mode['analysis']):
                                    settings = get_misc_setting(misc, use_eval=True)
                                    for key in analysis_dict:
                                        if key in settings:
                                            log.out(' - partial_name("{}", {}): {} -> {}'.format(method, key, settings[key], analysis_dict[key]))
                                            settings[key] = analysis_dict[key]
                                    new_misc = '; '.join( ['{}={}'.format(key, value) for key, value in settings.items()] )
                                    inner_mode['analysis'][n][3] = new_misc 
                        else:
                            log.out(' - partial_name() x != y: {} != {}'.format(xcompute, ycompute))
                            log.out(' - Skipping and hoping for the best (using supplied "xdata" and "ydata" columns)')
                    except: pass
                
                if xcompute.startswith('wild_match(') and ycompute.startswith('wild_match('):
                    try:
                        xwildcard_to_variable = wild_match(xcompute)
                        ywildcard_to_variable = wild_match(ycompute)
                        same_keys = set(xwildcard_to_variable.keys()) == set(ywildcard_to_variable.keys())
                        if same_keys:
                            wilcards_set = {('wildcards[{}]'.format(n), str(value)) for n, value in enumerate(wilcards_lst)}
                            wildcard_matched_key = None
                            for key in xwildcard_to_variable:
                                if key in wilcards_set:
                                    wildcard_matched_key = key
                                    break
    
                            # Set up variables
                            if wildcard_matched_key is not None:
                                x_variable, x_analysis_dict = xwildcard_to_variable[wildcard_matched_key]
                                y_variable, y_analysis_dict = ywildcard_to_variable[wildcard_matched_key]
                                log.out(' - wild_match()  : {} -> {}'.format(xcompute, x_variable))
                                log.out(' - wild_match()  : {} -> {}'.format(ycompute, y_variable))
                                
                                # Update the xdata, ydata, xcompute, and ycompute keys
                                compute_chars = ['${', '}']
                                if all(c in x_variable for c in compute_chars):
                                    inner_mode['xcompute'] = x_variable
                                else:
                                    inner_mode['xcompute'] = ''
                                    inner_mode['xdata'] = x_variable
                                    
                                if all(c in y_variable for c in compute_chars):
                                    inner_mode['ycompute'] = y_variable
                                else:
                                    inner_mode['ycompute'] = ''
                                    inner_mode['ydata'] = y_variable
                                    
                                # Update the analysis columns
                                analysis_dict = {**x_analysis_dict, **y_analysis_dict}
                                for n, (method, xlo, xhi, misc, name) in enumerate(inner_mode['analysis']):
                                    settings = get_misc_setting(misc, use_eval=True)
                                    for key in analysis_dict:
                                        if key in settings:
                                            log.out(' - partial_name("{}", {}): {} -> {}'.format(method, key, settings[key], analysis_dict[key]))
                                            settings[key] = analysis_dict[key]
                                    new_misc = '; '.join( ['{}={}'.format(key, value) for key, value in settings.items()] )
                                    inner_mode['analysis'][n][3] = new_misc 
                        else:
                            log.out(' - wild_match()   x != y: {} != {}'.format(xcompute, ycompute))
                            log.out(' - Skipping and hoping for the best (using supplied "xdata" and "ydata" columns)')
                    except: pass
                
                try: # we dont want crashes to exit this loop
                    analyzed = analysis(inner_mode, plot=True, savefig=savefig, dpi=dpi, log=log, log_clear=False)
                    
                    basename = analyzed.get_basename_and_builder_dirs()
                    basename =  io_functions.path_to_string(basename)
                    dirnames.append(os.path.dirname(basename))
                    processed_files.append(file)
                    
                    wildcards[file] = wilcards_lst
                    outputs[file] = analyzed.outputs
                    abouts[file] = analyzed.about
                    for fig_name, fig_instance in analyzed.figures.items():
                        if fig_name not in figures:
                            figures[fig_name] = [fig_instance]
                        else: figures[fig_name].append(fig_instance)
                    
                    log.out(' - analysis state: SUCCESS')
                    
                    # Clear all logs and data structs to keep
                    # RAM usage down and consistent logging 
                    #analyzed.figures.clear()
                    log.clear_all()
                    #del analyzed

                except: 
                    failed_files.append(file)
                    stack_trace_string = traceback.format_exc()
                    log.out(' - analysis state: FAILED')
                    log.out(' - traceback     :')
                    log.out(' {}'.format(stack_trace_string))
        
        # Finalize the array run
        array_file = mode.get('array_file', '')
        if outputs and dirnames and array_file:
            dirname = dirnames[-1]
            array_file = mode['array_file']
            output_basename = os.path.join(dirname, array_file)
            array_csv(outputs, abouts, wildcards, output_basename, log)
            for fig_name, figs in figures.items():
                try:
                    pdf_title = '{}_{}.pdf'.format(output_basename, fig_name)
                    pdf = PdfPages(pdf_title)
                    for n, fig in enumerate(figs):
                        try: filename = os.path.basename(processed_files[n])
                        except: filename = ''
                        try:
                            existing_title = fig.get_suptitle()
                            if existing_title: # If plot already has sup, newline the filename
                                filename = f'{existing_title}\n{filename}'
                            fig.suptitle(filename)
                            try: fig.tight_layout()
                            except: pass
                            pdf.savefig(fig)
                        except: 
                            stack_trace_string = traceback.format_exc()
                            log.out(f'\nAn error occurred: {stack_trace_string}')
                        finally:
                            # free memory now
                            if n + 1 < len(figs): plt.close(fig)
                    log.out(' - pdf_plots={}'.format(pdf_title))
                    pdf.close()
                except:
                    stack_trace_string = traceback.format_exc()
                    log.out(' - FAILED to write pdf_plots={}'.format(pdf_title))
                    log.out(f'\nAn error occurred: {stack_trace_string}')
                    
        if failed_files:
            log.out('\n\nFiles that failed to be analyzed for unknown reasons:')
            for file in failed_files:
                log.out(' - {}'.format(file))
        
        if ding: print('\a') # Alert
        #plt.show()
        #if '*' not in mode['logfile']: plt.show()
    else:
        analyzed = analysis(mode, plot=True, savefig=savefig, dpi=dpi, log=log, log_clear=log_clear)

    return analyzed


################################
# Generate main class analyzer #
################################
class analysis:
    def __init__(self, mode, plot=True, savefig=True, dpi=300, log=None, log_clear=True):
        #################
        # Set up logger #
        #################
        if log is None:
            log = io_functions.LUNAR_logger()
        log.configure(level='production')
        #log.configure(level='debug', print2console=False)
        if log_clear: log.clear_all()
        
        #################################
        # Setup some default attributes #
        #################################
        self.parent_directory = mode['parent_directory']
        self.logfile = mode['logfile']
        self.log = log
        self.columns = [] # list of columns from log file
        self.outputs = {} # {analysis-name:dict-of-outputs-based-on-method}
        self.about = {} # {analysis-name:dict-of-outputs-about-information}
        
        # Setup dict to store figures *NOTE: the 'figure_description', should be short as it
        # will be used for setting up a "{}_{}.pdf".format(array_file, figure_description)
        self.figures = {} # {'figure-description':fig_instance}
        
        
        #############################
        # Read log file to get data #
        #############################
        start_time = time.time()
        self.log.out('\n\n--------------------------------------------------------------------------------------')
        self.log.out('                              Inputs and Data2load sections')
        self.log.out('--------------------------------------------------------------------------------------')
        # Read logfile and get data
        try:
            logfile = self.logfile
            self.log.out('  logfile={}'.format(logfile))
            if os.path.isfile(logfile):
                keywords = mode['keywords']
                sections = mode['sections']
                
                # Sections may have additional settings behind the 1st ";" char so remove if needed
                if ';' in str(sections):
                    sections_only, additional_settings = str(sections).split(';', 1)
                    sections_only       = sections_only.strip()
                    additional_settings = additional_settings.strip()
                    sections = sections_only # update sections
                
                self.log.out('  keywords={}'.format(mode['keywords']))
                self.log.out('  sections={}'.format(sections))
                logged = read_log.file(logfile, keywords=keywords, pflag=self.log.print2console)
                data = logged.get_data(sections, remove_duplicates=True, pflag=self.log.print2console) # {column-name:[lst of data]}
                self.columns = list(data.keys())
            else: data = {}; self.log.GUI_error(f'ERROR lammps logfile {logfile} does not exist');
        except: 
            data = {}
            stack_trace_string = traceback.format_exc()
            self.log.out('Traceback:')
            self.log.out(stack_trace_string)
            self.log.GUI_error('ERROR failed to read logfile and load data')
        
        
        ##################################
        # If data exists, start analysis #
        ##################################
        if not data: 
            self.log.GUI_error('ERROR no data loaded, unable to analyze and plot none existant data')
        else:
            # Try to closing any currently open plots, before starting
            if plot:
                try: plt.close('all')
                except: pass
            
            # Get X- and Y-data
            xdata = mode['xdata']
            ydata = mode['ydata']
            
            self.log.out('  xdata="{}"'.format(xdata))
            self.log.out('  ydata="{}"'.format(ydata))
            #print( list(data.keys()) )
            
            # Check to see if xdata and ydata exists in data
            if xdata not in data and not mode['xcompute']:
                self.log.GUI_error(f'ERROR X-data column "{xdata}" does not exist in logfile. Please check that keywords, sections, and X-data is correct.')
            elif ydata not in data and not mode['ycompute']:
                self.log.GUI_error(f'ERROR Y-data column "{ydata}" does not exist in logfile. Please check that keywords, sections, and Y-data is correct.')
            else:
                # Perform compute on data
                xcompute = mode['xcompute']
                ycompute = mode['ycompute']
                indexes2remove = []
                if xcompute:
                    x, message = vectorized_string_compute.thermo_data(data, xcompute)
                    self.log.out('  xcompute={}'.format(xcompute))
                    self.log.out('  {}'.format(message))
                    xdata = 'xcompute'
                else: x = data[xdata] 
                if ycompute:
                    y, message = vectorized_string_compute.thermo_data(data, ycompute)
                    self.log.out('  ycompute={}'.format(ycompute))
                    self.log.out('  {}'.format(message))
                    ydata = 'ycompute'
                else: y = data[ydata] 
                
                # Remove paired x/y values if either value is NaN or infinite
                if xcompute or ycompute:
                    x_values = np.asarray(x, dtype=float)
                    y_values = np.asarray(y, dtype=float)
                    valid = np.isfinite(x_values) & np.isfinite(y_values)
                    indexes2remove = np.flatnonzero(~valid).tolist()
                    if indexes2remove:
                        self.log.out('  Removing non-finite x/y values at indexes: {}'.format(indexes2remove))
                        x = x_values[valid].tolist()
                        y = y_values[valid].tolist()
                
                x = list(x)
                y = list(y)
                self.log.out('  len(xdata)={}'.format(len(x)))
                self.log.out('  len(ydata)={}'.format(len(y)))
                
                # Check for zero lengths of data
                if len(x) == 0:
                    self.log.GUI_error('ERROR X-data has zero points. Please check that X-compute, keywords, sections, and X-data is correct.')
                if len(y) == 0:
                    self.log.GUI_error('ERROR Y-data has zero points. Please check that Y-compute, keywords, sections, and Y-data is correct.')
                
                # Grant access to Raw-data outside of this class
                about = {'xdata': 'Raw X-data that was the basis of the analysis.',
                         'ydata': 'Raw Y-data that was the basis of the analysis.'}
                output = {'xdata': x,
                          'ydata': y}
                self.outputs['Raw-data'] = output
                self.about['Raw-data'] = about
                
                # Log Labels
                self.log.out('  xlabel={}'.format(mode['xlabel']))
                self.log.out('  xlabel={}'.format(mode['ylabel']))
            
                # Apply nevery to x and y data
                try: nevery = int(mode['nevery'])
                except: nevery = 1
                self.log.out('  nevery={}'.format(nevery))
                if nevery > 1:
                    x = x[::nevery]
                    y = y[::nevery]
                    
                # Save xtmp and ytmp in-case they are needed due to LAMMPS data cleaning
                xlmp = x; ylmp = y;
                
                # Generate cite string to call later on for work to cite
                cite_string = '-'.join(10*['cite'])
                
                # Get any anaylsis that users may want
                analysis = []; write_plotted_data = False; rm_lmp_data = False; 
                LAMMPS_cleaning = {} # {'Method':[xlo, xhi, misc, name]}
                for n, (method, xlo, xhi, misc, name)  in enumerate(mode['analysis']):
                    try: 
                        # Get basic settings like xlo, xhi, and skips
                        if method == 'skip': continue
                        try: xlo = float(xlo)
                        except: xlo = min(x); self.log.GUI_error('  xlo was not provided for "{}", using minimum ({})'.format(method, xlo))
                        try: xhi = float(xhi)
                        except: xhi = max(x); self.log.GUI_error('  xhi was not provided for "{}", using maximum ({})'.format(method, xhi))
                        
                        # Set name if not provided ('LAMMPS data' methods will have "in-built" names if not provided)
                        if name == '' and 'LAMMPS data' not in method:
                            name = 'analysis-{}'.format(n)
                            self.log.warn(f'WARNING name was left empty, imposing {name}')
                            
                        # Get LAMMPS cleaning inputs
                        if 'LAMMPS data' in method:
                            if method == 'LAMMPS data (remove from plot)': rm_lmp_data = True
                            else: LAMMPS_cleaning[method] = [xlo, xhi, misc, name]
                            continue
                        
                        # Set flag for writing plotted data to .csv
                        if method == 'write plotted data to csv file':
                            write_plotted_data = True; continue
                        analysis.append([method, xlo, xhi, misc, name])
                    except: pass
                    
                # Save LAMMPS data to plot
                data2plot = [self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label='LAMMPS data', shiftable=True)]
                
                
                #####################################################
                # Apply different "cleaning" methods to LAMMPS data #
                #####################################################
                # Apply 'LAMMPS data (X-sort)' if present (this MUST Be done first before anything else occurs)
                if 'LAMMPS data (X-sort)' in LAMMPS_cleaning:
                    xlo, xhi, misc, name = LAMMPS_cleaning['LAMMPS data (X-sort)']
                    rx, ry = misc_funcs.reduce_data(x, y, xlo, xhi)
                    x, y, index_backmap = misc_funcs.increasing_independant_sort(x, y)
                    #x, y, removed_indexes = misc_funcs.increasing_independant_cleaning(x, y)
                    # for n, i in enumerate(index_backmap):
                    #     self.log.out('{}, {}'.format(n, i))
                
                # Apply "normal" LAMMPS data cleaning
                LAMMPS_data_xlos = set(); LAMMPS_data_xhis = set(); cleaning_settings = {} # {method:[lst of parameters]}
                for method in LAMMPS_cleaning:
                    xlo, xhi, misc, name = LAMMPS_cleaning[method]
                    LAMMPS_data_xlos.add(xlo); LAMMPS_data_xhis.add(xhi)
                    
                    #-------------------------------------------#
                    # Apply linear interpolation to LAMMPS data #
                    #-------------------------------------------#
                    if method == 'LAMMPS data (apply linear-interp)':
                        if name == '': name = 'LAMMPS data w/linear-interp'
                        setting = self.get_misc_setting(misc)
                        
                        if 'dx' in misc:
                            dx = setting['dx']
                        else: 
                            dx = None
                            misc = ' default-dx=None'
                            
                        if 'n_points' in misc:
                            n_points = setting['n_points']
                        else: 
                            n_points = None
                            misc = ' default-n_points=None'
                        
                        # Perform linear interp "cleaning"
                        self.log.out('  Implementing: LAMMPS data (apply linear-interp) with {} settings'.format(misc))
                        x, y = self.interp_constant_dx(x, y, xlo, xhi, dx, n_points)
                        data2plot.append(self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        self.log.out('  len(xdata)={}'.format(len(x)))
                        self.log.out('  len(ydata)={}'.format(len(y)))
                        
                        # Log info from method for use in other analysis such as RFR stres-strain analysis for transverse strain "cleaning"
                        cleaning_settings[method] = [dx, n_points]
                    
                    #-------------------------------------#
                    # Apply moving average to LAMMPS data #
                    #-------------------------------------#
                    if method == 'LAMMPS data (apply bin average)':
                        if name == '': name = 'LAMMPS data w/bin average'
                        setting = self.get_misc_setting(misc)
                        
                        if 'dx' in misc:
                            dx = setting['dx']
                        else: 
                            dx = 1
                            misc = ' default-dx=1'

                        if 'start' in misc:
                            start = setting['start']
                        else: 
                            start = None
                            misc = ' default-start=None'
                            
                        if 'stop' in misc:
                            stop = setting['stop']
                        else: 
                            stop = None
                            misc = ' default-stop=None'
                        
                        
                        # Perform moving average "cleaning"
                        self.log.out('  Implementing: LAMMPS data (apply bin average) with {} settings'.format(misc))
                        x, y = self.bin_average(x, y, xlo, xhi, dx, start, stop)
                        data2plot.append(self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        self.log.out('  len(xdata)={}'.format(len(x)))
                        self.log.out('  len(ydata)={}'.format(len(y)))
                        
                        # Log info from method for use in other analysis such as RFR stres-strain analysis for transverse strain "cleaning"
                        cleaning_settings[method] = [dx, start, stop]
                    
                    #-------------------------------------#
                    # Apply moving average to LAMMPS data #
                    #-------------------------------------#
                    if method == 'LAMMPS data (apply moving average)':
                        if name == '': name = 'LAMMPS data w/moving average'
                        setting = self.get_misc_setting(misc)
                        if 'window' in misc:
                            LAMMPS_window = setting['window']
                        else: 
                            LAMMPS_window = 100;
                            misc = ' default-window=100';
                        
                        # Perform moving average "cleaning"
                        self.log.out('  Implementing: LAMMPS data (apply moving average) with {} settings'.format(misc))
                        x, y = self.moving_average(x, y, xlo, xhi, LAMMPS_window)
                        data2plot.append(self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        
                        # Log info from method for use in other analysis such as RFR stres-strain analysis for transverse strain "cleaning"
                        cleaning_settings[method] = [LAMMPS_window]
                    
                        # Grant access to moving-average-data outside of this class
                        about = {'xdata': 'List of moving average X-data w/window={}.'.format(LAMMPS_window),
                                 'ydata': 'List of moving average Y-data w/window={}.'.format(LAMMPS_window)}
                        output = {'xdata': x,
                                  'ydata': y}
                        self.outputs[name] = output
                        self.about[name] = about
                    
                    #-----------------------------------------#
                    # Apply Butterworth filter to LAMMPS data #
                    #-----------------------------------------#
                    if method == 'LAMMPS data (apply Butterworth filter)':
                        if name == '': name = 'LAMMPS data w/butterworth filter'
                        setting = self.get_misc_setting(misc)
                        if 'order' in setting:
                            LAMMPS_order = setting['order']
                        else: LAMMPS_order=2; misc += ' default-order=2;';
                        
                        if 'wn' in setting:
                            LAMMPS_wn = setting['wn']
                        else: LAMMPS_wn=0.01; misc += ' default-wn=0.01;';
                        
                        if 'psd' in setting:
                            LAMMPS_psd = setting['psd']
                        else: LAMMPS_psd=False; misc = 'default-psd=False';
                        
                        if 'qm' in setting:
                            LAMMPS_qm = setting['qm']
                        else: LAMMPS_qm='1,1'; misc = 'default-qm=1,1-p';
                        
                        if 'csv' in setting:
                            LAMMPS_write_data = setting['csv']
                        else: LAMMPS_write_data = False; misc += ' default-csv=False';
                        
                        if 'savefig' in setting:
                            LAMMPS_savefig = str(setting['savefig'])
                        else:
                            LAMMPS_savefig = 'all'; misc += ' default-savefig=all';
                        
                        # Perform Butterworth filtering
                        self.log.out('  Implementing: LAMMPS data (apply butterworth filter) with {} settings'.format(misc))
                        figname = '{}_X={}_Y={}'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
                        x, y, wn, qm = self.butterworth_lowpass(x, y, xlo, xhi, LAMMPS_wn, LAMMPS_order, LAMMPS_qm, LAMMPS_psd, LAMMPS_write_data, LAMMPS_savefig, figname, dpi)
                        data2plot.append(self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        
                        # Log info from method for use in other analysis such as RFR stres-strain analysis for transverse strain "cleaning"
                        cleaning_settings[method] = [LAMMPS_wn, LAMMPS_order, LAMMPS_qm, LAMMPS_psd, LAMMPS_write_data, LAMMPS_savefig]
                        
                        # Grant access to Butterworth-data outside of this class
                        about = {'xdata': 'List of Butterworth X-data w/order={}; wn={}; qm={}.'.format(LAMMPS_order, LAMMPS_wn, LAMMPS_qm),
                                 'ydata': 'List of Butterworth Y-data w/order={}; wn={}; qm={}.'.format(LAMMPS_order, LAMMPS_wn, LAMMPS_qm),
                                 'wn': 'Value for the wn (normalized cutoff frequency used for filtering.',
                                 'qm': 'Value for the qm (quadrant mirroring used for filtering.'}
                        output = {'xdata': x,
                                  'ydata': y,
                                  'wn': wn,
                                  'qm': qm}
                        self.outputs[name] = output
                        self.about[name] = about
                    
                    #------------------------------------------------#
                    # Apply Whittaker-Eilers smoother to LAMMPS data #
                    #------------------------------------------------#
                    if method == 'LAMMPS data (apply Whittaker-Eilers)':
                        if name == '': name = 'LAMMPS data w/Whittaker-Eilers smoothing'
                        setting = self.get_misc_setting(misc)
                        if 'order' in setting:
                            LAMMPS_order = setting['order']
                        else: LAMMPS_order=2; misc += ' default-order=2;';
                        
                        if 'lambda' in setting:
                            LAMMPS_lmbda = setting['lambda']
                        else: LAMMPS_lmbda='op'; misc += ' default-lambda=op;';
                        
                        # Perform Whittaker-Eilers smoothing
                        self.log.out('  Implementing: LAMMPS data (Whittaker-Eilers smoothing) with {} settings'.format(misc))
                        figname = '{}_X={}_Y={}_Whittaker-Eilers'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
                        x, y = self.whittaker_eilers_smoothing(x, y, xlo, xhi, LAMMPS_order, LAMMPS_lmbda, mode['xlabel'], mode['ylabel'], figname)
                        data2plot.append(self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        
                        # Log info from method for use in other analysis such as RFR stres-strain analysis for transverse strain "cleaning"
                        cleaning_settings[method] = [LAMMPS_order, LAMMPS_lmbda]
                        
                        # Grant access to Whittaker-Eilers-data outside of this class
                        about = {'xdata': 'List of Whittaker-Eilers X-data w/order={}; lambda={}.'.format(LAMMPS_order, LAMMPS_lmbda),
                                 'ydata': 'List of Whittaker-Eilers Y-data w/order={}; lambda={}.'.format(LAMMPS_order, LAMMPS_lmbda)}
                        output = {'xdata': x,
                                  'ydata': y}
                        self.outputs[name] = output
                        self.about[name] = about
                        
                    #-------------------------------#
                    # Fit polynomial to LAMMPS data #
                    #-------------------------------#
                    if method == 'LAMMPS data (fit polynomial)':
                        if name == '': name = 'LAMMPS data w/fit polynomial'
                        setting = self.get_misc_setting(misc)
                        if 'deg' in setting:
                            LAMMPS_deg = setting['deg']
                        else: LAMMPS_deg=2; misc += ' default-deg=2;';
                        
                        # Fit polynomial
                        self.log.out('  Implementing: LAMMPS data (fit polynomial) with {} settings'.format(misc))
                        x, y = self.polynomial_fit(x, y, xlo, xhi, LAMMPS_deg)
                        data2plot.append(self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        
                        # Log info from method for use in other analysis such as RFR stres-strain analysis for transverse strain "cleaning"
                        cleaning_settings[method] = [LAMMPS_deg]
                        
                        # Grant access to Fit polynomial-data outside of this class
                        about = {'xdata': 'List of fit polynomial X-data w/deg={}.'.format(LAMMPS_deg),
                                 'ydata': 'List of fit polynomial Y-data w/deg={}.'.format(LAMMPS_deg)}
                        output = {'xdata': x,
                                  'ydata': y}
                        self.outputs[name] = output
                        self.about[name] = about
                        
                    #-----------------------------#
                    # Apply LOWESS to LAMMPS data #
                    #-----------------------------#
                    if method == 'LAMMPS data (LOWESS)':
                        if name == '': name = 'LAMMPS data w/LOWESS smoothing'
                        setting = self.get_misc_setting(misc)
                        if 'fraction' in setting:
                            LAMMPS_fraction = setting['fraction']
                        else: LAMMPS_fraction=0.2; misc += ' default-fraction=0.2;';
                        if 'max_iter' in setting:
                            LAMMPS_max_iter = setting['max_iter']
                        else: LAMMPS_max_iter=10; misc += ' default-max_iter=10;';
                        
                        # Smooth with LOWESS
                        self.log.out('  Implementing: LAMMPS data (LOWESS) with {} settings'.format(misc))
                        x, y = self.lowess(x, y, xlo, xhi, LAMMPS_fraction, LAMMPS_max_iter)
                        data2plot.append(self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        
                        # Log info from method for use in other analysis such as RFR stres-strain analysis for transverse strain "cleaning"
                        cleaning_settings[method] = [LAMMPS_fraction, LAMMPS_max_iter]
                        
                        # Grant access to Whittaker-Eilers-data outside of this class
                        about = {'xdata': 'List of LOWESS X-data w/fraction={}; max_iter={}.'.format(LAMMPS_fraction, LAMMPS_max_iter),
                                 'ydata': 'List of LOWESS Y-data w/fraction={}; max_iter={}.'.format(LAMMPS_fraction, LAMMPS_max_iter)}
                        output = {'xdata': x,
                                  'ydata': y}
                        self.outputs[name] = output
                        self.about[name] = about
    
    
                    #----------------------------------#
                    # Apply iFFT filter to LAMMPS data #
                    #----------------------------------#
                    if method == 'LAMMPS data (apply iFFT filter)':
                        if name == '': name = 'LAMMPS data (apply iFFT filter)'
                        setting = self.get_misc_setting(misc)
                        if 'threshold' in setting:
                            threshold = setting['threshold']
                        else: threshold='mean'; misc += ' default-threshold=mean;';
                        
                        if 'qm' in setting:
                            qm = setting['qm']
                        else: qm='1,1'; misc = 'default-qm=1,1-p';
                        
                        if 'savefig' in setting:
                            LAMMPS_savefig = str(setting['savefig'])
                        else:
                            LAMMPS_savefig = 'all'; misc += ' default-savefig=all';
                        
                        # Filter with inverse FFT
                        figname = '{}_X={}_Y={}_iFFT'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
                        self.log.out('  Implementing: LAMMPS data (apply iFFT filter) with {} settings'.format(misc))
                        y, qm_iFFT = self.iFFT_filter(x, y, xlo, xhi, threshold, qm, LAMMPS_savefig, figname, dpi)
                        data2plot.append(self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        
                        # Log info from method for use in other analysis such as RFR stres-strain analysis for transverse strain "cleaning"
                        cleaning_settings[method] = [threshold, qm, LAMMPS_savefig]
                        
                        # Grant access to iFFT-data outside of this class
                        about = {'xdata': 'List of iFFT X-data w/threshold={}; qm={}.'.format(threshold, qm_iFFT),
                                 'ydata': 'List of iFFT Y-data w/threshold={}; qm={}.'.format(threshold, qm_iFFT),
                                 'qm': 'Value for the qm (quadrant mirroring used for filtering.',
                                 'threshold': 'Threshold used for PSD cutoff'}
                        output = {'xdata': x,
                                  'ydata': y,
                                  'qm': qm_iFFT,
                                  'threshold': threshold}
                        self.outputs[name] = output
                        self.about[name] = about
                        
                    #----------------------------------------------#
                    # Apply Savitzky-Golay smoother to LAMMPS data #
                    #----------------------------------------------#
                    if method == 'LAMMPS data (Savitzky-Golay)':
                        if name == '': name = 'LAMMPS data (Savitzky-Golay)'
                        setting = self.get_misc_setting(misc)
                        if 'window' in setting:
                            window = setting['window']
                        else: threshold=100; misc += ' default-window=100;';
                        
                        if 'order' in setting:
                            order = setting['order']
                        else: order=1; misc = 'default-order=1';
                        
                        self.log.out('  Implementing: LAMMPS data (Savitzky-Golay) with {} settings'.format(misc))
                        x, y = self.SavitzkyGolay(x, y, xlo, xhi, window, order)
                        data2plot.append(self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        
                        # Log info from method for use in other analysis such as RFR stres-strain analysis for transverse strain "cleaning"
                        cleaning_settings[method] = [window, order]
                        
                        # Grant access to Whittaker-Eilers-data outside of this class
                        about = {'xdata': 'List of Savitzky-Golay X-data w/window={}; max_iter={}.'.format(window, order),
                                 'ydata': 'List of Savitzky-Golay Y-data w/window={}; max_iter={}.'.format(window, order)}
                        output = {'xdata': x,
                                  'ydata': y}
                        self.outputs[name] = output
                        self.about[name] = about

    
    
                # Set "Global" LAMMPS_data_xlo (maximize) and LAMMPS_data_xhi (minimize) values
                try: LAMMPS_data_xlo = max(LAMMPS_data_xlos)
                except: LAMMPS_data_xlo = min(x)
                try: LAMMPS_data_xhi = min(LAMMPS_data_xhis)
                except: LAMMPS_data_xhi = max(x)
                if len(LAMMPS_data_xlos) > 1: 
                    self.log.warn(f'WARNING different "LAMMPS cleaning methods" applied different xlo contraints. Using maximum xlo value ({LAMMPS_data_xlo}) for other analysis.')
                if len(LAMMPS_data_xlos) > 1: 
                    self.log.warn(f'WARNING different "LAMMPS cleaning methods" applied different xhi contraints. Using minimum xhi value ({LAMMPS_data_xhi}) for other analysis.')
                
                
                ############################################
                # Start performing any analysis if present #
                ############################################
                shift = False; shift_amount = 0; # Defaults
                for method, xlo, xhi, misc, name in analysis:                
                    #-----------------------------------------#
                    # Apply moving average to raw LAMMPS data #
                    #-----------------------------------------#
                    if method == 'moving average':
                        setting = self.get_misc_setting(misc)
                        if 'window' in setting:
                            window = setting['window']
                        else: 
                            window = 100;
                            misc += ' default-window=100';
                        xout, yout = self.moving_average(x, y, xlo, xhi, window)
                        label = '{} ({})'.format(name, misc)
                        self.log.out(self.format_analysis(method, xlo, xhi, misc, name))
                        self.log.out('  {}'.format(label))
                        movavgdata = self.plot_parms(x=xout, y=yout, style='point', marker='.', line='-', size=4, label=label, shiftable=True)
                        data2plot.append(movavgdata)
                        
                        # Grant access to moving-average-data outside of this class
                        about = {'xdata': 'List of moving average X-data w/window={}.'.format(window),
                                 'ydata': 'List of moving average Y-data w/window={}.'.format(window)}
                        output = {'xdata': xout,
                                  'ydata': yout}
                        self.outputs[name] = output
                        self.about[name] = about
                    
                    
                    #---------------------------------------------------------------------------------------------------------#
                    # Apply butterworth low pass filter for "seeing trends" in data (nothing will be done on raw LAMMPS data) #
                    #---------------------------------------------------------------------------------------------------------#
                    if method == 'Butterworth (low pass)':
                        setting = self.get_misc_setting(misc)
                        if 'order' in setting:
                            order = setting['order']
                        else: order=2; misc += ' default-order=2;';
                        
                        if 'wn' in setting:
                            wn = setting['wn']
                        else: wn=0.01; misc += ' default-wn=0.01;';
                        
                        if 'psd' in setting:
                            psd = setting['psd']
                        else: psd=False; misc = 'default-psd=False';
                        
                        if 'qm' in setting:
                            qm = setting['qm']
                        else: qm='1,1'; misc = 'default-qm=1,1-p';
                        
                        if 'csv' in setting:
                            write_data = setting['csv']
                        else: write_data = False; misc += ' default-csv=False';
                        
                        if 'savefig' in setting:
                            BW_savefig = str(setting['savefig'])
                        else:
                            BW_savefig = 'all'; misc += ' default-savefig=all';
    
                        label = '{} order = {}; wn = {}'.format(name, order, wn)
                        self.log.out('  {}'.format(label))
                        
                        figname = '{}_X={}_Y={}'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
                        xfilter, yfilter, wn, qm = self.butterworth_lowpass(x, y, xlo, xhi, wn, order, qm, psd, write_data, BW_savefig, figname, dpi)
                        self.log.out(self.format_analysis(method, xlo, xhi, misc, name))
                        bwfdata = self.plot_parms(x=xfilter, y=yfilter, style='point', marker='.', line='-', size=4, label=label, shiftable=True)
                        data2plot.append(bwfdata)
                        
                        # Grant access to Butterworth-data outside of this class
                        about = {'xdata': 'List of Butterworth X-data w/order={}; wn={}; qm={}.'.format(order, wn, qm),
                                 'ydata': 'List of Butterworth Y-data w/order={}; wn={}; qm={}.'.format(order, wn, qm),
                                 'wn': 'Value for the wn (normalized cutoff frequency used for filtering.',
                                 'qm': 'Value for the qm (quadrant mirroring used for filtering.'}
                        output = {'xdata': xfilter,
                                  'ydata': yfilter,
                                  'wn': wn,
                                  'qm': qm}
                        self.outputs[name] = output
                        self.about[name] = about
                        
                    #-------------------------------------------------------------#
                    # Apply iFFT filter (nothing will be done on raw LAMMPS data) #
                    #-------------------------------------------------------------#
                    if method == 'iFFT filter':
                        setting = self.get_misc_setting(misc)
                        if 'threshold' in setting:
                            threshold = setting['threshold']
                        else: threshold='mean'; misc += ' default-threshold=mean;';
                        
                        if 'qm' in setting:
                            qm = setting['qm']
                        else: qm='1,1'; misc = 'default-qm=1,1-p';
                        
                        if 'savefig' in setting:
                            savefig = str(setting['savefig'])
                        else:
                            savefig = 'all'; misc += ' default-savefig=all';
                        
                        # Filter with inverse FFT
                        figname = '{}_X={}_Y={}_iFFT'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
                        self.log.out('  Implementing: LAMMPS data (apply iFFT filter) with {} settings'.format(misc))
                        y_iFFT, qm_iFFT = self.iFFT_filter(x, y, xlo, xhi, threshold, qm, savefig, figname, dpi)
                        data2plot.append(self.plot_parms(x=x, y=y_iFFT, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        
                        # Grant access to iFFT-data outside of this class
                        about = {'xdata': 'List of iFFT X-data w/threshold={}; qm={}.'.format(threshold, qm_iFFT),
                                 'ydata': 'List of iFFT Y-data w/threshold={}; qm={}.'.format(threshold, qm_iFFT),
                                 'qm': 'Value for the qm (quadrant mirroring used for filtering.',
                                 'threshold': 'Threshold used for PSD cutoff'}
                        output = {'xdata': x,
                                  'ydata': y,
                                  'qm': qm_iFFT,
                                  'threshold': threshold}
                        self.outputs[name] = output
                        self.about[name] = about
                        
                    #---------------------------------#
                    # Apply Whittaker-Eilers smoother #
                    #---------------------------------#
                    if method == 'Whittaker-Eilers':
                        setting = self.get_misc_setting(misc)
                        if 'order' in setting:
                            order = setting['order']
                        else: order=2; misc += ' default-order=2;';
                        
                        if 'lambda' in setting:
                            lmbda = setting['lambda']
                        else: lmbda='op'; misc += ' default-lambda=op;';
                        
                        # Perform Whittaker-Eilers smoothing
                        figname = '{}_X={}_Y={}_Whittaker-Eilers-Non-LAMMPS'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
                        x, y = self.whittaker_eilers_smoothing(x, y, xlo, xhi, order, lmbda, mode['xlabel'], mode['ylabel'], figname)
                        data2plot.append(self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        
                        # Grant access to Whittaker-Eilers-data outside of this class
                        about = {'xdata': 'List of Whittaker-Eilers X-data w/order={}; lambda={}.'.format(order, lmbda),
                                 'ydata': 'List of Whittaker-Eilers Y-data w/order={}; lambda={}.'.format(order, lmbda)}
                        output = {'xdata': x,
                                  'ydata': y}
                        self.outputs[name] = output
                        self.about[name] = about
                        
                    #-----------------------#
                    # Apply LOWESS cleaning #
                    #-----------------------#
                    if method == 'LOWESS':
                        setting = self.get_misc_setting(misc)
                        if 'fraction' in setting:
                            fraction = setting['fraction']
                        else: fraction=0.2; misc += ' default-fraction=0.2;';
                        if 'max_iter' in setting:
                            max_iter = setting['max_iter']
                        else: max_iter=10; misc += ' default-max_iter=10;';
                        
                        # Smooth with LOWESS
                        self.log.out('  Implementing: LOWESS with {} settings'.format(misc))
                        x, y = self.lowess(x, y, xlo, xhi, fraction, max_iter)
                        data2plot.append(self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        
                        # Grant access to Whittaker-Eilers-data outside of this class
                        about = {'xdata': 'List of LOWESS X-data w/fraction={}; max_iter={}.'.format(LAMMPS_fraction, LAMMPS_max_iter),
                                 'ydata': 'List of LOWESS Y-data w/fraction={}; max_iter={}.'.format(LAMMPS_fraction, LAMMPS_max_iter)}
                        output = {'xdata': x,
                                  'ydata': y}
                        self.outputs[name] = output
                        self.about[name] = about
                    
                    
                    #-------------------------------#
                    # Fit polynomial to LAMMPS data #
                    #-------------------------------#
                    if method == 'fit polynomial':
                        setting = self.get_misc_setting(misc)
                        if 'deg' in setting:
                            deg = setting['deg']
                        else: deg=2; misc += ' default-deg=2;';
                        
                        # Fit polynomial
                        self.log.out('  Implementing: fit polynomial with {} settings'.format(misc))
                        x, y = self.polynomial_fit(x, y, xlo, xhi, deg)
                        data2plot.append(self.plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=name, shiftable=True))
                        
                        # Grant access to Fit polynomial-data outside of this class
                        about = {'xdata': 'List of fit polynomial X-data w/deg={}.'.format(LAMMPS_deg),
                                 'ydata': 'List of fit polynomial Y-data w/deg={}.'.format(LAMMPS_deg)}
                        output = {'xdata': x,
                                  'ydata': y}
                        self.outputs[name] = output
                        self.about[name] = about
                    
                    #-------------------------------------------------------------------#
                    # Apply a hyperbola fit to the data and return hyperbola parameters #
                    #-------------------------------------------------------------------#
                    if method == 'hyperbola':
                        minimum_convergence = None
                        initial_guess = False
                        slope_method = 'asymptote'
                        setting = self.get_misc_setting(misc)
                        if 'p' in setting:
                            minimum_convergence = setting['p']
                        if 'initial_guess' in setting:
                            initial_guess = setting['initial_guess']
                        if 'slope_method' in setting:
                            slope_method = setting['slope_method']
                        xout, yout, params, center, slopes, transition, tangent_intersection, tangents = self.fit_hyperbola(x, y, xlo, xhi, slope_method, minimum_convergence=minimum_convergence, initial_guess=initial_guess, maxiter=10**6)
                        self.log.out(self.format_analysis(method, xlo, xhi, misc, name))
                        
                        label = '{} slopes (lower-slope={};\nupper-slope={})'.format(name, slopes[0], slopes[1])
                        self.log.out('  {} - slopes (lower-slope={}; upper-slope={})'.format(name, slopes[0], slopes[1]))
                        hyperboladata = self.plot_parms(x=xout, y=yout, style='point', marker='.', line='-', size=2, label=label, shiftable=False)
                        data2plot.append(hyperboladata)
                        
                        label = '{} center (x={:.4f}; y={:.4f})'.format(name, center[0], center[1])
                        self.log.out('  {}'.format(label))
                        centerdata = self.plot_parms(x=center[0], y=center[1], style='point', marker='.', line='-', size=10, label=label, shiftable=True)
                        data2plot.append(centerdata)
                        
                        label = '{} tangents'.format(name)
                        tangentsdata = self.plot_parms(x=tangents[0], y=tangents[1], style='line', marker='.', line='--', size=2, label=label, shiftable=True)
                        data2plot.append(tangentsdata)
                        
                        label = '{} tangent intersection (x={:.4f}; y={:.4f})'.format(name, tangent_intersection[0], tangent_intersection[1])
                        self.log.out('  {}'.format(label))
                        tangentdata = self.plot_parms(x=tangent_intersection[0], y=tangent_intersection[1], style='point', marker='.', line='-', size=10, label=label, shiftable=True)
                        data2plot.append(tangentdata)
                        
                        if minimum_convergence is not None:
                            label = '{} transition region (P={}; xlo={:.4f}; xhi={:.4f})'.format(name, minimum_convergence, transition[0], transition[1])
                            self.log.out('  {}'.format(label))
                            miny = [min(yout) for _ in transition]
                            transitiondata = self.plot_parms(x=transition, y=miny, style='both', marker='|', line='-', size=10, label=label, shiftable=False)
                            data2plot.append(transitiondata)
                            
                        self.log.out(f'\n  {cite_string}')
                        self.log.out('  This method implements the work found at: https://doi.org/10.1016/j.polymer.2016.01.074')
                        self.log.out('  which should be cited if used for analyzing your LAMMPS data.')
                        self.log.out(f'  {cite_string}\n')
                            
                        # Grant access to Hyperbola-data outside of this class.
                        about = {'xdata': 'List of Hyperbola X-data w/p={}; initial_guess={}.'.format(minimum_convergence, initial_guess),
                                 'ydata': 'List of Hyperbola Y-data w/p={}; initial_guess={}.'.format(minimum_convergence, initial_guess),
                                 'params': 'List of Hyberbola parameters ordered [t0, v0, a, b, c]. Note: using volume data, not density data.',
                                 'center': 'List of center of hyberbola ordered [xc, yc].',
                                 'slopes': 'List of slopes at ends of hyberbola ordered [below-center-slope, above-center-slope]',
                                 'transition': 'List of range of transitioning data ordered [lo-end, hi-end]',
                                 'tan-inter': 'List of tangent intersections [xt, yt].'}                     
                        output = {'xdata': xout,
                                  'ydata': yout,
                                  'params': params,
                                  'center': center,
                                  'slopes': slopes,
                                  'transition': transition,
                                  'tan-inter': tangent_intersection}
                        self.outputs[name] = output
                        self.about[name] = about
                        
                    #------------------------------------------------------------------------#
                    # Apply a piecewise-regression to the LAMMPS data and return the outputs #
                    #------------------------------------------------------------------------#
                    if method == 'piecewise-regression':
                        setting = self.get_misc_setting(misc)
                        if 'n' in setting:
                            n = setting['n']
                        else: n = 1; misc += ' default-n=1';
                        if 'shift' in setting:
                            shift = setting['shift']
                        else: shift = False; misc += ' default-shift=False';
                            
                        xout, yout, xbreaks, ybreaks, slopes = self.piecewise_regression(x, y, xlo, xhi, n)
                        self.log.out(self.format_analysis(method, xlo, xhi, misc, name))
                        
                        shift_amount = 0
                        if shift:
                            shift_amount = ybreaks[0]
                            self.log.out(f'  Shifting all Y-data by {shift_amount}')
                        
                        label = '{} (n-breakpoints={})'.format(name, n)
                        self.log.out(f'  {label}')
                        piecewisedata = self.plot_parms(x=xout, y=yout, style='line', marker='.', line='-', size=3, label=label, shiftable=True)
                        data2plot.append(piecewisedata)
                        
                        label = '{} slopes and points: \n'.format(name)
                        for i, j in slopes:
                            slope = slopes[(i, j)]
                            label += '   slope between (p{} - p{}) = {}\n'.format(i, j, slope)
                            label += '    p{} = ({} - {})\n'.format(i, xbreaks[i], ybreaks[i]-shift_amount)
                            label += '    p{} = ({} - {})\n'.format(j, xbreaks[j], ybreaks[j]-shift_amount)
                        self.log.out(f'  {label}')
                        label = '{}'.format(name)
                        breakdata = self.plot_parms(x=xbreaks, y=ybreaks, style='point', marker='.', line='-', size=10, label=label, shiftable=True)
                        data2plot.append(breakdata)
                        
                        # Grant access to piecewise-regression-data outside of this class
                        about = {'xdata': 'List of piecewise-regression X-data w/n={}.'.format(n),
                                 'ydata': 'List of piecewise-regression Y-data w/n={}.'.format(n),
                                 'xbreaks': 'List of "break points" found from peicewise-regression in X-direction',
                                 'ybreaks': 'List of "break points" found from peicewise-regression in Y-direction',
                                 'slopes': 'Dictionary of slopes between break points format = {(i, i+1):slope}, where i and i+1 are indexes of break point in xbreaks and ybreaks'}
                        output = {'xdata': xout,
                                  'ydata': yout}
                                  # 'xbreaks': xbreaks,
                                  # 'ybreaks': ybreaks,
                                  # 'slopes': slopes}
                        for n, (xb, yb) in enumerate( zip(xbreaks, ybreaks) ):
                            xb_label = 'xbreaks[{}]'.format(n)
                            yb_label = 'ybreaks[{}]'.format(n)
                            output[xb_label] = xb
                            output[yb_label] = yb
                        for key in slopes:
                            s_label = 'slope[{}-{}]'.format(key[0], key[1])
                            slope = slopes[key]
                            output[s_label] = slope
                        self.outputs[name] = output
                        self.about[name] = about
                    
                    #-----------------------------------------------#
                    # Apply a spline integration to the LAMMPS data #
                    #-----------------------------------------------#
                    if method == 'spline-integration':
                        setting = self.get_misc_setting(misc)
                        if 'shift' in setting:
                            shift = setting['shift']
                        else: shift=False; misc += '  default-shift=False';
                            
                        self.log.out(self.format_analysis(method, xlo, xhi, misc, name))
                        xout, yout, area, shift_amount = self.spline_integration(x, y, xlo, xhi, shift)
                        label = '{} (area: {})'.format(name, area)
                        self.log.out(f'  {label}')
                        inetgrateddata = self.plot_parms(x=xout, y=yout, style='point', marker='.', line='-', size=6, label=label, shiftable=False)
                        data2plot.append(inetgrateddata)
                        
                        # Grant access to spline integration-data outside of this class
                        about = {'xdata': 'List of spline integration X-data w/shift={}.'.format(shift),
                                 'ydata': 'List of spline integration Y-data w/shift={}.'.format(shift),
                                 'area': 'Float of area under the curve in the given X-range ({}-{})'.format(xlo, xhi)}
                        output = {'xdata': xout,
                                  'ydata': yout,
                                  'area': area}
                        self.outputs[name] = output
                        self.about[name] = about
                    
                    #---------------------------------#
                    # Set a cursor to add to the plot #
                    #---------------------------------#
                    if method == 'cursor':
                        setting = self.get_misc_setting(misc)
                        if 'x' in setting or 'y' in setting:
                            if 'x' in setting:
                                xvalue = setting['x']
                            else: xvalue = None
                            
                            if 'y' in setting:
                                yvalue = setting['y']
                            else: yvalue = None
                            
                            # Set style and marker based on if x or y or x & y are set
                            plot_and_log = True
                            if xvalue is not None and yvalue is not None:
                                cursor_style = 'point'
                                cursor_marker = '+'
                            elif xvalue is not None:
                                cursor_style = 'vertical'
                                cursor_marker = ':'
                            elif yvalue is not None:
                                cursor_style = 'horizontal'
                                cursor_marker = ':'
                            else: plot_and_log = False
                            
                            if plot_and_log:
                                label = '{} ({})'.format(name, misc)
                                self.log.out(self.format_analysis(method, xlo, xhi, misc, name))
                                self.log.out('  {}'.format(label))
                                cursordata = self.plot_parms(x=xvalue, y=yvalue, style=cursor_style, marker=cursor_marker, line='-', size=4, label=label, shiftable=False)
                                data2plot.append(cursordata)
                                
                            # Grant access to cursor-data outside of this class
                            about = {'xvalue': 'Float of X-value of cursor.',
                                     'yvalue': 'Float of Y-value of cursor.'}
                            output = {'xvalue': xvalue,
                                      'yvalue': yvalue}
                            self.outputs[name] = output
                            self.about[name] = about
                
                    #-------------------------------------------#
                    # Find average of Y-data in a given X-range #
                    #-------------------------------------------#
                    if method == 'average':
                        # Perform analysis
                        reduced_x, reduced_y, average_y = self.average(x, y, xlo, xhi)
                        label = '{} (average = {:.6f})'.format(name, average_y)
                        self.log.out(self.format_analysis(method, xlo, xhi, misc, name))
                        self.log.out('  {}'.format(label))
                        avgdata = self.plot_parms(x=reduced_x, y=reduced_y, style='point', marker='.', line='-', size=4, label=label, shiftable=True)
                        data2plot.append(avgdata)
                        
                        # Grant access to average-data outside of this class
                        about = {'xdata': 'List of xdata used to average over.',
                                 'ydata': 'List of ydata used to average over.',
                                 'average': 'Float average Y-data in the given X-range {}={}.'.format(xlo, xhi)
                                 }
                        outputs = {'xdata': reduced_x,
                                   'ydata': reduced_y,
                                   'average': average_y}
                        self.outputs[name] = outputs
                        self.about[name] = about
                        
                    #-------------------------------------------------#
                    # Find linear regression within a certain X-range #
                    #-------------------------------------------------#
                    if method == 'linear regression':
                        setting = self.get_misc_setting(misc)
                        if 'shift' in setting:
                            shift = setting['shift']
                        else: shift=False; misc = 'default-shift=False';
    
                        if 'ci' in setting:
                            confidence_interval = setting['ci']
                        else: confidence_interval = False; misc += ' default-ci=0';
                            
                        b0, b1, xreg, yreg, r2, ci = self.linear_regression(x, y, xlo, xhi, confidence_interval)
                        self.log.out(self.format_analysis(method, xlo, xhi, misc, name))
                        if shift: 
                            shift_amount = b0
                            b0 = 0
                        label = '{} y = {:.10f}x + {:.10f} (shifted by {:.4f}; r2={:.4f})'.format(name, b1, b0, shift_amount, r2)
                        self.log.out('  {}'.format(label))
                        if 'extend' in setting:
                            xreg = [xlo, xhi]
                            extend = setting['extend']
                            if extend > 0:
                                xreg.insert(2, xreg[1]+extend)
                            else: xreg.insert(0, xreg[0]+extend)
                            yreg = [i * b1 + b0 for i in xreg]
                            regdata = self.plot_parms(x=xreg, y=yreg, style='both', marker='o', line='-', size=3, label=label, shiftable=True)
                        else:
                            regdata = self.plot_parms(x=xreg, y=yreg, style='line', marker='.', line='-', size=3, label=label, shiftable=True)
                        data2plot.append(regdata)
                        
                        # Grant access to linear-regression-data outside of this class
                        about = {'xdata': 'List of X-data defining the linear-regresion model, order [xlo, xhi].',
                                 'ydata': 'List of Y-data defining the linear-regresion model, order [ylo, yhi].',
                                 'b1': 'Float value of b1 parameter (slope) for general equation of line y = b1*x + b0.',
                                 'b0': 'Float value of b0 parameter (Y-intercept) for general equation of line y = b1*x + b0.'}
                        outputs = {'xdata': xreg,
                                   'ydata': yreg,
                                   'b1': b1,
                                   'b0': b0}
                        self.outputs[name] = outputs
                        self.about[name] = about
                    
                    #------------------------------------------------#
                    # Find minimum of maximum data within an X-range #
                    #------------------------------------------------#
                    if method in ['minimum', 'maximum']:
                        setting = self.get_misc_setting(misc)
                        xout, yout = misc_funcs.reduce_data(x, y, xlo, xhi)
                        if method == 'minimum': yvalue = min(yout)
                        if method == 'maximum': yvalue = max(yout)
                        xvalue = xout[yout.index(yvalue)]
                        label = '{} (x={}; y={})'.format(name, xvalue, yvalue)
                        minmaxdata = self.plot_parms(x=xvalue, y=yvalue, style='point', marker='o', line='-', size=8, label=label, shiftable=True)
                        data2plot.append(minmaxdata)
                        self.log.out('  {}'.format(label))
                        
                        # Grant access to minimum-or-maximum-data outside of this class
                        about = {'xvalue': 'Float of X-value {}.'.format(method),
                                 'yvalue': 'Float of Y-value {}.'.format(method)}
                        output = {'xvalue': xvalue,
                                  'yvalue': yvalue}
                        self.outputs[name] = output
                        self.about[name] = about
                        
                    #-------------------------------------------#
                    # Differentiate the data and produce a plot #
                    #-------------------------------------------#
                    if method == 'Calculus: Differentiate Data':
                        figname = '{}_X={}_Y={}_derivatives'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
                        setting = self.get_misc_setting(misc)
                        if 'order' in setting:
                            order = setting['order']
                        else: order='1,2'; misc += ' default-order=1,2;';
                        
                        if 'csv' in setting:
                            write_data = setting['csv']
                        else: write_data = False; misc += ' default-csv=False';
                        
                        if 'savefig' in setting:
                            c_savefig = str(setting['savefig'])
                        else:
                            c_savefig = 'all'; misc += ' default-savefig=all';
                        
                        self.log.out(self.format_analysis(method, xlo, xhi, misc, name))
                        critical_x, critical_y, inflection_x, inflection_y = self.differentiate_data(x, y, xlo, xhi, order, write_data, mode['xlabel'], mode['ylabel'], c_savefig, figname, dpi)
                        self.log.out('  +------------------------+')
                        self.log.out('  | Critical Points (x, y) |')
                        self.log.out('  +------------------------+')
                        for i, j in zip(critical_x, critical_y):
                            self.log.out('  {} {}'.format(i, j))
                        self.log.out('')
                        self.log.out('  +--------------------------+')
                        self.log.out('  | Inflection Points (x, y) |')
                        self.log.out('  +--------------------------+')
                        for i, j in zip(inflection_x, inflection_y):
                            self.log.out('  {} {}'.format(i, j))
                        self.log.out('\n')
                        
                        # Grant access to linear-regression-data outside of this class
                        about = {'xcritical': 'List of X-critical points.',
                                 'ycritical': 'List of Y-critical points.',
                                 'xinflection': 'List of X-inflection points.',
                                 'yinflection': 'List of Y-inflection points.'}

                        outputs = {'xcritical': critical_x,
                                   'ycritical': critical_y,
                                   'xinflection': inflection_x,
                                   'yinflection': inflection_y}
                        self.outputs[name] = outputs
                        self.about[name] = about
                        
                    #---------------------------------------#
                    # Integrate the data and produce a plot #
                    #---------------------------------------#
                    if method == 'Calculus: Integrate Data':
                        figname = '{}_X={}_Y={}_integrals'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
                        setting = self.get_misc_setting(misc)
                        if 'order' in setting:
                            order = setting['order']
                        else: order='1,2'; misc += ' default-order=1,2;';
                        
                        if 'csv' in setting:
                            write_data = setting['csv']
                        else: write_data = False; misc += ' default-csv=False';
                        
                        if 'savefig' in setting:
                            c_savefig = str(setting['savefig'])
                        else:
                            c_savefig = 'all'; misc += ' default-savefig=all';
                        
                        self.log.out(self.format_analysis(method, xlo, xhi, misc, name))
                        self.integrate_data(x, y, xlo, xhi, order, write_data, mode['xlabel'], mode['ylabel'], c_savefig, figname, dpi)
                        
                    #--------------------------------------------------------------------------#
                    # Find the elastic constants using the "Regression Fringe Response" Method #
                    #--------------------------------------------------------------------------#
                    if method == 'Regression Fringe Response Thermal':
                        setting = self.get_misc_setting(misc)
                        savefig = setting.get('savefig', True)
                        if savefig:
                            figname = '{}_X={}_Y={}_RFR_thermal'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
                        else: figname = ''

                        
                        if 'minspan' in setting:
                            min_span = setting['minspan']
                        else: min_span=25; misc += ' default-minspace=25';
                        
                        if 'rtemp' in setting:
                            rtemp = setting['rtemp']
                        else: rtemp=tuple([]); misc += ' default-rtemp=()';
                        
                        # Get X- and Y-units
                        if '(' in mode['xlabel'] and ')' in mode['xlabel']:
                            tmp1 = mode['xlabel'].rpartition('(')
                            tmp2 = tmp1[2].rpartition(')')
                            x_units = '({})'.format(tmp2[0])
                            x_label = tmp1[0]
                        else:
                            x_units = ''
                            x_label = ''
                        if '(' in mode['ylabel'] and ')' in mode['ylabel']:
                            tmp1 = mode['ylabel'].rpartition('(')
                            tmp2 = tmp1[2].rpartition(')')
                            y_units = '({})'.format(tmp2[0])
                            y_label = tmp1[0]
                        else:
                            y_units = ''
                            y_label = ''
                        
                        self.log.out(self.format_analysis(method, xlo, xhi, misc, name))
                        rfr_thermal_outputs, fig = self.rfr_thermal(x, y, xlo, xhi, min_span, rtemp, setting, x_label, y_label, x_units, y_units, figname, dpi)
                        if fig is not None: self.figures['RFR_thermal'] = fig

                        # plot and log lo-CTE Results                        
                        x_reg = rfr_thermal_outputs['glassy_reg_x']
                        y_reg = rfr_thermal_outputs['glassy_reg_y']
                        b0, b1, r2 = rfr_thermal_outputs['glassy_line']
                        label = '{} y = {:.6f}x + {:.6f} ("Glassy-CTE"; $r^2$={:.6f})'.format(name, b1, b0, r2)
                        self.log.out('  {}'.format(label))
                        regdata = self.plot_parms(x=x_reg, y=y_reg, style='both', marker='|', line='--', size=12, label=label, shiftable=False)
                        data2plot.append(regdata)
                        
                        # plot and log hi-CTE Results                        
                        x_reg = rfr_thermal_outputs['rubbery_reg_x']
                        y_reg = rfr_thermal_outputs['rubbery_reg_y']
                        b0, b1, r2 = rfr_thermal_outputs['rubbery_line']
                        label = '{} y = {:.6f}x + {:.6f} ("Rubbery-CTE"; $r^2$={:.6f})'.format(name, b1, b0, r2)
                        self.log.out('  {}'.format(label))
                        regdata = self.plot_parms(x=x_reg, y=y_reg, style='both', marker='|', line='--', size=12, label=label, shiftable=False)
                        data2plot.append(regdata)
                        
                        # plot and log Tg results
                        intersection = rfr_thermal_outputs['tg']
                        label = '{} CTE intersection (Tg; x={:.4f}; y={:.4f})'.format(name, intersection[0], intersection[1])
                        self.log.out('  {}'.format(label))
                        tg_data = self.plot_parms(x=intersection[0], y=intersection[1], style='point', marker='p', line='-', size=12, label=label, shiftable=False)

                        # plot "linear regression extensions"
                        x_extension = [rfr_thermal_outputs['glassy_reg_x'][-1], intersection[0], rfr_thermal_outputs['rubbery_reg_x'][0]]
                        y_extension = [rfr_thermal_outputs['glassy_reg_y'][-1], intersection[1], rfr_thermal_outputs['rubbery_reg_y'][0]]
                        label = '{} CTE intersection extensions'.format(name)
                        ext_data = self.plot_parms(x=x_extension, y=y_extension, style='both', marker='|', line='--', size=12, label=label, shiftable=False)
                        data2plot.append(ext_data)
                        data2plot.append(tg_data) # Append this here so it is "in front"
                    
                        # plot and log transition data
                        transition = rfr_thermal_outputs['transition']
                        label = '{} transition region (xlo={:.4f}; xhi={:.4f})'.format(name, transition[0], transition[1])
                        self.log.out('  {}'.format(label))
                        miny = [min(y) for _ in transition]
                        transitiondata = self.plot_parms(x=transition, y=miny, style='both', marker='|', line='-', size=12, label=label, shiftable=False)
                        data2plot.append(transitiondata)
                        
                        
                        # Create nice table of outputs for easy copy and pasting
                        self.log.out('\n{:>2}-------------------Table of outputs-------------------'.format(''))
                        self.log.out('{:>4}{:<28}: {}'.format('', 'Tg', rfr_thermal_outputs['tg'][0]))
                        self.log.out('{:>4}{:<28}: {}'.format('', 'Glassy-CTE', rfr_thermal_outputs['glassy_cte']))
                        self.log.out('{:>4}{:<28}: {}'.format('', 'Rubbery-CTE', rfr_thermal_outputs['rubbery_cte']))
                        self.log.out('{:>4}{:<28}: {}'.format('', 'Glassy-slope', rfr_thermal_outputs['glassy_line'][1]))
                        self.log.out('{:>4}{:<28}: {}'.format('', 'Rubbery-slope', rfr_thermal_outputs['rubbery_line'][1]))
                        self.log.out('{:>4}{:<28}: {}'.format('', 'Glassy Volume', rfr_thermal_outputs['glassy_volume']))
                        self.log.out('{:>4}{:<28}: {}'.format('', 'Rubbery Volume', rfr_thermal_outputs['rubbery_volume']))
                        self.log.out('{:>4}{:<28}: {}'.format('', 'Glassy-transition', min(transition)))
                        self.log.out('{:>4}{:<28}: {}'.format('', 'Rubbery-transition', max(transition)))
                        self.log.out('{:>4}{:<28}: {}'.format('', 'Transition span',  max(transition)-min(transition)))
                        self.log.out('{:>2}{}\n'.format('', '------------------------------------------------------'))
                        
                        
                        # Grant access to Muzzy-modulus-data outside of this class
                        about = {'glassy_bounds': 'Lower bounds for "cold CTE"',
                                 'rubbery_bounds': 'Upper bounds for "hot CTE"',
                                 'glassy_line': 'Linear regression parameters for lo-CTE [b0, b1, r^2]',
                                 'rubbery_line': 'Linear regression parameters for hi-CTE [b0, b1, r^2]',
                                 'tg': 'Glass transition point [xpoint, ypoint]',                                 
                                 # 'glassy_slope':   'slope of glassy region',
                                 # 'rubbery_slope':  'slope of rubbery region',
                                 'glassy_volume':  'if rtemp provide this is the computed reference volume for the glassy-CTE as a float; else None',
                                 'rubbery_volume': 'if rtemp provide this is the computed reference volume for the rubbery-CTE as a float; else None',
                                 'glassy_cte':     'if rtemp provided Glassy-CTE as a float; else None',
                                 'rubbery_cte':    'if rtemp provided Rubbery-CTE as a float; else None',
                                 'transition': 'Transition region [xlo, xhi]',
                                 'glassy_reg_x':  'Lower linear regression X-data [x1, x2]', 
                                 'glassy_reg_y':  'Lower linear regression Y-data [y1, y2]', 
                                 'rubbery_reg_x': 'Upper linear regression X-data [x1, x2]', 
                                 'rubbery_reg_y': 'Upper linear regression Y-data [y1, y2]'}
                        
                        self.outputs[name] = rfr_thermal_outputs
                        self.about[name] = about

                     
                    #--------------------------------------------------------------------------#
                    # Find the elastic constants using the "Regression Fringe Response" Method #
                    #--------------------------------------------------------------------------#
                    # Check for 'Kemppainen-Muzzy Modulus' for backwards compatability. This will eventually be depricated
                    if method == 'Regression Fringe Response Mechanical' or method == 'Regression Fringe Response Modulus' or method == 'Kemppainen-Muzzy Modulus': 
                        figname = '{}_X={}_Y={}_RFR_Mechanical'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
                        setting = self.get_misc_setting(misc)
                        
                        # Raw is a hidden flag due to Tristan. He wants people to see
                        # the difference between "raw" fit or "cleaned" fit. This has
                        # caused to much confussion, so it is now hidden for developers
                        # ONLY!                         
                        if 'raw' in setting:
                            raw = setting['raw']
                        else: raw=False
                        
                        if 'shift' in setting:
                            shift = setting['shift']
                        else: shift='no'; misc += ' default-shift=no';
                        
                        if 'minxhi' in setting:
                            minxhi = setting['minxhi']
                        else: minxhi = 0.01; misc += ' default-minxhi=0.01';
                        
                        if 'maxxhi' in setting:
                            maxxhi = setting['maxxhi']
                        else: maxxhi = 0; misc += ' default-maxxhi=0';
                        
                        if 'yp' in setting:
                            yp = setting['yp']
                        else: yp = 0; misc += ' default-yp=0';
                        
                        if 'up' in setting:
                            up = setting['up']
                        else: up = 0; misc += ' default-up=0';
                        
                        if 'offset' in setting:
                            offset = setting['offset']
                        else: offset = 0; misc += ' default-offset=0';
                        
                        if 'xlo' in setting:
                            xlo_method = setting['xlo']
                        else: xlo_method = 'rfs'; misc += ' default-xlo=rfs';
                        
                        if 'csv' in setting:
                            write_data = setting['csv']
                        else: write_data = False; misc += ' default-csv=False';
                        
                        if 'ci' in setting:
                            confidence_interval = setting['ci']
                        else: confidence_interval = False; misc += ' default-ci=0';
                        
                        if 'ds' in setting:
                            derivative_span = setting['ds']
                        else: derivative_span = 0; misc += ' default-ds=0';
                        
                        if 'dd' in setting:
                            derivative_degree = setting['dd']
                        else: derivative_degree = 3; misc += ' default-dd=3';
                        
                        if 'grid' in setting:
                            grid = setting['grid']
                        else: grid = 'off'; misc += ' default-grid=off';
                        
                        if 'savefig' in setting:
                            rfr_savefig = str(setting['savefig'])
                        else:
                            rfr_savefig = 'all'; misc += ' default-savefig=all';
                            
                        if 't12_avg' in setting:
                            t12_avg = setting['t12_avg']
                        else:
                            t12_avg = False

                        
                        # get settings for poisson's ratio
                        t1 = []; t2 = [];
                        if 't1' in setting:
                            try: 
                                column = setting['t1']
                                if column in data: 
                                    t1 = [data[column]]
                                elif all(char in column for char in '${}'):
                                    log.out('  Attempting to use compute option for t1:')
                                    t1_data, message = vectorized_string_compute.thermo_data(data, column)
                                    t1 = [t1_data]
                                    log.out('    t1_compute={}'.format(column))
                                    log.out('    {}'.format(message))
                            except: t1 = []; self.log.GUI_error(f'ERROR t1 column {setting["t1"]} does not exist in logfile. Skipping poissons ratio calculation from -d(trans)/d(axial). Please check that t1, keywords, and sections is correct.')
                        else: t1 = []; misc += ' default-t1=[]'
                        
                        
                        if 't2' in setting:
                            try: 
                                column = setting['t2']
                                if column in data: 
                                    t2 = [data[column]]
                                elif all(char in column for char in '${}'):
                                    log.out('  Attempting to use compute option for t2:')
                                    t2_data, message = vectorized_string_compute.thermo_data(data, column)
                                    t2 = [t2_data]
                                    log.out('    t1_compute={}'.format(column))
                                    log.out('    {}'.format(message))
                            except: t2 = []; self.log.GUI_error(f'ERROR t2 column {setting["t2"]} does not exist in logfile. Skipping poissons ratio calculation from -d(trans)/d(axial). Please check that t1, keywords, and sections is correct.')
                        else: t2 = []; misc += ' default-t2=[]'
                        
                        
                        #------------------------------------------------------------------------------------------------------------#
                        # This entire function assumes strain increases with time (i.e. tensile or shear tests). Therefore if strain #
                        # decreases with time for compression or unloading tests, we need to "reverse" the stress and strain lists.  #
                        #------------------------------------------------------------------------------------------------------------#
                        flip_data = False
                        if x[-1] < x[0]:
                            flip_data = True
                            x.reverse()
                            y.reverse()
                            xlmp.reverse()
                            ylmp.reverse()
    
    
                        # For t1 and t2 we are going to have to process the data to make sure it is consistent with the axial stress (really only the 
                        # Moving average is required as that changes the data length. Overs are commented out, but can be uncommented at anytime).
                        if t1 and t2:                            
                            if nevery > 1:
                                t1 = [t1[0][::nevery]]
                                t2 = [t2[0][::nevery]]
                            t1_tmp = t1[0]; t2_tmp = t2[0]
                            if flip_data:
                                t1_tmp.reverse()
                                t2_tmp.reverse()
                            for lmp_cleaning in cleaning_settings:
                                if lmp_cleaning == 'LAMMPS data (apply linear-interp)':
                                    dx, n_points = cleaning_settings[lmp_cleaning]
                                    if t1_tmp: 
                                        self.log.out('\n  Applying linear-interp to "t1" transverse strain direcion.')
                                        xdummy, t1_tmp = self.interp_constant_dx(xlmp, t1_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, dx, n_points)
                                        t1 = [t1_tmp] # Bin average can not be used with raw data, so reset t1 completely
                                    if t2_tmp:
                                        self.log.out('\n  Applying linear-interp to "t2" transverse strain direcion.')
                                        xdummy, t2_tmp = self.interp_constant_dx(xlmp, t2_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, dx, n_points)
                                        t2 = [t2_tmp] # Bin average can not be used with raw data, so reset t1 completely
                                
                                if lmp_cleaning == 'LAMMPS data (apply bin average)':
                                    dx, start, stop = cleaning_settings[lmp_cleaning]
                                    if t1_tmp: 
                                        self.log.out('\n  Applying bin average to "t1" transverse strain direcion.')
                                        xdummy, t1_tmp = self.bin_average(xlmp, t1_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, dx, start, stop)
                                        t1 = [t1_tmp] # Bin average can not be used with raw data, so reset t1 completely
                                    if t2_tmp:
                                        self.log.out('\n  Applying bin average to "t2" transverse strain direcion.')
                                        xdummy, t2_tmp = self.bin_average(xlmp, t2_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, dx, start, stop)
                                        t2 = [t2_tmp] # Bin average can not be used with raw data, so reset t1 completely
                                
                                
                                if lmp_cleaning == 'LAMMPS data (apply moving average)':
                                    LAMMPS_window = cleaning_settings[lmp_cleaning]
                                    if t1_tmp: 
                                        self.log.out('\n  Applying moving average to "t1" transverse strain direcion.')
                                        xdummy, t1_tmp = self.moving_average(xlmp, t1_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, LAMMPS_window)
                                        t1 = [t1_tmp] # Moving average can not be used with raw data, so reset t1 completely
                                    if t2_tmp:
                                        self.log.out('\n  Applying moving average to "t2" transverse strain direcion.')
                                        xdummy, t2_tmp = self.moving_average(xlmp, t2_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, LAMMPS_window)
                                        t2 = [t2_tmp] # Moving average can not be used with raw data, so reset t1 completely
                                
                                if lmp_cleaning == 'LAMMPS data (apply Whittaker-Eilers)':
                                    LAMMPS_order, LAMMPS_lmbda = cleaning_settings[lmp_cleaning]
                                    if t1_tmp: 
                                        self.log.out('\n  Applying Whittaker-Eilers to "t1" transverse strain direcion.')
                                        xdummy, t1_tmp = self.whittaker_eilers_smoothing(xlmp, t1_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, LAMMPS_order, LAMMPS_lmbda, 'Axial strain', 'Transverse strain (t1)', figname+'_t1')
                                        t1.append(t1_tmp)
                                    if t2_tmp:
                                        self.log.out('\n  Applying Whittaker-Eilers to "t2" transverse strain direcion.')
                                        xdummy, t2_tmp = self.whittaker_eilers_smoothing(xlmp, t2_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, LAMMPS_order, LAMMPS_lmbda, 'Axial strain', 'Transverse strain (t2)', figname+'_t2') 
                                        t2.append(t2_tmp)
                                
                                if lmp_cleaning == 'LAMMPS data (apply Butterworth filter)':
                                    LAMMPS_wn, LAMMPS_order, LAMMPS_qm, LAMMPS_psd, LAMMPS_write_data, LAMMPS_savefig = cleaning_settings[lmp_cleaning]
                                    if t1_tmp:
                                        self.log.out('\n  Applying Butterworth filter to "t1" transverse strain direcion.')
                                        xdummy, t1_tmp, t1_wn, t1_qm = self.butterworth_lowpass(xlmp, t1_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, LAMMPS_wn, LAMMPS_order, LAMMPS_qm, LAMMPS_psd, LAMMPS_write_data, LAMMPS_savefig, figname+'_t1', dpi)
                                        t1.append(t1_tmp)
                                    if t2_tmp:
                                        self.log.out('\n  Applying Butterworth filter to "t2" transverse strain direcion.')
                                        xdummy, t2_tmp, t2_wn, t2_qm = self.butterworth_lowpass(xlmp, t2_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, LAMMPS_wn, LAMMPS_order, LAMMPS_qm, LAMMPS_psd, LAMMPS_write_data, LAMMPS_savefig, figname+'_t2', dpi)
                                        t2.append(t2_tmp)
                                
                                if lmp_cleaning == 'LAMMPS data (fit polynomial)':
                                    LAMMPS_deg = cleaning_settings[lmp_cleaning]
                                    if t1_tmp:
                                        self.log.out('\n  Fitting polynomial to "t1" transverse strain direcion.')
                                        xdummy, t1_tmp = self.polynomial_fit(xlmp, t1_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, LAMMPS_deg)
                                        t1.append(t1_tmp)
                                    if t2_tmp:
                                        self.log.out('\n  Fitting polynomial to "t2" transverse strain direcion.')
                                        xdummy, t2_tmp = self.polynomial_fit(xlmp, t2_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, LAMMPS_deg)
                                        t2.append(t2_tmp)
                                        
                                if lmp_cleaning == 'LAMMPS data (LOWESS)':
                                    LAMMPS_fraction, LAMMPS_max_iter = cleaning_settings[lmp_cleaning]
                                    if t1_tmp:
                                        self.log.out('\n  Applying LOWESS to "t1" transverse strain direcion.')
                                        xdummy, t1_tmp = self.lowess(xlmp, t1_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, LAMMPS_fraction, LAMMPS_max_iter)
                                        t1.append(t1_tmp)
                                    if t2_tmp:
                                        self.log.out('\n  Applying LOWESS to to "t2" transverse strain direcion.')
                                        xdummy, t2_tmp = self.lowess(xlmp, t2_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, LAMMPS_fraction, LAMMPS_max_iter)
                                        t2.append(t2_tmp)
                                        
                                if lmp_cleaning == 'LAMMPS data (apply iFFT filter)':
                                    threshold, LAMMPS_qm, LAMMPS_savefig = cleaning_settings[lmp_cleaning]
                                    if t1_tmp:
                                        self.log.out('\n  Applying iFFT filter to "t1" transverse strain direcion.')
                                        t1_tmp, qm_iFFT = self.iFFT_filter(xlmp, t1_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, threshold, LAMMPS_qm, LAMMPS_savefig,  figname+'_t1', dpi)
                                        t1.append(t1_tmp)
                                    if t2_tmp:
                                        self.log.out('\n  Applying iFFT filter to "t2" transverse strain direcion.')
                                        t2_tmp, qm_iFFT = self.iFFT_filter(xlmp, t2_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, threshold, LAMMPS_qm, LAMMPS_savefig,  figname+'_t2', dpi)
                                        t2.append(t2_tmp)
                                        
                                if lmp_cleaning == 'LAMMPS data (Savitzky-Golay)':
                                    window, order = cleaning_settings[lmp_cleaning]
                                    if t1_tmp:
                                        self.log.out('\n  Applying Savitzky-Golay smoother to "t1" transverse strain direcion.')
                                        x_tmp, t1_tmp = self.SavitzkyGolay(xlmp, t1_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, window, order)
                                        t1.append(t1_tmp)
                                    if t2_tmp:
                                        self.log.out('\n  Applying Savitzky-Golay smoother to "t2" transverse strain direcion.')
                                        x_tmp, t2_tmp = self.SavitzkyGolay(xlmp, t2_tmp, LAMMPS_data_xlo, LAMMPS_data_xhi, window, order)
                                        t2.append(t2_tmp)
                        
                        # Get stress units
                        if '(' in mode['ylabel'] and ')' in mode['ylabel']:
                            tmp1 = mode['ylabel'].rpartition('(')
                            tmp2 = tmp1[2].rpartition(')')
                            stress_units = '({})'.format(tmp2[0])
                        else: stress_units = ''
                        
                        # Perform the RFR-analysis
                        self.log.out(self.format_analysis(method, xlo, xhi, misc, name))
                        xlo_out, xhi_out, yield_point_derivative, yield_point_offset, nu1, nu2, nu12_avg, ultimate_point, fig1, fig2 = self.rfr_mechanical(x, y, xlo, xhi, minxhi, maxxhi, xlo_method, yp, up,
                                                                                                                                                           offset, t1, t2, stress_units, write_data, derivative_span,
                                                                                                                                                           derivative_degree, grid, t12_avg, rfr_savefig, figname, dpi)
                        
                        # Fit linear regression model to x/y and xlmp/ylmp (if no "LAMMPS data (i)" is used x/y and xlmp/ylmp will be the same )
                        b0_clean, b1_clean, xreg_clean, yreg_clean, r2_clean, ci_clean = self.linear_regression(x, y, xlo_out, xhi_out, confidence_interval)
                        b0_raw, b1_raw, xreg_raw, yreg_raw, r2_raw, ci_raw = self.linear_regression(xlmp, ylmp, xlo_out, xhi_out, confidence_interval)
                        b0 = (b0_clean + b0_raw)/2 # use an average of b0's for shifting and determining apparent offset
                        b1 = (b1_clean + b1_raw)/2 # use an average of b1's for shifting and determining apparent offset 
                        
                        if fig1 is not None: self.figures['RFR_mechanical_stress_strain'] = fig1
                        if fig2 is not None: self.figures['RFR_mechanical_strain_strain'] = fig2
                        
                        # If shift equals 'yint', shift data by the yintercept
                        if shift == 'yint': 
                            shift_amount = b0
                            self.log.out(f'  Shifting all Y-data by the y-intercept={shift_amount} of the linear regression model')
                            b0 = 0 
                            b0_raw -= shift_amount
                            b0_clean -= shift_amount
                        
                        # If shift equals 'ymin', shift data by the minimum of y
                        if shift == 'ymin':
                            max_index = y.index(max(y))
                            ymin = min(y[:max_index])
                            shift_amount = ymin
                            self.log.out(f'  Shifting all Y-data by min(y-data)={shift_amount}')
                            b0 -= shift_amount
                            b0_raw -= shift_amount
                            b0_clean -= shift_amount
                        
                        # Apply shift to yield points
                        yp_derivative = yield_point_derivative
                        yp_offset = yield_point_offset
                        if yp_derivative:
                            yp_derivative[1] -= shift_amount
                        if yp_offset:
                            yp_offset[1] -= shift_amount
                        if ultimate_point:
                            ultimate_point[1] -= shift_amount

                        # Clean the data based on standing stress wave
                        # if yp_derivative:
                        #     rfr_mechanical.check_yp(x, y, xlmp, ylmp, yp_derivative)
                            
                        # Plot "raw fit"
                        if raw:
                            label = '{} y = {:.6f}x + {:.6f} ("raw fit"; shifted by {:.4f}; $r^2$={:.6f})'.format(name, b1_raw, b0_raw, shift_amount, r2_raw)
                            self.log.out('  RFR equation of line {}'.format(label))
                            regdata_raw = self.plot_parms(x=xreg_raw, y=yreg_raw, style='line', marker='.', line='-', size=3, label=label, shiftable=True)
                            data2plot.append(regdata_raw)
                        
                        # Plot "clean fit"
                        if raw:
                            label = '{} y = {:.6f}x + {:.6f} ("clean fit"; shifted by {:.4f}; $r^2$={:.6f})'.format(name, b1_clean, b0_clean, shift_amount, r2_clean)
                        else:
                            label = '{} y = {:.6f}x + {:.6f} (shifted by {:.4f}; $r^2$={:.6f})'.format(name, b1_clean, b0_clean, shift_amount, r2_clean)
                        self.log.out('  RFR equation of line {}'.format(label))
                        regdata_clean = self.plot_parms(x=xreg_clean, y=yreg_clean, style='line', marker='.', line='-', size=3, label=label, shiftable=True)
                        data2plot.append(regdata_clean)
                                            
                        # Plot yield point calculations
                        if yp_derivative:
                            label = '{} yield point from derivative (x={}; y={})'.format(name, yp_derivative[0], yp_derivative[1])
                            yielding = self.plot_parms(x=yp_derivative[0], y=yp_derivative[1], style='point', marker='o', line='-', size=8, label=label, shiftable=False)
                            data2plot.append(yielding)
                            
                            # Determine apparent offset based on yield strength from yield_point_derivative
                            yield_yvalue = yp_derivative[1]
                            slope_xvalue = (yield_yvalue - b0)/b1
                            apparent_offset = abs(slope_xvalue - yp_derivative[0])
                            self.log.out('  RFR Modulus found the "apparent offset" to compute the yield strength to be {:.6f}'.format(apparent_offset))
                        if yp_offset:
                            label = '{} yield point from offset (x={}; y={})'.format(name, yp_offset[0], yp_offset[1])
                            yielding = self.plot_parms(x=yp_offset[0], y=yp_offset[1], style='point', marker='o', line='-', size=8, label=label, shiftable=False)
                            data2plot.append(yielding)
                            
                        # Plot ultimate point calculations
                        if ultimate_point:
                            label = '{} ultimate point from derivative (x={}; y={})'.format(name, ultimate_point[0], ultimate_point[1])
                            ultimate = self.plot_parms(x=ultimate_point[0], y=ultimate_point[1], style='point', marker='o', line='-', size=8, label=label, shiftable=False)
                            data2plot.append( ultimate)
                        
                        # Log some outputs
                        self.log.out(f'  RFR Modulus found xlo={xlo_out} and xhi={xhi_out}')
                        if yp_derivative:
                            self.log.out('  RFR Modulus found the yield point from derivative as ({:.6f}, {:.6f})'.format(yp_derivative[0], yp_derivative[1]))
                        if yp_offset:
                            self.log.out('  RFR Modulus found the yield point from offset as ({:.6f}, {:.6f})'.format(yp_offset[0], yp_offset[1]))
                        if nu1 is not None and nu2 is not None:
                            nu_avg = (nu1 + nu2)/2
                            self.log.out("  RFR Modulus found the poissons ratio's to be nu_1={:.6f}, nu_2={:.6f}, nu_avg={:.6f}".format(nu1, nu2, nu_avg))
                        else: nu_avg = None
    
                        # Create nice table of outputs for easy copy and pasting
                        self.log.out('\n{:>2}-------------------Table of outputs-------------------'.format(''))
                        if raw:
                            self.log.out('{:>4}{:<28}: {}'.format('', 'Modulus "raw"', b1_raw))
                            self.log.out('{:>4}{:<28}: {}'.format('', '  - proportional limit (x)', xreg_raw[-1]))
                            self.log.out('{:>4}{:<28}: {}'.format('', '  - proportional limit (y)', yreg_raw[-1]))
                            self.log.out('{:>4}{:<28}: {}'.format('', '  - Y-intercept', b0_raw))
                            if ci_raw is not None:
                                b0_name = '{:<28}'.format('  - Y-intercept (C.I. '+str(confidence_interval)+')')
                                b1_name = '{:<28}'.format('  - Slope (C.I. '+str(confidence_interval)+')')
                                ci_b0 = ci_raw[0]
                                ci_b1 = ci_raw[1]
                                self.log.out('{:>4}{:<28}: {} to {}'.format('', b0_name, ci_b0[0], ci_b0[1]))
                                self.log.out('{:>4}{:<28}: {} to {}'.format('', b1_name, ci_b1[0], ci_b1[1]))
                            self.log.out('{:>4}{:<28}: {}'.format('', '  - r^2', r2_raw))
                            self.log.out('')
                        
                        if raw:
                            self.log.out('{:>4}{:<28}: {}'.format('', 'Modulus "clean"', b1_clean))
                        else:
                            self.log.out('{:>4}{:<28}: {}'.format('', 'Modulus', b1_clean))
                        self.log.out('{:>4}{:<28}: {}'.format('', '  - proportional limit (x)', xreg_clean[-1]))
                        self.log.out('{:>4}{:<28}: {}'.format('', '  - proportional limit (y)', yreg_clean[-1]))
                        self.log.out('{:>4}{:<28}: {}'.format('', '  - Y-intercept', b0_clean))
                        if ci_clean is not None:
                            b0_name = '{:<28}'.format('  - Y-intercept (C.I. '+str(confidence_interval)+')')
                            b1_name = '{:<28}'.format('  - Slope (C.I. '+str(confidence_interval)+')')
                            ci_b0 = ci_clean[0]
                            ci_b1 = ci_clean[1]
                            self.log.out('{:>4}{:<28}: {} to {}'.format('', b0_name, ci_b0[0], ci_b0[1]))
                            self.log.out('{:>4}{:<28}: {} to {}'.format('', b1_name, ci_b1[0], ci_b1[1]))
                        self.log.out('{:>4}{:<28}: {}'.format('', '  - r^2', r2_clean))
                        if yp_derivative:
                            self.log.out('')
                            self.log.out('{:>4}{:<28}: {}'.format('', 'Yield point "derivative"', yp_derivative[1]))
                            self.log.out('{:>4}{:<28}: {}'.format('', '  - "apparent offset"', apparent_offset))
                        if yp_offset:
                            self.log.out('')
                            self.log.out('{:>4}{:<28}: {}'.format('', 'Yield point "offset"', yp_offset[1]))
                        if nu1 is not None and nu2 is not None:
                            self.log.out('')
                            self.log.out('{:>4}{:<28}: {}'.format('', 'Poissons ratio "nu_1"', nu1))
                            self.log.out('{:>4}{:<28}: {}'.format('', 'Poissons ratio "nu_2"', nu2))
                            self.log.out('{:>4}{:<28}: {}'.format('', 'Poissons ratio "nu_avg"', nu_avg))
                            if t12_avg: self.log.out('{:>4}{:<28}: {}'.format('', 'Poissons ratio "nu12_avg"', nu12_avg))
                        self.log.out('{:>2}{}\n'.format('', '----------------------------------------------------'))
                        
                        self.log.out(f'\n  {cite_string}')
                        self.log.out('  This method implements the work found at: https://doi.org/10.26434/chemrxiv-2025-fk935')
                        self.log.out('  which should be cited if used for analyzing your LAMMPS data.')
                        self.log.out(f'  {cite_string}\n')
                        
                        # Provide a statement about "raw" vs "clean" Modulus
                        if raw:
                            self.log.out('  NOTE: Modulus "raw" means fit a linear regression to "raw" LAMMPS data and Modulus "clean" means fit')
                            self.log.out('  linear regression to cleaned data. The data could be cleaned using any of the following methods:')
                            self.log.out('    - LAMMPS data (apply moving average)')
                            self.log.out('    - LAMMPS data (apply Butterworth filter)')
                            self.log.out('    - LAMMPS data (apply Whittaker-Eilers)')
                            self.log.out('    - LAMMPS data (fit polynomial)')
                            self.log.out('  If none of the above methods are used Modulus "clean" will be fit to the "raw" LAMMPS data and will')
                            self.log.out('  be the same as Modulus "raw".\n\n')
                
                        # Grant access to Muzzy-modulus-data outside of this class
                        about = {'xdata': 'List of X-data defining the linear-regresion model fit to "cleaned MD data", order [xlo, xhi].',
                                 'ydata': 'List of Y-data defining the linear-regresion model fit to "cleaned MD data", order [ylo, yhi].',
                                 'b1': 'Float value of b1 parameter (slope) for general equation of line y = b1*x + b0, when linear regression is fit to "cleaned MD data".',
                                 'b0': 'Float value of b0 parameter (Y-intercept) for general equation of line y = b1*x + b0, when linear regression is fit to "cleaned MD data".',
                                 'b1-ci': 'b1 confidence_interval of {} or None if not computed'.format(confidence_interval),
                                 'b0-ci': 'b0 confidence_interval of {} or None if not computed'.format(confidence_interval),
                                 'xlo': 'Float value of X-lo value of the found linear-region of the stress-strain curve',
                                 'xhi': 'Float value of X-hi value of the found linear-region of the stress-strain curve',
                                 'yield_point_derivative': 'List of floats found from yield point determination, using derivative methods. Order [X-yp, Y-yp] or [None, None] if not using yp. NOTE: SHIFTED based on shift method',
                                 'yield_point_offset': 'List of floats found from yield point determination, using offset methods. Order [X-yp, Y-yp] or [None, None] if not using offset. NOTE: SHIFTED based on shift method',
                                 'nu1': 'Float defining Poissons ratio from tranverse1 (t1) and axial strain (X-data) linear-relation. If not using t1, returns None',
                                 'nu2': 'Float defining Poissons ratio from tranverse1 (t2) and axial strain (X-data) linear-relation. If not using t2, returns None',
                                 'nu_avg': 'Float defining Poissons ratio found be averaging nu1 and nu2 together. If not using t1 and t2, returns None',
                                 'ultimate_point': 'List of floats found from ultimate point determination, using derivative methods. Order [X-up, Y-up] or [None, None] if not using up. NOTE: SHIFTED based on shift method',}
                        if not yp_derivative:
                            yp_derivative = [None, None]
                        if not yp_offset:
                            yp_offset = [None, None]
                        if not ultimate_point:
                            ultimate_point = [None, None]
                        outputs = {'xdata': xreg_clean,
                                   'ydata': yreg_clean,
                                   'b1': b1_clean,
                                   'b0': b0_clean,
                                   'xlo': xlo_out,
                                   'xhi': xhi_out,
                                   'yield_point_derivative': yp_derivative,
                                   'yield_point_offset': yp_offset,
                                   'nu1': nu1,
                                   'nu2': nu2,
                                   'nu_avg': nu_avg,
                                   'ultimate_point': ultimate_point}
                        
                        if t12_avg:
                            about['nu12_avg'] = 'Float defining Poissons ratio from average t1 and t2 and axial strain (X-data) linear-relation. If not using t1 and t2, returns None'
                            outputs['nu12_avg'] = nu12_avg
                        if ci_clean is not None:
                            outputs['b0-ci'] = ci_clean[0]
                            outputs['b1-ci'] = ci_clean[1]
                        else:
                            outputs['b0-ci'] = None
                            outputs['b1-ci'] = None
                            
                        if ci_raw is not None:
                            outputs['b0-ci-raw'] = ci_raw[0]
                            outputs['b1-ci-raw'] = ci_raw[1]
                        else:
                            outputs['b0-ci-raw'] = None
                            outputs['b1-ci-raw'] = None
                            
                        # Update based on raw setting
                        if raw:
                            about['xdata-raw'] = 'List of X-data defining the linear-regresion model fit to "raw MD data", order [xlo, xhi].'
                            about['ydata-raw'] = 'List of Y-data defining the linear-regresion model fit to "raw MD data", order [ylo, yhi].'
                            about['b1-raw'] = 'Float value of b1 parameter (slope) for general equation of line y = b1*x + b0, when linear regression is fit to "raw MD data".'
                            about['b0-raw'] = 'Float value of b0 parameter (Y-intercept) for general equation of line y = b1*x + b0, when linear regression is fit to "raw MD data".'
                            about['b1-raw-ci'] = 'b1 confidence_interval of {} or None if not computed'.format(confidence_interval)
                            about['b0-raw-ci'] = 'b0 confidence_interval of {} or None if not computed'.format(confidence_interval)
                            about['b0-raw-ci'] = 'b0 confidence_interval of {} or None if not computed'.format(confidence_interval)
                            about['b1-raw-ci'] = 'b1 confidence_interval of {} or None if not computed'.format(confidence_interval)
                            outputs['xdata-raw'] = xreg_raw
                            outputs['ydata-raw'] = yreg_raw
                            outputs['b1-raw'] = b1_raw
                            outputs['b0-raw'] = b0_raw  
                            
                        self.outputs[name] = outputs
                        self.about[name] = about
                
                #----------------------------#
                # Generate plot if asked for #
                #----------------------------#
                if plot: 
                    fig, ax = plt.subplots()
                    color_index = 0
                
                maxpoints = 5_000 # Maximum number of points to plot and save to .csv
                csv_data  = {} # { label: [[xdata], [ydata]], ... }
                for data in data2plot:
                    xx = data['x']
                    yy = data['y']
                    label = data['label'].replace('\n', ' ')
                    style = data['style']
                    shiftable = data['shiftable']
                    marker = data['marker']
                    line = data['line']
                    size = data['size']
                    if plot: 
                        color = colors[color_index]
                        
                    # Limit the max points plotted
                    if isinstance(xx, (list, tuple, np.ndarray)) and style in ['line', 'point', 'both']:
                        if len(xx) > maxpoints:
                            self.log.warn('WARNING during plotting: label="{}" has npoints={}, for performance reasons data is being sampled to npoints={}'.format(label, len(xx), maxpoints))
                            xx, yy = self.limit_points(xx, yy, maxpoints)
                    
                    # Skip lammps data if rm_lmp_data and label is 'LAMMPS data'
                    if rm_lmp_data and label == 'LAMMPS data': continue
                    
                    # If shift, shift all data
                    if shift and shiftable:
                        try: yy = [i - shift_amount for i in yy]
                        except: 
                            try: yy -= shift_amount
                            except: pass
                    
                    # Plot data with or without line width size
                    if style == 'line':
                        if plot: plt.plot(xx, yy, line, color=color, lw=size, label=label)
                        csv_data[label] = [list(xx), list(yy)]
                    elif style == 'point':
                        if plot: plt.plot(xx, yy, marker, color=color, ms=size, label=label)
                        try: tmpx = list(xx); tmpy = list(yy)
                        except: tmpx = [xx]; tmpy = [yy]
                        csv_data[label] = [tmpx, tmpy]
                    elif style == 'both':
                        if plot: plt.plot(xx, yy, marker, color=color, ls=line, lw=size/4, ms=2*size, label=label)
                        csv_data[label] = [list(xx), list(yy)]
                    elif style == 'vertical':
                        if plot: plt.axvline(xx, color=color, ls=marker, lw=size, label=label)
                        csv_data[label] = [[xx], ['cursor-x']]
                    elif style == 'horizontal':
                        if plot: plt.axhline(yy, color=color, ls=marker, lw=size, label=label)
                        csv_data[label] = [['cursor-x'], [yy]]
                    else: raise Exception(f'ERROR {style} not supported')
                    
                    # increment color index and rest to zero if to large
                    if plot:
                        color_index += 1
                        if color_index + 1 > len(colors):
                            color_index = 0
                        
                # Set labels
                if plot:
                    if mode['xlabel']:
                        ax.set_xlabel(mode['xlabel'])
                    else: ax.set_xlabel(xdata)
                    if mode['ylabel']:
                        ax.set_ylabel(mode['ylabel'])
                    else: ax.set_ylabel(ydata)
                    
                    # Set legend and size
                    ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fancybox=True, ncol=1)
                    fig.set_size_inches(12, 4, forward=True)
                    fig.tight_layout()
                    
                    self.figures['main'] = fig
                    
                    
                    # Save current image
                    if savefig:
                        figname = '{}_X={}_Y={}.jpeg'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
                        fig.savefig(figname, dpi=dpi)
                
                #----------------------------------#
                # Save data from plot to csv files #
                #----------------------------------#
                if write_plotted_data: self.savedata_to_csv(xdata, ydata, csv_data)
            
                #------------------#
                # Analysis wrap-up #
                #------------------#
                # Script run time
                execution_time = (time.time() - start_time)
                self.log.out('Execution time in seconds: ' + str(execution_time))
                
                # write log
                logname = '{}_X={}_Y={}.log.lunar'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
                self.log.write_logged(logname)
        
        
    ########################################################
    # Start of methods that are used throughout this class #
    ########################################################
    #----------------------------------------------#
    # Method to limit the number of points plotted #
    #----------------------------------------------#
    def limit_points(self, x, y, maxpoints):
        if len(x) != len(y):
            raise ValueError("x and y must have the same length")
        if maxpoints < 1:
            raise ValueError("maxpoints must be at least 1")
        if len(x) <= maxpoints:
            return list(x), list(y)
    
        n = len(x)
    
        # Choose evenly spaced indices. When keeping multiple points,
        # this includes the first and last points.
        if maxpoints == 1:
            indices = [0]
        else:
            indices = [round(i * (n - 1) / (maxpoints - 1)) for i in range(maxpoints)]
    
        return [x[i] for i in indices], [y[i] for i in indices]
    
    #------------------------------------------#
    # Method to write plotted data to csv file #
    #------------------------------------------#
    def savedata_to_csv(self, xdata, ydata, csv_data):
        csvname = '{}_X={}_Y={}.csv'.format(self.get_basename_and_builder_dirs(), xdata, ydata)
        self.log.out(f'\n\nWriting {csvname} ...')
        misc_funcs.savedata_to_csv(csv_data, csvname)
        return

    #----------------------------------------------------------------#
    # Method to get basename of file and build directories if needed #
    #----------------------------------------------------------------#
    def get_basename_and_builder_dirs(self):
        # Get currently defined parent_directory variable
        logfile = self.logfile
        parent_directory = self.parent_directory
        
        # Find present working directory
        pwd = os.getcwd()
        
        # Find/create paths to store code results
        path = os.path.join(pwd, parent_directory)
        
        # Going to use io_functions get_dir_from_topofile() function which uses 'topofile'
        # string not 'logfile' string. So convert 'logfile' to 'topofile' if applicable
        if 'logfile' in parent_directory:
            parent_directory = parent_directory.replace('logfile', 'topofile')
            self.log.out('  Using path from logfile to set parent_directory ...')
            path = io_functions.get_dir_from_topofile(logfile, parent_directory)
            
        # Check if path exists. IF not create.
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
            
        # Set basename 
        root = os.path.basename(logfile)
        basename = os.path.join(path, root)
        return basename
                
    #-------------------------------------------#
    # Method to get misc dict of key/value pair #
    #-------------------------------------------#
    def get_misc_setting(self, misc, use_eval=True):
        # Setup the globals namespace to limit scope of what eval() can do
        allowed_builtins = ['min','max','sum','abs','len','map','range','reversed']
        copied_builtins = globals()['__builtins__'].copy()
        globals_dict = {}
        globals_dict['__builtins__'] = {key:copied_builtins[key] for key in allowed_builtins}
        
        # Parse misc string
        setting = {} # {keyword:float or int or Boolean}
        tmp1 = misc.split(';')
        for tmp2 in tmp1:
            tmp3 = tmp2.split('=')
            if len(tmp3) >= 2:
                i = tmp3[0].strip()
                if use_eval:
                    try: j = eval(tmp3[1], globals_dict)
                    except: j = str(tmp3[1])
                else: j = str(tmp3[1])
                setting[i] = j
        return setting
                
    #----------------------------------------#
    # Method to format analysis for log file #
    #----------------------------------------#
    def format_analysis(self, method, xlo, xhi, misc, name):
        text = '\n\n--------------------------------------------------------------------------------------\n'
        txt = 'method={}; xlo={}; xhi={}; misc={}; name="{}"'.format(method, xlo, xhi, misc, name)
        chunks = len(text)
        nchunks = math.ceil(len(txt)/chunks)
        chunk_size = math.ceil(len(txt)/nchunks)
        tmp = [ txt[(i-1)*chunk_size:i*chunk_size] for i in range(1, nchunks+1) ]
        for i in tmp:
            text += '{}\n'.format(i)
        text += '--------------------------------------------------------------------------------------'
        return text
    
    #-------------------------------------------------------------------------#
    # Method to build data2plot dict. The following meanings:                 #
    #    x = list/array of xdata to plot                                      #
    #    y = list/array of xdata to plot                                      #
    #    style = 'marker' or 'line' or 'both' or 'horizontal' or 'vertical'   #
    #    marker = marker style                                                #
    #    line = line style                                                    #
    #    size = int (for 'point' or 'line' style)                             #
    #    label = 'to put in legend'                                           #
    #    shiftable = Boolean, whether the data is shiftable or not (lin-reg   #
    #                for stress-v-strain)                                     #
    #-------------------------------------------------------------------------#
    def plot_parms(self, x=[], y=[], style='point', marker='.', line='-', size=4, label='default', shiftable=False):
        plot_setup = {}
        plot_setup['x'] = x
        plot_setup['y'] = y
        plot_setup['size'] = size
        plot_setup['style'] = style
        plot_setup['marker'] = marker
        plot_setup['line'] = line
        plot_setup['label'] = label
        plot_setup['shiftable'] = shiftable
        return plot_setup
    
    #---------------------------------------------------------#
    # Method to linear interpolate data onto constant dx grid #
    #---------------------------------------------------------#
    def interp_constant_dx(self, x, y, xlo, xhi, dx, n_points):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            xinterp, yinterp = misc_funcs.interp_constant_dx(reduced_x, reduced_y, dx, n_points)
        else:
            xinterp, yinterp = x, y
            self.log.GUI_error(f'ERROR (bin average) no LAMMPS data in xrange {xlo} - {xhi}. Returning uninterpolated data')
        return list(xinterp), list(yinterp)
    
    
    #---------------------------------------#
    # Method implementing a box car average #
    #---------------------------------------#
    def bin_average(self, x, y, xlo, xhi, dx, start, stop):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            x_binavg, y_binavg = misc_funcs.bin_average(reduced_x, reduced_y, dx, start=start, stop=stop)
        else:
            x_binavg = [0, 1]; y_binavg = [0, 1]
            self.log.GUI_error(f'ERROR (bin average) no LAMMPS data in xrange {xlo} - {xhi}')
        return list(x_binavg), list(y_binavg)
    
    #-----------------------------------------------#
    # Method implementing a moving average function #
    #-----------------------------------------------#
    def moving_average(self, x, y, xlo, xhi, window):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            x_movavg = misc_funcs.moving_average(reduced_x, window)
            y_movavg = misc_funcs.moving_average(reduced_y, window)
        else:
            x_movavg = [0, 1]; y_movavg = [0, 1]
            self.log.GUI_error(f'ERROR (moving average) no LAMMPS data in xrange {xlo} - {xhi}')
        return list(x_movavg), list(y_movavg)
    
    #-----------------------------------------------#
    # Method implementing a LOWESS function #
    #-----------------------------------------------#
    def lowess(self, x, y, xlo, xhi, fraction, max_iter):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            xout, yout = misc_funcs.lowess(reduced_x, reduced_y, fraction=fraction, max_iter=max_iter)
        else:
            xout = [0, 1]; yout = [0, 1]
            self.log.GUI_error(f'ERROR (LOWESS) no LAMMPS data in xrange {xlo} - {xhi}')
        return list(xout), list(yout)
    
    #--------------------------------------------------#
    # Method implementing a lowpass butterworth filter #
    #--------------------------------------------------#
    def butterworth_lowpass(self, x, y, xlo, xhi, wn, order, qm, psd, write_data, savefig, figname, dpi):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            # Convert to numpy arrays before passing to the butterworth filter
            reduced_x = np.array(reduced_x); reduced_y = np.array(reduced_y)
            
            # Setup the different savefig options
            wn_savefig = False; qm_savefig = False
            if '0' not in savefig and '1' in savefig or 'all' in savefig:
                wn_savefig = True
            if '0' not in savefig and '2' in savefig or 'all' in savefig:
                qm_savefig = True

            # Optimized wn if user desires and the log value            
            wn_method = str(wn)
            if wn_method.startswith('op') or wn_method.startswith('oa'):
                wn = signal_processing.butter_optimize_wn_with_FFT(reduced_x, reduced_y, wn_method, write_data, wn_savefig, figname+'_FFT_PSD_wn', dpi)
            if wn_method.startswith('or'):
                wn = signal_processing.butter_optimize_wn_with_residuals(reduced_x, reduced_y, order, qm, wn_method, wn_savefig, figname+'_Residuals', dpi, delta=0.001)
            if wn_method.startswith('cve'):
                wn_error, delta, stop, n = False, 0.001, 1.0, 1
                if '<' in wn_method and '>' in wn_method:
                    try:
                        values = [float(i) for i in wn_method.split('<')[-1].split('>')[0].split(',')]
                        if len(values) == 3: delta, stop, n = values
                        else: wn_error = True
                    except: wn_error = True
                    if wn_error:
                        self.log.GUI_error(f'ERROR could not parse wn={wn_method} string. Correct syntax: wn=cve<delta, stop, n>-p')
                        self.log.GUI_error('Using defaults of delta={}, stop={}, n={} for cross-validation calculations'.format(delta, stop, n))
                wn = signal_processing.butter_optimize_wn_with_CVE(reduced_x, reduced_y, n, wn_method, order, delta, stop, qm='1,1', savefig=wn_savefig, figname=figname+'_LOOCVE', dpi=dpi)

            # Plot PSD
            if psd and isinstance(wn, float) and not wn_method.startswith('op'):
                signal_processing.generate_and_plot_fft_wn(reduced_x, reduced_y, wn, write_data, wn_savefig, figname, dpi)

            # Apply filter
            xfilter = reduced_x
            yfilter, qm = signal_processing.butter_lowpass_filter(reduced_x, reduced_y, wn, order, qm, write_data, qm_savefig, figname+'_QM', dpi, pflag=True)
        else:
            xfilter = reduced_x
            yfilter = reduced_y
            self.log.GUI_error(f'ERROR (butterworth low pass) no LAMMPS data in xrange {xlo} - {xhi}')
        return list(xfilter), list(yfilter), wn, qm
    
    
    #--------------------------------------#
    # Method implementing a Fourier filter #
    #--------------------------------------#
    def iFFT_filter(self, x, y, xlo, xhi, threshold, qm, savefig, figname, dpi):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            iFFT_y, qm = signal_processing.iFFT_filter(reduced_x, reduced_y, threshold, qm, savefig, figname, dpi)
        else:
            iFFT_y = y
        return list(iFFT_y), qm
    

    #----------------------------------------------#
    # Method implementing a SavitzkyGolay smoother #
    #----------------------------------------------#
    def SavitzkyGolay(self, x, y, xlo, xhi, window, order):
        from scipy.signal import savgol_filter
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            xout = reduced_x.copy()
            yout = savgol_filter(np.array(reduced_y), window, order)
        else:
            xout = x
            yout = y
        return list(xout), list(yout)
    
    #--------------------------------------------------#
    # Method implementing a Whittaker-Eilers smoothing #
    #--------------------------------------------------#
    def whittaker_eilers_smoothing(self, x, y, xlo, xhi, order, lmbda, xlabel, ylabel, figname):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            smoothed, optimal_lambda, fig = whittaker_smoother.smooth(reduced_x, reduced_y, order, lmbda, xlabel=xlabel, ylabel=ylabel, basename=figname)
            if optimal_lambda is not None:
                self.log.out(f'    Optimized lambda for Whittaker-Eilers smoothing was found to be {optimal_lambda}.')
        else:
            reduced_x = x.copy(); smoothed = y.copy()
            self.log.GUI_error(f'ERROR (linear regression) no LAMMPS data in xrange {xlo} - {xhi}')
        return list(reduced_x), list(smoothed)

    #---------------------------------------------------------#
    # Method implementing a Y-data average in a given X-range #
    #---------------------------------------------------------#
    def average(self, x, y, xlo, xhi):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            average_y = misc_funcs.avg(reduced_y)
        else:
            average_y = 0
            self.log.GUI_error(f'ERROR (average) no LAMMPS data in xrange {xlo} - {xhi}')
        return reduced_x, reduced_y, average_y
    
    #------------------------------------------------------------------#
    # Method implementing a linear regression model in a given X-range #
    #------------------------------------------------------------------#
    def linear_regression(self, x, y, xlo, xhi, confidence_interval):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            linear_regression = misc_funcs.linear_regression(reduced_x, reduced_y)
            b1 = linear_regression.b1; b0 = linear_regression.b0; xreg = [xlo, xhi];
            r2 = linear_regression.r2
            yreg = [i * b1 + b0 for i in xreg]
            
            # Compute confidence interval if value is not zero
            if 0 < confidence_interval < 1:
                alpha = 1 - confidence_interval
                b0_ci = linear_regression.confidence_interval_b0(alpha)
                b1_ci = linear_regression.confidence_interval_b1(alpha)
                ci = [b0_ci, b1_ci]
            else:
                ci = None
        else:
            b0 = 0; b1 = 1; xreg = [0, 1]; yreg = [0, 1]; r2 = 0; ci = None
            self.log.GUI_error(f'ERROR (linear regression) no LAMMPS data in xrange {xlo} - {xhi}')
        return b0, b1, xreg, yreg, r2, ci
    
    #-----------------------------------------------------------------#
    # Method implementing a "piecewise regression" in a given X-range #
    #-----------------------------------------------------------------#
    def piecewise_regression(self, x, y, xlo, xhi, n):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            xout, yout, xbreaks, ybreaks, slopes = misc_funcs.piecewise_regression(reduced_x, reduced_y, xlo, xhi, n)
        else:
            xout = [0, 1]; yout = [0, 1]; slopes = {(0,1):1}
            xbreaks = [0, 1]; ybreaks = [0, 1];
            self.log.GUI_error(f'ERROR (peicewise-regression) no LAMMPS data in xrange {xlo} - {xhi}')
        return xout, yout, xbreaks, ybreaks, slopes
    
    #------------------------------------------#
    # Method implementing a spline-integration #
    #------------------------------------------#
    def spline_integration(self, x, y, xlo, xhi, shift):
        from scipy.interpolate import InterpolatedUnivariateSpline
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            xs = np.array(reduced_x)
            ys = np.array(reduced_y)
            if shift:
                shift_amount = ys[0]
                ys = ys - shift_amount               
            else: shift_amount = 0
            
            spl = InterpolatedUnivariateSpline(xs, ys, k=1)  # k=1 gives linear interpolation
            area = spl.integral(xlo, xhi)
            xout = list(xs)
            yout = list(spl(xs))
        else:
            xout = [0, 1]; yout = [0, 1]; area = 0;
            self.log.GUI_error(f'ERROR (spline-integration) no LAMMPS data in xrange {xlo} - {xhi}')
        return xout, yout, area, shift_amount
    
    #-------------------------------------------------------------------------------------#
    # Method implenting a hyperbola fit for Tg/CTE calculations from the following paper: #
    # Uncertainty quantification in molecular dynamics studies of the glass transition    #
    # temperature - Paul N. Patrone, Andrew Dienstfrey, ... - Polymer Volume 87 - 2016    #
    #-------------------------------------------------------------------------------------#
    def fit_hyperbola(self, x, y, xlo, xhi, slope_method, minimum_convergence=None, initial_guess=False, maxiter=10**4):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:        
            xout, yout, params, center, slopes, transition, tangent_intersection, tangents = misc_funcs.fit_hyperbola(reduced_x, reduced_y, xlo, xhi, minimum_convergence, initial_guess, maxiter, slope_method=slope_method)
        else:
            xout = [0, 1]; yout = [0, 1]; slopes = [0, 1]
            center = [0, 0]; params = [0, 0, 0, 0, 0];
            transition = []; tangent_intersection = [0, 0]
            tangents = [(0, 0, 0), (0, 0, 0)] # {(x-points), (ypoints)}
            self.log.GUI_error(f'ERROR no (hyperbola) LAMMPS data in xrange {xlo} - {xhi}')
        return xout, yout, params, center, slopes, transition, tangent_intersection, tangents
    
    #--------------------------------------#
    # Method implementing a polynomial fit #
    #--------------------------------------#
    def polynomial_fit(self, x, y, xlo, xhi, degree):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            # Convert to numpy arrays before passing to np_curve_fit
            reduced_x = np.array(reduced_x); reduced_y = np.array(reduced_y)
            xin, yin, xfit, yfit, param = misc_funcs.np_curve_fit(reduced_x, reduced_y, degree=degree, domain=None)
        else:
            xfit = x.copy(); yfit = y.copy()
            self.log.GUI_error(f'ERROR (linear regression) no LAMMPS data in xrange {xlo} - {xhi}')
        return list(xfit), list(yfit)
    
    #------------------------------------------------------------------------------#
    # Method implementing the Kemppainen-Muzzy stress-strain analysis calculations #
    #------------------------------------------------------------------------------#
    def rfr_mechanical(self, x, y, xlo, xhi, minxhi, maxxhi, xlo_method, yp, up, offset, t1, t2,stress_units, write_data, derivative_span, derivative_degree, grid, t_12_avg, savefig, figname, dpi):
        strain, stress = misc_funcs.reduce_data(x, y, xlo, xhi)
        t1_reduced = []; t2_reduced = []
        fig1, fig2 = None, None
        for trans1 in t1:
            xdummy, t1_tmp = misc_funcs.reduce_data(x, trans1, xlo, xhi)
            t1_reduced.append(t1_tmp)
        for trans2 in t2:
            xdummy, t2_tmp = misc_funcs.reduce_data(x, trans2, xlo, xhi)
            t2_reduced.append(t2_tmp)
        if strain and stress:
            out = rfr_mechanical.compute(strain, stress, minxhi, maxxhi, xlo_method, yp, up, offset, t1_reduced, t2_reduced, stress_units, write_data,
                                         derivative_span, derivative_degree, grid, t_12_avg, savefig, figname, dpi)
            xlo, xhi, yield_point_derivative, yield_point_offset, nu1, nu2, nu12, ultimate_point, fig1, fig2 = out
        else:
            xlo = min(x); xhi = max(x);  yield_point_derivative = []; yield_point_offset = []; nu1 = None; nu2 = None; ultimate_point = []
            self.log.GUI_error(f'ERROR (Regression Fringe Response Modulus) no LAMMPS data in xrange {xlo} - {xhi}')
        return xlo, xhi, yield_point_derivative, yield_point_offset, nu1, nu2, nu12, ultimate_point, fig1, fig2
    
    #--------------------------------------------------#
    # Method implementing the RFR thermal calculations #
    #--------------------------------------------------#
    def rfr_thermal(self, x, y, xlo, xhi, min_span, rtemp, settings, x_label, y_label, x_units, y_units, figname, dpi):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        fig = None
        if reduced_x and reduced_y:
            outputs, fig = rfr_thermal.compute(reduced_x, reduced_y, min_span, rtemp, settings, x_label, y_label, x_units, y_units, figname, dpi)
        else:
            outputs = {'lo_bounds': [min(x), max(x)],
                       'hi_bounds': [min(x), max(x)],
                       'lo_line': 3*[0],
                       'hi_line': 3*[0],
                       'tg': [0.0, 0.0],
                       'lo_reg_x': [min(x), max(x)], 
                       'lo_reg_y': [min(y), max(y)], 
                       'hi_reg_x': [min(x), max(x)], 
                       'hi_reg_y': [min(y), max(y)]}
        return outputs, fig
    
    #---------------------------------------------#
    # Method implementing derivative calculations #
    #---------------------------------------------#
    def differentiate_data(self, x, y, xlo, xhi, order, write_data, xlabel, ylabel, savefig, figname, dpi):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            dxn, dy1, dy2 = misc_funcs.compute_derivative(reduced_x, reduced_y)
            if '-zc' in str(order):
                critical_x, critical_y = misc_funcs.value_crossing(dxn, dy1, yvalue=0, style='low')
                inflection_x, inflection_y = misc_funcs.value_crossing(dxn, dy2, yvalue=0, style='low')
                crit_y = [reduced_y[reduced_x.index(x)] for x in critical_x]
                inflec_y = [reduced_y[reduced_x.index(x)] for x in inflection_x]
            else:
                critical_x = [];  critical_y = []; crit_y = []
                inflection_x = []; inflection_y = []; inflec_y = []
            if write_data:
                csv_data = {} # { 'name' : [[xdata], [ydata]], ... }
                csv_data['data'] = [reduced_x, reduced_y]
                csv_data['d1y/dx1'] = [dxn, dy1]
                csv_data['d2y/dx2'] = [dxn, dy2]
                csv_data['critical points'] = [critical_x, critical_y]
                csv_data['inflection points'] = [inflection_x, inflection_y]
                misc_funcs.savedata_to_csv(csv_data, figname+'.csv')
            
            # Set plot settings
            d1_color = 'tab:blue'    # 1st derivative color
            d2_color = 'tab:red'     # 2nd derivative color
            z1_color = 'tab:cyan'    # Zero crossing values for 1st derivatives (critical points)
            z2_color = 'tab:orange'  # Zero crossing values for 2nd derivatives (inflection points)
            fontsize = 12            # Font size for axis and labels
            ymin = min(reduced_y)    # Lower bound for veritical lines
            ymax = max(reduced_y)    # Lower bound for veritical lines
            
            if '1' in str(order) and '2' in str(order):
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
                ax3 = ax2.twinx()
                ax1.plot(reduced_x, reduced_y)
                ax2.plot(dxn, dy1, color=d1_color)
                ax3.plot(dxn, dy2, color=d2_color)
                if critical_x:
                    ax1.vlines(critical_x, ymin, ymax, color=z1_color, label='Critical Points')
                    ax2.plot(critical_x, critical_y, 'o', ms=8, color=z1_color, label='Critical Points')
                if inflection_x:
                    ax1.vlines(inflection_x, ymin, ymax, color=z2_color, label='Inflection Points')
                    ax3.plot(inflection_x, inflection_y, 'o', ms=8, color=z2_color, label='Inflection Points')

                
                #ax1.set_xlabel(xlabel, fontsize=fontsize)
                ax1.set_ylabel(ylabel, fontsize=fontsize)
                ax2.set_xlabel(xlabel, fontsize=fontsize)
                ax2.set_ylabel('$d^1y/dx^1$', fontsize=fontsize, color=d1_color)
                ax3.set_ylabel('$d^2y/dx^2$', fontsize=fontsize, color=d2_color)
                
                ax1.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
                ax2.legend(loc='lower left', bbox_to_anchor=(0, 0), fancybox=True, ncol=1, fontsize=8)
                ax3.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
                
                ax2.tick_params(axis='y', colors=d1_color)
                ax3.tick_params(axis='y', colors=d2_color)
                ax2.spines['left'].set_color(d2_color)
                ax3.spines['right'].set_color(d2_color)
                
                fig.tight_layout()
                if '0' not in savefig and '1' in savefig or 'all' in savefig:
                    fig.savefig(figname+'_1,2.jpeg', dpi=dpi)
            
            elif '1' in str(order):
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
                fontsize = 12
                
                ax1.plot(reduced_x, reduced_y)
                ax2.plot(dxn, dy1, color=d1_color)
                if critical_x:
                    ax1.vlines(critical_x, ymin, ymax, color=z1_color, label='Critical Points')
                    ax2.plot(critical_x, critical_y, 'o', ms=8, color=z1_color, label='Critical Points')
                
                #ax1.set_xlabel(xlabel, fontsize=fontsize)
                ax1.set_ylabel(ylabel, fontsize=fontsize)
                ax2.set_xlabel(xlabel, fontsize=fontsize)
                ax2.set_ylabel('$d^1y/dx^1$', fontsize=fontsize, color=d1_color)
                
                ax1.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
                ax2.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
                
                ax2.tick_params(axis='y', colors=d1_color)
                ax2.spines['left'].set_color(d1_color)
                
                fig.tight_layout()
                if '0' not in savefig and '1' in savefig or 'all' in savefig:
                    fig.savefig(figname+'_1.jpeg', dpi=dpi)
                
            elif '2' in str(order):
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
                fontsize = 12
                
                ax1.plot(reduced_x, reduced_y)
                ax2.plot(dxn, dy2, color=d2_color)
                if inflection_x:
                    ax1.vlines(inflection_x, ymin, ymax, color=z2_color, label='Inflection Points')
                    ax2.plot(inflection_x, inflection_y, 'o', ms=8, color=z2_color, label='Inflection Points')
                
                #ax1.set_xlabel(xlabel, fontsize=fontsize)
                ax1.set_ylabel(ylabel, fontsize=fontsize)
                ax2.set_xlabel(xlabel, fontsize=fontsize)
                ax2.set_ylabel('$d^2y/dx^2$', fontsize=fontsize, color=d2_color)
                
                ax1.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
                ax2.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
                
                ax2.tick_params(axis='y', colors=d2_color)
                ax2.spines['left'].set_color(d2_color)
                
                fig.tight_layout()
                if '0' not in savefig and '1' in savefig or 'all' in savefig:
                    fig.savefig(figname+'_2.jpeg', dpi=dpi)
            else:
                self.log.GUI_error(f'ERROR (Calculus: Differentiate Data) order can only be "1" or "2" or "1,2" or "2,1". Inputted order = {order}')
        else:
            self.log.GUI_error(f'ERROR (Calculus: Differentiate Data) no LAMMPS data in xrange {xlo} - {xhi}')
        return critical_x, crit_y, inflection_x, inflec_y
    
    #-------------------------------------------#
    # Method implementing integral calculations #
    #-------------------------------------------#
    def integrate_data(self, x, y, xlo, xhi, order, write_data, xlabel, ylabel, savefig, figname, dpi):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            ix1, iy1 = misc_funcs.compute_integral(reduced_x, reduced_y)
            ix2, iy2 = misc_funcs.compute_integral(ix1, iy1)
            if write_data:
                csv_data = {} # { 'name' : [[xdata], [ydata]], ... }
                csv_data['data'] = [reduced_x, reduced_y]
                csv_data['1st-integral'] = [ix1, iy1]
                csv_data['2nd-integral'] = [ix2, iy2]
                misc_funcs.savedata_to_csv(csv_data, figname+'.csv')
            if '1' in str(order) and '2' in str(order):
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
                ax2_color = 'tab:blue'
                ax3_color = 'tab:red'
                fontsize = 12
                
                ax3 = ax2.twinx()
                ax1.plot(reduced_x, reduced_y)
                ax2.plot(ix1, iy1, color=ax2_color)
                ax3.plot(ix2, iy2, color=ax3_color)
                
                #ax1.set_xlabel(xlabel, fontsize=fontsize)
                ax1.set_ylabel(ylabel, fontsize=fontsize)
                ax2.set_xlabel(xlabel, fontsize=fontsize)
                ax2.set_ylabel('1st order integral', fontsize=fontsize, color=ax2_color)
                ax3.set_ylabel('2nd order integral', fontsize=fontsize, color=ax3_color)
                
                ax2.tick_params(axis='y', colors=ax2_color)
                ax3.tick_params(axis='y', colors=ax3_color)
                ax2.spines['left'].set_color(ax2_color)
                ax3.spines['right'].set_color(ax3_color)
                
                fig.tight_layout()
                if '0' not in savefig and '1' in savefig or 'all' in savefig:
                    fig.savefig(figname+'_1,2.jpeg', dpi=dpi)
    
            elif '1' in str(order):
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
                ax2_color = 'tab:blue'
                fontsize = 12
                
                ax1.plot(reduced_x, reduced_y)
                ax2.plot(ix1, iy1, color=ax2_color)
                
                #ax1.set_xlabel(xlabel, fontsize=fontsize)
                ax1.set_ylabel(ylabel, fontsize=fontsize)
                ax2.set_xlabel(xlabel, fontsize=fontsize)
                ax2.set_ylabel('1st order integral', fontsize=fontsize, color=ax2_color)
                
                ax2.tick_params(axis='y', colors=ax2_color)
                ax2.spines['left'].set_color(ax2_color)
                
                fig.tight_layout()
                if '0' not in savefig and '1' in savefig or 'all' in savefig:
                    fig.savefig(figname+'_1.jpeg', dpi=dpi)
                
            elif '2' in str(order):
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
                ax2_color = 'tab:red'
                fontsize = 12
                
                ax1.plot(reduced_x, reduced_y)
                ax2.plot(ix2, iy2, color=ax2_color)
                
                #ax1.set_xlabel(xlabel, fontsize=fontsize)
                ax1.set_ylabel(ylabel, fontsize=fontsize)
                ax2.set_xlabel(xlabel, fontsize=fontsize)
                ax2.set_ylabel('2nd order integral', fontsize=fontsize, color=ax2_color)
                
                ax2.tick_params(axis='y', colors=ax2_color)
                ax2.spines['left'].set_color(ax2_color)
                
                fig.tight_layout()
                if '0' not in savefig and '1' in savefig or 'all' in savefig:
                    fig.savefig(figname+'_2.jpeg', dpi=dpi)
                    
                else:
                    self.log.GUI_error(f'ERROR (Calculus: Integrate Data) order can only be "1" or "2" or "1,2" or "2,1". Inputted order = {order}')
        else:
            self.log.GUI_error(f'ERROR (Calculus: Integrate Data) no LAMMPS data in xrange {xlo} - {xhi}')
        
        return