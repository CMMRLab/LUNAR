# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
March 1st, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


This script implements an automatic method for calling the 
RFR modes for tensile simulations. Which will automatically
log the results to a CSV file, set by the basename variable.

If working from Anaconda's Spyder IDE, you can:
    - Enable plots to pop-up via: %matplotlib qt
    - Stop plots from popping up via: %matplotlib inline
"""
#########################################
### Import Necessary Python Libraries ###
#########################################
import math
import glob
import time
import os

##############
### Inputs ###
##############
# Set path of LUNAR folder. For example if full path:
# fullpath = 'C:/Users/USER/Desktop/LUNAR'
# path2lunar = 'C:/Users/USER/Desktop/LUNAR'
path2lunar = 'C:/Users/jdkem/Desktop/LUNAR'


# Set the relative directory from this script, where LAMMPS logfiles are stored. 
logfile = '../../**EXAMPLES/log_analysis/RFR_array/*/tensile_*_PBZ_pxld_*_replicate_*_FF_PCFF.log.lammps'
logfile = 'logfiles/Tensile/*.log.lammps'


# Set the direction specific settings for log_analysis RFR method
directions = {'x':{'sections': '1',
                   'xdata': 'v_etruex',
                   'ydata': 'f_sxx_ave',
                   't1': 'v_etruey',
                   't2': 'v_etruez'},
              
              'y':{'sections': '1',
                   'xdata': 'v_etruey',
                   'ydata': 'f_syy_ave',
                   't1': 'v_etruex',
                   't2': 'v_etruez'},
              
              'z':{'sections': '1',
                   'xdata': 'v_etruez',
                   'ydata': 'f_szz_ave',
                   't1': 'v_etruex',
                   't2': 'v_etruey'}}

# Set the output csv file name (without the '.csv' extension) to automatically write logged results to.
basename = 'RFR_default_array'


# Plotting options for how log_analysis handles things internally
savefig = False  # True or False to write figure to a file.
dpi = 300        # Integer to set dots per inch of saved file (a good option is around 300 DPI).


#########################################################
### Start analysis using LUNAR/log_analysis scripting ###
#########################################################
if __name__ == "__main__": 
    # Get present working directory
    pwd = os.getcwd()
    
    # Move to LUNAR, import what is needed and then move back to pwd
    os.chdir(path2lunar)
    import src.glob_wildcards as glob_wildcards
    import src.io_functions as io_functions
    import src.log_analysis.main as main
    os.chdir(pwd)
    
    # Set up logging dictionary
    start_time = time.time()
    logger = {'Filename': [],
              'Yield Shift': [],
              'Yield Strain': [],
              'Yield Strength': [],
              'Elastic modulus': [],
              'xlo': [],
              'xhi': [],
              'nu1': [],
              'nu2': [],
              'nu_avg': [],
              'wildcards': []}
    
    # Loop through files and start analysis
    files = sorted(glob.glob(logfile))
    for n, file in enumerate(files, 1):
        rootname = os.path.basename(file)
        wildcards = glob_wildcards.get_glob_wildcards(logfile, file)
        
        # Set direction information
        if 'tensile_1' in rootname:
            direction = directions['x']
        if 'tensile_2' in rootname:
            direction = directions['y']
        if 'tensile_3' in rootname:
            direction = directions['z']
            
        # Construct the analysis list for log_analysis
        t1 = direction['t1'] 
        t2 = direction['t2'] 
        analysis = [['LAMMPS data (apply Butterworth filter)', '', '', 'qm=msr; psd=False; csv=False; savefig=all; order=2; wn=op', 'Butterworth'],
                    ['Regression Fringe Response Modulus', '', '', 
                     f'shift=ymin; minxhi=0.0025; maxxhi=0.0; xlo=rfs; yp=1; offset=0.0; t1={t1}; t2={t2}; csv=False; savefig=1', 'RFR']]
        
        # Construct the log_analysis mode dictionary
        mode = {}
        mode['parent_directory'] = 'logfile/v7'
        mode['logfile'] = file
        mode['keywords'] = ['Step']
        mode['sections'] = direction['sections']
        mode['xdata'] = direction['xdata']
        mode['ydata'] = direction['ydata']
        mode['nevery'] = 1
        mode['xcompute'] = None
        mode['ycompute'] = None
        mode['xlabel'] = 'Strain'
        mode['ylabel'] = 'Stress (MPa)'
        mode['analysis'] = analysis
        
        # Analylze files
        try:
            log = io_functions.LUNAR_logger(level='production', print2console=False, write2log=False)
            analysis = main.analysis(mode, plot=True, savefig=savefig, dpi=dpi, log=log, log_clear=True)
            
            # Get Butterworth RFR outputs
            butterworth = analysis.outputs['Butterworth']
            RFR = analysis.outputs['RFR']

            
            # Log results
            logger['Filename'].append(rootname)
            logger['Yield Shift'].append(RFR['b0-clean'])
            logger['Yield Strain'].append(RFR['yield_point_derivative'][0])
            logger['Yield Strength'].append(RFR['yield_point_derivative'][1])
            logger['Elastic modulus'].append(RFR['b1-clean'])
            logger['xlo'].append(RFR['xlo'])
            logger['xhi'].append(RFR['xhi'])
            logger['nu1'].append(RFR['nu1'])
            logger['nu2'].append(RFR['nu2'])
            logger['nu_avg'].append(RFR['nu_avg'])
            logger['wildcards'].append(' '.join([str(i) for i in wildcards]))
            
            print(f'Completed {n} of {len(files)}:', file)
            print(f'Completed {n} of {len(files)}:', file)
        except: print('Failed:', file)
        
    # Go through and find averages and standard deviations for each logged value
    statistics = {} # {'logger-key':(mean, std, min, max)}
    for key in logger:
        values = logger[key]
        if all(isinstance(x, (int, float)) for x in values) and values:
            mean = sum(values)/len(values)
            deviations = [(x - mean)**2 for x in values]
            variance = sum(deviations)/len(deviations)
            std = math.sqrt(variance)
            min_value = min(values)
            max_value = max(values)
        else:
            mean = None
            std = None
            min_value = None
            max_value = None
        statistics[key] = (mean, std, min_value, max_value)
        
    # Invert data and store in a matrix to write to csv file
    ncolumns = len(logger); nrows = max(map(len, list(logger.values())))
    matrix = [[str(None)]*ncolumns for n in range(nrows)]; titles = list(logger.keys())
    for i in range(nrows):
        for j, name in enumerate(titles):
            matrix[i][j] = str(logger[name][i])
            
    # Add statistics to matrix
    means = []; stds = []; spaces = []
    for key in statistics:
        mean, std, min_value, max_value = statistics[key]
        space = ''.join(['-' for _ in range(len(key))])
        spaces.append(space)
        if mean is not None:
            means.append(str(mean))
            stds.append(str(std))
        else:
            means.append('Mean')
            stds.append('STD')
    matrix.append(spaces)
    matrix.append(means)
    matrix.append(stds)
            
    # Write data to csv
    try:
        print(f'\n\nWriting {basename}')
        with open(basename+'.csv', 'w') as f:        
            f.write('{}\n'.format(', '.join(titles)))
            for row in matrix:
                f.write('{}\n'.format(', '.join(row)))
    except: print('ERROR could not write csv file. Likely opened by another software.')
        
    # Write out the statistics
    print('\n\nStatistics:')
    for key in statistics:
        mean, std, min_value, max_value = statistics[key]
        if mean is not None:
            print(' - {}'.format(key))
            print('    mean: {}'.format(mean))
            print('    std: {}'.format(std))
            print('    min: {}'.format(min_value))
            print('    max: {}\n'.format(max_value))
    
    # Wrap up script run
    execution_time = (time.time() - start_time)
    print()
    print('Execution time in seconds: ' + str(execution_time))
    