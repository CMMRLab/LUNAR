# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
December 11, 2025
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
import matplotlib.pyplot as plt
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
#
# This example is using a rel-path so this file can be instantly
# called from inside LUNAR/EXAMPLES/log_analysis scripting, but
# if abs-path is provide, this file can be placed anywhere on 
# your C-drive and LUNAR can still be found
path2lunar = '../../'


# Set the relative directory from this script, where LAMMPS logfiles are stored. 
logfile = 'logfiles/bmi_wideTg/c44M-*_*_TgAnneal_BMI_3comp.log'
logfile = 'logfiles/pbz_tg/PBZ_*_R*.log.lammps'


# Set the output csv file name (without the '.csv' extension) to automatically write logged results to.
basename = 'RFR_thermal_default_array'


# Plotting options for how log_analysis handles things internally
savefig = True   # True or False to write figure to a file.
dpi = 300        # Integer to set dots per inch of saved file (a good option is around 300 DPI).


# Material/LAMMPS sim specific logfile sections
log_section = 'all'
if 'PBZ' in logfile:
    log_section = '-1' #'1'
if 'BMI' in logfile:
    log_section = '-1' #'3'


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
    logger = {'filename': [],
              'lo-slope': [],
              'hi-slope': [],
              'lo-CTE': [],
              'hi-CTE': [],
              'tg': [],
              'room-temp-vol':[],
              'wildcards': []}
    
    # Loop through files and start analysis
    files = sorted(glob.glob(logfile))
    for n, file in enumerate(files, 1):
        rootname = os.path.basename(file)
        wildcards = glob_wildcards.get_glob_wildcards(logfile, file)
        if n > 5: break

        # Construct the analysis list for log_analysis
        analysis = [['LAMMPS data (X-sort)', '', '', '', 'Keeps X-data increasing'],
                    ['LAMMPS data (apply Whittaker-Eilers)', '', '', 'order=2; lambda=op<1e-2, 1e12, 10>', 'LAMMPS Whittaker-Eilers'],
                    ['average', 290, 310, '', 'room temp volume average (A$^3$)'],
                    ['skip', '', '', 'p=0.9; initial_guess=True', 'Hyperbola fit'],
                    ['Regression Fringe Response Thermal', '', '', 'minspan=25; rtemp=290,310', 'RFR-thermal']]
        
        # Construct the log_analysis mode dictionary
        mode = {}
        mode['parent_directory'] = 'logfile/v0'
        mode['logfile'] = file
        mode['keywords'] = ['Step', 'Temp', 'Volume', 'v_Time']
        mode['sections'] = log_section
        mode['xdata'] = 'Temp'
        mode['ydata'] = 'Volume'
        mode['nevery'] = 1
        mode['xcompute'] = None
        mode['ycompute'] = None
        mode['xlabel'] = 'Temperature (K)'
        mode['ylabel'] = 'Volume (A$^3$)'
        mode['analysis'] = analysis
        
        # Analylze files
        successful = False
        try:
            print('  Analyzing {:>4} of {:<4} File={}'.format(n, len(files), file))
            log = io_functions.LUNAR_logger(level='production', print2console=False, write2log=True)
            analysis = main.analysis(mode, plot=True, savefig=savefig, dpi=dpi, log=log, log_clear=True)
            
            # Get Butterworth RFR outputs
            RFR = analysis.outputs['RFR-thermal']
            room_temp_avg = analysis.outputs['room temp volume average (A$^3$)']
    
            # Log results
            room_temp_vol = room_temp_avg['average']
            lo_slope = RFR['lo_line'][1]
            hi_slope = RFR['hi_line'][1]
            lo_CTE = RFR['lo_cte']
            hi_CTE = RFR['hi_cte']
            tg = RFR['tg'][0]
            successful = True
            print(f'Completed {n} of {len(files)}:', file)
        except: 
            print('  FAILED {:>4} of {:<4} File={}'.format(n, len(files), file))
            room_temp_vol = None
            lo_slope = None
            hi_slope = None
            lo_CTE = None
            hi_CTE = None
            tg = None
            
        # Log desired results into logger
        print('\n\n')
        if successful:
            logger['filename'].append( rootname )
            logger['lo-slope'].append( lo_slope )
            logger['hi-slope'].append( hi_slope )
            logger['lo-CTE'].append( lo_CTE )
            logger['hi-CTE'].append( hi_CTE )
            logger['tg'].append( tg )
            logger['room-temp-vol'].append( room_temp_vol )
            logger['wildcards'].append(' '.join([str(i) for i in wildcards]))
        
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
    
    
    # Optional plot (adjust as needed)
    conversion = [float(i.split()[0]) for i in logger['wildcards']] # convert strings to floats (index-0 is conversion)
    tgs = logger['tg']
    cte_lo = logger['lo-CTE']
    cte_hi = logger['hi-CTE']
    
    # Tristan specific madness
    if 'PBZ' in logfile:
        cte_lo = [CTE*10**6/3*(1/spec_volume) for CTE, spec_volume in zip(logger['lo-slope'], logger['room-temp-vol'])]
        cte_hi = [CTE*10**6/3*(1/spec_volume) for CTE, spec_volume in zip(logger['hi-slope'], logger['room-temp-vol'])]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    
    ax1.plot(conversion, tgs, '*', ms=12, color='tab:blue')
    ax1.set_xlabel('Conversion (%)')
    ax1.set_ylabel('Tg (K)')
    
    
    ax2.plot(conversion, cte_lo,'*', ms=12, color='tab:blue', label='Lo-CTE')
    ax2.plot(conversion, cte_hi, '*', ms=12, color='tab:red', label='Hi-CTE')
    ax2.set_xlabel('Conversion (%)')
    ax2.set_ylabel('CTE')
    ax2.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=2)
    
    fig.tight_layout()
    fig.savefig(basename+'.jpeg', dpi=dpi)
    plt.show()
    

    