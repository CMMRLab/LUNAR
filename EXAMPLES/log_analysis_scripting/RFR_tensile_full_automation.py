# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
December 11, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


This script implements an automatic method for calling the 
Kemppainen-Muzzy modes for tensile simulations. Which will
automatically log the results to a CSV file, set by the
csvname variable.

If working from Anaconda's Spyder IDE, you can:
    - Enable plots to pop-up via: %matplotlib qt
    - Stop plots from popping up via: %matplotlib inline
"""
#########################################
### Import Necessary Python Libraries ###
#########################################
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
files_directory = 'logfiles/Tensile'


# Set relative paths from LUNAR to desired mode files for each direction of strain
xfile = 'src/log_analysis/modes/default_RFR_Butterworth_X_modulus.py'
yfile = 'src/log_analysis/modes/default_RFR_Butterworth_Y_modulus.py'
zfile = 'src/log_analysis/modes/default_RFR_Butterworth_Z_modulus.py'


# Set a "global" parent_dictory option. This is not required, but can be useful to keep inputs/outputs organized. Remember the following:
#  'logfile' will write outputs in the same directory of the analyzed file.
#  'logfile/PATH' will write outputs in a PATH relative to logfile. For example:
#                 'logfile/../' will write results "one directory up"
#                 'logfile/../NEW' will write results "one directory up" in a folder called NEW
#                 'logfile/NEW/TEST' will write results "two directories deeper" with relative path 'NEW/TEST'
parent_directory = 'logfile/test'


# Set the output csv file name (without the '.csv' extension) to automatically write logged results to.
csvname = 'Automatic_tensile_testing'


# Plotting options for how log_analysis handles things internally
savefig = True  # True or False to write figure to a file.
dpi = 300       # Integer to set dots per inch of saved file (a good option is around 300 DPI).




#########################################################
### Start analysis using LUNAR/log_analysis scripting ###
#########################################################
if __name__ == "__main__": 
    # Get present working directory
    pwd = os.getcwd()
    
    # Move to LUNAR, import what is needed and then move back to pwd
    os.chdir(path2lunar)
    import src.io_functions as io_functions
    import src.log_analysis.main as main
    os.chdir(pwd)
    
    # Load the modes from a combined path
    xmode = main.import_file(os.path.join(path2lunar, xfile)).mode
    ymode = main.import_file(os.path.join(path2lunar, yfile)).mode
    zmode = main.import_file(os.path.join(path2lunar, zfile)).mode
    
    # Set path to logfiles and get all logfiles from the path that have the '.log.lammps' ending
    path = os.path.join(pwd, files_directory);
    logfiles = sorted([file for file in os.listdir(path) if file.endswith('.log.lammps')])
    
    # Setup a custom logging object
    log = io_functions.LUNAR_logger(level='production', print2console=False, write2log=True)
    
    # Setup logging dictionary to add values to for writing to .csv file later on
    logger = {'Filename': [],
              'Elastic modulus': [],
              'Yield Strength': [],
              'Poissons ratio - nu1': [],
              'Poissons ratio - nu2': [],
              'Poissons ratio - nu_avg': []}
    
    # Start automated analysis
    start_time = time.time()
    print('\n\nStarting automatic RFR calculations ...')
    for n, logfile in enumerate(logfiles, 1):
        successful = False
        try:
            print('  Analyzing {:>4} of {:<4} File={}'.format(n, len(logfiles), logfile))
            
            # Set mode based on file naming (this will have to be updated for each
            # user, depending on the naming scheme that is used for their project).
            if   'tensile_1' in logfile: mode = xmode
            elif 'tensile_2' in logfile: mode = ymode
            elif 'tensile_3' in logfile: mode = zmode
            else: raise Exception(f'ERROR can not determine mode based on file naming convention for file {logfile}.')
            
            # Update mode['logfile'] and mode['parent_directory'] for loaded mode
            mode['logfile'] = os.path.join(path, logfile)
            mode['parent_directory'] = parent_directory
            
            # Run analysis and log results
            analyzed = main.analysis(mode, plot=True, savefig=savefig, dpi=dpi, log=log)
            #print( list(analyzed.outputs.keys()) )
            
            # Access outputs from log analysis
            results = analyzed.outputs['RFR-mechanical'] # name of analysis 
    
            # Log desired results into logger
            modulus = results['b1']
            strength = results['yield_point_derivative'][1]
            poissons_1 = results['nu1']
            poissons_2 = results['nu2']
            poissons_avg = results['nu_avg']
            successful = True
        except:
            print('  FAILED {:>4} of {:<4} File={}'.format(n, len(logfiles), logfile))
            modulus = None
            strength = None
            poissons_1 = None
            poissons_2 = None
            poissons_avg = None
            
        # Log desired results into logger
        print('\n\n')
        if successful:
            logger['Filename'].append(logfile)
            logger['Elastic modulus'].append(modulus)
            logger['Yield Strength'].append(strength)
            logger['Poissons ratio - nu1'].append(poissons_1)
            logger['Poissons ratio - nu2'].append(poissons_2)
            logger['Poissons ratio - nu_avg'].append(poissons_avg)
        
    # Invert data and store in a matrix to write to csv file
    ncolumns = len(logger); nrows = max(map(len, list(logger.values())))
    matrix = [[str(None)]*ncolumns for n in range(nrows)]; titles = list(logger.keys())
    for i in range(nrows):
        for j, name in enumerate(titles):
            matrix[i][j] = str(logger[name][i])
    
    # Write data to csv
    print(f'\n\nWriting {csvname}')
    with open(csvname+'.csv', 'w') as f:        
        f.write('{}\n'.format(', '.join(titles)))
        for row in matrix:
            f.write('{}\n'.format(', '.join(row)))
    
    # Print the total analysis execution time
    execution_time = (time.time() - start_time)
    print('\n\nExecution time in seconds: ' + str(execution_time))