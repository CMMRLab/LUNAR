# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
October 23rd, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


This script implements an automatic method for calling the 
Kemppainen-Muzzy modes for shear simulations and computing
only the elastic responses. Which will automatically log
the results to a CSV file, set by the csvname variable.

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
path2lunar = 'C:/Users/jdkem/Desktop/LUNAR'


# Set the relative directory from this script, where LAMMPS logfiles are stored. 
files_directory = 'logfiles/Shear'


# Set relative paths from LUNAR to desired mode files for each direction of strain
xyfile = 'src/log_analysis/modes/default_RFR_Butterworth_XY_modulus.py'
xzfile = 'src/log_analysis/modes/default_RFR_Butterworth_XZ_modulus.py'
yzfile = 'src/log_analysis/modes/default_RFR_Butterworth_YZ_modulus.py'


# Set a "global" parent_dictory option. This is not required, but can be useful to keep inputs/outputs organized. Remember the following:
#  'logfile' will write outputs in the same directory of the analyzed file.
#  'logfile/PATH' will write outputs in a PATH relative to logfile. For example:
#                 'logfile/../' will write results "one directory up"
#                 'logfile/../NEW' will write results "one directory up" in a folder called NEW
#                 'logfile/NEW/TEST' will write results "two directories deeper" with relative path 'NEW/TEST'
parent_directory = 'logfile/elastic_test'


# Set the output csv file name (without the '.csv' extension) to automatically write logged results to.
csvname = 'Automatic_shear_elastic_testing'


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
    xymode = main.import_file(os.path.join(path2lunar, xyfile)).mode
    xzmode = main.import_file(os.path.join(path2lunar, xzfile)).mode
    yzmode = main.import_file(os.path.join(path2lunar, yzfile)).mode
    
    # Set path to logfiles and get all logfiles from the path that have the '.log.lammps' ending
    path = os.path.join(pwd, files_directory);
    logfiles = sorted([file for file in os.listdir(path) if file.endswith('.log.lammps')])
    
    # Setup a custom logging object
    log = io_functions.LUNAR_logger(level='production', print2console=False, write2log=True)
    
    # Setup logging dictionary to add values to for writing to .csv file later on
    logger = {'Filename': [],
              'Shear modulus': []}
    
    # Start automated analysis
    start_time = time.time()
    print('\n\nStarting automatic Kemppainen-Muzzy calculations ...')
    for n, logfile in enumerate(logfiles, 1):
        try:
            print('  Analyzing {:>4} of {:<4} File={}'.format(n, len(logfiles), logfile))
            
            # Set mode based on file naming (this will have to be updated for each
            # user, depending on the naming scheme that is used for their project).
            if   'shear_1' in logfile: mode = xymode
            elif 'shear_2' in logfile: mode = xzmode
            elif 'shear_3' in logfile: mode = yzmode
            else: raise Exception(f'ERROR can not determine mode based on file naming convention for file {logfile}.')
            
            # Update mode['logfile'] and mode['parent_directory'] for loaded mode
            mode['logfile'] = os.path.join(path, logfile)
            mode['parent_directory'] = parent_directory
            
            # Run analysis and log results
            analyzed = main.analysis(mode, plot=True, savefig=savefig, dpi=dpi, log=log)
            
            # Log desired results into logger (uncomment print statement to see all available keys)
            results = analyzed.outputs['Modulus'] # name of analysis 
            #print(results.keys())
            logger['Filename'].append(logfile)
            logger['Shear modulus'].append(results['b1-clean'])
            #print(analyzed.outputs['LAMMPS Butterworth Filter'])
        except:
            print('  FAILED {:>4} of {:<4} File={}'.format(n, len(logfiles), logfile))
    
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