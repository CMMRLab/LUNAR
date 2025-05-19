# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
October 21st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


This script will show the basics of how you can script
and run log_analysis to get desired values. These values
could be automatically logged to a file based on how the
scripter wants to write methods to dumping data to files.

This script will load already built modes from the 
LUNAR/src/log_analysis/modes directory. Then the script
will show how to customize LUNAR's logging class to stop
the printing and writing of files, as it is likely desired
that if you are script, you want to only see your printouts.

Additionally, it will be shown how to pass the "plot" 
parameter as False, to stop log_analysis from plotting
as if you are scripting, you will likely want control
over not seeing the plot for automation purposes or
to generate your own plots. Please see get_started_4.py
for seeing how to plot data the log_analysis is returning
from the outputs attribute.
"""
#########################################
### Import Necessary Python Libraries ###
#########################################
import os


##############
### Inputs ###
##############
# Set path of LUNAR folder. For example if full path:
# fullpath = 'C:/Users/USER/Desktop/LUNAR'
# path2lunar = 'C:/Users/USER/Desktop/LUNAR'
path2lunar = 'C:/Users/jdkem/Desktop/LUNAR'


# In "get_started_1.py", the mode dictionary was "hand built". You can also just open
# a mode file from within LUNAR by setting its relative path from LUNAR. *NOTE: this
# will load everything from that mode file, including the logfile, so we will have to
# update the logfile to the one we want to analyze
modefile = 'src/log_analysis/modes/default_density_average.py'


# Set the logfile we want to load and analyze. We will manually have to adjust from the
# loaded mode.
logfile = 'logfiles/property=density_ts=0.5.log.lammps'


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
    
    
    # Load the mode from a combined path and update 'logfile' value
    path = os.path.join(path2lunar, modefile)
    mode = main.import_file(path).mode
    mode['logfile'] = logfile
    
    
    # Setup a custom logging object to pass into main.analysis(). The following procedure meanings
    #   1.) Generate an logging object by setting all parameters explicitly. Meaning of each:
    #        - level=<string>, where the string can be 'production' or 'debug'. The 'production'
    #                          option only prints and writes to the log file the develop of LUNAR
    #                          deems is important for users to see. The 'debug' option prints and 
    #                          writes everything the 'production' option does, but additional 
    #                          information, the developer of LUNAR uses for debugging purposes.
    #        - write2log=<Boolean>, where boolean is True to write all printed information to a
    #                               logfile and False is to not write all printed information to 
    #                               a log file.
    #        - print2console=<Boolean>, where boolean is True to print things to the console and
    #                                   False to not print things to the console.
    #
    #   2.) Configure the log class to ones desires. *NOTE: if you setup the log object how you
    #       liked in the first place, this step is not needed. However understanding how to re-
    #       configure the logger can be useful for quickly adjusting the logger in different use
    #       cases.
    #
    # For this example, we will intialize the logging object with some settings and then change
    # them to show the configure method.
    log = io_functions.LUNAR_logger(level='production', print2console=True, write2log=True)
    log.configure(level='production', print2console=False, write2log=False)
    
    
    # Perform a log_analysis in a script, passing in the log object as a parameter and
    # setting the plot parameter as False. The plot=<Boolean>, where True plots the data
    # and False does not plot the data. The default if not provided is to plot the data.
    analyzed = main.analysis(mode, plot=False, savefig=savefig, dpi=dpi, log=log)
    
    
    # Get access to the logged results in Python data structures. Each analysis
    # method will have a name, with logged results in an "outputs" attribute and
    # information about the outputs in an "about" attribute. Both attributes are
    # Python dictionaries. To see all outputs and abouts, simply print the keys
    # of each.
    print('------------------------------------')
    print('LUNAR log_analysis scripting methods')
    print('------------------------------------')
    print('Available outputs: ', list(analyzed.outputs.keys()))
    print('Available about: ', list(analyzed.about.keys()))
    
    
    # Understanding how to "dive" deeper into the structures
    print('\n\nAvailable per analysis data lists, floats, or ints that can be logged:')
    for name in analyzed.outputs:
        print(name, list(analyzed.outputs[name].keys()))
        
    
    # Lets get the final average density now
    name = 'Density average'
    average_density = analyzed.outputs[name]['average']
    about_average_density = analyzed.about[name]['average']
    print('\n\nThe final averaged density is: ', average_density)
    print('The about information on this value is: ', about_average_density)