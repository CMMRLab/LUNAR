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
LUNAR/src/log_analysis/modes directory. See
get_started_3.py for enabling control over LUNAR's
logging class, to stop print outs and files to be
written and to stop the plot from automatically
being generated.
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
    import src.log_analysis.main as main
    os.chdir(pwd)
    
    
    # Load the mode from a combined path and update 'logfile' value
    path = os.path.join(path2lunar, modefile)
    mode = main.import_file(path).mode
    mode['logfile'] = logfile
    
    
    # Perform a log_analysis in a script
    analyzed = main.analysis(mode, savefig=savefig, dpi=dpi)
    
    
    # Get access to the logged results in Python data structures. Each analysis
    # method will have a name, with logged results in an "outputs" attribute and
    # information about the outputs in an "about" attribute. Both attributes are
    # Python dictionaries. To see all outputs and abouts, simply print the keys
    # of each.
    print('\n\n\n------------------------------------')
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