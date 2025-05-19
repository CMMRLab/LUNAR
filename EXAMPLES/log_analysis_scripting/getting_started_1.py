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

This script will "hand build" a mode dictionary that
log_analysis uses for storing settings. See get_started_2.py
for loading already built mode files.
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
# path2lunar = 'C:/Users/USER/Desktop'
path2lunar = 'C:/Users/jdkem/Desktop'


# Build a log_analysis mode dictionary. Somethings will be documented hear, but
# once you get used to these mode dictionaries you can likely generate them without
# proper documentation or you can "build" your analysis in LUNAR/log_analysis GUI
# and click "save settings as mode", which will generate the mode dictionary for 
# you from GUI settings.
#
# Define the analysis list (NOTE the add more analysis, just have many sublists)
#             method        X-lo            X-hi          Misc         Name
#            <string>, <int or float>, <int or float>,  <string>,    <string>
analysis = [['average',    1500,           2100,         '',       'Density average']]

# Define the mode dictionary
mode = {'logfile': 'logfiles/property=density_ts=0.5.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1-2',
        'xdata': 'Step',
        'ydata': 'Density',
        'xlabel': 'Time (ps)',
        'ylabel': 'Density ($g/cm^3$)',
        'xcompute': '${Step}*(1/2000)',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }


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
    os.chdir(os.path.join(path2lunar, 'LUNAR'))
    import src.log_analysis.main as main
    os.chdir(pwd)
    
    
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