# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 20, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


    ****************************************************************
    * Requirements:                                                *
    *   python 3.7+                                                *
    *                                                              *
    * Dependencies:                                                *
    *   python matplotlib module:                                  *
    *    - pip3 install matplotlib (if pip manager is installed)   *
    *                                                              *
    * Run methods:                                                 *
    *   - GUI (manipulate variables and run from. Default          *
    *          settings set from this script)                      *
    *                                                              *
    * Notes for Anaconda Spyder IDE users:                         *
    *   - If running this from the Anaconda Spyder IDE, before     *
    *     running you will have to turn on the interactive plot    *
    *     if you want the interactive plot to pop-up. You may do   *
    *     this by executing the following command in the console:  *
    *           %matplotlib qt                                     *
    *     Then to turn back on the inline plotting, where the plots*
    *     will appear in the Plots tab execute the following       *
    *     command in the console:                                  *
    *        %matplotlib inline                                    *
    ****************************************************************
"""

################
# GUI controls #
##################################################################################################################
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account for the  #
# different ways individuals have their OS display settings setup. A GUI_zoom = 100, means that default settings #
# are used, whereas a GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80 means decrease GUI by 20%. #
# Examples:                                                                                                      #
#   GUI_zoom = 100 # use default GUI size                                                                        #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                                    #
#                                                                                                                #
# Update use_GUI and GUI_zoom as desired.                                                                        #
##################################################################################################################
GUI_zoom = 100
use_GUI = True


##################################################################################################################
# Define basic settings, such as the 'mode' to load with, the predefined modes, and the figure image quality.    #
#                                                                                                                #
# Update settings as desired.                                                                                    #
##################################################################################################################
settings = {'mode': 'src/atom_typing/log_processor/modes/default_pyro_CH.py', # Default run mode to load GUI with.
            'modes-dir': 'src/atom_typing/log_processor/modes',         # Set directory where all modes are saved.
            'replace_logfile_when_loading_mode': False,                 # True or False default for replacing logfile when loading a mode.
            }

# Josh's quick way for changing default load mode
#settings['mode'] = 'src/atom_typing/log_processor/modes/test_2.py'
settings['mode'] = 'src/atom_typing/log_processor/modes/test_9_multiple_wildcards.py'


################
# IDE controls #
##################################################################################################################
# These IDE controls are mainly to demonstrate how to use this code in an IDE setting, particularly for the use  #
# of generating a script for personally processing. If you look at the main function call when not using the GUI #
# you will see how main() can be called and what the return values are:                                          # 
#     header, values, data = main(mode)                                                                          # 
#                                                                                                                #
# If you want to use main() in a script, please use the "Quick help" page accessible from the GUI to see how to  #
# setup a mode dictionary and look at the modes in 'LUNAR/src/atom_typing/log_processor/modes/' directory to see #
# examples. The header, values, data return values are briefly described below:                                  #
#   header = [list of logger keywords]                                                                           #
#   values = [[File1: value that maps onto index in headers list], [File2], [FileN] ]                            #
#   data = {keyword : [list of values ordered from how files where iterated through]}                            #
#                                                                                                                #
# Update settings as desired.                                                                                    #
##################################################################################################################
# Define logger in a doc string
logger = """
filename,    numeric,    wildcards,    rings[6]['%Ring'],    rings[6]['C']['%natoms'],    hybridizations['Sp1-C']['%natoms']
hybridizations['Sp2-C']['%natoms'],    hybridizations['Sp3-C']['%natoms'],    hybridizations['all-H']['%natoms']
"""

# Define mode for IDE controls
mode = {'independent_variable': 'numeric',  # if this is an empty string a plot will not appear
        'dependent_variables': 'logger',    # if this is an empty string a plot will not appear
        'parent_directory': 'logfile',
        'sorting_direction': 'ascending',
        'sorting_method': 'sort',
        'logfile': 'EXAMPLES/array_processing/atom_typing_logfile_processor/poly_tracking_replicate_1_time_*ps_typed.log.lunar',
        'newfile': 'atom_typing_logfile_processor',
        'logger': logger}



###################################
### Import needed files and run ###
###################################
if __name__ == "__main__":  
    from src.atom_typing.log_processor.main import main
    from src.atom_typing.log_processor.GUI import GUI
    
    if use_GUI:
        GUI(settings, GUI_zoom)
    else:
        header, values, data = main(mode)

