# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 27th, 2024
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
    *        %matplotlib qt                                        *
    *     Then to turn back on the inline plotting, where the plots*
    *     will appear in the Plots tab execute the following       *
    *     command in the console:                                  *
    *        %matplotlib inline                                    *
    ****************************************************************
"""
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


##################################################################################################################
# Define basic settings, such as the 'mode' to load with, the predefined modes, and the figure image quality.    #
#                                                                                                                #
# Update settings as desired.                                                                                    #
##################################################################################################################
settings = {'mode': 'src/log_analysis/modes/default_blank.py', # Default run mode to load GUI with
            'modes-dir': 'src/log_analysis/modes',     # Set directory where all modes are saved
            'save-fig' : True,                         # True or False option to save figure are .jpeg
            'image-dpi': 300,                          # Set the image quality of the generated image (NOTE higher dpi's will result in slower code)
           }


###################################
### Import needed files and run ###
###################################
if __name__ == "__main__":  
    import src.log_analysis.GUI as GUI
    GUI.GUI(settings, GUI_zoom)