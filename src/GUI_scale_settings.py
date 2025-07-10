# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
June 1, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import platform


# Adjust the below screen_width variable to fill up the "correct" amount of the screen width:
#   1080p is 1920x1080 (width = 1920) 
#   4k    is 3840x2160 (width = 3840)
screen_width = 1920 #3840


# Function adjust screen scale (development resolution = 1080p, (1920px wide))
def get_screen_ratio(screen_width):
    screen_ratio = (screen_width/1920)
    return screen_ratio


# The values in this dictionary can be set to None, at which point the GUI's will load with the default
# settings of Tkinter or set to user defined int/string to adjust default settings. Each OS-name can be
# setting to different settings - to eliminate issues of working across different OS's (such as using 
# Windows and WSL using the same LUNAR version).
os_name = platform.system()

#------------------#
# Windows Specific #
#------------------#
if os_name == 'Windows':
    font_size = 10
    # default settings
    font_settings   = {'size': font_size, # None or int or scaled_font_size (12 is a good option)
                       'type': None,      # None or string ('Segoe UI' is a good option)
                       }
    
    # 4k settings
    if screen_width > 1920:
        ratio = get_screen_ratio(screen_width)
        screen_settings = {'scaling_factor' : ratio}
    else: # 1080 settings
        screen_settings = {'scaling_factor' : None}
    
    
#----------------#
# Linux Specific #
#----------------#  
elif os_name == 'Linux':
    font_size = 10
    # default settings
    font_settings   = {'size'             : font_size,   # None or int or scaled_font_size (12 is a good option)
                       'type'             : None,        # None or string ('Segoe UI' is a good option)
                       }
    
    # 4k settings
    if screen_width > 1920:
        font_size = 9
        ratio = get_screen_ratio(screen_width)
        scaled_font_size = int(font_size*ratio)
        font_settings   = {'size'           : scaled_font_size,  # None or int or scaled_font_size (12 is a good option)
                           'type'           : None,              # None or string ('Segoe UI' is a good option)
                           'dialog_size'    : scaled_font_size,  # None or int or scaled_font_size (12 is a good option) 
                           'dialog_type'    : None,              # None or string ('Segoe UI' is a good option) 
                           }
        
        screen_settings = {'scaling_factor' : ratio}
    else: # 1080 settings
        screen_settings = {'scaling_factor' : None}
    

#----------------#
# MacOS Specific #
#----------------#
elif os_name == 'Darwin':
    font_settings   = {'size'             : None,    # None or int or scaled_font_size (12 is a good option)
                       'type'             : None,    # None or string ('Segoe UI' is a good option)
                       }
    
    screen_settings = {'scaling_factor'   : None}
    
    
else:
    font_settings   = {'size'             : None,   # None or int or scaled_font_size (12 is a good option)
                       'type'             : None,   # None or string ('Segoe UI' is a good option)
                       }
    
    screen_settings = {'scaling_factor'   : None}