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


# Find OS name
os_name = platform.system()


# The values in this dictionary can be set to None, at which point the GUI's will
# load with the default settings of Tkinter or set to user defined int/string to
# adjust default settings. Each OS-name can be set to different settings - to
# eliminate issues of working across different OS's (such as using Windows and WSL
# using the same LUNAR version).
if os_name == 'Windows':
    font_settings = {'size': None,           # None or int (12 is a good option)
                     'type': None,           # None or string ('Segoe UI' is a good option)
                     }
    
elif os_name == 'Linux':
    font_settings = {'size': None,           # None or int (12 is a good option)
                     'type': None,           # None or string ('Segoe UI' is a good option)
                     }
    
elif os_name == 'Darwin':
    font_settings = {'size': None,           # None or int (12 is a good option)
                     'type': None,           # None or string ('Segoe UI' is a good option)
                     }
    
else:
    font_settings = {'size': None,           # None or int (12 is a good option)
                     'type': None,           # None or string ('Segoe UI' is a good option)
                     }