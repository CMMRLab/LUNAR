# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 05/13/2024
Author: Josh Kemppainen
Purpose: Blank mode to load GUI with

"""
# analysis list
analysis = []

# loadable mode
mode = {'logfile': 'UPDATE-ME',
        'keywords': ['Step', 'Temp'],
        'sections': 'all',
        'xdata': '',
        'ydata': '',
        'xlabel': '',
        'ylabel': '',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }

