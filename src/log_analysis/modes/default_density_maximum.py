# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 05/13/2024
Author: Josh Kemppainen
Purpose: Mode to find largest density in a given
         X-range of data

"""
# analysis list
analysis = [['LAMMPS data (apply moving average)', '', '', 'windw=100', ''],
            ['maximum', '', '', '', 'Density maximum']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/property=density_ts=0.5.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1,2',
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

