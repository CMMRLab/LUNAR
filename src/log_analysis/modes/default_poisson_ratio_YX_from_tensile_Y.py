# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 05/13/2024
Author: Josh Kemppainen
Purpose: Mode to find YX-poissions strained in
         Y-dir

"""
# analysis list
analysis = [['moving average', '', '', 'window=50', 'moving average'],
            ['linear regression', '', 0.02, '', '$\\nu_{yx}$  ']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/property=tensile_modulus_y_strain_rate=2e8.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1',
        'xdata': 'v_etruey',
        'ydata': 'v_etruex',
        'xlabel': 'True Strain (in Y)',
        'ylabel': 'True Strain (in X)',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }

