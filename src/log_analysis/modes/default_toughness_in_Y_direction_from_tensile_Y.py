# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 05/13/2024
Author: Josh Kemppainen
Purpose: Mode to find toughness when strained in
         Y-dir

"""
# analysis list
analysis = [['spline-integration', 0, 0.1, 'window=100; shift=True;', 'Toughness']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/property=tensile_modulus_y_strain_rate=2e8.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1',
        'xdata': 'v_etruey',
        'ydata': 'f_syy_ave',
        'xlabel': 'True Strain',
        'ylabel': 'True Stress (MPa)',
        'xscale': '',
        'yscale': '',
        'analysis': analysis
        }
