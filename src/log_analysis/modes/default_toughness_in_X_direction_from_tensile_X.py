# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 05/13/2024
Author: Josh Kemppainen
Purpose: Mode to find toughness when strained in
         X-dir
         using linear regressions


"""
# analysis list
analysis = [['LAMMPS data (apply moving average)', '', '', 'windw=100', ''],
            ['spline-integration', '', '', 'shift=True;', 'Toughness']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/property=tensile_modulus_x_strain_rate=2e8.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1',
        'xdata': 'v_etruex',
        'ydata': 'f_sxx_ave',
        'xlabel': 'True Strain',
        'ylabel': 'True Stress (MPa)',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }

