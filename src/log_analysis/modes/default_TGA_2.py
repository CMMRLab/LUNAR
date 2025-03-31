# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: XX/YY/ZZZZ
Author: John Doe
Purpose: Mode to find xxx

"""
# analysis list
analysis = [['LAMMPS data (fit polynomial)', '', '', 'deg=10', 'Polynomial'],
            ['Calculus: Differentiate Data', '', '', 'order=1,2-zc; csv=False; savefig=1', 'Derivatives']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/heating_trimer_ReaxFF_rep_1_del_2-50.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': 'all',
        'xdata': 'v_ramp_temp',
        'ydata': 'v_char_yield',
        'xlabel': 'Temperature (K)',
        'ylabel': 'Mass (%)',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }

