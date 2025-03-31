# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 05/13/2024
Author: Josh Kemppainen
Purpose: Mode to find v0 and v1 for bulk modulus
         calulation

"""
# analysis list
analysis = [['average', 0, 400, '', 'v0 average (A$^3$)'],
            ['average', 1000, 2100, '', 'v1 average (A$^3$)']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/property=bulk_modulus_v0_v1_calc_ts=0.5_p0=1atm_p1=5000atm.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1-3',
        'xdata': 'Step',
        'ydata': 'Volume',
        'xlabel': 'Time (ps)',
        'ylabel': 'Volume (A$^3$)',
        'xcompute': '${Step}*(1/2000)',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }

