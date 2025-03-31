# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 05/13/2024
Author: Josh Kemppainen
Purpose: Mode to find Tg and CTE using piecewise
         regression method


"""
# analysis list
analysis = [['average', 290, 310, '', 'room temp volume average (A$^3$)'],
            ['moving average', '', '', 'window=100', 'moving average'],
            ['piecewise-regression', '', '', 'n=1', 'Segmented piecewise-regression']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/properties=Tg_and_CTE_heating_rate=50_k_ns.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': 'all',
        'xdata': 'Temp',
        'ydata': 'Volume',
        'xlabel': 'Temperature (K)',
        'ylabel': 'Volume (A$^3$)',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }

