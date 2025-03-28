# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 05/13/2024
Author: Josh Kemppainen
Purpose: Mode to find Tg and CTE using hyperbola
         method

"""
# analysis list
analysis = [['average', 290, 310, '', 'room temp volume average (A$^3$)'],
            ['moving average', 100, 800, 'window=100', 'moving average'],
            ['hyperbola', 100, 800, 'p=0.9', 'Hyperbola fit']]

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
        'analysis': analysis
        }

