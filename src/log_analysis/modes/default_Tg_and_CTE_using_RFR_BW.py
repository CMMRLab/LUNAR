# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 05/13/2024
Author: Josh Kemppainen
Purpose: Mode to find Tg and CTE using hyperbola
         method

Tg ~= 450 K

"""
# analysis list
analysis = [['LAMMPS data (apply Butterworth filter)', '', '', 'qm=msr; psd=False; csv=False; savefig=all; order=2; wn=op', 'LAMMPS Butterworth'],
            ['Regression Fringe Response Thermal', '', '', 'rtemp=290,310; minspan=50; tg=PDF(prob, 0); weights=bw*rw*rs', 'RFR-thermal'],
            ['hyperbola', '', '', 'p=0.9; initial_guess=True', 'Hyperbola fit'],
            ['write plotted data to csv file', '', '', '', '']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/properties=Tg_and_CTE_heating_rate=50_k_ns.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '-1; wavg=False',
        'xdata': 'Temp',
        'ydata': 'Volume',
        'xlabel': 'Temperature (K)',
        'ylabel': 'Volume (A$^3$)',
        'xcompute': 'np.linspace(min(${Temp}), max(${Temp}), len(${Temp}))',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile/RFR_thermal',
        'array_file': 'array_logging_results',
        }