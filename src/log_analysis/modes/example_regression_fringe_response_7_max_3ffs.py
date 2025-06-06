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
analysis = [['LAMMPS data (LOWESS)', '', '', 'fraction=0.015; max_iter=30', 'LOWESS'],
            ['LAMMPS data (apply Whittaker-Eilers)', '', '', 'order=2; lambda=100', 'Whittaker-Eilers'],
            ['LAMMPS data (apply moving average)', '', '', 'window=10', 'Moving average'],
            ['Regression Fringe Response Modulus', '', '0.25', 'shift=ymin; minxhi=0.005; maxxhi=0.07; xlo=rfs; yp=max-3ffs; offset=0.0; t1=trans; t2=trans; csv=True; ci=0.95', 'Modulus'],
            ['write plotted data to csv file', '', '', '', '']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/regression_fringe_response_example_7.txt',
        'keywords': ['strain', 'stress'],
        'sections': '1',
        'xdata': 'strain',
        'ydata': 'stress',
        'xlabel': 'True Strain',
        'ylabel': 'True Stress (MPa)',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile/JK_TW_TM_figs_no_legends_nu'
        }

