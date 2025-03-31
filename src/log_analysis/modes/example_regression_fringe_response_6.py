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
analysis = [['LAMMPS data (LOWESS)', '', '', 'fraction=0.01; max_iter=30', 'LOWESS'],
            ['LAMMPS data (apply Whittaker-Eilers)', '', '', 'order=2; lambda=300', 'Whittaker-Eilers'],
            ['Regression Fringe Response Modulus', '', 0.19, 'shift=ymin; minxhi=0.000; maxxhi=0.0; xlo=rfs; yp=1; offset=0.0; t1=trans; t2=trans; csv=True; ci=0.95', 'Modulus'],
            ['write plotted data to csv file', '', '', '', '']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/regression_fringe_response_example_6.txt',
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

