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
analysis = [['LAMMPS data (apply Whittaker-Eilers)', '', '', 'order=2; lambda=1_000', 'Whittaker-Eilers'],
            ['Regression Fringe Response Modulus', '', 0.25, 'shift=ymin; minxhi=0.000; maxxhi=0.0; xlo=rfs; yp=1; offset=0.0; t1=v_etruey; t2=v_etruez; csv=False; ci=0.95', 'Modulus']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/regression_fringe_response_example_2.txt',
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
        'parent_directory': 'logfile/Published'
        }

