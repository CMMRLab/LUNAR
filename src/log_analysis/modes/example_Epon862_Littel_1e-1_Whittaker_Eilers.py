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
analysis = [['LAMMPS data (apply Whittaker-Eilers)', '', '', 'order=2; lambda=op<1e-2, 1e12, 50>-p', 'LAMMPS Whittaker-Eilers'],
            ['Kemppainen-Muzzy Modulus', '', '', 'shift=no; minxhi=0.005; maxxhi=0.0; xlo=rfs; yp=1; offset=0.01; csv=False', 'Modulus']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/Figure_12_1e-1.log.Littell',
        'keywords': ['strain', 'stress'],
        'sections': '1',
        'xdata': 'strain',
        'ydata': 'stress',
        'xlabel': 'Strain',
        'ylabel': 'Stress (MPa)',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }

