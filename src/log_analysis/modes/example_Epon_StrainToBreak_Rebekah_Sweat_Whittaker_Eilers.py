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
            ['Kemppainen-Muzzy Modulus', '', '', 'shift=no; minxhi=0.005; maxxhi=0.0; xlo=rfs; yp=1; offset=0.0; csv=False', 'Modulus']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/Epon_StrainToBreak.Rebekah.Sweat',
        'keywords': ['Stress', 'Eyy-DIC'],
        'sections': '1',
        'xdata': 'Eyy-DIC',
        'ydata': 'Stress',
        'xlabel': 'Strain',
        'ylabel': 'Stress ($N/mm^2$)',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }

