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
            ['Regression Fringe Response Modulus', '', '', 'shift=ymin; minxhi=0.01; maxxhi=0.05; xlo=rfs; yp=0; csv=False; savefig=all', 'Modulus']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/shear_2_EPON_862_pxld_88.2_replicate_1_FF_PCFF-class2xe.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1',
        'xdata': 'v_etruexz',
        'ydata': 'f_sxz_ave',
        'xlabel': 'True Strain',
        'ylabel': 'True Stress (MPa)',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }

