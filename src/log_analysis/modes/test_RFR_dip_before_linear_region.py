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
analysis = [['LAMMPS data (apply Butterworth filter)', '', '', 'qm=msr; psd=True; csv=False; order=2; wn=op', 'LAMMPS Butterworth Filter'],
            ['Regression Fringe Response Modulus', '', '', 'shift=ymin; minxhi=0.005; maxxhi=0.0; xlo=rfs; yp=1; offset=0.0; t1=v_etruex; t2=v_etruey; csv=False', 'Modulus']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/tensile_Z_dip_before_linear_region.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1',
        'xdata': 'v_etruez',
        'ydata': 'f_szz_ave',
        'xlabel': 'True Strain',
        'ylabel': 'True Stress (MPa)',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }

