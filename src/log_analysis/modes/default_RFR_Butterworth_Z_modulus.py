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
analysis = [['LAMMPS data (apply Butterworth filter)', '', '', 'qm=msr; psd=False; csv=False; savefig=all; order=2; wn=op', 'Butterworth'],
            ['Regression Fringe Response Modulus', '', '', 'shift=ymin; minxhi=0.0025; maxxhi=0.0; xlo=rfs; yp=1; up=-1; offset=0.0; t1=v_etruex; t2=v_etruey; csv=False; savefig=all', 'RFR-mechanical']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/tensile_3_EPON_862_pxld_88.2_replicate_1_FF_PCFF-class2xe.log.lammps',
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

