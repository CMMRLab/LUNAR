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
            ['Regression Fringe Response Modulus', '', '', 'shift=ymin; minxhi=0.01; maxxhi=0.0; xlo=rfs; yp=1; csv=False; savefig=all', 'RFR-mechanical'],
            ['maximum', '', '', '', 'UTS']]

# loadable mode

mode = {'logfile': 'EXAMPLES/log_analysis/shear_*_*_pxld_*_replicate_*_FF_*.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '-1',
        'xdata': 'v_etruexy',
        'ydata': 'f_sxy_ave',
        'xlabel': 'True Strain',
        'ylabel': 'True Stress (MPa)',
        'xcompute': 'wild_match( wildcards[0]=1=v_etruexy,   wildcards[0]=2=v_etruexz,   wildcards[0]=3=v_etrueyz  )',
        'ycompute': 'wild_match( wildcards[0]=1=f_sxy_ave,   wildcards[0]=2=f_sxz_ave,   wildcards[0]=3=f_syz_ave )',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile/array_wild_match_RFR_shear',
        'array_file': 'array_logging_results',
        }