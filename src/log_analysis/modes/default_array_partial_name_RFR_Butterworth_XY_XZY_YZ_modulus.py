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
        'xcompute': 'partial_name( shear_1=v_etruexy,   shear_2=v_etruexz,   shear_3=v_etrueyz  )',
        'ycompute': 'partial_name( shear_1=f_sxy_ave,   shear_2=f_sxz_ave,   shear_3=f_syz_ave )',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile/array_partial_name_RFR_shear',
        'array_file': 'array_logging_results',
        }