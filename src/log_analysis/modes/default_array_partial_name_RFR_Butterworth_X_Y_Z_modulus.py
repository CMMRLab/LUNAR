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
            ['Regression Fringe Response Modulus', '', '', 'shift=ymin; minxhi=0.005; maxxhi=0.0; xlo=rfs; yp=1; offset=0.0; up=-1; t1=v_etruey; t2=v_etruez; csv=False; savefig=all', 'RFR-mechanical'],
            ['maximum', '', '', '', 'UTS']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/tensile_*_*_pxld_*_replicate_*_FF_*.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '-1',
        'xdata': 'v_etruex',
        'ydata': 'f_sxx_ave',
        'xlabel': 'True Strain',
        'ylabel': 'True Stress (MPa)',
        'xcompute': 'partial_name( tensile_1=v_etruex; t1=v_etruey; t2=v_etruez,   tensile_2=v_etruey; t1=v_etruex; t2=v_etruez,   tensile_3=v_etruez; t1=v_etruex; t2=v_etruey  )',
        'ycompute': 'partial_name( tensile_1=f_sxx_ave,                            tensile_2=f_syy_ave,                            tensile_3=f_szz_ave )',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile/array_partial_name_RFR_tensile',
        'array_file': 'array_logging_results',
        }

