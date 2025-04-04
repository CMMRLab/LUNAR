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
analysis = [['iFFT filter', '', '', 'qm=3,2-p; threshold=mean:1.0; savefig=all', 'iFFT Filter'],
            ['Butterworth (low pass)', '', '', 'qm=msr; psd=True; csv=False; savefig=all; order=2; wn=op-p', 'Butterworth Filter']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/tensile_1_EPON_862_pxld_88.2_replicate_1_FF_PCFF-class2xe.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1',
        'xdata': 'v_etruex',
        'ydata': 'f_sxx_ave',
        'xlabel': 'True Strain',
        'ylabel': 'True Stress (MPa)',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }

