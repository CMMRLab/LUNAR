# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 05/13/2024
Author: Josh Kemppainen
Purpose: Mode to find shear modulus in YZ-dir
         using linear regressions


"""
# analysis list
analysis = [['moving average', 0, 0.1, 'window=50', 'moving average'],
            ['linear regression', 0, 0.02, 'shift=True', 'Modulus'],
            ['cursor', '', '', 'x=0.02; y=5;', 'Yield strength']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/property=shear_modulus_yz_strain_rate=2e8.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1',
        'xdata': 'v_etrueyz',
        'ydata': 'f_syz_ave',
        'xlabel': 'True Strain',
        'ylabel': 'True Stress (MPa)',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis
        }

