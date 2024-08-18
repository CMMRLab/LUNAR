# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 05/13/2024
Author: Josh Kemppainen
Purpose: Mode to find XZ-poissions strained in
         X-dir


"""
# analysis list
analysis = [['moving average', 0, 0.1, 'window=50', 'moving average'],
            ['linear regression', 0, 0.02, '', '$\\nu_{xz}$  ']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/property=tensile_modulus_x_strain_rate=2e8.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1',
        'xdata': 'v_etruex',
        'ydata': 'v_etruez',
        'xlabel': 'True Strain (in X)',
        'ylabel': 'True Strain (in Z)',
        'xcompute': '',
        'ycompute': '',
        'analysis': analysis
        }

