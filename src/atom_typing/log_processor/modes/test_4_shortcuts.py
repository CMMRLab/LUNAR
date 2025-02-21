# -*- coding: utf-8 -*-
"""
Created by atom_typing_log_processor.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: XX/YY/ZZZZ
Author: John Doe
Purpose: Mode to find xxx

"""

# dependant and independant variable strings
independent_variable = "numeric"
dependent_variables = """
numeric,    %Ring-Count:6,    %Ring-Count:C6,    hybridizations['all-H']['%natoms']
"""

# logger string
logger = """
filename,    numeric,    wildcards,    %Ring-Count:6,    %Ring-Count:C6,    %Hybrid-Mass:Sp1-C,
%Hybrid-Mass:Sp2-C,   %Hybrid-Mass:Sp3-C,   %Hybrid-Count:all-C,    hybridizations['all-H']['%natoms']
"""

# loadable mode
mode = {'independent_variable': independent_variable,
        'dependent_variables': dependent_variables,
        'parent_directory': 'logfile/tests',
        'sorting_direction': 'ascending',
        'sorting_method': 'natsort',
        'logfile': 'EXAMPLES/array_processing/atom_typing_logfile_processor/poly_tracking_replicate_1_time_*',
        'newfile': 'test_4_shortcuts',
        'logger': logger
        }

