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
rings[6]['%Ring'],    rings[6]['C']['%natoms']
"""

# logger string
logger = """
filename,    wildcards,    numeric,    rings[6]['%Ring']
"""

# loadable mode
mode = {'independent_variable': independent_variable,
        'dependent_variables': dependent_variables,
        'parent_directory': 'logfile/tests',
        'sorting_direction': 'descending',
        'sorting_method': 'natsort',
        'logfile': 'EXAMPLES/array_processing/atom_typing_logfile_processor/poly_tracking_replicate_1_time_*',
        'newfile': 'test_2',
        'logger': logger
        }

