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
independent_variable = "wildcards[0]"
dependent_variables = """
logger
"""

# logger string
logger = """
filename,    wildcards,    numeric,    wildcards[0],    wildcards[1],    wildcards[-1],    wildcards[-2],    wildcards[-3]

rings[5]['Count'],   rings[5]['%Ring'],    rings[5]['C']['%Mass'],
"""

# loadable mode
mode = {'independent_variable': independent_variable,
        'dependent_variables': dependent_variables,
        'parent_directory': 'logfile/tests',
        'sorting_direction': 'ascending',
        'sorting_method': 'wildcards[0]',
        'logfile': 'EXAMPLES/array_processing/atom_typing_logfile_processor/heating_temp_*K_time_*ps_bonds_typed.log.lunar',
        'newfile': 'test_9_multiple_wildcards',
        'logger': logger
        }

