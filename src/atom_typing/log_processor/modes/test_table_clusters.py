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
logger
"""

# logger string
logger = """
filename,    numeric,    wildcards,    clusters[1]['Molecule Size'],  clusters[2]['Molecule Size']
clusters[3]['Molecule Size'],    clusters[1]['Molecule Formula']
"""

# loadable mode
mode = {'independent_variable': independent_variable,
        'dependent_variables': dependent_variables,
        'parent_directory': 'logfile',
        'sorting_direction': 'ascending',
        'sorting_method': 'sort',
        'logfile': 'EXAMPLES/array_processing/atom_typing_logfile_processor/heating_temp_*K_time_*ps_bonds_typed.log.lunar',
        'newfile': 'atom_typing_logfile_processor',
        'logger': logger
        }

