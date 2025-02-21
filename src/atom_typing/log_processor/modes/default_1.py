# -*- coding: utf-8 -*-
"""
Created by atom_typing_log_processor.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 02/18/2025
Author: Josh Kemppainen
Purpose: Mode to find show basic setup for
         atom_typing_log_processor.py.

"""

# dependant and independant variable strings
independent_variable = "numeric"
dependent_variables = """
rings[6]['%Ring'],    rings[6]['C']['%natoms'],    hybridizations['Sp1-C']['%natoms']
hybridizations['Sp2-C']['%natoms'],    hybridizations['Sp3-C']['%natoms'],    hybridizations['all-H']['%natoms']
"""

# logger string
logger = """
filename,    numeric,    wildcards,    rings[6]['%Ring'],    rings[6]['C']['%natoms'],    hybridizations['Sp1-C']['%natoms']
hybridizations['Sp2-C']['%natoms'],    hybridizations['Sp3-C']['%natoms'],    hybridizations['all-H']['%natoms']
"""

# loadable mode
mode = {'independent_variable': independent_variable,
        'dependent_variables': dependent_variables,
        'parent_directory': 'logfile',
        'sorting_direction': 'ascending',
        'sorting_method': 'sort',
        'logfile': 'EXAMPLES/array_processing/atom_typing_logfile_processor/poly_tracking_replicate_1_time_*ps_typed.log.lunar',
        'newfile': 'atom_typing_logfile_processor',
        'logger': logger
        }

