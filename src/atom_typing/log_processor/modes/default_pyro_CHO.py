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
logger
"""

# logger string
logger = """
filename,    numeric,    wildcards,    mass,    CharYield:78829.55,

rings[5]['Count'],   rings[5]['%Ring'],    rings[5]['C']['%Mass'],
rings[6]['Count'],   rings[6]['%Ring'],    rings[6]['C']['%Mass'],
rings[7]['Count'],   rings[7]['%Ring'],    rings[7]['C']['%Mass'],

hybridizations['Sp1-C']['%Mass'],  hybridizations['Sp2-C']['%Mass'],  hybridizations['Sp3-C']['%Mass'],
hybridizations['Sp1-O']['%Mass'],  hybridizations['Sp2-O']['%Mass'],  hybridizations['Sp3-O']['%Mass'],

hybridizations['all-C']['%Mass'], hybridizations['all-H']['%Mass'],   hybridizations['all-O']['%Mass'],
"""

# loadable mode
mode = {'independent_variable': independent_variable,
        'dependent_variables': dependent_variables,
        'parent_directory': 'logfile',
        'sorting_direction': 'ascending',
        'sorting_method': 'natsort',
        'logfile': 'EXAMPLES/array_processing/atom_typing_logfile_processor/poly_tracking_replicate_1_time_*ps_typed.log.lunar',
        'newfile': 'atom_typing_logfile_processor',
        'logger': logger
        }

