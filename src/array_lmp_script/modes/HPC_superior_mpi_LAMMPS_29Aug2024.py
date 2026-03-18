# -*- coding: utf-8 -*-
"""
Created by array_lmp_script.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within array_lmp_script.py.

Date: 10/09/2025
Author: Josh Kemppainen
Purpose: Example mode

"""
# LAMMPS script find and replace. Rules:
#    find = string (recommend using $<string>)
#    replace = iterable that eval() understands
lmp_array = {}


# Batch script find and replace. Rules:
#    find = string (recommend using $<string>)
#    replace = float, int, or string that eval() understands
batch_array = {'$<ncores>': '1',
               '$<inv_ncores>': '1/$<ncores>',
               '$<que>': '',
               '$<email>': '',
               '$<input>': '$<build_script>',
               '$<output>': '$<build_script>.out'}


# List of files to copy into the build directories (find
# and replace operations work on these file strings).
copy_files = []


# The full mode dictionary
mode = {'lmp_script': 'superior_mpi_LAMMPS_29Aug2024',
        'batch_script': '$<lmp_script>.sh',
        'auto_script': 'auto_submit_$<lmp_script>',
        'build_script': '$<lmp_script>',
        'build_directory': '.',
        'lmp_script_find': '$<lmp_script>',
        'build_script_find': '$<build_script>',
        'auto_submit': 'qsub',
        'lmp_array': lmp_array, 
        'batch_array': batch_array,
        'copy_files': copy_files}