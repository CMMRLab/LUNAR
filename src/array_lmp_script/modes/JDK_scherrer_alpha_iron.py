# -*- coding: utf-8 -*-
"""
Created by array_lmp_script.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within array_lmp_script.py.

Date: XX/YY/ZZZZ
Author: John Doe
Purpose: Mode to find xxx

"""
# LAMMPS script find and replace. Rules:
#    find = string (recommend using $<string>)
#    replace = iterable that eval() understands
lmp_array = {'$<replicate>' : '[1, 2, 3]',
             '$<density>' : '[0.25, 0.5, 1.0]'}


# Batch script find and replace. Rules:
#    find = string (recommend using $<string>)
#    replace = float, int, or string that eval() understands
batch_array = {'$<ncores>' : '16',
               '$<inv_ncores>' : '1/$<ncores>',
               '$<que>' : 'long.q',
               '$<email>' : 'jdkemppa@mtu.edu',
               '$<input>' : '$<build_script>',
               '$<output>' : '$<build_script>.out'}


# List of files to copy into the build directories (find
# and replace operations work on these file strings).
copy_files = ['poly_crystalline_iron_25x25x25_$<density>gcc_rep_$<replicate>.data']


# loadable mode
mode = {'lmp_script': '../XRD_paper/XRD_Scherrer_eqn/Iron_models/LAMMPS_XRD/Poly_crystal_25x25x25/PolyIron_XRD_PCFF.template',
        'batch_script': '$<lmp_script>.sh',
        'auto_script': 'auto_submit_$<lmp_script>',
        'build_script': '$<lmp_script>.script',
        'build_directory': 'replicate_$<replicate>_$<density>gcc',
        'lmp_script_find': '$<lmp_script>',
        'build_script_find': '$<build_script>',
        'auto_submit': 'qsub',
        'lmp_array': lmp_array,
        'batch_array': batch_array,
        'copy_files': copy_files}
