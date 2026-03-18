# -*- coding: utf-8 -*-
"""
Created by array_lmp_script.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within array_lmp_script.py.

"""
# LAMMPS script find and replace. Rules:
#    find = string (recommend using $<string>)
#    replace = iterable that eval() understands
lmp_array = {'$<d2theta>' : '[8, 4, 2, 1, 0.5, 0.25]'}


# Batch script find and replace. Rules:
#    find = string (recommend using $<string>)
#    replace = float, int, or string that eval() understands
batch_array = {'$<ncores>' : '32',
               '$<inv_ncores>' : '0.03125',
               '$<que>' : 'medium.q',
               '$<email>' : 'jdkemppa@mtu.edu',
               '$<input>' : '$<build_script>',
               '$<output>' : '$<build_script>.out'}


# List of files to copy into the build directories (find
# and replace operations work on these file strings).
copy_files = ['Fe_typed_PCFF.data',
              'Fe_typed_PCFF.ecoeffs',
              'Fe_typed_PCFF.lmpmol']


# loadable mode
mode = {'lmp_script': '../XRD_paper/XRD_pure_crystal/Iron_models/d2theta_testing/FE_crystal_manual.template',
        'batch_script': '$<lmp_script>.sh',
        'auto_script': 'auto_submit_$<lmp_script>',
        'build_script': '$<lmp_script>.script',
        'build_directory': 'out/d2theta_$<d2theta>',
        'lmp_script_find': '$<lmp_script>',
        'build_script_find': '$<build_script>',
        'auto_submit': 'qsub',
        'lmp_array': lmp_array,
        'batch_array': batch_array,
        'copy_files': copy_files}
