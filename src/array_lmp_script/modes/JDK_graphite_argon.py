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
             '$<atomic_frac>' : '[0.0, 0.25, 0.5, 0.75, 0.80, 0.85, 0.9, 0.95, 1.0]',
             '$<amorphous_density>' : '[0.125, 0.25, 0.50, 0.75, 1.0, 1.5, 2.26]'}


# Batch script find and replace. Rules:
#    find = string (recommend using $<string>)
#    replace = float, int, or string that eval() understands
batch_array = {'$<ncores>' : '1',
               '$<inv_ncores>' : '1/$<ncores>',
               '$<que>' : 'short.q',
               '$<email>' : 'jdkemppa@mtu.edu',
               '$<input>' : '$<build_script>',
               '$<output>' : '$<build_script>.out'}


# List of files to copy into the build directories (find
# and replace operations work on these file strings).
copy_files = ['Ar-atom.lmpmol',
              'graphite_AB_25x25x25_PCFF.data']


# loadable mode
mode = {'lmp_script': '../XRD_paper/Graphite_argon/atomic_fraction/AF_Graphite_argon.template',
        'batch_script': '$<lmp_script>.sh',
        'auto_script': 'auto_submit_$<lmp_script>',
        'build_script': '$<lmp_script>_$<replicate>',
        'build_directory': 'Graphite_Argon/Replicate_$<replicate>/AmorDensity_$<amorphous_density>_AtomicFrac_$<atomic_frac>',
        'lmp_script_find': '$<lmp_script>',
        'build_script_find': '$<build_script>',
        'auto_submit': 'qsub',
        'lmp_array': lmp_array,
        'batch_array': batch_array,
        'copy_files': copy_files}
