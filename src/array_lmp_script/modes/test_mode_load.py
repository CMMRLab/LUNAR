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
lmp_array = {'$<replicate>': '[1, 2]',
             '$<atomic_frac>': [0.0, 0.25],
             '$<amorphous_density>': [1.5, 2.5],
             '$<crystal_density>': [2.26, 2.27]}


# Batch script find and replace. Rules:
#    find = string (recommend using $<string>)
#    replace = float, int, or string that eval() understands
batch_array = {'$<ncores>' : '1.TEST',
               '$<inv_ncores>' : '1/$<ncores>.TEST',
               '$<que>' : 'meem3.q.TEST',
               '$<email>' : 'jdkemppa@mtu.edu.TEST',
               '$<input>' : '$<build_script>.TEST'}


# List of files to copy into the build directories (find
# and replace operations work on these file strings).
copy_files = ['Ar-atom.lmpmol', 'graphite_AB_25x25x25_PCFF.data']


# loadable mode
mode = {'lmp_script': 'EXAMPLES/array_lmp_script/xrd_composite/atomic_fraction_XRD.script.TEST',
        'batch_script': 'atomic_fraction_XRD.script.TEST.sh',
        'auto_script': 'auto_submit_$<lmp_script>.TEST',
        'build_script': '$<lmp_script>_$<replicate>.TEST',
        'build_directory': 'Graphite_Argon/Replicate_$<replicate>/AF_$<atomic_frac>_Density_$<amorphous_density>.TEST',
        'lmp_script_find': '$<lmp_script>',
        'build_script_find': '$<build_script>',
        'auto_submit': 'qsub.TEST',
        'lmp_array': lmp_array,
        'batch_array': batch_array,
        'copy_files': copy_files}
