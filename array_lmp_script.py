# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
March 18, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############
### Inputs ###
##############
##################################################################################################################
# Python boolean variable to use GUI or not (True or False). All GUI settings are intialized from this script so #
# adjusting this script will set the default settings that the GUI will load with. Please NOTE that when using   #
# the GUI to check the console or terminal print outs every time a system is run through the code because that   #
# is where the import information will be displayed. Examples:                                                   #
#   use_GUI = True  # Will launch GUI to manipulate variables and run from GUI                                   #
#   use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have to be from the IDE     #
#                                                                                                                #
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account for the  #
# different ways individuals have their OS display settings setup. A GUI_zoom = 100, means that default settings #
# are used, whereas a GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80 means decrease GUI by 20%. #
# Examples:                                                                                                      #
#   GUI_zoom = 100 # use default GUI size                                                                        #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                                    #
#                                                                                                                #
# Update use_GUI and GUI_zoom as desired.                                                                        #
##################################################################################################################
use_GUI = True
GUI_zoom = 100



# LAMMPS script find and replace. Rules:
#    find = string (recommend using $<string>)
#    replace = iterable that eval() understands
lmp_array = {'$<replicate>': '[1, 2]',
             '$<atomic_frac>': [0.0, 0.25],
             '$<amorphous_density>': [1.5, 2.5],
             '$<elements>': 'lmp_data.elements_from_masses(graphite_AB_replicate_$<replicate>.data, Ar)'
             }

# Batch script find and replace. Rules:
#    find = string (recommend using $<string>)
#    replace = float, int, or string that eval() understands
batch_array = {'$<ncores>': '1',
               '$<inv_ncores>': '1/$<ncores>',
               '$<que>': 'meem3.q',
               '$<email>': 'jdkemppa@mtu.edu',
               '$<input>': '$<build_script>',
               '$<output>': '$<build_script>.out'
               }

# List of files to copy into the build directories (find
# and replace operations work on these file strings).
copy_files = ['*-atom.lmpmol', 'graphite_AB_replicate_$<replicate>.data']
#copy_files = []


# The full mode dictionary
mode = {'lmp_script': 'EXAMPLES/array_lmp_script/xrd_composite/example_xrd_composite.script',
        'batch_script': '$<lmp_script>.sh',
        'auto_script': 'auto_submit_$<lmp_script>',
        'build_script': '$<lmp_script>_$<replicate>.script',
        'build_directory': 'Graphite_Argon/Replicate_$<replicate>/AF_$<atomic_frac>_Density_$<amorphous_density>',
        'lmp_script_find': '$<lmp_script>',
        'build_script_find': '$<build_script>',
        'auto_submit': 'qsub',
        'lmp_array': lmp_array, 
        'batch_array': batch_array,
        'copy_files': copy_files 
        }


#########################
# GUI run mode settings #
#########################
settings = {'mode': 'src/array_lmp_script/modes/example_xrd_composite.py', # Default run mode to load GUI with.
            'modes-dir': 'src/array_lmp_script/modes',       # Set directory where all modes are saved.
            'replace_lmp_script_when_loading_mode': False}   # True or False default for replacing lmp script when loading a mode.



#####################################################
### Import main from src.array_lmp_script and run ###
#####################################################
if __name__ == "__main__":  
    if use_GUI:
        from src.array_lmp_script.GUI import GUI
        GUI(settings, GUI_zoom)
    else:
        from src.array_lmp_script.main import main
        main(mode, get_mode_from_script=True)