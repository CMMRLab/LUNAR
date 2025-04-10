# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.5
April 10th, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


    **************************************************************
    * Requirements:                                              *
    *   python 3.7+                                              *
    *                                                            *
    * Run methods:                                               *
    *   - IDE (manipulate variables and run from IDE)            *
    *   - GUI (manipulate variables and run from. Default        *
    *          settings set from this script)                    *
    *   - command line (python3 cluster_analysis.py -man to get  *
    *                   command line override options and        *
    *                   examples)                                *
    *                                                            *
    **************************************************************
"""
##############
### Inputs ###
##############
#########################################################################################################
# Python boolean variable to use GUI or not (True or False). All GUI settings are intialized from this  #
# script so adjusting this script will set the default settings that the GUI will load with. Please     #
# NOTE that when using the GUI to check the console or terminal print outs every time a system is run   #
# through the code because that is where the import information will be displayed. Examples:            #
#   use_GUI = True  # Will launch GUI to manipulate variables and run from GUI                          #
#   use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have to be from the#
#                     IDE or command line                                                               #
#                                                                                                       #
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account #
# for the different ways individuals have their OS display settings setup. A GUI_zoom = 100, means that #
# default settings are used, whereas a GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80  #
# means decrease GUI by 20%. Examples:                                                                  #
#   GUI_zoom = 100 # use default GUI size                                                               #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                           #
#                                                                                                       #
# Update use_GUI and GUI_zoom as desired.                                                               #
#########################################################################################################
use_GUI = True
GUI_zoom = 100


#########################################################################################################
# LAMMPS topofile with bonds in it (THIS IS A MUST!). The LAMMPS datafile can have new LAMMPS "type     #
# labels" in it. The only supported atom styles are full, charge, or molecular. More atom styles can be #
# coded if needed.                                                                                      #
#########################################################################################################
#topofile = 'EXAMPLES/auto_cluster_analysis/EPON_862_pxld_60_typed.data'
#topofile = 'EXAMPLES/periodicity/PAN_crystal_PCFF.data'
#topofile = 'EXAMPLES/periodicity/CNT_PCFF.data'
#topofile = 'EXAMPLES/periodicity/EPON_862_pxld_60_typed.data'
#topofile = 'EXAMPLES/periodicity/dgebf_typed_IFF.data'
topofile = 'EXAMPLES/periodicity/Cellulose-supercell_PCFF.data'


#########################################################################################################
# Inputs for computing degree of polymerization (Xn) and Extent of reaction (p). Equations and info can #
# be found in: Polymers  Chemistry and Physics of Modern Materials, Third Edition by Cowie page 42      #
#                                                                                                       #
# Extent of reaction:                                                                                   #
#    p = (2*(N0 - N))/(N0*fav)                                                                          #
#                                                                                                       #
# Critical Extent of reaction (gel point):                                                              #
#   pg = 2/fav                                                                                          #
#                                                                                                       #
# Degree of polymerization:                                                                             #
#    Xn = 2/(2 - p*fav)                                                                                 #
#                                                                                                       #
# Where:                                                                                                #
#    p   = Extent of reaction                                                                           #
#    pg  = Critical Extent of reaction (gel point)                                                      #
#    N0  = Number of intial molecules before polymerization                                             #
#    N   = Number of molecules at any stage or polymerization                                           #
#    fav = Average number of functional groups present per monomer unit                                 #
#    Xn  = Degree of polymerization                                                                     #
#                                                                                                       #
# Determining fav for EPON 862 example:                                                                 #
#    1 DETDA molecule (4-functional groups)                                                             #
#    2 DGEBF molecule (2-fuctional groups)                                                              #
#    mix = 1*DETDA + 2*DGEBF (2 DETDA to every 1 DGEBF)                                                 #
#    fav = (1-DETDA*4-functional-groups + 2-DGEBF*2-functional-groups)/3-molecules                      #
#    fav = (1*4 + 2*2)/3 = 2.66                                                                         #
#########################################################################################################
N0 = 150
fav = 2.66


#########################################################################################################
# Option to write a file of BASENAME.txt file of all connectivty information found (True or False).     #
#########################################################################################################
txtfile = True


###################################
### Import needed files and run ###
###################################
if __name__ == "__main__":  
    ##############################
    # Import Necessary Libraries #
    ##############################
    import src.command_line as cl
    import src.clusters as clusters
    import sys

    ####################################################################################################
    # Set IDE/command line mode option:                                                                #
    # - if run in IDE commandline_inputs lst will have zero entries. Which means use hard coded inputs #
    #   from inputs section.                                                                           #
    #                                                                                                  #
    # - if run in command line, but commandline_inputs still has zero entries, this means use hard     #
    #   coded inputs from inputs section.                                                              #
    #                                                                                                  #
    # - else run in command line and commandline_inputs has entries determine what entries will be     #
    #   overrided from hard coded inputs from inputs section. This overide will occur in in the inputs #
    #   section called Command Line Override. *NOTE: if user does not specify certain command line     #
    #   options, the default settings from inputs section will be enforced. A message will tell user   #
    #   if inputs from inputs section have been enforced. The user can find all command line options   #
    #   by running:                                                                                    #
    #       python3 cluster_analysis.py -opt                                                           #
    #              or                                                                                  #
    #       python3 cluster_analysis.py -man                                                           #
    #   in the terminal and the code will print out all command line option and terminate before any   #
    #   further analysis is done                                                                       #
    ####################################################################################################
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: 
        topofile, N0, txtfile, fav = cl.clusters_single(topofile, N0, txtfile, fav, commandline_inputs)
        use_GUI = False

    if use_GUI or '-gui' in commandline_inputs:
        print('\n\n\ncluster_analysis is currently running in GUI mode, where all GUI inputs are intialized from cluster_analysis.\n\n\n')
        from src.clusters_GUI import GUI 
        GUI(topofile, N0, txtfile, fav, GUI_zoom, pflag=True)
    else: molecules = clusters.analysis(topofile, N0, txtfile, fav)