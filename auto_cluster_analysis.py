# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
November 4th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    **************************************************************
    * Requirements:                                              *
    *   python 3.7+                                              *
    *                                                            *
    * Dependencies:                                              *
    *   python natsort module:                                   *
    *    - pip3 install natsort (if pip manager is installed)    *
    *    - if not installed files will be sorted using the       *
    *      standard python sort function (which may not sort     *
    *      files in a human logical manner).                     *
    *                                                            *
    *   ALL read in topofiles from the files_directory MUST have *
    *   some NUMBER OR LETTER that represents the "progression"  *
    *   or "evolution" sequence of the topofiles. Usually this   *
    *   will be a "time stamp" being used to write the LAMMPS    *
    *   datafile using the nevery command in LAMMPS, but also    *
    *   could be a crosslink density value as well that provides *
    *   a logical sequencing of the sorted topofiles.            *
    *                                                            *
    * Run methods:                                               *
    *   - IDE (manipulate variables and run from IDE)            *
    *   - GUI (manipulate variables and run from. Default        *
    *          settings set from this script)                    *
    *   - command line (python3 auto_cluster_analysis.py -man to *
    *                   get command line override options and    *
    *                   examples)                                *
    *                                                            *
    **************************************************************
"""


##############
### Inputs ###
##############
#################################################################################################################
# Python boolean variable to use GUI or not (True or False). All GUI settings are intialized from this script   #
# so adjusting this script will set the default settings that the GUI will load with. Please NOTE that when     #
# using the GUI to check the console or terminal print outs every time a system is run through the code         #
# because that is where the import information will be displayed. Examples:                                     #
#   use_GUI = True  # Will launch GUI to manipulate variables and run from GUI                                  #
#   use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have to be from the        #
#                     IDE or command line                                                                       #
#                                                                                                               #
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account for the #
# different ways individuals have their OS display settings setup. A GUI_zoom = 100, means that default settings#
# are used, whereas a GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80 means decrease GUI by 20%.#
# Examples:                                                                                                     #
#   GUI_zoom = 100 # use default GUI size                                                                       #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                                   #
#                                                                                                               #
# Update use_GUI and GUI_zoom as desired.                                                                       #
#################################################################################################################
use_GUI = True
GUI_zoom = 100


#################################################################################################################
# Python string variable type to set directory where LAMMPS datafiles are stored. The code will change to this  #
# directory to read ALL files in this directory with the *.data extension. This means ONLY have LAMMPS *.data   #
# files in this directory you want to read into this code to automate the molecule logging process. Examples:   #
#     files_directory = 'EXMAPLES/auto_cluster_analysis' # Will read all *.data files in this directory         #
#                                                                                                               #
# Update files_directory as desired.                                                                            #
#################################################################################################################
files_directory = 'EXAMPLES/auto_cluster_analysis'


#################################################################################################################
# Python string variable type to set new file name(s). Example for newfile = 'auto_logging_molecule_data'       #
#     filename = auto_logging_molecule_data.csv  # *NOTE info in file will be dependant on info in logger dict  #
#                                                                                                               #
# Update newfile as desired.                                                                                    #
#################################################################################################################
newfile = 'auto_logging_molecule_data'


#################################################################################################################
# Inputs for computing degree of polymerization (Xn) and Extent of reaction (p). Equations and info can be      #
# found in: Polymers  Chemistry and Physics of Modern Materials, Third Edition by Cowie page 42                 #
#                                                                                                               #
# Extent of reaction:                                                                                           #
#    p = (2*(N0 - N))/(N0*fav)                                                                                  #
#                                                                                                               #
# Critical Extent of reaction (gel point):                                                                      #
#   pg = 2/fav                                                                                                  #
#                                                                                                               #
# Degree of polymerization:                                                                                     #
#    Xn = 2/(2 - p*fav)                                                                                         #
#                                                                                                               #
# Where:                                                                                                        #
#    p   = Extent of reaction                                                                                   #
#    pg  = Critical Extent of reaction (gel point)                                                              #
#    N0  = Number of intial molecules before polymerization                                                     #
#    N   = Number of molecules at any stage or polymerization                                                   #
#    fav = Average number of functional groups present per monomer unit                                         #
#    Xn  = Degree of polymerization                                                                             #
#                                                                                                               #
# Determining fav for EPON 862 example:                                                                         #
#    1 DETDA molecule (4-functional groups)                                                                     #
#    2 DGEBF molecule (2-fuctional groups)                                                                      #
#    mix = 1*DETDA + 2*DGEBF (2 DETDA to every 1 DGEBF)                                                         #
#    fav = (1-DETDA*4-functional-groups + 2-DGEBF*2-functional-groups)/3-molecules                              #
#    fav = (1*4 + 2*2)/3 = 2.66                                                                                 #
#################################################################################################################
N0 = 150
fav = 2.66


#################################################################################################################
# Option to write a .txt file of BASENAME.txt (topofile) of all connectivty information found (True or False).  #
#################################################################################################################
txtfile = False


###########################################
### Main auto_cluster_analysis function ###
###########################################
import src.io_functions as io_functions
def main(files_directory, N0, fav, txtfile, newfile, log=io_functions.LUNAR_logger()):
    import src.clusters as clusters
    import time
    import os
    
    # Get time
    start_time = time.time()
    
    # Configure log (default is level='production', switch to 'debug' if debuging)
    log.configure(level='production')
    #log.configure(level='debug')
    
    ###########################################################################
    # Setting up directories and where to write final files and results to    #
    ###########################################################################
    # Find present working directory and Find/create paths to store code results
    pwd = os.getcwd(); 
    path = os.path.join(pwd, files_directory);
    
    # Check if path exists. IF not create.
    if not os.path.isdir(path): 
        os.makedirs(path, exist_ok=True)
        
    # Change the current working directory to path to write all new files to outputs directory
    os.chdir(path)
    log.out('Moved to {} directory to find LAMMPS datafiles to iterate through'.format(path))
    
    # Check for newfile, if already exists delete
    if os.path.exists(os.path.join(path, newfile+'.csv')):
        os.remove(os.path.join(path, newfile+'.csv'))
    
    ################################################
    # Function to write .csv file from logger dict #
    ################################################
    def write_csv(logger, newfile):
        # Function to write to .csv file with mode='w' or mode='a'
        def write_from_logger(logger, mode):
            with open(newfile, mode) as f:
                # Find titles and data to write
                titles = []; data = [];
                for i in logger:
                    titles.append('{}'.format(i)); data.append('{}'.format(logger[i]))
                # Join with comma's/Write title and data if file does not exists 
                if mode == 'w': f.write('{}\n'.format(', '.join(titles))) # only write headers in 'w' mode
                f.write('{}\n'.format(', '.join(data))) # always write data in either 'w' or 'a' mode
            return
        
        # Check if path+newfile exists, if not write title and data. else append to file
        if not os.path.exists(os.path.join(path, newfile)): write_from_logger(logger, mode='w')                
        else: write_from_logger(logger, mode='a')  
        return

    ###############################################################################################################
    # Find molecules and store them in the molecules instance. Info available from the instance:                  #
    # molecules.ATTRIBUTE with the following attributes:                                                          #
    #                                                                                                             #
    #  - molecules.clusters = [ {molID 1 atomIDs}, {molID2 atomIDs}, NmolIDs, ... ]                               #
    #  - molecules.info[molID].ATTRIBUTE where the following attributes are available:                            #
    #      * .atoms = set(atomIDs in molID)                                                                       #
    #      * .size = number of atoms in molID                                                                     #
    #      * .mass = mass of molID                                                                                #
    #      * .psize = % number of atoms in molID                                                                  #
    #      * .pmass = % mass of molID fragment                                                                    #
    #      * Examples:                                                                                            #
    #          molecules.info[1].atoms (gets set of atomIDs in molID 1)                                           #
    #          molecules.info[2].mass (gets mass of molID 2)                                                      #
    #          molecules.info[3].size (gets natoms of molID 3)                                                    #
    #      * NOTES:                                                                                               #
    #          molIDs are set based on cluster size that is 1st sorted by number of atoms.can ONLY access molIDs  #
    #          that exists (the 1st 10 molIDs are intialized w/ ZERO's to help access molIDs that may not exist). #
    #  - molecules.mass_total = total mass of system in mass units of datafile                                    #
    #  - molecules.size_total = total number of atoms in the read in datafile                                     #
    #  - molecules.p = Extent of reaction (chemistry/eqn specific. May need to be changed)                        #
    #  - molecules.Xn = Degree of polymerization (chemistry/eqn specific. May need to be changed)                 #
    #  - molecules.Mw = Weight-average molar mass  (Mw)                                                           #
    #  - molecules.Mn = Number-average molar mass  (Mn)                                                           #
    #  - molecules.Mz = Higher-average molar mass  (Mz)                                                           #
    #  - molecules.Mz1 = Higher-average molar mass (Mz+1)                                                         #
    #  - molecules.RMW = Weight-averge reduced molecular weight (RMW)                                             #
    ###############################################################################################################
    try:
        from natsort import natsorted
        topofiles = natsorted([file for file in os.listdir(path) if file.endswith('.data')])
        log.out('Read in topofiles were sorted using natsort, thus they should be iterated through')
        log.out('in a sequence that represents evolution or progression and logged in that order.')
    except:
        topofiles = sorted([file for file in os.listdir(path) if file.endswith('.data')])
        log.warn('WARNING natsort was not installed and read in topofiles were sorted using pythons')
        log.out('sorted function. The file iteration sequence may be unorder from a evolution or')
        log.out('progression stand point. Depending on the names given to the read in files. You may')
        log.out('intall natsort with:  pip3 install natsort   (if pip manager is installed)')
    log.out('\n\nIterating through topofiles to find molecule information ....')
    for iteration, topofile in enumerate(topofiles, 1):
        #---------------------------------------------------------------------------------------#
        # Find molecules and store data in a class w/multiple layers of attributes listed above #
        #---------------------------------------------------------------------------------------#
        molecules = clusters.analysis(topofile, N0, txtfile, fav, pflag=False, log=log)
        
        #---------------------------------------------------------------------------------#
        # Set up logger dict to hold the info to write to the .csv file. Update as desired#
        #         header name in file:  written value                                     #
        #---------------------------------------------------------------------------------#
        logger = {'Iteration':          iteration,
                  '1st cluster Mass':   molecules.info[1].mass,
                  '2nd cluster Mass':   molecules.info[2].mass,
                  '1st cluster %Mass':  molecules.info[1].pmass,
                  '2nd cluster %Mass':  molecules.info[2].pmass,
                  'RMW':                molecules.RMW,
                  'filename':           molecules.filename,
                  }
 
        #-----------------------------------#
        # Write logger results to .csv file #
        #-----------------------------------#
        write_csv(logger, newfile+'.csv')
        
        
    ####################################
    # Cleanup and finalization of code #
    ####################################
    # Change back to the intial directory to keep directory free for deletion
    os.chdir(pwd)
    
    # Print file locations and completion of code
    log.out(f'\n\nAll outputs can be found in {path} directory')
    log.out('\n\nNormal program termination\n\n')
    
    # Script run time
    execution_time = (time.time() - start_time)
    log.out('Execution time in seconds: ' + str(execution_time))
    
    # Show number of warnings and errors
    log.out_warnings_and_errors()
    return

###################################
### Import needed files and run ###
###################################
if __name__ == "__main__": 
    import src.command_line as cl
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
    #       python3 auto_cluster_analysis.py -opt                                                      #
    #              or                                                                                  #
    #       python3 auto_cluster_analysis.py -man                                                      #
    #   in the terminal and the code will print out all command line option and terminate before any   #
    #   further analysis is done                                                                       #
    ####################################################################################################
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: 
        files_directory, N0, txtfile, fav, newfile = cl.clusters_auto(files_directory, N0, txtfile, fav, newfile, commandline_inputs)
        use_GUI = False

    if use_GUI or '-gui' in commandline_inputs:
        print('\n\n\nauto_cluster_analysis is currently running in GUI mode, where all GUI inputs are intialized from auto_cluster_analysis.\n\n\n')
        from src.auto_clusters_GUI import GUI 
        GUI(files_directory, N0, fav, txtfile, newfile, GUI_zoom)
    else: main(files_directory, N0, fav, txtfile, newfile)