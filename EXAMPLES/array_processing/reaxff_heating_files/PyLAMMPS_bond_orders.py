# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 20, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##################################
### Import Necessary Libraries ###
##################################
from lammps import lammps
from mpi4py import MPI
import difflib
import glob
import sys
import os




##############
### Inputs ###
##############
# Setup topofile and newfile (exactly like a LUNAR setup)
topofile = 'heating_dimer_ReaxFF_CHON-17_Weak_rep_1_del_2-50_temp_*K_time_*ps.data'
newfile = ':_bonds'

# Setup temp (K) and press (atm) settings. If press is set to 0, NVT will be
# used and if press is greater then 0, NPT will be used. It is advised to use
# the same settings as what was used during a pyrolysis run.
press = 0
temp = 300

# Boolean to check to see if the bond order file has already been generated. This
# will significantly speed up reprocessing a directory of partially filled with
# processed bond orders. Set as False to reprocess all files.
check_existances_of_bond_order_file = True



########################
### Helper functions ###
########################
# Function to get basename (using LUNAR's newfile method)
def get_basename(topofile, newfile, character=':'):
    root = topofile[:topofile.rfind('.')]
    character_check = newfile.strip() # avoid whitespace errors
    if character_check.startswith(character): # prefix
        basename = '{}{}'.format(root, newfile[1:])
    elif character_check.endswith(character): # suffix
        basename = '{}{}'.format(newfile[:-1], root)
    elif newfile in ['', ' ', '   ']: # single, double, or triple spacing
        basename = '{}'.format(root)
    else: # reset completely
        basename = '{}'.format(newfile)
    return basename

# Function to convert string to float or int and will default to string if needed
def string2digit(string):
    digit = string
    try: 
        digit = float(string)
        if digit.is_integer():
            digit = int(digit)
    except: pass
    return digit

# Function to get wildcards between a glob_string and a string glob matched with.
def get_glob_wildcards(glob_string, matched_string):
    # Force both strings to follow the same format rules
    glob_string = os.path.normpath(glob_string)
    matched_string = os.path.normpath(matched_string)
    
    # Find all strings that diff sets as adds "+"
    wildcards = []; add = False; tmp = ''
    for i, s in enumerate(difflib.ndiff(glob_string, matched_string)):
        if s[0] == ' ':
            add = False
        if s[0] == '-':
            if tmp != '': wildcards.append(tmp)
            add = False; tmp = ''
        if s[0] == '+':
            add = True 
        if add:
            tmp += s[-1]
    if tmp != '': wildcards.append(tmp)
            
    # Check to see if any wildcards show up in glob_string, if they do, remove them
    # from wildcards list. Additionally try converting from a string to float or int
    wildcards = [string2digit(i) for i in wildcards if i not in glob_string]
    return wildcards




#############################################
### Run PyLAMMPS to iterate through files ###
#############################################
if __name__ == "__main__": 

    # Generate a LAMMPS instance as lmp
    lmp = lammps() 
    
    # Start looping through files that glob() matched with
    lmp.command('print "Using array input option:"')
    files = sorted(glob.glob(topofile))
    for file in files:
        # Perform some checks to see if files need to be processed or not
        root = file[:file.rfind('.')]
        basename = get_basename(file, newfile)
        if newfile.endswith(':') and root.startswith(newfile[:-1]): # prefix
            lmp.command('print " - WARNING matched file {} already has newfile {} extension and was skipped"'.format(file, newfile))
            lmp.command('clear')
            continue
        elif newfile.startswith(':') and root.endswith(newfile[1:]): # suffix
            lmp.command('print " - WARNING matched file {} already has newfile {} extension and was skipped"'.format(file, newfile))
            lmp.command('clear')
            continue
        elif os.path.exists('{}.data'.format(basename)) and check_existances_of_bond_order_file:
            lmp.command('print " - WARNING basename {} already exists and was skipped. To avoid reprocessing"'.format(basename))
            lmp.command('clear')
            continue
        
        # Show how to access wildcards from match (incase temp needs to be generated from a difference
        # in wildcards to set temp based on the wildcard - which is useful for temperature ramping)
        wildcards = get_glob_wildcards(topofile, file)
        
        # We will just print this to the logfile to show the wildcards. So we need to convert everything to a string
        wildcard_string = ', '.join([str(i) for i in wildcards])
        lmp.command('print " - Matched wildards: [{}]"'.format(wildcard_string))
        
        
        #############################################################################################
        # I will demonstrate how to access the temp from the wildcards. Assume the following:       #
        #    topofile = 'heating_dimer_ReaxFF_CHON-17_Weak_rep_1_del_2-50_temp_*K_time_*ps.data'    #
        #    matched =  'heating_dimer_ReaxFF_CHON-17_Weak_rep_1_del_2-50_temp_300K_time_0ps.data'  #
        #                                                                                           #
        # The returned wildcards would be (NOTE that if the wildcard could be converted to an int   #
        # or a float it would have been via the string2digit() function call):                      #
        #    wildcards = [300, 0]                                                                   #
        #                                                                                           #
        # We know that the first wild card is the temp, so we could simply just generated a temp    # 
        # variable below as:                                                                        #
        #    temp = wildcards[0]                                                                    #
        #                                                                                           #
        # Then the temp will be dynamically set based on information from the string. This could be #
        # done for pressure as well or whatever else is needed.                                     #
        #############################################################################################
        
        
        # Standard LAMMPS script.
        #----------------------Variables----------------------
        lmp.command('variable        file string {}'.format(file))
        lmp.command('variable        myid string {}'.format(basename))
        
        #----------------------Initialization----------------------
        lmp.command('units           real')
        lmp.command('dimension       3')
        lmp.command('boundary        p p p')
        lmp.command('atom_style      charge')
        
        #----------------------ForceField----------------------
        lmp.command('read_data       ${file}')
        lmp.command('pair_style	     reaxff NULL  safezone 5.0 mincap 100')
        lmp.command('pair_coeff      * * CHON_2017_weak.reax C H')
        lmp.command('fix             charges all qeq/reaxff 1 0.0 10.0 1.0e-6 reax/c')
                   
        #----------------------Reset atomids-----------------
        lmp.command('reset_atoms id  sort yes')
        
        #----------------------Settings----------------------
        lmp.command('timestep        0.1')
        lmp.command('thermo          100')
        lmp.command('thermo_style    custom step temp press etotal ke pe ebond eangle edihed evdwl density lx ly lz')
        lmp.command('log             ${myid}.log.lammps')
    
        #----------------------Get Bonding info----------------------
        if press <= 0:
            lmp.command('fix         1 all nvt temp {} {} $(100*dt)'.format(temp, temp))
        else:
            lmp.command('fix          1 all npt temp {} {} $(100*dt) aniso {} {} $(1000*dt)'.format(temp, temp, press, press)) 
        lmp.command('fix             2 all reaxff/bonds 100 ${myid}.reaxff')
        lmp.command('run             1000 # 0.1 ps') 
        lmp.command('write_data      ${myid}.data')
        
        # Clear lammps for next iteration
        lmp.command('clear')
        

    
