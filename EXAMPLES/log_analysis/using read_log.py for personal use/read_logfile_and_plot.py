# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 28, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    ****************************************************************
    * Requirements:                                                *
    *   python 3.7+                                                *
    *                                                              *
    * Dependencies:                                                *
    *   python matplotlib module:                                  *
    *    - pip3 install matplotlib (if pip manager is installed)   *
    ****************************************************************
"""
##############################
# Import Necessary Libraries #
##############################
import read_log as read_log
import matplotlib.pyplot as plt


##########
# Inputs #
##########
# logfile to read in (using file from a directory backwards)
logfile = '../property=density_ts=0.5.log.lammps'

# Set some column keywords to find sections in logfile with thermo data.
# *NOTE: You minimally need one keyword where 'Step' is a good choice.*
keywords = ['Step']

# Set the sections of the logfile to get data from. Sections are counted
# starting from 1 and end at N-sections. The following options are available
# in example format:
#    'all'    get data for all sections from logfile
#    '1'      get data for only section 1 from log file
#    '1,3'    get data for sections 1 and 3 from log file
#    '1-3'    get data for sections 1, 2, and 3 from log file
#    '1-3,5'  get data for sections 1, 3, 3, and 5 from log file
# *NOTE if reading multiple sections from the log file and a column
#  if missing in a certain section, the get_data function will create
#  zeros.*
sections = 'all'


####################################################
# Read logfile, show how to access some attributes #
# of logfile, and get data from logfile to plot    #
####################################################
# Read logfile
log = read_log.file(logfile, keywords=keywords)
print('\n\nAccess some logfile attributes:')
print(' - log file sections:          ', log.sections)
print(' - number of lines in log file:', log.nlines)
print(' - Show how to access data directly from log class:')
print('     log.data[sectionID][column-name] to access section')
print('     and data from section. *NOTE: it is easier to use')
print('     the get data function shown below.*')
print('         log.data[1]["Step"]    = ', log.data[1]['Step'])
print('         log.data[1]["Density"] = ', log.data[1]['Density'])

# get data based on log class and sections string. pflag is True
# or False and will print or not print sections being loaded.
print('\n\nGet data from log class via get_data function:')
data = log.get_data(sections, remove_duplicates=True, pflag=True) # {column-name:[lst of data]}
print('\nAccess some data values via keys:')
print('   - truncated data["Step"] = ', data['Step'][:6])
print('   - truncated data["Temp"] = ', data['Temp'][:6])


##########################
# Generate a simple plot #
##########################
# Create new plot
fig, ax = plt.subplots()
plt.plot(data['Step'], data['Density'], '.', mfc='white', ms=5, markeredgecolor='tab:blue', lw=0.01, label='LAMMPS density data')
ax.set_xlabel('Step', fontsize=12)
ax.set_ylabel('Density ($g/cm^3$)', fontsize=12)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, ncol=3, fontsize=12)