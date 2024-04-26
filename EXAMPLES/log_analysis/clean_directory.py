# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 28, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

This script is for Josh Kemppainen ONLY! It will allow
him to quickly delete all *.txt *.csv and *.jpeg files
from this directory to reduce the size of LUNAR zip.
"""
##############################
# Import Necessary Libraries #
##############################
import os


##################
# Start cleaning #
##################
# get present working directory
pwd = os.getcwd() 
print('pwd =', pwd)

# Find files to delete
lunar_files = [file for file in os.listdir(pwd) if file.endswith('.lunar') and 'README' not in file]
csv_files = [file for file in os.listdir(pwd) if file.endswith('.csv')]
jpeg_files = [file for file in os.listdir(pwd) if file.endswith('.jpeg')]


# Start deleting files
print('\n\nStart deleting files')
files2delete = lunar_files +csv_files + jpeg_files
for file in files2delete:
    print('deleting : ', file)
    os.remove(file)