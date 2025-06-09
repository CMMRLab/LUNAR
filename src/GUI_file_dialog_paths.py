# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
June 4, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

## These are the defaults that NEED to be checked before each LUNAR send out
# xrdfile = os.path.normpath(os.path.join(script_path, '../'))
# logfile = os.path.normpath(os.path.join(script_path, '../'))
# topofile = os.path.normpath(os.path.join(script_path, '../'))
# nta_file = os.path.normpath(os.path.join(script_path, '../'))
"""
##############################
# Import Necessary Libraries #
##############################
import os

# This scripts path and name
script_path = os.path.dirname(os.path.realpath(__file__))
script_name = os.path.basename(__file__)


# LUNAR related file paths and names. The default will place the dialog
# box in LUNAR's "main" directory by using the following sequence:
#  file_path = os.path.normpath(os.path.join(script_path, '../'))
xrdfile = os.path.normpath(os.path.join(script_path, '../'))
logfile = os.path.normpath(os.path.join(script_path, '../'))
topofile = os.path.normpath(os.path.join(script_path, '../'))



