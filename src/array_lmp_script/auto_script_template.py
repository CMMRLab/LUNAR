# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
October 8, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.io_functions as io_functions
import os


# Set the first part of the bash script
part1 = """#!/bin/bash
#
# This script is meant for automating the batch 
# script submission process and has been written
# by: LUNAR/array_lmp_script.py $<version>
#
# Usage:
#   bash $<build_script>
#
# Set current working directory
pwd=$PWD


# Set the root directory where this script will be run from. Here is a few examples:
#   "."   -> Run this script relative to where array_lmp_script.py  wrote this script
#   "../" -> Run this script one directory backwards relative to where array_lmp_script.py wrote this script
#   "New" -> Run this script one directory forwards ("New") relative to where array_lmp_script.py wrote this script
root_dir="."


# Set the directories where batch scripts are and the batch scripts. NOTE: This will be done using the appending
# method for easier deleting/commenting if the user desires.
dirs2submit=()
sims2submit=()
"""

part2 = """

# Start iterating through batch dirs and scripts to submit			
for i in "${!dirs2submit[@]}"; do
  dir=${dirs2submit[$i]}
  sim=${sims2submit[$i]}
  
  # Check if user supplied a root directory with
  # an ending backslach and adjust accordingly
  if [[ $root_dir == */ ]]; then
    dir="${root_dir}${dir}"
  else
    dir="${root_dir}/${dir}"
  fi

  echo ""
  echo "Iteration: $i"
  echo "Directory: ${dir}"
  echo "Batch script: ${sim}"
  echo ""
  
  # Change to sim-dir and submit; then change
  # back to pwd for the next iteration
  cd "$dir"
  $<submit> ${sim}
  cd "$pwd"  
  
done

echo ""
echo "All done!"
"""

# Write the bash script
def write(filename, auto_submit, version, submit='qsub'):    
    with open(filename, 'w', encoding='us-ascii', newline='\n') as f:
        parta = part1.replace('$<build_script>', os.path.basename(filename))
        parta = parta.replace('$<version>', version)
        f.write(parta)
        
        for submit_dir, submit_sim in auto_submit:
            f.write('\n')
            f.write('sims2submit+=("{}")\n'.format(io_functions.path_to_string(submit_sim)))
            f.write('dirs2submit+=("{}")\n'.format(io_functions.path_to_string(submit_dir)))
        
        partb = part2.replace('$<submit>', submit)
        f.write(partb)

    return