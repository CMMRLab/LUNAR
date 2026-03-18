#!/bin/bash
#
# This script is meant for automating the batch 
# script submission process and has been written
# by: LUNAR/array_lmp_script.py v1.2 / 18 March 2026
#
# Usage:
#   bash auto_submit_atomic_fraction_XRD.sh
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

sims2submit+=("atomic_fraction_XRD_1.sh")
dirs2submit+=("Graphite_Argon/Replicate_1/AF_0.0_Density_0.125")

sims2submit+=("atomic_fraction_XRD_1.sh")
dirs2submit+=("Graphite_Argon/Replicate_1/AF_0.0_Density_0.25")

sims2submit+=("atomic_fraction_XRD_1.sh")
dirs2submit+=("Graphite_Argon/Replicate_1/AF_0.5_Density_0.125")

sims2submit+=("atomic_fraction_XRD_1.sh")
dirs2submit+=("Graphite_Argon/Replicate_1/AF_0.5_Density_0.25")

sims2submit+=("atomic_fraction_XRD_1.sh")
dirs2submit+=("Graphite_Argon/Replicate_1/AF_1.0_Density_0.125")

sims2submit+=("atomic_fraction_XRD_1.sh")
dirs2submit+=("Graphite_Argon/Replicate_1/AF_1.0_Density_0.25")

sims2submit+=("atomic_fraction_XRD_2.sh")
dirs2submit+=("Graphite_Argon/Replicate_2/AF_0.0_Density_0.125")

sims2submit+=("atomic_fraction_XRD_2.sh")
dirs2submit+=("Graphite_Argon/Replicate_2/AF_0.0_Density_0.25")

sims2submit+=("atomic_fraction_XRD_2.sh")
dirs2submit+=("Graphite_Argon/Replicate_2/AF_0.5_Density_0.125")

sims2submit+=("atomic_fraction_XRD_2.sh")
dirs2submit+=("Graphite_Argon/Replicate_2/AF_0.5_Density_0.25")

sims2submit+=("atomic_fraction_XRD_2.sh")
dirs2submit+=("Graphite_Argon/Replicate_2/AF_1.0_Density_0.125")

sims2submit+=("atomic_fraction_XRD_2.sh")
dirs2submit+=("Graphite_Argon/Replicate_2/AF_1.0_Density_0.25")


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
  qsub ${sim}
  cd "$pwd"  
  
done

echo ""
echo "All done!"
