#! /bin/bash
# 
# This is a custom batch script to usage with LUNAR/array_lmp_script.py for
# automatic array generation of LAMMPS scripts and batch scripts. Currently
# there are the following find strings set in this script that should have
# a match in the batch_array section:
#   $<ncores>
#   $<inv_ncores>
#   $<que>
#   $<email>
#   $<input>
#   $<output>
# In LUNAR/src/array_lmp_script/modes the template batch script for these
# find and replace strings is: HPC_superior_PyLAMMPS_parallel_30June2021.py
# 
# This file contains critical settings for appropriate submission and execution
# of the simulation. Editing this file [and/or its bulk (re-)generation]
# without explicit permission of (or discussion with) the administrators can
# lead to improper use of resources, extended wait times in the queue, etc.,
# and will be grounds for removing your account from the HPC infrastructure.
#
# Refer to HPC 101 Training Camp v2.0 (https://mtu.instructure.com/courses/1374830)
# for additional information.
 
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q cmmr.q
#$ -pe mpichg 16
# Not an array simulation
#$ -M jdkemppa@mtu.edu
#$ -m abes
# No dependent simulation
#$ -hard -l mem_free=2G
#$ -hard -l lammps_lic=.0625000000
# Uses traditional CPU
#
#$ -notify

# Load and list modules
source ${HOME}/.bashrc
module load mpi/impi/2018.3.222-iccifort-2018.3.222-GCC-7.3.0-2.30
module list

# Input/Output files
INPUT_FOLDER="${PWD}"
INPUT_FILE="$<input>"
OUTPUT_FILE="$<output>"
ARRAY_TASK_ID=""
ADDITIONAL_OPTIONS=""

# Prepare to run the simulation
cd /research/${USER}/lammps/build
source myenv/bin/activate

# Run the simulation
cd ${INPUT_FOLDER}
mpirun -np ${NSLOTS} python3 ${INPUT_FILE}${ARRAY_TASK_ID}

# List modules
module list
