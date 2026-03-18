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
# find and replace strings is: HPC_superior_mpi_LAMMPS_29Oct2020.py
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
#$ -q $<que>
#$ -pe mpichg $<ncores>
# Not an array simulation
#$ -M $<email>
#$ -m abes
# No dependent simulation
#$ -hard -l mem_free=2G
#$ -hard -l lammps_lic=$<inv_ncores>
# Uses traditional CPU
#
#$ -notify

# Load and list modules
source ${HOME}/.bashrc
module load mpi/impi/2018.3.222-iccifort-2018.3.222-GCC-7.3.0-2.30
# module load lammps/2020.10.29-CPU-stable-python
module list

# Input/Output files
INPUT_FOLDER="${PWD}"
INPUT_FILE="$<input>"
OUTPUT_FILE="$<output>"
ARRAY_TASK_ID=""
ADDITIONAL_OPTIONS=""

# Prepare to run the simulation
LAMMPS="/research/${USER}/lammps-29Oct20" 
export LD_LIBRARY_PATH="${LAMMPS}/src:${LD_LIBRARY_PATH}"
cd ${LAMMPS}/build
source myenv/bin/activate
cd ${INPUT_FOLDER}

# Run the simulation
mpirun -n ${NSLOTS} -machine ${TMP}/machines ${LAMMPS}/src/lmp_mpi -log ${INPUT_FOLDER}/${OUTPUT_FILE}${ARRAY_TASK_ID} -in ${INPUT_FOLDER}/${INPUT_FILE}${ARRAY_TASK_ID}

# List modules
module list
