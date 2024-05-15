#! /bin/bash
#
# USAGE: bash run_all.sh
#
# ARGUMENTS:
# None

# Set current directory
DIR="automatic_tests/outputs"
SECONDS=0


# print intialization
echo "*******************************"
echo "*******************************"
echo "***  LUNAR automated tests  ***"
echo "*******************************"
echo "*******************************"
echo "DIR = $DIR"


# run LUNAR/all2lmp.py tests
bash all2lmp.sh $DIR/"all2lmp"


# Set run time of script
ELAPSED="Elapsed time during tests: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo "$ELAPSED"