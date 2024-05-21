#! /bin/bash
#
# USAGE: bash run_all.sh
#
# ARGUMENTS:
# None

# Set current directory
DIR_FROM_LUNAR="automatic_tests"
DIR_IN_TESTS="outputs"
DIR=$DIR_FROM_LUNAR/$DIR_IN_TESTS
SECONDS=0


# print intialization
echo "*******************************"
echo "*******************************"
echo "***  LUNAR automated tests  ***"
echo "*******************************"
echo "*******************************"
echo "DIR = $DIR"
mkdir $DIR_IN_TESTS


# run LUNAR/all2lmp.py tests
bash all2lmp.sh $DIR/all2lmp | tee $DIR_IN_TESTS/all2lmp_TEST_LOG.txt

# run LUNAR/atom_typing.py tests
bash atom_typing.sh $DIR/atom_typing | tee $DIR_IN_TESTS/atom_typing_TEST_LOG.txt


# Set run time of script
ELAPSED="Elapsed time during tests: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo ""
echo "$ELAPSED"