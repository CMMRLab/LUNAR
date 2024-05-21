#! /bin/bash
#
# USAGE: bash all2lmp.sh <DIR>
#
# ARGUMENTS:
# 1. <DIR> 
#      LUNAR's parent directory variable to
#      write all files to that location

# Check if correct number of arguments
if [ "$#" -lt 1 ]; then
  echo "Incorrect number of arguments"
  echo "USAGE: bash atom_typing.sh <DIR>"
  exit 1
fi

# Set current directory
PWD_DIR=$PWD
DIR=$1

# Set python call (i.e. "python3 filename.py" or "python filename.py" or "py filename.py" ...)
PYTHON="python3"

# print intialization
echo ""
echo "****************************************"
echo "* LUNAR/atom_typing.py automated tests *"
echo "****************************************"

# change to LUNAR directory
cd ../

############################
# Run atom_typing.py tests #
############################
echo "Test1: PCFF-IFF, .mol file, and all CLI options"
FF="frc_files/pcff.frc"
NEWFILE=":_PCFF-IFF_test"
FF="PCFF-IFF"
TOPO="EXAMPLES/EPON_862/atom_typing_Inputs/detda.mol"
BONDFILE="n.u."
CHARGEFILE="frc_files/Gasteiger_parameters.txt"
$PYTHON atom_typing.py -dir $DIR -topo $TOPO -bond BONDFILE -charge-file $CHARGEFILE -newfile $NEWFILE -ff $FF \
                       -nta-comments F -reset-charges F -vdw-scale 1.1 -pdb skip -boundary f-f-f -bond-reset F
					   
echo "Test2: DREIDING, .smiles file, Gasteiger charging, and all CLI options"
FF="DREIDING"
NEWFILE="DETDA_SMILES_DREIDING_test"
FF="PCFF-IFF"
TOPO="CCC1=CC(=C(C(=C1N)C)N)CC.smiles"
BONDFILE="n.u."
CHARGEFILE="frc_files/Gasteiger_parameters.txt"
$PYTHON atom_typing.py -dir $DIR -topo $TOPO -bond BONDFILE -charge-file $CHARGEFILE -newfile $NEWFILE -ff $FF \
                       -nta-comments F -reset-charges T -vdw-scale 1.1 -pdb skip -boundary f-f-f -bond-reset F
					   
echo "Test3: PCFF, .mol2 file, bond finding, and all CLI options"
FF="PCFF"
NEWFILE=":_PCFF_test"
FF="PCFF"
TOPO="EXAMPLES/EPON_862/atom_typing_Inputs/dgebf.mol2"
BONDFILE="n.u."
CHARGEFILE="frc_files/Gasteiger_parameters.txt"
$PYTHON atom_typing.py -dir $DIR -topo $TOPO -bond BONDFILE -charge-file $CHARGEFILE -newfile $NEWFILE -ff $FF \
                       -nta-comments F -reset-charges F -vdw-scale 1.1 -pdb skip -boundary f-f-f -bond-reset T
					   
echo "Test4: PCFF, .pdb file, and all CLI options"
FF="PCFF"
NEWFILE=":_COMPASS_test"
FF="compass"
TOPO="EXAMPLES/EPON_862/atom_typing_Inputs/dgebf.pdb"
BONDFILE="n.u."
CHARGEFILE="frc_files/Gasteiger_parameters.txt"
$PYTHON atom_typing.py -dir $DIR -topo $TOPO -bond BONDFILE -charge-file $CHARGEFILE -newfile $NEWFILE -ff $FF \
                       -nta-comments F -reset-charges F -vdw-scale 1.1 -pdb skip -boundary f-f-f -bond-reset F
					   
echo "Test5: CVFF-IFF, .data file, and all CLI options"
FF="PCFF"
NEWFILE=":_CVFF-IFF_test"
FF="CVFF-IFF"
TOPO="EXAMPLES/EPON_862/atom_typing_Inputs/dgebf.data"
BONDFILE="n.u."
CHARGEFILE="frc_files/Gasteiger_parameters.txt"
$PYTHON atom_typing.py -dir $DIR -topo $TOPO -bond BONDFILE -charge-file $CHARGEFILE -newfile $NEWFILE -ff $FF \
                       -nta-comments F -reset-charges F -vdw-scale 1.1 -pdb skip -boundary f-f-f -bond-reset F
					   
echo "Test6: CVFF, .data file, writing .pdb file and all CLI options"
FF="PCFF"
NEWFILE=":_CVFF_test"
FF="CVFF"
TOPO="EXAMPLES/EPON_862/atom_typing_Inputs/dgebf.data"
BONDFILE="n.u."
CHARGEFILE="frc_files/Gasteiger_parameters.txt"
$PYTHON atom_typing.py -dir $DIR -topo $TOPO -bond BONDFILE -charge-file $CHARGEFILE -newfile $NEWFILE -ff $FF \
                       -nta-comments F -reset-charges F -vdw-scale 1.1 -pdb types -boundary f-f-f -bond-reset F
					   
echo "Test6: Clay-FF, .data file, writing .pdb file, bond finding, and all CLI options"
FF="PCFF"
NEWFILE=":_Clay-FF_test"
FF="Clay-FF"
TOPO="EXAMPLES/EPON_862/atom_typing_Inputs/dgebf.data"
BONDFILE="n.u."
CHARGEFILE="frc_files/Gasteiger_parameters.txt"
$PYTHON atom_typing.py -dir $DIR -topo $TOPO -bond BONDFILE -charge-file $CHARGEFILE -newfile $NEWFILE -ff $FF \
                       -nta-comments F -reset-charges F -vdw-scale 1.1 -pdb skip -boundary f-f-f -bond-reset T
					   
echo "Test7: OPLS-AA, .data file, writing .pdb file, and all CLI options"
FF="PCFF"
NEWFILE=":_OPLS-AA_test"
FF="OPLS-AA"
TOPO="EXAMPLES/EPON_862/atom_typing_Inputs/dgebf.data"
BONDFILE="n.u."
CHARGEFILE="frc_files/Gasteiger_parameters.txt"
$PYTHON atom_typing.py -dir $DIR -topo $TOPO -bond BONDFILE -charge-file $CHARGEFILE -newfile $NEWFILE -ff $FF \
                       -nta-comments F -reset-charges F -vdw-scale 1.1 -pdb skip -boundary f-f-f -bond-reset F
					   
echo "Test8: general:2, .data file, writing .pdb file, and all CLI options"
FF="PCFF"
NEWFILE=":_general2_test"
FF="general:2"
TOPO="EXAMPLES/EPON_862/atom_typing_Inputs/dgebf.data"
BONDFILE="n.u."
CHARGEFILE="frc_files/Gasteiger_parameters.txt"
$PYTHON atom_typing.py -dir $DIR -topo $TOPO -bond BONDFILE -charge-file $CHARGEFILE -newfile $NEWFILE -ff $FF \
                       -nta-comments F -reset-charges F -vdw-scale 1.1 -pdb skip -boundary f-f-f -bond-reset F



# change back to testing directory
cd $PWD_DIR