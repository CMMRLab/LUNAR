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
  echo "USAGE: bash all2lmp.sh <DIR>"
  exit 1
fi

# Set current directory
PWD_DIR=$PWD
DIR=$1

# print intialization
echo ""
echo "************************************"
echo "* LUNAR/all2lmp.py automated tests *"
echo "************************************"

# change to LUNAR directory
cd ../

########################
# Run all2lmp.py tests #
########################
echo "Test1: class2 w/PCFF and all CLI options"
FF="frc_files/pcff.frc"
NEWFILE=":_PCFF_test"
ASM="frc_files/general_assumed_equivs.coeffs"
TOPO="EXAMPLES/EPON_862/atom_typing_Outputs/detda_typed.data"
NTA="EXAMPLES/EPON_862/atom_typing_Outputs/detda_typed.nta"
python3 all2lmp.py -dir $DIR -class 2 -frc $FF -newfile $NEWFILE -topo $TOPO -nta $NTA -ignore T \
                   -type-labels T -reset-charges T -atomstyle full -auto-equivs T -asm $ASM -assumed-equivs F \
				   -reset-molids T -write-comments T -write-bond-react T -morse-bond F -rx 45 -ry 45 -rz 45 \
				   -sx 10 -sy 20 -sz 30 -add2box 5

echo "Test2: class1 w/CVFF, MSI outputs, triclinic cell, and all CLI options"
FF="frc_files/cvff.frc"
NEWFILE=":_CVFF_test"
ASM="frc_files/general_assumed_equivs.coeffs"
TOPO="EXAMPLES/all2lmp_car_mdf/cnt-hexagonal-class1.mdf"
NTA="EXAMPLES/all2lmp_car_mdf/cnt-hexagonal-class1.car"
python3 all2lmp.py -dir $DIR -class 1 -frc $FF -newfile $NEWFILE -topo $TOPO -nta $NTA -ignore T \
                   -type-labels T -reset-charges T -atomstyle full -auto-equivs T -asm $ASM -assumed-equivs F \
				   -reset-molids T -write-comments T -write-bond-react F -morse-bond F -rx 0 -ry 0 -rz 0 \
				   -sx 0 -sy 0 -sz 0 -add2box 0
				   
echo "Test3: class0 w/OPLS-AA, MSI outputs, orthangol cell, and all CLI options"
FF="frc_files/oplsaa.frc"
NEWFILE=":_OPLS-AA_test"
ASM="frc_files/general_assumed_equivs.coeffs"
TOPO="EXAMPLES/all2lmp_car_mdf/ethane-oplsaa.mdf"
NTA="EXAMPLES/all2lmp_car_mdf/ethane-oplsaa.car"
python3 all2lmp.py -dir $DIR -class 0 -frc $FF -newfile $NEWFILE -topo $TOPO -nta $NTA -ignore T \
                   -type-labels T -reset-charges T -atomstyle full -auto-equivs T -asm $ASM -assumed-equivs F \
				   -reset-molids T -write-comments T -write-bond-react F -morse-bond F -rx 0 -ry 0 -rz 0 \
				   -sx 0 -sy 0 -sz 0 -add2box 0
				   
echo "Test4: classr for ReaxFF and all CLI options"
FF="frc_files/all2lmp_reaxff.frc"
NEWFILE=":_ReaxFF_test"
ASM="frc_files/general_assumed_equivs.coeffs"
TOPO="EXAMPLES/EPON_862/atom_typing_Inputs/detda.mol"
NTA="empty.nta"
python3 all2lmp.py -dir $DIR -class r -frc $FF -newfile $NEWFILE -topo $TOPO -nta $NTA -ignore T \
                   -type-labels F -reset-charges F -atomstyle charge -auto-equivs T -asm $ASM -assumed-equivs F \
				   -reset-molids F -write-comments T -write-bond-react T -morse-bond F -rx 0 -ry 0 -rz 0 \
				   -sx 0 -sy 0 -sz 0 -add2box 0
				   
echo "Test5: classd for DREIDING-q, calling atom_typing, and all CLI options"
FF="frc_files/all2lmp_dreiding.frc"
NEWFILE=":_DREIDING_test"
ASM="frc_files/general_assumed_equivs.coeffs"
TOPO="EXAMPLES/EPON_862/atom_typing_Inputs/detda.mol"
NTA="DREIDING-q"
python3 all2lmp.py -dir $DIR -class d -frc $FF -newfile $NEWFILE -topo $TOPO -nta $NTA -ignore T \
                   -type-labels T -reset-charges F -atomstyle full -auto-equivs T -asm $ASM -assumed-equivs F \
				   -reset-molids F -write-comments T -write-bond-react T -morse-bond F -rx 0 -ry 0 -rz 0 \
				   -sx 0 -sy 0 -sz 0 -add2box 0
				   
echo "Test6: classd for DREIDING, calling atom_typing, and all CLI options"
FF="frc_files/all2lmp_dreiding.frc"
NEWFILE=":_DREIDING_test"
ASM="frc_files/general_assumed_equivs.coeffs"
TOPO="EXAMPLES/EPON_862/atom_typing_Inputs/dgebf.mol"
NTA="DREIDING"
python3 all2lmp.py -dir $DIR -class d -frc $FF -newfile $NEWFILE -topo $TOPO -nta $NTA -ignore T \
                   -type-labels F -reset-charges T -atomstyle full -auto-equivs T -asm $ASM -assumed-equivs F \
				   -reset-molids F -write-comments T -write-bond-react T -morse-bond F -rx 0 -ry 0 -rz 0 \
				   -sx 0 -sy 0 -sz 0 -add2box 0
				   
echo "Test7: classs1 for Skeleton1 and all CLI options"
FF="empty.frc"
NEWFILE=":_Skeleton1_test"
ASM="frc_files/general_assumed_equivs.coeffs"
TOPO="EXAMPLES/EPON_862/atom_typing_Outputs/detda_typed.data"
NTA="EXAMPLES/EPON_862/atom_typing_Outputs/detda_typed.nta"
python3 all2lmp.py -dir $DIR -class s1 -frc $FF -newfile $NEWFILE -topo $TOPO -nta $NTA -ignore T \
                   -type-labels F -reset-charges F -atomstyle charge -auto-equivs T -asm $ASM -assumed-equivs F \
				   -reset-molids F -write-comments T -write-bond-react T -morse-bond F -rx 0 -ry 0 -rz 0 \
				   -sx 0 -sy 0 -sz 0 -add2box 0
				   
echo "Test8: classs2 for Skeleton2 and all CLI options"
FF="empty.frc"
NEWFILE=":_Skeleton2_test"
ASM="frc_files/general_assumed_equivs.coeffs"
TOPO="EXAMPLES/EPON_862/atom_typing_Outputs/detda_typed.data"
NTA="EXAMPLES/EPON_862/atom_typing_Outputs/detda_typed.nta"
python3 all2lmp.py -dir $DIR -class s2 -frc $FF -newfile $NEWFILE -topo $TOPO -nta $NTA -ignore T \
                   -type-labels F -reset-charges F -atomstyle charge -auto-equivs T -asm $ASM -assumed-equivs F \
				   -reset-molids F -write-comments T -write-bond-react T -morse-bond F -rx 0 -ry 0 -rz 0 \
				   -sx 0 -sy 0 -sz 0 -add2box 0
				   
echo "Test9: class2 w/PCFF and equivs options and all CLI options"
FF="frc_files/pcff.frc"
NEWFILE=":_FF_applied"
ASM="frc_files/general_assumed_equivs.coeffs"
TOPO="EXAMPLES/atom_typing_general/detda_general.data"
NTA="EXAMPLES/atom_typing_general/detda_general_filled_in.nta"
python3 all2lmp.py -dir $DIR -class 2 -frc $FF -newfile $NEWFILE -topo $TOPO -nta $NTA -ignore T \
                   -type-labels T -reset-charges T -atomstyle full -auto-equivs T -asm $ASM -assumed-equivs F \
				   -reset-molids T -write-comments T -write-bond-react T -morse-bond F -rx 45 -ry 45 -rz 45 \
				   -sx 10 -sy 20 -sz 30 -add2box 5

# change back to testing directory
cd $PWD_DIR