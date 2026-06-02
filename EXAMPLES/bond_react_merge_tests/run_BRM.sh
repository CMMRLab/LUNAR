#! /bin/bash
# 
# Path2LUNAR
CURRENT=$PWD
LUNAR="../../"

# Change to LUNAR
cd $LUNAR
pwd

# Run all BRM tests
start=1
end=17
for ((i=start; i<=end; i++)); do
    infile="EXAMPLES/bond_react_merge_tests/test$i/merge_files.txt"
    echo "Running BRM: $infile"
	python3 bond_react_merge.py  -files  infile:$infile -atomstyle full -map T -wrd T -wrm T -wmf T -tl T
done

# Wrap-up
cd $CURRENT

# Run the check mapping functions
python3 check_maps.py

# Clean-up run directories
echo ""
echo "Cleaning run directories:"
for ((i=start; i<=end; i++)); do
    rm_dir="test$i/run"
    echo "cleaning: $rm_dir"
	rm -rf $rm_dir
done
echo
echo "All done!"