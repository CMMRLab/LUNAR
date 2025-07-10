# -*- coding: utf-8 -*-
"""
This script is for Josh to be able
to have a few quick tests that he 
can perform by copying and pasting
the variables defined below into
LUNAR/sheet_builder.py, so that
he can test the development of the
grafting_files = '' option 
"""


#use_GUI = False
run_mode = 'sheet'
run_mode = 'sheet'
seed = 13

length_in_edgetype = 12.0
length_in_perpendicular = 12.0


sheet_nlayers = 3
types = {1: 'C', 2: 'C', 3: 'C', 4: 'C'}
functional_atoms = 'C<25,-,1>|O;   C<25,+,3>|H;'



sheet_nlayers = 1
types = {1: 'C', 2: 'C', 3: 'C', 4: 'C'}
functional_atoms = 'C<30,+,1>|O'


types = {1: 'C', 2: 'C', 3: 'C', 4: 'C'}
length_in_edgetype = 30.0
length_in_perpendicular = 30.0
sheet_nlayers = 3
grafting_files = 'C<7,+,3><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol; C<7,-,1><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol'


# types = {1:'cg1|cge', 2:'cg1|cge', 3:'cg1|cge', 4:'cg1|cge'}
# length_in_edgetype = 30.0
# length_in_perpendicular = 30.0
# sheet_nlayers = 3
# grafting_files = 'cg1<7,+,3><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol; cg1<7,-,1><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol'


# Test: 1
use_GUI = False
run_mode = 'sheet'
length_in_edgetype = 30.0
length_in_perpendicular = 30.0
sheet_nlayers = 1
types = {1: 'C', 2: 'C', 3: 'C', 4: 'C'}
grafting_files = ''
functional_atoms = 'C<30,+,1>|O|'


# Test: 2
use_GUI = False
run_mode = 'sheet'
length_in_edgetype = 30.0
length_in_perpendicular = 30.0
sheet_nlayers = 3
grafting_files = ''
types = {1:'cg1|cge', 2:'cg1|cge', 3:'cg1|cge', 4:'cg1|cge'}
functional_atoms = 'cg1<25,-,1>|O|;   cg1<25,+,3>|O|H;'


# Test: 3
use_GUI = False
run_mode = 'sheet'
# plane = 'xy'
# plane = 'yz'
length_in_edgetype = 30.0
length_in_perpendicular = 30.0
sheet_nlayers = 3
types = {1: 'C', 2: 'C', 3: 'C', 4: 'C'}
grafting_files = 'C<1,+,3><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol; C<1,-,1><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol'
minimum_distance = 4.0
functional_atoms = 'C<30,+,2>|O|'


# Test: 4
use_GUI = False
run_mode = 'sheet'
length_in_edgetype = 30.0
length_in_perpendicular = 30.0
sheet_nlayers = 3
types = {1:'cg1|cge', 2:'cg1|cge', 3:'cg1|cge', 4:'cg1|cge'}
grafting_files = 'cg1<7,+,3><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol; cg1<7,-,1><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol'
minimum_distance = 4.0
functional_atoms = ''


# Test 5: Ringed test
use_GUI = False
run_mode = 'sheet'
plane = 'xy'
length_in_edgetype = 30.0
length_in_perpendicular = 30.0
sheet_nlayers = 3
types = {1: 'C1', 2: 'C1', 3: 'C1', 4: 'C1'}
grafting_files = 'C1<5,-,1><3,4>|EXAMPLES/sheet_builder/ring_graft.3.4.mol; C1<5,+,3><3,4>|EXAMPLES/sheet_builder/ring_graft.3.4.mol'
minimum_distance = 8.0 #'cylinder'
functional_atoms = ''


# Test 6: CNT
use_GUI = False
symmetric_tube_axis = 'z'
run_mode = 'symmetric-tube'
symmetric_length = 25.0
diameter = 12.0
symmetric_ntubes = 3
types = {1: 'C1', 2: 'C1', 3: 'C1', 4: 'C1'}
grafting_files = 'C1<7,+,3><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol'
grafting_files = 'C1<5,+,3><3,4>|EXAMPLES/sheet_builder/ring_graft.3.4.mol'
distance = 8.0
functional_atoms = ''


# Test: 7 hydroxyl
use_GUI = False
run_mode = 'sheet'
seed = 13
plane = 'xy'
#plane = 'yz'
#plane = 'xz'
length_in_edgetype = 30.0
length_in_perpendicular = 30.0
sheet_nlayers = 3
types = {1: 'C', 2: 'C', 3: 'C', 4: 'C'}
grafting_files = 'C<7,+,3><1>|EXAMPLES/sheet_builder/hydroxyl.1.mol2; C<7,-,1><1>|EXAMPLES/sheet_builder/hydroxyl.1.mol2'
minimum_distance = 8.0
functional_atoms = ''


# # Test: 8 lmp data w/o type labels
use_GUI = True
run_mode = 'sheet'
seed = 13
plane = 'xy'
#plane = 'yz'
#plane = 'xz'
length_in_edgetype = 30.0
length_in_perpendicular = 30.0
sheet_nlayers = 3
types = {1: 'C', 2: 'C', 3: 'C', 4: 'C'}
grafting_files = 'C<7,+,3><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol; C<7,-,1><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol'
grafting_files = 'C<7,+,3><15>|EXAMPLES/sheet_builder/PBZ_graft_w_type_labels.15.data; C<7,-,1><15>|EXAMPLES/sheet_builder/PBZ_graft_wo_type_labels.15.data'
minimum_distance = 8.0
functional_atoms = ''


# Test 9: Chiral CNT
use_GUI = False
n = 5
m = 10
chiral_length = 20.0
chiral_tube_axis = 'z'
run_mode = 'chiral-tube'
types = {1: 'C1', 2: 'C1', 3: 'C1', 4: 'C1'}
grafting_files = 'C1<7,+,1><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol'
grafting_files = 'C1<5,+,1><3,4>|EXAMPLES/sheet_builder/ring_graft.3.4.mol'
distance = 0.0

# grafting_files = ''
# functional_atoms = 'C<30,+,*>|O|'

# grafting_files = ''
# functional_atoms = ''



# Sagars setup
run_mode = 'sheet'
plane = 'xy'
#plane = 'yz'
#plane = 'xz'
length_in_edgetype = 30.0
length_in_perpendicular = 30.0
sheet_nlayers = 3
types = {1: 'C', 2: 'C', 3: 'C', 4: 'C'}
grafting_files = 'C<7,+,3><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol; C<7,-,1><15>|EXAMPLES/sheet_builder/PBZ_graft.15.mol'
minimum_distance = 8.0
functional_atoms = ''