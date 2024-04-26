This src file folder contains all modularized codes to run the codes found ../ (a directory backwards in path).
The location of these files MATTER AND SHOULD NOT BE MOVED! Each code in the ../ directory has a corresponding
folder in /src (without the .py ext) were those codes will call from. Some folders such as atom_typing and 
auto_morse_bond have similar codes with the same name, but THEY HAVE SUBTLE DIFFERENCES in them depending on
what the code is meant for. Inside /src directory there are some python files that are called on by (one directory
up):
  - all2lmp.py
  - atom_typing.py
  - atom_removal.py
  - auto_cluster_analysis.py
  - auto_morse_bond.py
  - bond_react_merge.py
  - bond_react_merge_prep.py
  - cell_builder.py
  - cluster_analysis.py
  - convert2graphite.py
  - lmp2SYBYLmol2.py
  - free_volume.py
The files that are in /src have been generalized enough to be called by all code codes listed above. Future 
work will involve generalizing more modules and sticking them in the /src directory. 