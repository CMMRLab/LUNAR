This folder contains an example for .pbd/packmol integration into the LUNAR codes. File describitions:
  - water.inp, simple packmol script to read in water.pdb
  - water.pdb, generated using ChemDraw 3D and atom name(s) where switched to PCFF atom types
  - water_packmol_atom_name_as_atom_type.pdb, output file from packmol, where packmol maintains the atom types
  - water_packmol_atom_name_as_atom_type_typed.data file where atom_typing.py was used to retype every atom in the .pdb file
  - water_packmol_atom_name_as_atom_type_typed.nta  file where atom_typing.py set the new atom types for every atom in the .pdb file
  - water_packmol_atom_name_as_atom_type_typed_IFF.data
  
atom_typing.py settings used:
	topofile = 'EXAMPLES/packmol_pdb_methods/atom_typing_from_pdb/water_packmol_atom_name_as_atom_type.pdb'
	ff_name = 'PCFF-IFF'
	newfile = 'typed'
	parent_directory = 'topofile'
	
all2lmp.py settings used:
	topofile = 'EXAMPLES/packmol_pdb_methods/atom_typing_from_pdb/water_packmol_atom_name_as_atom_type_typed.data'
	nta_file = 'EXAMPLES/packmol_pdb_methods/atom_typing_from_pdb/water_packmol_atom_name_as_atom_type_typed.nta'
	frc_file = 'frc_files/pcff_iff_v1_5_CNT_poly_solv_Hpan_mod.frc'
	assumed = 'frc_files/general_assumed_equivs.coeffs'
	ff_class = 2
	parent_directory = 'topofile'
	
This method is useful because if atom_typing.py has automatic atom typing set up for the force field you wish to
work with, because typically if you are using packmol you will be generating large systems, where manual atom
typing is impossible. Then you can rely on atom_typing.py to type your entire system from a .pdb file. Changing
atom types can be a bit difficult here if you want to try different atom types then what atom_typing.py can provide.