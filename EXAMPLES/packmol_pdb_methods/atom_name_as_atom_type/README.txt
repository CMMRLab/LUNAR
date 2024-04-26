This folder contains an example for .pbd/packmol integration into the LUNAR codes. File describitions:
  - water.inp, simple packmol script to read in water.pdb
  - water.pdb, generated using ChemDraw 3D and atom name(s) where switched to PCFF atom types
  - water_packmol_atom_name_as_atom_type.pdb, output file from packmol, where packmol maintains the atom types
  - water_packmol_atom_name_as_atom_type_IFF.data, output from all2lmp.py to generate a LAMMPS datafile using the
    atom types that are set in the atom name column from the water_packmol_atom_name_as_atom_type.pdb.
	
all2lmp.py settings used:
	topofile = 'EXAMPLES/packmol_pdb_methods/atom_name_as_atom_type/water_packmol_atom_name_as_atom_type.pdb'
	nta_file = 'types_from_pdb.nta'
	frc_file = 'frc_files/pcff_iff_v1_5_CNT_poly_solv_Hpan_mod.frc'
	assumed = 'frc_files/general_assumed_equivs.coeffs'
	ff_class = 2
	parent_directory = 'topofile'
	
	where, nta_file = 'types_from_pdb.nta' will envoke all2lmp to using the atom name from the .pdb
	file to set the atom types to use for force field assignment.
	
This method is useful because it allows users to atom type a small number of atom's and then can use
packmol to increase the number of atoms, where packmol will maintain the atom name/atom type in the 
files.