This folder contains an example for .pbd/packmol integration into the LUNAR codes. File describitions:
  - water.inp, simple packmol script to read in water.pdb
  - water.pdb, generated using ChemDraw 3D and atom name(s) where switched to atomTypeIDs (1=o* and 2=hw)
  - water.nta, manually generated setting the atomTypeID map for the atomTypeIDs manually set in water.pdb (1=o* and 2=hw)
  - water_packmol_atom_name_as_atom_typeID.pdb, output file from packmol, where packmol maintains the atom types
  - water_packmol_atom_name_as_atom_typeID_IFF.data, output from all2lmp.py to generate a LAMMPS datafile using the
    atomtypeIDs that are set in the atom name column from the water_packmol_atom_name_as_atom_typeID.pdb in
	combination w/ the atomTypeIDs in water_packmol_atom_name_as_atom_typeID.pdb.
	
all2lmp.py settings used:
	topofile = 'EXAMPLES/packmol_pdb_methods/atom_name_as_atom_type_using_nta_file/water_packmol_atom_name_as_atom_typeID.pdb'
	nta_file = 'EXAMPLES/packmol_pdb_methods/atom_name_as_atom_type_using_nta_file/water.nta'
	frc_file = 'frc_files/pcff_iff_v1_5_CNT_poly_solv_Hpan_mod.frc'
	assumed = 'frc_files/general_assumed_equivs.coeffs'
	ff_class = 2
	parent_directory = 'topofile'
	
	where, in the nta_file the "style type" was used to map to the manually set atomTypeIDs in the
	water.pdb that packmol maintained in the water_packmol_atom_name_as_atom_typeID.pdb output.
	
This method is useful because it allows users to atom type a small number of atom's and then can use
packmol to increase the number of atoms, where packmol will maintain the atom name/atomTypeID in the 
files. This method also allows the user to add other optional commands to the nta file and also allows
the user to easily change atom types even after packmol has ran. Josh recommends this method as the 
primary method if users want to manually set atom types and wants to integrate w/ packmol with the
maximum possible control.