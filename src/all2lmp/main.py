# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.20
October 4th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

##############################
# Import Necessary Libraries #
##############################
import src.all2lmp.interatomic_force_fields as interatomic_force_fields
import src.all2lmp.bond_react_templates as bond_react_templates
import src.all2lmp.fill_in_parameters as fill_in_parameters
import src.all2lmp.assumed_auto_fill as assumed_auto_fill
import src.all2lmp.charge_neutralize as charge_neutralize
import src.all2lmp.remove_things as remove_things
import src.all2lmp.ff_functions as ff_functions
import src.all2lmp.check_inputs as check_inputs
import src.all2lmp.command_line as command_line
import src.all2lmp.write_lmp as write_lmp
import src.all2lmp.find_BADI as find_BADI
import src.all2lmp.read_frc as read_frc
import src.all2lmp.read_nta as read_nta
import src.mol2SYBYL2lmp as mol2SYBYL2lmp
import src.io_functions as io_functions
import src.read_lmp as read_lmp
import src.mol2lmp as mol2lmp
import src.pdb2lmp as pdb2lmp
import src.msi2lmp as msi2lmp
import time
import sys
import os





####################################################################
### Main function to perform all topology and FF assigning tasks ###
####################################################################
def main(topofile, nta_file, frc_file, assumed, parent_directory, newfile, atom_style, ff_class, use_auto_equivalence, use_assumed_auto_fill, reset_molids,
         reset_charges, write_txt_comments, write_bond_react, print_options, use_morse_bonds, include_type_labels, add2box, ignore_missing_parameters,
         shift, rotate, commandline_inputs, log=io_functions.LUNAR_logger()):
    start_time = time.time()
    
    # Configure log (default is level='production', switch to 'debug' if debuging)
    log.configure(level='production')
    #log.configure(level='debug', print2console=False)

    #########################
    # Command Line Override #
    #########################
    # if -opt or -man option is in commandline_inputs print options and stop code execution
    if '-opt' in commandline_inputs or  '-man' in commandline_inputs or print_options:
        # call man page and exit if '-opt' or '-man' is provided at the command line
        command_line.print_man_page(topofile, nta_file, frc_file, assumed, parent_directory, newfile, atom_style, ff_class, use_auto_equivalence, use_assumed_auto_fill, reset_molids,
                                    reset_charges, write_txt_comments, write_bond_react, print_options, use_morse_bonds, include_type_labels, add2box, ignore_missing_parameters,
                                    shift, rotate)
        sys.exit()
    
    ###################################################################################
    # if -opt option is NOT in commandline_inputs start assiging over ride parameters #
    ###################################################################################
    if '-opt' not in commandline_inputs and commandline_inputs:    
        # call inputs for commandline over rides
        over_rides = command_line.inputs(commandline_inputs, topofile, nta_file, frc_file, assumed, parent_directory, newfile, atom_style, ff_class, use_auto_equivalence,
                                         use_assumed_auto_fill, reset_molids, reset_charges, write_txt_comments, write_bond_react, print_options, use_morse_bonds,
                                         include_type_labels, add2box, ignore_missing_parameters, shift, rotate)
        
        # Set new inputs from over_rides class
        topofile = over_rides.topofile; nta_file = over_rides.nta_file; add2box = over_rides.add2box;
        frc_file = over_rides.frc_file; assumed = over_rides.assumed; ignore_missing_parameters = over_rides.ignore_missing_parameters
        parent_directory = over_rides.parent_dir;  newfile = over_rides.newfile; shift = over_rides.shift; rotate = over_rides.rotate
        atom_style = over_rides.atom_style; ff_class = over_rides.ff_class;
        use_auto_equivalence = over_rides.use_auto_equivalence; use_assumed_auto_fill = over_rides.use_assumed_auto_fill;
        reset_molids = over_rides.reset_molids; reset_charges = over_rides.reset_charges;
        write_txt_comments = over_rides.write_txt_comments; write_bond_react = over_rides.write_bond_react;
        print_options = over_rides.print_options; use_morse_bonds = over_rides.use_morse_bonds; include_type_labels = over_rides.include_type_labels;
    
    ######################################################
    # Check for valid inputs; if not valid exit the code #
    ######################################################
    if not check_inputs.safety(topofile, nta_file, frc_file, assumed, parent_directory, newfile, atom_style, ff_class, use_auto_equivalence, use_assumed_auto_fill, reset_molids,
                               reset_charges, write_txt_comments, write_bond_react, print_options, use_morse_bonds, include_type_labels, log):
        log.error('ERROR safety check failed')
    
    ###########################################
    # Initialize some preliminary information #
    ###########################################
    # set version and print starting information to screen
    version = 'v1.20 / 4 October 2024'
    log.out(f'\n\nRunning all2lmp {version}')
    log.out(f'Using Python version {sys.version}')
    log.out('Trying Atom Equivalences if needed')
    if use_auto_equivalence: log.out('Trying Atom Auto-Equivalences if needed')
    if use_assumed_auto_fill: log.out('Trying assumed auto fill if needed')
    if reset_charges: log.out('Trying to find atom charge based on bond-incs')
    if reset_molids: log.out('Resetting molecule-ids based on cluster analysis')
    log.out('\n')
    
    #####################################################################################
    # if nta_file is a force field name that atom_typing.py supports, use the nta_file  #
    # name to set the force field atom types. Currently supported ff_name's are:        #
    #    - PCFF-IFF                                                                     #
    #    - PCFF                                                                         #
    #    - compass                                                                      #
    #    - CVFF                                                                         #
    #    - CVFF-IFF                                                                     #
    #    - Clay-FF                                                                      #
    #    - DREIDING                                                                     #
    #    - OPLS-AA                                                                      #
    # if a '-q' is added to the name such as DREIDING-q, then use the Gasteiger charges #
    # from atom_typing.py and reset the reset_charge boolean of all2lmp.py to False.    #
    #####################################################################################
    # Setup atom_typing charging methods
    gasteiger_charging = False
    if '-q' in nta_file:
        check_ff = nta_file.replace('-q', '')
        gasteiger_charging = True
    else: check_ff = nta_file
    
    # see if check_ff is a supported type, if so find new topofile and nta_file with ff specifc atom types
    check_ff = os.path.basename(check_ff)
    if check_ff in ['PCFF-IFF', 'PCFF', 'compass', 'CVFF-IFF', 'CVFF', 'Clay-FF', 'DREIDING', 'OPLS-AA']:
        from src.atom_typing.main import main as atom_typing
        import atom_typing as at
        
        # Set atom_typing defaults from atom_typing.py file
        log.out(f'Calling atom_typing.py from all2lmp.py, because nta_file = {nta_file}.')
        log.out('   The following defaults are loaded from the variables in atom_typing.py file:')
        log.out(f'      delete_atoms = {at.delete_atoms}')
        key1 = list(at.mass_map.keys())[0]
        log.out(f'      mass_map = {key1}:{at.mass_map[key1]}, ...')
        key1 = list(at.bondorder.keys())[0]
        log.out(f'      bondorder = {key1}:{at.bondorder[key1]}, ...')
        key1 = list(at.maxbonded.keys())[0]
        log.out(f'      maxbonded = {key1}:{at.maxbonded[key1]}, ...')
        log.out(f'      boundary = {at.boundary}')
        log.out(f'      vdw_radius_scale = {at.vdw_radius_scale}')
        log.out(f'      bonds_via_distance_override = {at.bonds_via_distance_override}')
        log.out(f'      pdb_file = {at.pdb_file}')
        log.out(f'      chargefile = {at.chargefile}')
        log.out(f'      include_comments_nta = {at.include_comments_nta}')
        
        # Set all2lmp defaults for when calling atom_typing
        at_bondfile = 'n.u.'
        at_newfile = ':_typed_from_within_all2lmp'
        at_ff_name = check_ff
        log.out('   The following variables are being adjusted based on all2lmp.py:')
        log.out(f'      topofile = {topofile}')
        log.out(f'      bondfile = {at_bondfile}')
        log.out(f'      newfile = {at_newfile}')
        log.out(f'      ff_name = {at_ff_name}')
        log.out(f'      reset_charges = {gasteiger_charging} (for gasteiger charges from atom_typing.py)')
        log.out(f'      parent_directory = {parent_directory}')
        
        # Run atom_typing.py
        log.out('   Running atom_typing.py from all2lmp.py based on the above atom_typing.py inputs')
        atom_typing(topofile, at_bondfile, parent_directory, at_newfile, at_ff_name, at.delete_atoms, at.mass_map, at.bondorder, at.maxbonded, 
                    at.boundary, at.vdw_radius_scale, gasteiger_charging, False, [], at.bonds_via_distance_override, at.pdb_file, at.chargefile,
                    at.include_comments_nta, log=log)
        
        # Update the topofile, nta_file and reset_charges variables based on the atom_typing.py run
        if 'topofile' in parent_directory:
            dummy_parent_directory = ''
            if parent_directory == 'topofile':
                pwd = os.path.dirname(os.path.abspath(topofile))
            else:
                pwd = io_functions.get_dir_from_topofile(topofile, parent_directory)
                parent_directory = 'topofile'
        else: 
            pwd = os.getcwd()
            dummy_parent_directory = parent_directory
        basename_topofile = os.path.basename(topofile)
        basename_topofile = os.path.join(pwd, dummy_parent_directory, basename_topofile)
        basename_from_atom_typing = '{}{}'.format(basename_topofile[:basename_topofile.rfind('.')], at_newfile.replace(':', ''))
        topofile = '{}.data'.format(basename_from_atom_typing)
        nta_file = '{}.nta'.format(basename_from_atom_typing)
        log.out('   Updated the following all2lmp.py variables to use the files generated by atom_typing.py:')
        log.out(f'      topofile = {topofile}')
        log.out(f'      nta_file = {nta_file}')
        if gasteiger_charging:
            reset_charges = False
            log.out(f'      reset_charges = {reset_charges}')

        
        # Check for consistency between check_ff and ff_class
        ff2class = {'PCFF-IFF':2, 'PCFF':2, 'compass':2, 'CVFF-IFF':1,
                    'CVFF':1, 'Clay-FF':1, 'DREIDING':'d', 'OPLS-AA':0}
        if ff_class != ff2class[check_ff]:
            ff_class = ff2class[check_ff]
            log.out(f'      WARNING inconsistency between ff_name and ff_class. ff_name = {check_ff}.')
            log.out(f'      Interanlly updating ff_class to be consistent with ff_name. Updated ff_class = {ff_class}.')
        
        
    ######################################################
    # Read frc file and check frc file based on FF class #
    ######################################################
    # Find info of ff_class for class 1 or 2
    if ff_class in [0, 1, 2, 'd']:
        if os.path.isfile(frc_file):
            frc = read_frc.forcefield_file(frc_file, log)
            log.out(f'Read in {frc_file} forcefeild file')
        else: log.error(f'ERROR .frc file: {frc_file} does not exist')
        
    # If ff_class is 'i' or 'ilmp' for interatomic force fields like ReaxFF, SNAP, REBO, ... build nta dictionary from m class 
    if ff_class in ['i', 'ilmp']:
        if os.path.isfile(frc_file):
            frc = interatomic_force_fields.forcefield_file(frc_file)
            log.out(f'Read in {frc_file} forcefield file')
        else: log.error(f'ERROR .frc file: {frc_file} does not exist')
            
    # If ff_class is s1 or s2 build empty frc class
    if ff_class in ['s1', 's2']:
        log.out(f'Generating and empty frc class to pass as a variable to use the skeleton option: {ff_class}')
        class Empty_frc:  pass
        frc = Empty_frc()
        frc.type = 'skeleton' # add frc.type attribute as skeleton
    
    ########################################################################
    # Warn and exit if .frc file is inconsitent with data in the .frc file #
    ########################################################################
    # Check frc type for fix-bond FF's (read in with read_frc.forcefield_file)
    if frc.type == 'fix-bond':
        # Warn if FF class is 1 by looking if 12-6 is empty
        if ff_class == 0 and len(frc.torsion_1_opls) == 0:
            log.error('\n\nERROR Inconsistent data in force field file for Class0 FF: Zero 12-6 pair coeffs or Zero torsion_1 opls - Use a class0 .frc file\n')
    
        # Warn if FF class is 1 by looking if 12-6 is empty
        if ff_class == 1 and len(frc.pair_coeffs_12_6) == 0:
            log.error('\n\nERROR Inconsistent data in force field file for Class1 FF: Zero 12-6 pair coeffs - Use a class1 .frc file\n')
        elif ff_class == 1 and len(frc.torsion_1_opls) > 0:
            log.error('\n\nERROR Inconsistent data in force field file for Class1 FF: Non-Zero torsion_1 opls - Use a class1 .frc file\n')
        
        # Warn if FF class is 2 by looking if pair coeffs 9-6 is empty
        if ff_class == 2 and len(frc.pair_coeffs_9_6) == 0:
            log.error('\n\nERROR Inconsistent data in force field file for Class2 FF: Zero 9-6 pair coeffs - Use a class2 .frc file\n')
            
        # Warn if FF class is 1 by looking if 12-6 and out_of_plane_DREIDING is empty
        if ff_class == 'd' and len(frc.pair_coeffs_12_6) == 0 and len(frc.out_of_plane_DREIDING) == 0:
            log.error('\n\nERROR Inconsistent data in force field file for Classd FF: Zero 12-6 pair coeffs or Zero DREIDING OOP parms - Use a classd .frc file\n')
            
        # Warn if reaxff file was read in as a fix bond file
        if len(frc.equivalences) == 0 and len(frc.auto_equivalences) == 0:
            log.error('\n\nERROR Inconsistent data in force field file: Zero equivalence or Zero auto equivalences - Use a class 0, 1, or 2 .frc file, likely read in the reaxFF file.\n')
            
    # Check frc type for interatomic (read in with interatomic_force_fields.forcefield_file - Currently no error logging or exiting conditions known, thus pass)
    elif frc.type == 'interatomic': pass
    
    # Check if frc type is skeleton (Currently no error logging or exiting conditions known, thus pass)
    elif frc.type == 'skeleton': pass
    
    # Check if frc.type was flagged as 'Inconsistent force field file' (happens when trying to read a fix-bond frc file with ff_class = 'r')
    # (read in with reaxff.forcefield_file, but the file provided was actually a fix bond FF file)
    elif frc.type == 'Inconsistent force field file':
        log.error('\n\nERROR Inconsistent data in force field file for reaxFF: Use all2lmp reaxff specific .frc file\n')
    
    # else exit with general message of can not resole
    else: log.error('\n\nERROR Inconsistent force field file - can not resolve inconsistency\n')
        
    #########################################################################
    # Read .data file or .mol or .mol2 or .car/.mdf file based on extension #
    #########################################################################
    # Read lammps data file
    if topofile.endswith('data') or topofile.endswith('dat') or topofile.endswith('data.gz') or topofile.endswith('dat.gz'):
        if os.path.isfile(topofile):
            m = read_lmp.Molecule_File(topofile, method='forward', sections = ['Atoms', 'Bonds', 'Velocities'])
            log.out(f'Read in {m.filename} LAMMPS datafile')
        else: log.error(f'ERROR lammps datafile: {topofile} does not exist'); sys.exit();
        
    # Read .mol file
    elif topofile.endswith('mol') or topofile.endswith('sdf'):
        if os.path.isfile(topofile):
            m = mol2lmp.Molecule_File(topofile)
            log.out(f'Read in {m.filename} chemdraw .mol or .sdf file')
        else: log.error(f'ERROR .mol or .sdf file: {topofile} does not exist')
        
    # Read .mol2 file (VMD MDL MOL2 file)
    elif topofile.endswith('mol2'):
        if os.path.isfile(topofile):
            m = mol2SYBYL2lmp.Molecule_File(topofile)
            log.out(f'Read in {m.filename} SYBYL MOL2 file')
        else: log.error(f'ERROR .mol2 file: {topofile} does not exist')
        
    # Read .mol2 file (pdb file)
    elif topofile.endswith('pdb'):
        if os.path.isfile(topofile):
            m = pdb2lmp.Molecule_File(topofile)
            log.out(f'Read in {m.filename} pdb file')
        else: log.error(f'ERROR .pdb file: {topofile} does not exist')
            
    # Read .car and .mdf files. topofile=.mdf and nta_file=.car
    elif topofile.endswith('mdf') and nta_file.endswith('car') or nta_file == 'topofile':
        # update nta_file name if nta_file == 'topofile'
        if nta_file == 'topofile':
            base = topofile[:topofile.rfind('.')]; ext = 'car'
            nta_file = '{}.{}'.format(base, ext)
            log.out(f'Using path and topfile name from topofile to set nta_file = {nta_file}')
    
        # Read msi files
        if os.path.isfile(topofile) and os.path.isfile(nta_file):
            # wrap_atoms set to True will tell msi2lmp to wrap atoms and False will not and let LAMMPS wrap the atoms (support for orthogonal or triclinc
            #  False). no_center set to True will tell msi2lmp to not center atoms; False will keep default of centering atoms except for triclinic cells
            m = msi2lmp.Molecule_File(nta_file, topofile, wrap_atoms=True, no_center=False)
            log.out(f'Read in {m.filename} material studio file')
        else: log.error(f'ERROR .mdf file or .car file: {topofile} or {nta_file} does not exist')
    
    # else file is not supported
    else: log.error(f'ERROR unsupported topofile extension {topofile}')
    
    # Read nta file if topofile had .data or .data or .mol or .mol2 or .pdb ending (nta = {atomid:atomtype}, edge={atomid:[lst of edge extend atomtypes]}) and class1 or class2 FF's
    if ff_class in [0, 1, 2, 'd', 's1', 's2']:        
        # Read .nta file for .data, .dat, .mol, .sdf or .mol2  or .pdb file formats. Pass in m class for checking and to support style id or style type .nta files
        if topofile.endswith('data') or topofile.endswith('dat') or topofile.endswith('mol') or topofile.endswith('sdf') or topofile.endswith('mol2') or topofile.endswith('pdb') or topofile.endswith('data.gz') or topofile.endswith('dat.gz'):
            # update nta_file name if nta_file == 'topofile'
            if nta_file == 'topofile':
                base = topofile[:topofile.rfind('.')]; ext = 'nta'
                nta_file = '{}.{}'.format(base, ext)
                log.out(f'Using path and topofile name from topofile to set nta_file = {nta_file}')
                
            # get nta dict from atom names in .pdb if topofile is .pdb and nta_file string is 'types_from_pdb.nta'
            if topofile.endswith('pdb') and nta_file == 'types_from_pdb.nta':
                log.out('Using atom name (columns 13-16 in .pdb for assigning atom types)')
                nta = {}; name = {}; edge = {}; charges = {}; neutralize = {'all': False, 'bond-inc': False, 'user-defined': False, 'zero': False}
                remove = {'angle-nta':[], 'dihedral-nta':[], 'improper-nta':[], 'angle-ID':[], 'dihedral-ID':[], 'improper-ID':[], 'zero':{'angle':False, 'dihedral':False, 'improper':False}}
                for i in m.atoms:
                    atom = m.atoms[i]
                    nta[i] = atom.atom_name
                    name[i] = atom.atom_name
                    
            # Read nta_file
            else:
                if os.path.isfile(nta_file):
                    nta, name, edge, charges, neutralize, remove = read_nta.atoms(nta_file, m, log)
                    log.out(f'Read in {nta_file} new-type-assignment file')
                else: log.error(f'ERROR .nta file: {nta_file} does not exist')
                
        # elif set nta diectionary from .car and .mdf files and set edge as an empty dict
        elif topofile.endswith('mdf') and nta_file.endswith('car'):
            # Get nta dict from m and set edge dict and charges as empty since and set all neutralization methods as False along with
            # a null remove dict. Since these option are not available for msi files, since the .nta info comes from the .car file.
            nta = m.nta; name = m.name; edge = {}; charges = {}; neutralize = {'all': False, 'bond-inc': False, 'user-defined': False, 'zero': False}
            remove = {'angle-nta':[], 'dihedral-nta':[], 'improper-nta':[], 'angle-ID':[], 'dihedral-ID':[], 'improper-ID':[], 'zero':{'angle':False, 'dihedral':False, 'improper':False}}
            log.out(f'nta dictionary came from msi {nta_file} file')
            
    # If ff_class is 'i' or 'ilmp' for interatomic force fields like reaxFF, REBO, AIREBO, SNAP,... build nta dictionary from m class 
    if ff_class in ['i', 'ilmp']:
        nta, name, edge, charges, neutralize, remove = interatomic_force_fields.generate_nta(m, frc, reset_charges, log)
    
    # Read asssumed auto fill coeffs (aafc) file
    if use_assumed_auto_fill:
        if os.path.isfile(assumed):
            aafc = assumed_auto_fill.read_assumed_coeffs(assumed)
            log.out(f'Read in {assumed} assumed-auto-fill file')
        else: log.error(f'ERROR assumed file: {assumed} does not exist')
    else: aafc = {}
        
    #########################################################
    # Check that assumed coeffs exist in the read .frc file #
    #########################################################
    if use_assumed_auto_fill:
        aafc = assumed_auto_fill.check_assumed_auto_fill_types(frc, aafc, frc_file, log)
    
    ####################
    # Finding topology #
    ####################
    # Find bonds, angles, dihedrals, impropers, atom types, bond types, angle types, dihedral types, and improper types
    BADI = find_BADI.BADI(m, nta, name, ff_class, log)
    
    # Fill in parameters for all coeff types in read in system and assign charge to each atom via summing of bond increments of bonded atoms
    parameters = fill_in_parameters.get(nta, name, edge, frc, BADI, m, use_auto_equivalence, use_assumed_auto_fill, aafc, reset_charges, reset_molids,
                                        ff_class, use_morse_bonds, add2box, ignore_missing_parameters, shift, rotate, log)
    
    # If anything from remove dict use remove_things.post_processor to remove certain user-defined topology/parameters from the parameters class
    remove_booleans = [True for i in remove if remove[i] and i in ['angle-nta', 'angle-ID', 'dihedral-nta', 'dihedral-ID', 'improper-nta', 'improper-ID']]
    remove_booleans.extend([remove['zero'][i] for i in remove['zero']])
    if any(remove_booleans) and ff_class in [0, 1, 2, 'd']:
        parameters = remove_things.post_processor(parameters, remove, ff_class, BADI, log)
    
    ##############################################
    # Update charges or neutralize system charge #
    ##############################################
    if neutralize['all'] or neutralize['bond-inc'] or neutralize['user-defined'] or neutralize['zero'] or charges:
        parameters = charge_neutralize.update(neutralize, parameters, nta_file, charges, log)
    
    ####################################################################################################
    # Compute density and volume along system quantities. Generate percetage table and use ignore flag #
    ####################################################################################################
    ff_functions.compute_mass_volume_density(parameters, BADI, ff_class, remove_booleans, reset_charges, ignore_missing_parameters, frc_file, log)
    
    ########################################################################
    # Setting up directories and where to write final files and results to #
    ########################################################################
    # Find present working directory
    pwd = os.getcwd()
    
    # Find/create paths to store code results
    path = os.path.join(pwd, parent_directory)
    
    # If parent_directory == 'topofile' use topofile path as path
    if 'topofile' in parent_directory:
        log.out('\n\nUsing path from topofile to set parent_directory ...')
        path = io_functions.get_dir_from_topofile(topofile, parent_directory)
        
    # Check if path exits. IF not create
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
            
    # Change the current working directory to path so all files get written to that directory
    os.chdir(path)

    ##############################
    # Write new lammps data file #
    ##############################
    # Set new file name and write new datafile and version number
    basename = io_functions.get_basename(topofile, newfile=newfile, character=':', pflag=True)
       
    # write lammps datafile
    log.out('\n\nWriting LAMMPS datafile')
    if newfile != ':':
        write_lmp.datafile(basename, atom_style, parameters, ff_class, version, include_type_labels, log)
    
    # write comments file if user desires this file
    if write_txt_comments and newfile != ':':
        log.out('Writing verbose comment file')
        write_lmp.comments(basename, atom_style, parameters, ff_class, version, log)
        
    # write fix bond/react template files if user desires these files
    if write_bond_react:
        bond_react_templates.write_moltemp(basename, parameters, ff_class, include_type_labels, version)
        log.out('\n\n ********************************************************* ')
        log.out(' *** The molecule template file successfully created!  *** ')
        
        bond_react_templates.write_ecoeffs(basename, parameters, ff_class, version, include_type_labels)
        log.out(' *** The extra coefficients file successfully created! *** ')
        log.out(' *** Note - Addition of new atom types (post reaction) *** ')
        log.out(' *** need to be added manually in the Massess and Pair *** ')
        log.out(' *** Coeffs sections of the *.ecoeffs file and in the  *** ')
        log.out(' *** Types section of the .moltemp file.               *** ')
        log.out(' ********************************************************* ')
    
    # Print file locations
    log.out(f'\n\nAll outputs can be found in {path} directory')
    
    # Print style hints info
    log.out('\n\nEach energy coeff Header in the written datafile contains style hint flags for')
    log.out('applying proper settings in LAMMPS to be consistent with the energy definitions')
    log.out('and ordering of the energy parameters.')
    
    # Print completion of code
    log.out('\n\nNormal program termination\n\n')
    
    # Script run time
    execution_time = (time.time() - start_time)
    log.out('Execution time in seconds: ' + str(execution_time))
    
    # Show number of warnings and errors
    log.out_warnings_and_errors()
    
    # write log
    if newfile != ':':
        log.write_logged(basename+'.log.lunar')
            
    # Change back to the intial directory to keep directory free for deletion
    os.chdir(pwd)
    return parameters