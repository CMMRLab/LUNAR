# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.8
October 4th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

fill_in_parameters is long script with 1 class called get which will find all info that is
needed to write a class2 FF data file. The get class calls many functions to find each set
of info needed to be found. Underneath the get class there are seperate functions to find 
the following data for directly writing to the new lammps data file:
    
    - Structure/Topology:
        - Box dimensions (from data file, not currenlt searching atom positions to build box)
        - natoms, nbonds, nangles, ndihedrals, nimpropers
        - natom types, nbond types, nangle types, ndihedral types, nimproper types
        - atoms class with type, molid, charge, x pos, y pos, z pos, image flags: x, y, z (from data file, not currently finding new ones)
        - bonds class with type and bond atom ids
        - angles class with type and angle atom ids
        - dihedrals class with type and dihedral atom ids
        - impropers class with type and improper atom ids
        
    - Potential Parameters/Coeffs
        - masses class with type, mass, and equivalences used
        - pair coeffs class with type, pair coeffs, and equivalences used
        - bond coeffs class with type, bond coeffs, equivalence used, comments
        - angle coeffs class with type, angle coeffs, equivalence used, comments
        - dihedral coeffs class with type, dihedral coeffs, equivalence used, comments
        - improper coeffs class with type, improper coeffs, equivalence used, comments
        - bondbond coeffs class with type, bondbond coeffs, equivalence used, comments
        - bondangle coeffs class with type, bondangle coeffs, equivalence used, comments   
        - angleangle coeffs class with type, angleangle coeffs, equivalence used, comments   
        - angleangletorsion coeffs class with type, angleangletorsion coeffs, equivalence used, comments   
        - endbondtorsion coeffs class with type, endbondtorsion coeffs, equivalence used, comments  
        - middlebondtorsion coeffs class with type, middlebondtorsion coeffs, equivalence used, comments  
        - bondbond13 coeffs class with type, bonddbond13 coeffs, equivalence used, comments  
        - angletorsion coeffs class with type, angletorsion coeffs, equivalence used, comments  
        
This script is long and complex, but hopefully written well enough to take apart and reuse/debug
if necessary. Each function only serves one purpose and that is to find each set of info that it
tries to find. This should for allow quick browsing and understanding of the structure of this script.
This script is perhaps the most complex script in the code, but also has to prefrom some of the most
complex tasks of finding all info needed to write the final data file. Some parts of this script uses
read in values from the inputed lammps data file such as (box dimensions and image flags), to build this
set of scripts to work on other topology files, you will need to examine read_PCFF and see how it reads,
organizes,and stores data in a class and then build similar scripts to read files such as .mol, .pdb, ...
if you desire to use this set of scripts as a starting point to convert other files to lammps data files.
mol2lmp as serves as a good example of what could be used to read in other molecule file formats.
"""


################################
# Import Necessary Libraries   #
################################
import src.all2lmp.ff_functions as ff_functions
import sys

class Box: pass # .xlo .xhi .ylo .yhi .zlo .zhi .lx .ly .lz
class Atom: pass # .type .molid .charge .x .y .z .iflag_x .iflag_y .iflag_z .symbol
class Bond: pass  # .type .atomids = [atom1id, atom2id] .symbol
class Angle: pass  # .type .atomids = [atom1id, atom2id, atom3id] .symbol
class Dihedral: pass  # .type .atomids = [atom1id, atom2id, atom3id, atom4id] .symbol
class Improper: pass  # .type .atomids = [atom1,atom2,atom3,atom4] .symbol
class Type: pass # .coeffs .type .equivalent .comment
class get:
    def __init__(self, nta, name, edge, frc, BADI, m, use_auto_equivalence, use_assumed_auto_fill, aafc, reset_charges, reset_molids, ff_class,
                 use_morse_bonds, add2box, ignore_missing_parameters, shift, rotate, log):
        
        # Perform an atom types check to warn user if atom-type does not exists if .frc file
        log.out('\n\n\nGetting force field parameters for the read in system')
        if ff_class in [0, 1, 2, 'r', 'd']: atom_types_check(nta, frc, ignore_missing_parameters, log)
        
        # Option to skip print outs. For development only!!! (True to skip False to print out)
        skip_printouts = False 
        
        # If ff_class == 'd' reset some options thad dont exists in a DREIND force feild
        if ff_class == 'd':
            use_auto_equivalence = False; reset_charges = False; # No bond-incs or not auto-equivalences
            log.out('\nUsing class d for DREIDING force field, resetting use_auto_equivalence and reset_charges to False')
            log.out('since DREIDING does not currently offer this level of support. Also the all2lmp_dreiding.frc almost fully')
            log.out('supports all attributes of the DREIDING force field. However the H-bond parameters are not in the file and')
            log.out('consequently any H-bond LJ parameters will have to be dealt with manually if your systems has H-bonds.')
        
        # Access info through init method        
        self.header = m.header
        
        # shift and rotate atoms (if user desires)
        phi = rotate['x']; theta = rotate['y']; psi = rotate['z'];
        if shift['x'] != 0 or shift['y'] != 0 or shift['z'] != 0:
            m = ff_functions.shift_system(m, shift, log)
            log.out(f'Shifted system by {shift["x"]} in x-axis, by {shift["y"]} in y-axis, and by {shift["z"]} in z-axis.')
        if phi != 0 or theta != 0 or psi != 0:
            sum_iflags = 0; x = []; y = []; z = [];
            for i in m.atoms:
                atom = m.atoms[i]
                sum_iflags += abs(atom.ix)
                sum_iflags += abs(atom.iy)
                sum_iflags += abs(atom.iz)
                x.append(atom.x)
                y.append(atom.y)
                z.append(atom.z)
            if sum_iflags == 0:
                xc = sum(x)/len(x); yc = sum(y)/len(y); zc = sum(z)/len(z)
                xs = (0 - xc); ys = (0 - yc); zs = (0 - zc);
                m = ff_functions.shift_system(m, {'x':xs, 'y':ys, 'z':zs}, log)
                m = ff_functions.rotate_system(m, phi, theta, psi)
                m = ff_functions.shift_system(m, {'x':xc, 'y':yc, 'z':zc}, log)
                log.out(f'Rotated system about x-axis by {phi}, y-axis by {theta}, and z-axis by {psi} degrees.')
            else: log.warn('WARNING could not rotate system as system is periodically bonded (non-zero image flags).')

            
        
        
        # Find box dimensions which is done from the read in box dimensions of the file or reset_box_dims=True will search all atom coordinates and find min/max values to set box dimensions
        self.box = get_box(m, add2box, log)
        
        # If ff_class is 'i' or 'ilmp' for interatomic force fields like reaxFF, REBO, AIREBO, SNAP,... find atom types only
        if ff_class in ['i', 'ilmp']:
            # Set total number of atoms/bonds/angles/dihedrals/impropers
            self.natoms = len(m.atoms); self.nbonds = 0;
            self.nangles = 0; self.ndihedrals = 0; self.nimpropers = 0;
            
            # Set total number of atoms/bonds/angles/dihedrals/impropers types
            self.natomtypes = len(BADI.atom_types_lst); self.nbondtypes = 0;
            self.nangletypes = 0; self.ndihedraltypes = 0; self.nimpropertypes = 0;
            
            # Find atoms info such as setting up atoms charge, molid, x-pos, y-pos, z-pos, and image flags along with the masses and pair_coeffs
            self.atoms, self.velocities = find_atoms(nta, name, edge, frc, BADI, m, reset_charges, reset_molids, use_auto_equivalence, ff_class, skip_printouts, log)    
            self.masses, self.mass_comment = find_interatomic_atom_parameters(frc, BADI, m, ff_class, log)
            
        
        # Find info of ff_class for class 0 or 1 or 2 or 'd' or 's1' or 's2'
        if ff_class in [0, 1, 2, 'd', 's1', 's2']:
                
            # Set total number of atoms/bonds/angles/dihedrals/impropers
            self.natoms = len(m.atoms)
            self.nbonds = len(BADI.bonds)
            self.nangles = len(BADI.angles)
            self.ndihedrals = len(BADI.dihedrals)
            self.nimpropers = len(BADI.impropers)
            
            # Set total number of atoms/bonds/angles/dihedrals/impropers types
            self.natomtypes = len(BADI.atom_types_lst)
            self.nbondtypes = len(BADI.bond_types_lst)
            self.nangletypes = len(BADI.angle_types_lst)
            self.ndihedraltypes = len(BADI.dihedral_types_lst)
            self.nimpropertypes = len(BADI.improper_types_lst)
            
            # Option to reset and resort atomTypes in bond/angle/dihedral/improper and crossterm coeffs based on FF ordering and then reset bonds/angles/dihedrals/impropers atomID list orderings
            # to make the lists of atomIDs consistent with the atomTypes ordering of the coeff. THIS SHOULD ALWAYS BE TRUE for dihedrals, impropers and crossterms, but may be adjusted for all other
            # coeff and atomID lists, but might as well sort them as well.
            sort_remap_atomids = True
            
            # Find atoms info such as setting up atoms charge, molid, x-pos, y-pos, z-pos, and image flags along with the masses and pair_coeffs
            self.atoms, self.velocities = find_atoms(nta, name, edge, frc, BADI, m, reset_charges, reset_molids, use_auto_equivalence, ff_class, skip_printouts, log)        
            self.masses, self.pair_coeffs, self.mass_comment, self.pair_comment = find_atom_parameters(frc, BADI, use_auto_equivalence, ff_class, skip_printouts, log)
            
            # Find bond_coeffs and then find bonds and reorder atomids to match bond coeff using sort_remap_atomids (Reordering is Not Necessary - might as well though for cleanliness)
            self.bond_coeffs, self.bond_map, self.bond_comment = find_bond_parameters(frc, BADI, use_auto_equivalence, use_assumed_auto_fill, aafc, ff_class, use_morse_bonds, skip_printouts, log)
            self.bonds = find_bonds(nta, BADI, self.bond_map, sort_remap_atomids, log)
            
            # Find angle_coeffs and then find angles and reorder atomids to match angle coeff using sort_remap_atomids (Reordering is Not Necessary - might as well though for cleanliness)
            self.angle_coeffs, self.angle_map, self.angle_comment = find_angle_parameters(frc, BADI, use_auto_equivalence, use_assumed_auto_fill, aafc, ff_class, skip_printouts, log)
            self.angles = find_angles(nta, BADI, self.angle_map, sort_remap_atomids, ff_class, log)
            
            # Find dihedral_coeffs and then find dihedrals and reorder atomids to match dihedral coeff using sort_remap_atomids (Reordering is Not Necessary - might as well though for cleanliness)   
            self.dihedral_coeffs, self.dihedral_map, self.dihedral_comment = find_dihedral_parameters(frc, BADI, use_auto_equivalence, use_assumed_auto_fill, aafc, ff_class, skip_printouts, log)
            self.dihedrals = find_dihedrals(nta, BADI, self.dihedral_map, sort_remap_atomids, ff_class, log)
    
            # Find improper_coeffs and then find impropers and reorder atomids to match improper coeff using sort_remap_atomids (Reordering MAYBE Necessary - class2 has symmetry, but class1 does not)
            self.improper_coeffs, self.improper_map, self.improper_comment = find_improper_parameters(frc, BADI, use_auto_equivalence, use_assumed_auto_fill, aafc, ff_class, skip_printouts, log)
            self.impropers = find_impropers(nta, BADI, self.improper_map, m, sort_remap_atomids, ff_class, log)
            
            # Find cross terms based on ordering from *_map dictionaries if sort_remap_atomids: else the ordering will be from default ordering of intial pass through to find coeff types
            if ff_class in [2, 's2']:
                self.bondbond_coeffs = find_bondbond_parameters(frc, BADI, self.angle_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log)
                self.bondangle_coeffs = find_bondangle_parameters(frc, BADI, self.angle_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log)
                self.angleangletorsion_coeffs = find_angleangletorsion_parameters(frc, BADI, self.dihedral_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log)
                self.endbondtorsion_coeffs = find_endbondtorsion_parameters(frc, BADI, self.dihedral_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log)
                self.middlebondtorsion_coeffs = find_middlebondtorsion_parameters(frc, BADI, self.dihedral_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log)
                self.bondbond13_coeffs = find_bondbond13_parameters(frc, BADI, self.dihedral_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log)
                self.angletorsion_coeffs = find_angletorsion_parameters(frc, BADI, self.dihedral_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log)
                self.angleangle_coeffs = find_angleangle_parameters(frc, BADI, self.improper_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log)
        

################################################
# Function to perform atom-types check to warn #
# user if atom-type exists in .frc file or not #
################################################
def atom_types_check(nta, frc, ignore_missing_parameters, log):
    log.out('\nChecking for atom-type assignment consistency with atom-types in .frc file')
    unique_types = sorted({nta[i] for i in nta}); failed = [];
    for Type in unique_types:
        if Type not in frc.atom_types:
            failed.append(Type)
            if ignore_missing_parameters: log.warn(f'WARNING atom-type: {Type} not found in {frc.filename} file.')
            else: log.out(f'ERROR atom-type: {Type} not found in {frc.filename} file.')
    if failed:
        log.out(f'This indicates incorrect atom-types for usage with the {frc.filename} file.')
        if not ignore_missing_parameters: sys.exit()
    return


################################
# Function for retrieving box  #
# dimensions from orginal file #
################################    
def get_box(m, add2box, log):
    xline = m.xbox_line.split()
    yline = m.ybox_line.split()
    zline = m.zbox_line.split()
    xlo = float(xline[0]); xhi = float(xline[1]);
    ylo = float(yline[0]); yhi = float(yline[1]);
    zlo = float(zline[0]); zhi = float(zline[1]);
    lx = xhi - xlo; ly = yhi - ylo; lz = zhi - zlo;
    
    # add to box
    if add2box != 0:
        sum_iflags = 0;
        for i in m.atoms:
            atom = m.atoms[i]
            sum_iflags += abs(atom.ix)
            sum_iflags += abs(atom.iy)
            sum_iflags += abs(atom.iz)
        if sum_iflags == 0:
            log.out(f'Adjusting box on all 6-faces by {add2box}')
            xlo -= add2box
            xhi += add2box
            ylo -= add2box
            yhi += add2box
            zlo -= add2box
            zhi += add2box
        else: log.warn('WARNING add2box is none zero and read-in topofile has none-zero image flags. Box dimensions left alone.')
    
    # Set new box dimensions
    b = Box()
    b.xlo = xlo; b.xhi = xhi;
    b.ylo = ylo; b.yhi = yhi;
    b.zlo = zlo; b.zhi = zhi;
    b.lx = lx; b.ly = ly; b.lz = lz;
    b.xy = m.xy; b.xz = m.xz; b.yz = m.yz;
    box = b 
    return box
        

##############################
# Function for finding atoms #
##############################
def find_atoms(nta, name, edge, frc, BADI, m, reset_charges, reset_molids, use_auto_equivalence, ff_class, skip_printouts, log):
    
    # Dictionaries to add info to
    atoms = {}   # {atom number : atom object}
    velocities = {}  # {atom number : tuple of velocities}

    # If user wants molids compute them and reset them
    if reset_molids:
        log.out('\n\nFinding molecules to reset molids....')
        molecules, molids = ff_functions.clusters(m, log, pflag=True)
        log.out('\n')
    
    # Find new atoms info
    system_charge_sum = 0; molecules = set([]); system_periodicity = [];
    for id1 in m.atoms:
        
        # Find atom info and intialize other info
        atom = m.atoms[id1]  
        new_atom_type = nta[id1]
        additional_comments = ''
        try: velocity = m.velocities[id1]
        except: velocity = (0, 0, 0)
        velocities[id1] = velocity
        
        # Check if atom type exists in .frc file
        if ff_class in [0, 1, 2, 'd', 'r']:
            if new_atom_type in frc.atom_types:
                try: # Try finding connections and printing to user (not all .frc files have connections)
                    atom_connections = frc.atom_types[new_atom_type].connection
                
                    # Warn if inconsistent # of connects if .frc file has connections column (n.u. stands for not used - from ReaxFF nomenclature)
                    if atom_connections != 'n.u.' and len(BADI.graph[id1]) != atom_connections and new_atom_type != 'cg1': # reduce errors for cg1 with 2IFF if cge is used
                        if not skip_printouts and ff_class in [0, 1, 2, 'r', 'd']: log.warn('WARNING inconsistent # of connects on atom-id {} type {}'.format(id1, new_atom_type))
                except: pass
            
            # Else warn that atom type does not exists in force field file
            else:
                if not skip_printouts and ff_class in [0, 1, 2, 'r', 'd']:
                    log.warn('WARNING atom-id {} type {} does not exists in the force field file'.format(id1, new_atom_type))
        
        # Reset molids is user wants, else try getting from file, except set as 1
        if reset_molids in [True, False]:
            molid = molids[id1]
        else:
            try: molid = atom.molid
            except: molid = 1
            try: molid = int(reset_molids)
            except: pass
            
        # try to get image flags, if not set as zero (future development could include re-imaging the atoms)
        try: ix = atom.ix; iy = atom.iy; iz = atom.iz;
        except: ix = 0; iy = 0; iz = 0;
            
        # try getting element; except set element as atom type
        try: new_atom_element = frc.atom_types[new_atom_type].element
        except: 
            if ff_class in [0, 1, 2, 'r', 'd']: log.warn(f'WARNING element type for atomID {id1} type {new_atom_type} could not be found. This will create an error in bond_react_merge.py for map file generation')
            new_atom_element = new_atom_type # set as new_atom_type just in case
    
        # Reset charge via bond-increments if user wants, else use charge from file
        if reset_charges and ff_class in [0, 1, 2]:
            charge, bond_incs = ff_functions.charge_from_bond_increments(id1, BADI.graph, nta, frc, use_auto_equivalence, skip_printouts, log)
            
            # If atomid in edge dict loop through user defined edge types to sum charges, add extended_charge to charge and print usage status
            if id1 in edge:
                edges = edge[id1]; string_edges = '  '.join(edges);
                extended_charge = ff_functions.nta_edge_id_bond_incs_extend(id1, new_atom_type, edges, frc, use_auto_equivalence, log)
                charge = charge + extended_charge; additional_comments = '  [ edge id used: {} ]'.format(string_edges)
                if not skip_printouts: log.out(f'nta file had edge id section for atomid: {id1}. Extended edges for charge from bond-incs: {string_edges}')
        else: charge = atom.charge; bond_incs = {}; # set empty bond_incs dictionary
                
        # Sum charges and append molids to track system charge and periodicity
        system_charge_sum += charge;  molecules.add(molid); system_periodicity.append([ix, iy, iz])        

        # Save atom info into class
        a = Atom()
        a.type = BADI.atom_types_dict[name[id1]]
        a.symbol =  new_atom_type #name[id1] #new_atom_type
        a.name = name[id1]
        a.comments = additional_comments
        a.element = new_atom_element
        a.molid = molid
        a.charge = charge
        a.x = atom.x
        a.y = atom.y
        a.z = atom.z
        a.ix = ix
        a.iy = iy
        a.iz = iz
        a.bond_incs = bond_incs
        atoms[id1] = a
        
        
    # Function to check is there are none-zero image flags to tell user if system is periodic or not
    def check_images(system_periodicity):
        periodicity = 'has all image flags set to zero'
        for iflags in system_periodicity:
            if any(iflags) != 0: periodicity = 'has some non-zero image flags'; break
        return periodicity    
    
    # Print atoms info
    log.out('\n\nAtoms information has been found')
    log.out('  {} {}'.format(m.filename, check_images(system_periodicity)) )
    log.out('  There are {} atoms in {} molecules in this file'.format(len(atoms), len(molecules)) )
    if abs(system_charge_sum) == 0.0: system_charge_sum = 0.0 # python tallying sometimes produces -0.0000 so reset to zero in this case
    log.out('  The total charge in the system is {:.4f}'.format(system_charge_sum) )
    return atoms, velocities

##########################################
# Function for finding reaxff atom types #
##########################################
def find_interatomic_atom_parameters(frc, BADI, m, ff_class, log):
    """
    Function to fill in interatomic force fields like reaxFF, REBO, AIREBO, SNAP,... specific masses
    """
    # Add printing buffer
    log.out('')

    # Dictionaries to add info to    
    masses = {} # { atom type : mass }
    
    # set comments
    mass_comment = 'line1: pair_style reaxff NULL    line2: pair_coeff * * ffield.reax {}'.format(' '.join(BADI.atom_types_lst))
    
    # Find masses
    for number, i in enumerate(BADI.atom_types_lst, 1):
        
        # atom number
        #number = BADI.atom_types_dict[i]
        mass_equiv = ''; mass_coeff_comment = '';
        
        # strip any numbers from i 
        i = ''.join([j for j in i if j.isalpha()])

        # Find mass, try type 1st, then equiv, else set as zeros
        if i in frc.atom_types:
            mass_equiv = i
            mass = frc.atom_types[i].mass
            mass_coeff_comment = 'standard type used'    
        else:
            # Try getting mass from currently defined masses
            try: 
                mass = m.masses[number].coeffs[0]
                mass_equiv = 'N/A'
                mass_coeff_comment = 'Using mass from previously defined mass'
                log.warn('WARNING unable to find atom info for Masses {} {}, but using previously define mass for atomTypeID'.format(number, i))
            except:
                mass = 0.0; mass_equiv = 'N/A';
                mass_coeff_comment = 'UNABLE to find coeff parameters'
                log.warn('WARNING unable to find atom info for Masses {} {}'.format(number, i))
            
        # Save mass info into class
        t1 = Type()
        t1.type = i
        t1.equivalent = mass_equiv
        t1.coeffs = mass
        t1.comments = mass_coeff_comment
        masses[number] = t1
    return masses, mass_comment


###################################
# Function for finding atom types #
###################################
def find_atom_parameters(frc, BADI, use_auto_equivalence, ff_class, skip_printouts, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html
    https://docs.lammps.org/pair_class2.html
    https://docs.lammps.org/pair_lj_cut_coul.html
    
    Nonbond Lennard-Jones
    The style of nonbond potential is specified in the input command file.
    
    (1) lj/cutoff
    
      E = 4 epsilon [ (sigma/r)^12 - (sigma/r)^6 ]
    
      standard Lennard Jones potential
    
      r = distance (computed by LAMMPS)
    
      coeff1 = epsilon (energy)
      coeff2 = sigma (distance)
    
      2 coeffs are listed in data file or set in input script
      1 cutoff is set in input script
        
    (5) class2/cutoff
    
      E = epsilon [ 2 (sigma/r)^9 - 3 (sigma/r)^6 ]
    
      used with class2 bonded force field
    
      r = distance (computed by LAMMPS)
    
      coeff1 = epsilon (energy)
      coeff2 = sigma (distance)
    
      2 coeffs are listed in data file or set in input script
      1 cutoff is set in input script
    """
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # Dictionaries to add info to    
    masses = {} # { atom type : mass }
    pair_coeffs = {}  # {atom type : list of coeffs }

              
    # Set comments, equivs/auto coeff dicts and equivs/auto mapping dicts
    if ff_class == 0:
        mass_comment = 'class1' # class1 Header comment
        pair_comment = 'lj/cut/coul/long' # class1 Header comment
        equivalents = frc.equivalences # class1 equivalents dict
        auto_equivs = frc.auto_equivalences # class1 auto-equivalents dict
        pair_equivalent_coeffs = frc.pair_coeffs_12_6 # class1 equivalent coeffs
        pair_auto_equiv_coeffs = frc.pair_coeffs_12_6 # class1 auto-equivalent coeffs
    if ff_class == 1:
        mass_comment = 'class1' # class1 Header comment
        pair_comment = 'lj/cut/coul/long' # class1 Header comment
        equivalents = frc.equivalences # class1 equivalents dict
        auto_equivs = frc.auto_equivalences # class1 auto-equivalents dict
        pair_equivalent_coeffs = frc.pair_coeffs_12_6 # class1 equivalent coeffs
        pair_auto_equiv_coeffs = frc.pair_coeffs_12_6 # class1 auto-equivalent coeffs
    elif ff_class == 2:
        mass_comment = 'class2' # class2 Header comment
        pair_comment = 'lj/class2/coul/long' # class2 Header comment
        equivalents = frc.equivalences # class2 equivalents dict
        auto_equivs = frc.auto_equivalences # class2 auto-equivalents dict
        pair_equivalent_coeffs = frc.pair_coeffs_9_6 # class2 equivalent coeffs
        pair_auto_equiv_coeffs = frc.pair_coeffs_9_6 # class2 auto-equivalent coeffs
    if ff_class == 'd':
        mass_comment = 'class1' # class1 Header comment
        pair_comment = 'lj/cut/coul/long' # class1 Header comment
        equivalents = frc.equivalences # class1 equivalents dict
        auto_equivs = {} # class1 auto-equivalents dict
        pair_equivalent_coeffs = frc.pair_coeffs_12_6 # class1 equivalent coeffs
        pair_auto_equiv_coeffs = {} # class1 auto-equivalent coeffs
    if ff_class in ['s1', 's2']:
        mass_comment = 'Skeleton: {}'.format(ff_class) # skeleton Header comment
        pair_comment = 'Skeleton: {}'.format(ff_class) # skeleton Header comment
    
    # Find pair coeffs
    for i in BADI.atom_types_lst:
        atom_type = i # set as "fullname" 
        if ':' in atom_type:
            i = atom_type[:atom_type.rfind(':')] # strip any ':NAME ending'
        
        # atom number
        number = BADI.atom_types_dict[atom_type]
        mass_equiv = ''; pair_equiv = '';
        mass_coeff_comment = ''
        pair_coeff_comment = ''
        
        # if ff_class 's1' or 's2' log generic and continue to next iterations
        if ff_class in ['s1', 's2']:
            # Save mass info into class
            t1 = Type()
            t1.type = i
            t1.equivalent = 'N/A'
            t1.coeffs = 0
            t1.comments = 'N/A'
            masses[number] = t1
            
            # Save pair_coeff info into class
            t2 = Type()
            t2.type = i
            t2.equivalent = 'N/A'
            t2.coeffs = []
            t2.comments = 'N/A'
            pair_coeffs[number] = t2
            continue # continue to next iteration
            
        
        ##########################################################
        # Find mass, try type 1st, then equiv, else set as zeros #
        ##########################################################
        # try for standard type
        if i in frc.atom_types:
            mass_equiv = i
            mass = frc.atom_types[i].mass
            mass_coeff_comment = 'standard type used'
        # next try for equivalence
        elif i in equivalents:
            # Find equivalence
            mass_equiv = equivalents[i].nonb
            if mass_equiv in frc.atom_types:
                mass = frc.atom_types[mass_equiv].mass
                mass_coeff_comment = 'equivalent type used'
            else:
                mass = 0.0; mass_equiv = 'N/A';
                mass_coeff_comment = 'UNABLE to find coeff parameters'
                if not skip_printouts: log.warn('WARNING unable to find atom info for Masses {} {}'.format(number, i))
        # next try for auto equivalence
        elif i in auto_equivs:
            # Find auto equivalence
            mass_equiv = auto_equivs[i].nonb
            if mass_equiv in frc.atom_types:
                mass = frc.atom_types[mass_equiv].mass
                mass_coeff_comment = 'auto equivalent type used'
            else:
                mass = 0.0; mass_equiv = 'N/A';
                mass_coeff_comment = 'UNABLE to find coeff parameters'
                if not skip_printouts: log.warn('WARNING unable to find atom info for Masses {} {}'.format(number, i))
        # else log failure
        else:
            mass = 0.0; mass_equiv = 'N/A';
            mass_coeff_comment = 'UNABLE to find coeff parameters'
            if not skip_printouts: log.warn('WARNING unable to find atom info for Masses {} {}'.format(number, i))
            
        
        ################################################################
        # Find pair_coeff, try type 1st, then equiv, else set as zeros #
        ################################################################
        # try for standard type
        if i in pair_equivalent_coeffs:
            pair_equiv = i
            pair_coeff = pair_equivalent_coeffs[i]
            pair_coeff_comment = 'standard type used'
            
            
            if ff_class == 0:                
                A = pair_coeff.A; B = pair_coeff.B;  
                if A > 0: one = (B*B)/(4.0*A);
                else: one = 0.0
                if B > 0: two = (A/B)**(1.0/6.0);
                else: two = 0.0
            elif ff_class == 1:                
                A = pair_coeff.A; B = pair_coeff.B;  
                if A > 0: one = (B*B)/(4.0*A);
                else: one = 0.0 
                if B > 0: two = (A/B)**(1.0/6.0);
                else: two = 0.0
            elif ff_class == 2:
                one = pair_coeff.eps; two = pair_coeff.r;
            elif ff_class == 'd':                
                A = pair_coeff.A; B = pair_coeff.B;  
                if A > 0: one = (B*B)/(4.0*A);
                else: one = 0.0 
                if B > 0: two = (A/B)**(1.0/6.0);
                else: two = 0.0
                
        # next try for equivalence
        elif i in equivalents:
            # Find equivalence
            pair_equiv = equivalents[i].nonb
            if pair_equiv in pair_equivalent_coeffs:
                pair_coeff = pair_equivalent_coeffs[pair_equiv]
                pair_coeff_comment = 'equivalent type used'
                
                if ff_class == 0:                
                    A = pair_coeff.A; B = pair_coeff.B;  
                    if A > 0: one = (B*B)/(4.0*A);
                    else: one = 0.0
                    if B > 0: two = (A/B)**(1.0/6.0);
                    else: two = 0.0
                elif ff_class == 1:                
                    A = pair_coeff.A; B = pair_coeff.B;  
                    if A > 0: one = (B*B)/(4.0*A);
                    else: one = 0.0 
                    if B > 0: two = (A/B)**(1.0/6.0);
                    else: two = 0.0
                elif ff_class == 2:
                    one = pair_coeff.eps; two = pair_coeff.r;
                elif ff_class == 'd':                
                    A = pair_coeff.A; B = pair_coeff.B;  
                    if A > 0: one = (B*B)/(4.0*A);
                    else: one = 0.0 
                    if B > 0: two = (A/B)**(1.0/6.0);
                    else: two = 0.0
                    
            else:
                one = 0.0; two = 0.0; pair_equiv = 'N/A';
                pair_coeff_comment = 'UNABLE to find coeff parameters'
                if not skip_printouts: log.warn('WARNING unable to find atom info for Pair Coeff {} {}'.format(number, i))
                
        # next try for auto_equivalence
        elif i in auto_equivs and use_auto_equivalence:
            # Find auto equivalence
            pair_equiv = auto_equivs[i].nonb
            if pair_equiv in pair_auto_equiv_coeffs:
                pair_coeff = pair_auto_equiv_coeffs[pair_equiv]
                pair_coeff_comment = 'auto equivalent type used'
                
                if ff_class == 0:                
                    A = pair_coeff.A; B = pair_coeff.B;  
                    one = (B*B)/(4.0*A); two = (A/B)**(1.0/6.0);
                elif ff_class == 1:                
                    A = pair_coeff.A; B = pair_coeff.B;  
                    one = (B*B)/(4.0*A); two = (A/B)**(1.0/6.0);
                elif ff_class == 2:
                    one = pair_coeff.eps; two = pair_coeff.r;
                elif ff_class == 'd':                
                    A = pair_coeff.A; B = pair_coeff.B;  
                    one = (B*B)/(4.0*A); two = (A/B)**(1.0/6.0);
                    
            else:
                one = 0.0; two = 0.0; pair_equiv = 'N/A';
                pair_coeff_comment = 'UNABLE to find coeff parameters'
                if not skip_printouts: log.warn('WARNING unable to find atom info for Pair Coeff {} {}'.format(number, i))
                
        # else log failure
        else:
            one = 0.0; two = 0.0; pair_equiv = 'N/A';
            pair_coeff_comment = 'UNABLE to find coeff parameters'
            if not skip_printouts: log.warn('WARNING unable to find atom info for Pair Coeff {} {}'.format(number, i))
                        
        # Save mass info into class
        t1 = Type()
        t1.type = atom_type
        t1.equivalent = mass_equiv
        t1.coeffs = mass
        t1.comments = mass_coeff_comment
        masses[number] = t1
        
        # Save pair_coeff info into class
        t2 = Type()
        t2.type = atom_type
        t2.equivalent = pair_equiv
        t2.coeffs = [one, two]
        t2.comments = pair_coeff_comment
        pair_coeffs[number] = t2   
    return masses, pair_coeffs, mass_comment, pair_comment


###################################
# Function for finding bond types #
###################################
def find_bond_parameters(frc, BADI, use_auto_equivalence, use_assumed_auto_fill, aafc, ff_class, use_morse_bonds, skip_printouts, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html
    
    Bonds
    The style of bond potential is specified in the input command file.
    
    (1) harmonic
    
      E = K (r - r0)^2
    
      standard harmonic spring
    
      r = distance (computed by LAMMPS)
    
      coeff1 = K (energy/distance^2)  (the usual 1/2 is included in the K)
      coeff2 = r0 (distance)
    
      2 coeffs are listed in data file or set in input script
    
    (5) class2
    
      E = K2 (r - r0)^2  +  K3 (r - r0)^3  +  K4 (r - r0)^4
    
      r = distance (computed by LAMMPS)
    
      coeff1 = r0 (distance)
      coeff2 = K2 (energy/distance^2)
      coeff3 = K3 (energy/distance^3)
      coeff4 = K4 (energy/distance^4)
    
      4 coeffs are listed in data file - cannot be set in input script
    """
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # Bond coeffs dictionary
    bond_coeffs = {} # {bond type : list of coeffs}  {1: [340,1.5], 2: [450,1.2], ...}
    bond_map = {} # { Ordered tuple(type1, type2) : numeric id }
    
    # Set comments, equivs/auto coeff dicts and equivs/auto mapping dicts
    if ff_class == 0:
        bond_comment = 'harmonic' # class1 Header comment
        equivalents = frc.equivalences # class1 equivalents dict
        auto_equivs = frc.auto_equivalences # class1 auto-equivalents dict
        equivalent_coeffs = frc.quadratic_bonds # class1 equivalent coeffs (both equiv and auto-equiv types will be in this dict)
        auto_equiv_coeffs = frc.quadratic_bonds_auto # class1 auto-equivalent coeffs
    if ff_class == 1 or ff_class == 'd':
        bond_comment = 'harmonic' # class1 Header comment
        equivalents = frc.equivalences # class1 equivalents dict
        auto_equivs = frc.auto_equivalences # class1 auto-equivalents dict
        equivalent_coeffs = frc.quadratic_bonds # class1 equivalent coeffs (both equiv and auto-equiv types will be in this dict)
        auto_equiv_coeffs = frc.quadratic_bonds_auto # class1 auto-equivalent coeffs
        
        # Set morse bond options
        if use_morse_bonds:
            bond_comment = 'morse' # class1 Header comment
            equivalents = frc.equivalences # class1 equivalents dict
            auto_equivs = frc.auto_equivalences # class1 auto-equivalents dict
            equivalent_coeffs = frc.morse_bonds # class1 equivalent coeffs (both equiv and auto-equiv types will be in this dict)
            auto_equiv_coeffs = frc.morse_bonds_auto # class1 auto-equivalent coeffs
    elif ff_class == 2:
        bond_comment = 'class2' # class2 Header comment
        equivalents = frc.equivalences # class2 equivalents dict
        auto_equivs = frc.auto_equivalences # class2 auto-equivalents dict
        equivalent_coeffs = frc.quartic_bonds # class2 equivalent coeffs
        auto_equiv_coeffs = frc.quadratic_bonds # class2 auto-equivalent coeffs
        
    if ff_class in ['s1', 's2']:
        bond_comment = 'Skeleton: {}'.format(ff_class) # skeleton Header comment
    
    ####################
    # Find bond coeffs #
    ####################
    for i in BADI.bond_types_lst:
        
        # Flag things
        equivalent_flag = False
        auto_equivalent_flag = False
        add_onto_assumed_file_flag = False
        assumed_equivalent_flag = False
        assumed_auto_equivalent_flag = False
        
        # bonded atom id types
        type1 = i[0]; type2 = i[1];
        
        # Set order/equiv/match and update later if found
        order = (type1, type2)
        equiv = (type1, type2)
        match = (type1, type2)
        
        # bond number
        number = BADI.bond_types_dict[i]
        
        # if ff_class 's1' or 's2' log generic and continue to next iterations
        if ff_class in ['s1', 's2']:
            # set order as type1, type2
            order = (type1, type2)
            t = Type()
            t.type = order
            t.match = ['N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A']
            t.coeffs = []
            t.comments = 'Skeleton option used'
            bond_coeffs[number] = t
            bond_map[number] = order
            continue
        
        # Try matching equivalent bonds in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
        if ff_functions.match_2_body(type1, type2, log, equivalent_coeffs, equivalences=equivalents, form='equiv')[0]:
            boolean, match, order, equiv = ff_functions.match_2_body(type1, type2, log, equivalent_coeffs, equivalences=equivalents, form='equiv')
            bond_coeff = equivalent_coeffs[match]
            equivalent_flag = True
            
        # Try matching auto equivalent bonds in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
        elif ff_functions.match_2_body(type1, type2, log, auto_equiv_coeffs, equivalences=auto_equivs, form='auto-equiv')[0] and use_auto_equivalence:
            boolean, match, order, equiv = ff_functions.match_2_body(type1, type2, log, auto_equiv_coeffs, equivalences=auto_equivs, form='auto-equiv')
            bond_coeff = auto_equiv_coeffs[match]
            auto_equivalent_flag = True
            
        # Try matching with assumed types if user wants and quadratic and quartic fails
        elif use_assumed_auto_fill:
            # Find basic elements in coeff and see if they exists in the assumed auto fill types
            if type1 in frc.atom_types and type2 in frc.atom_types:
                element1 = frc.atom_types[type1].element
                element2 = frc.atom_types[type2].element
                
                # Find assumed types
                if (element1, element2) in aafc.bond_coeffs:
                    assumed1, assumed2 = aafc.bond_coeffs[(element1, element2)]
                elif (element2, element1) in aafc.bond_coeffs:
                    assumed2, assumed1  = aafc.bond_coeffs[(element2, element1)]
                else:
                    add_onto_assumed_file_flag = True
                
                # Try using assumed1 and assumed2 if they were in assumed types (set equivalences to False since they are already in equivalent form)
                if not add_onto_assumed_file_flag:
                    
                    # Function to remap assumed types back onto real types
                    def assumed_bondtypes_order_remap2_realtypes_order(assumed_order, type1, type2, log):
                        if assumed_order == (assumed1, assumed2): order = (type1, type2)
                        elif assumed_order == (assumed2, assumed1): order = (type2, type1)
                        else:
                            order = (type1, type2)
                            if not skip_printouts: log.warn('WARNING: assumed auto fill option used for Bond Coeff {} {} {}, but ordering of real types COULD NOT BE REMAPPED from assumed types. DEFAULT ORDERING USED.'.format(number, type1, type2))
                        return order
                        
                    # Try matching assumed type with equivalent bonds
                    if ff_functions.match_2_body(assumed1, assumed2, log, equivalent_coeffs, equivalences=False, form='equiv')[0]:
                        boolean, match, aa_order, equiv = ff_functions.match_2_body(assumed1, assumed2, log, equivalent_coeffs, equivalences=False, form='equiv')
                        bond_coeff = equivalent_coeffs[match]
                        assumed_equivalent_flag = True
                        
                        # Find correct order of real types based on assumed type ordering (aa_order)
                        order = assumed_bondtypes_order_remap2_realtypes_order(aa_order, type1, type2, log)

                    # Try matching assumed type with auto equivalent bonds  (set equivalences to False since they are already in equivalent form)
                    elif ff_functions.match_2_body(assumed1, assumed2, log, auto_equiv_coeffs, equivalences=False, form='auto-equiv')[0]:
                        boolean, match, aa_order, equiv = ff_functions.match_2_body(assumed1, assumed2, log, auto_equiv_coeffs, equivalences=False, form='auto-equiv')
                        bond_coeff = auto_equiv_coeffs[match]
                        assumed_auto_equivalent_flag = True

                        # Find correct order of real types based on assumed type ordering (aa_order)
                        order = assumed_bondtypes_order_remap2_realtypes_order(aa_order, type1, type2, log)
                            
                            
        # Save bond_coeff info into class for equivalent form
        if equivalent_flag and not auto_equivalent_flag:
            # Set coeff based on FF class
            if ff_class == 0:
                coeff = [bond_coeff.k2, bond_coeff.r0]
            elif ff_class == 1 or ff_class == 'd':
                if not use_morse_bonds:
                    coeff = [bond_coeff.k2, bond_coeff.r0]
                else:
                    coeff = [bond_coeff.d, bond_coeff.alpha, bond_coeff.r0]
            elif ff_class == 2:
                coeff = [bond_coeff.r0, bond_coeff.k2, bond_coeff.k3, bond_coeff.k4]
                    
            t = Type()
            t.type = order
            t.match = match
            t.equivalent = equiv
            t.coeffs = coeff
            t.comments = 'equivalent type used'
            bond_coeffs[number] = t
            bond_map[number] = order
        
        # Save bond_coeff info into class for auto equivalent form
        elif auto_equivalent_flag and not equivalent_flag:
            # Set coeff based on FF class
            if ff_class == 0:
                coeff = [bond_coeff.k2, bond_coeff.r0]
            elif ff_class == 1 or ff_class == 'd':
                if not use_morse_bonds:
                    coeff = [bond_coeff.k2, bond_coeff.r0]
                else:
                    coeff = [bond_coeff.d, bond_coeff.alpha, bond_coeff.r0]
            elif ff_class == 2:
                coeff = [bond_coeff.r0, bond_coeff.k2, 0.0, 0.0]
            
            t = Type()
            t.type = order
            t.match = match
            t.equivalent = equiv
            t.coeffs = coeff
            t.comments = 'auto equivalent type used'
            bond_coeffs[number] = t
            bond_map[number] = order
            
        # Save bond_coeff info into class for assumed auto fill option
        elif use_assumed_auto_fill and assumed_equivalent_flag or assumed_auto_equivalent_flag:
            if not skip_printouts: log.warn('WARNING unable to find bond info for Bond Coeff {} {} {}, but assumed auto fill option filled in with assumed type'.format(number, type1, type2))
        
            # If equivalent coeffs where used (set order as type1, type2)
            if assumed_equivalent_flag:
                
                # Set coeff based on FF class
                if ff_class == 0:
                    coeff = [bond_coeff.k2, bond_coeff.r0]
                elif ff_class == 1 or ff_class == 'd':
                    if not use_morse_bonds:
                        coeff = [bond_coeff.k2, bond_coeff.r0]
                    else:
                        coeff = [bond_coeff.d, bond_coeff.alpha, bond_coeff.r0]
                elif ff_class == 2:
                    coeff = [bond_coeff.r0, bond_coeff.k2, bond_coeff.k3, bond_coeff.k4]
                
                t = Type()
                t.type = order
                t.match = match
                t.equivalent = equiv
                t.coeffs = coeff
                t.comments = '{:^10} {:2} {:2}'.format('ASSUMED equivalent type used', assumed1, assumed2) 
                bond_coeffs[number] = t
                bond_map[number] = order
                
            # elif auto equivalent coeffs where used (set order as type1, type2)
            elif assumed_auto_equivalent_flag:
                # Set coeff based on FF class
                if ff_class == 0:
                    coeff = [bond_coeff.k2, bond_coeff.r0]
                elif ff_class == 1 or ff_class == 'd':
                    if not use_morse_bonds:
                        coeff = [bond_coeff.k2, bond_coeff.r0]
                    else:
                        coeff = [bond_coeff.d, bond_coeff.alpha, bond_coeff.r0]       
                elif ff_class == 2:
                    coeff = [bond_coeff.r0, bond_coeff.k2, 0.0, 0.0]
                
                t = Type()
                t.type = order
                t.match = match
                t.equivalent = equiv
                t.coeffs = coeff
                t.comments = '{:^10} {:2} {:2}'.format('ASSUMED auto equivalent type used', assumed1, assumed2)
                bond_coeffs[number] = t 
                bond_map[number] = order
            
        # Save bond_coeff info into class as zeros if all attempts fail
        else:
            # If this block runs due to assumed class being used, but missing parameters in file let user know
            if add_onto_assumed_file_flag:
                if not skip_printouts: log.warn('WARNING unable to find bond info for Bond Coeff {} {} {}, assumed_auto_fill was used but currently {} {} coeff is not in the file, please add to file'.format(number, type1, type2, element1, element2))
            else:
                if not skip_printouts: log.warn('WARNING unable to find bond info for Bond Coeff {} {} {}'.format(number, type1, type2))
                
            # Set coeff based on FF class
            if ff_class == 0:
                coeff = 2*[0.0]
            elif ff_class == 1 or ff_class == 'd':
                if not use_morse_bonds:
                    coeff = 2*[0.0]
                else:
                    coeff = 3*[0.0]
            elif ff_class == 2:
                coeff = 4*[0.0]
            
            # set order as type1, type2
            order = (type1, type2)
            t = Type()
            t.type = order
            t.match = ['N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A']
            t.coeffs = coeff
            t.comments = 'UNABLE to find coeff parameters'
            bond_coeffs[number] = t
            bond_map[number] = order 
    return bond_coeffs, bond_map, bond_comment


##############################
# Function for finding bonds #
##############################
def find_bonds(nta, BADI, bond_map, sort_remap_atomids, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html
    """
    
    # Dictionaries to add info to
    bonds = {}   # {bond number : bond object}
    
    # Find bonds info
    for n, (id1, id2) in enumerate(BADI.bonds):
        type1 = nta[id1]; type2 = nta[id2];
        
        # Try to find type in forward order 1st
        if (type1, type2) in BADI.bond_types_dict:
            number = BADI.bond_types_dict[(type1, type2)]
        elif (type2, type1) in BADI.bond_types_dict:
            number = BADI.bond_types_dict[(type2, type1)]
        else:
            log.error('ERROR Finding Bond Type Number ID based on atom-id ordering')
        
        # Find bond ordering and set bond_type and atomids
        if sort_remap_atomids:
            ordering = bond_map[number]
            if ordering == (type1, type2):
                bond_type = (type1, type2)
                atomids = [id1, id2]
            elif ordering == (type2, type1):
                bond_type = (type2, type1)
                atomids = [id2, id1]
            else:
                log.error('ERROR Finding Bond Type for bonds ordering')
        # if not sort_remap_atomids use orginal ordering
        else:
            bond_type = (type1, type2)
            atomids = [id1, id2]
        
        # Save bond info into class
        b = Bond()
        b.type = number
        b.symbol = bond_type
        b.atomids = atomids
        bonds[n+1] = b     
    return bonds


####################################
# Function for finding angle types #
####################################
def find_angle_parameters(frc, BADI, use_auto_equivalence, use_assumed_auto_fill, aafc, ff_class, skip_printouts, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html

    Angles
    The style of angle potential is specified in the input command file.
    
    (1) harmonic
    
      E = K (theta - theta0)^2
    
      theta = radians (computed by LAMMPS)
    
      coeff1 = K (energy/radian^2) (the usual 1/2 is included in the K)
      coeff2 = theta0 (degrees) (converted to radians within LAMMPS)
    
      2 coeffs are listed in data file
    
    (2) class2
    
      E = K2 (theta - theta0)^2 +  K3 (theta - theta0)^3 + 
           K4 (theta - theta0)^4
    
      theta = radians (computed by LAMMPS)
    
      coeff1 = theta0 (degrees) (converted to radians within LAMMPS)
      coeff2 = K2 (energy/radian^2)
      coeff3 = K3 (energy/radian^3)
      coeff4 = K4 (energy/radian^4)
    
      4 coeffs are listed in data file
    """
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # Angle coeffs dictionary
    angle_coeffs = {} # {angle type : list of coeffs}
    angle_map = {} # { Ordered tuple(type1, type2, type3) : numeric id }
        
    # Set comments, equivs/auto coeff dicts and equivs/auto mapping dicts
    if ff_class == 0:
        angle_comment = 'harmonic' # class1 Header comment
        equivalents = frc.equivalences # class1 equivalents dict
        auto_equivs = frc.auto_equivalences # class1 auto-equivalents dict
        equivalent_coeffs = frc.quadratic_angles # class1 equivalent coeffs (both equiv and auto-equiv types will be in this dict)
        auto_equiv_coeffs = frc.quadratic_angles_auto # class1 auto-equivalent coeffs
    elif ff_class == 1 or ff_class == 'd':
        angle_comment = 'harmonic' # class1 Header comment
        equivalents = frc.equivalences # class1 equivalents dict
        auto_equivs = frc.auto_equivalences # class1 auto-equivalents dict
        equivalent_coeffs = frc.quadratic_angles # class1 equivalent coeffs (both equiv and auto-equiv types will be in this dict)
        auto_equiv_coeffs = frc.quadratic_angles_auto # class1 auto-equivalent coeffs
    elif ff_class == 2:
        angle_comment = 'class2' # class2 Header comment
        equivalents = frc.equivalences # class2 equivalents dict
        auto_equivs = frc.auto_equivalences # class2 auto-equivalents dict
        equivalent_coeffs = frc.quartic_angles # class2 equivalent coeffs
        auto_equiv_coeffs = frc.quadratic_angles # class2 auto-equivalent coeffs
        
    if ff_class in ['s1', 's2']:
        angle_comment = 'Skeleton: {}'.format(ff_class) # skeleton Header comment
        
    #####################
    # Find angle coeffs #
    #####################
    for i in BADI.angle_types_lst:
        
        # Flag things
        equivalent_flag = False
        auto_equivalent_flag = False
        add_onto_assumed_file_flag = False
        assumed_equivalent_flag = False
        assumed_auto_equivalent_flag = False

        
        # bonded atom id types
        type1 = i[0]; type2 = i[1]; type3 = i[2];
        
        # Set order/equiv/match and update later if found
        order = (type1, type2, type3)
        equiv = (type1, type2, type3)
        match = (type1, type2, type3)
        
        # angle number
        number = BADI.angle_types_dict[i]
        
		# if ff_class 's1' or 's2' log generic and continue to next iterations
        if ff_class in ['s1', 's2']:
            # set order as type1, type2, type3
            order = (type1, type2, type3)
            t = Type()
            t.type = order
            t.match = ['N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A']
            t.coeffs = []
            t.comments = 'Skeleton option used'
            angle_coeffs[number] = t
            angle_map[number] = order
            continue
        
        # Try matching equivalent angles in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
        if ff_functions.match_3_body(type1, type2, type3, log, equivalent_coeffs, equivalences=equivalents, wildcard_search=True, form='equiv')[0]:
            boolean, match, order, equiv = ff_functions.match_3_body(type1, type2, type3, log, equivalent_coeffs, equivalences=equivalents, wildcard_search=True, form='equiv')
            angle_coeff = equivalent_coeffs[match]
            equivalent_flag = True
            
        # Try matching auto equivalent angles in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences, 3) with equivalences and wild cards
        elif ff_functions.match_3_body(type1, type2, type3, log, auto_equiv_coeffs, equivalences=auto_equivs, wildcard_search=True, form='auto-equiv')[0] and use_auto_equivalence:
            boolean, match, order, equiv = ff_functions.match_3_body(type1, type2, type3, log, auto_equiv_coeffs, equivalences=auto_equivs, wildcard_search=True, form='auto-equiv')
            angle_coeff = auto_equiv_coeffs[match]
            auto_equivalent_flag = True

        # Try matching with assumed types if user wants and equivalent and auto equivalent fails
        elif use_assumed_auto_fill:
            # Find basic elements in coeff and see if they exists in the assumed auto fill types
            if type1 in frc.atom_types and type2 in frc.atom_types and type3 in frc.atom_types:
                element1 = frc.atom_types[type1].element
                element2 = frc.atom_types[type2].element
                element3 = frc.atom_types[type3].element
                
                # Find assumed types
                if (element1, element2, element3) in aafc.angle_coeffs:
                    assumed1, assumed2, assumed3 = aafc.angle_coeffs[(element1, element2, element3)]
                elif (element3, element2, element1) in aafc.angle_coeffs:
                    assumed3, assumed2, assumed1  = aafc.angle_coeffs[(element3, element2, element1)]
                else:
                    add_onto_assumed_file_flag = True
                
                # Try using assumed1, assumed2, and assumed3 if they were in assumed types (set equivalences to False since they are already in equivalent form)
                if not add_onto_assumed_file_flag:
                    
                    # Function to remap assumed types back onto real types
                    def assumed_angletypes_order_remap2_realtypes_order(assumed_order, type1, type2, type3, log):
                        if assumed_order == (assumed1, assumed2, assumed3):
                            order = (type1, type2, type3)
                        elif assumed_order == (assumed3, assumed2, assumed1):
                            order = (type3, type2, type1)
                        else:
                            order = (type1, type2, type3)
                            if not skip_printouts: log.warn('WARNING: assumed auto fill option used for Angle Coeff {} {} {} {}, but ordering of real types COULD NOT BE REMAPPED from assumed types. DEFAULT ORDERING USED.'.format(number, type1, type2, type3))
                        return order
                    
                    # Try matching assumed type with auto equivalent angles
                    if ff_functions.match_3_body(assumed1, assumed2, assumed3, log, equivalent_coeffs, equivalences=False, wildcard_search=True, form='equiv')[0]:
                        boolean, match, aa_order, equiv = ff_functions.match_3_body(assumed1, assumed2, assumed3, log, equivalent_coeffs, equivalences=False, wildcard_search=True, form='equiv')
                        angle_coeff = equivalent_coeffs[match]
                        assumed_equivalent_flag = True
                        
                        # Find correct order of real types based on assumed type ordering (aa_order)
                        order = assumed_angletypes_order_remap2_realtypes_order(aa_order, type1, type2, type3, log)
                        
                    # Try matching assumed type with auto equivalent angles (set equivalences to False since they are already in equivalent form)
                    elif ff_functions.match_3_body(assumed1, assumed2, assumed3, log, auto_equiv_coeffs, equivalences=False, wildcard_search=True, form='auto-equiv')[0]:
                        boolean, match, aa_order, equiv = ff_functions.match_3_body(assumed1, assumed2, assumed3, log, auto_equiv_coeffs, equivalences=False, wildcard_search=True, form='auto-equiv')
                        angle_coeff = auto_equiv_coeffs[match]
                        assumed_auto_equivalent_flag = True
                        
                        # Find correct order of real types based on assumed type ordering (aa_order)
                        order = assumed_angletypes_order_remap2_realtypes_order(aa_order, type1, type2, type3, log)
                        
                        
        # Save angle_coeff info into class for equivalent form
        if equivalent_flag and not auto_equivalent_flag:
            
            # Set coeff based on FF class
            if ff_class == 0:
                coeff = [angle_coeff.k2, angle_coeff.theta0]
            elif ff_class == 1 or ff_class == 'd':
                coeff = [angle_coeff.k2, angle_coeff.theta0]
            elif ff_class == 2:
                coeff = [angle_coeff.theta0, angle_coeff.k2, angle_coeff.k3, angle_coeff.k4]
            
            t = Type()
            t.type = order
            t.match = match
            t.equivalent = equiv
            t.coeffs = coeff
            t.comments = 'equivalent type used'
            angle_coeffs[number] = t
            angle_map[number] = order
        
        # Save angle_coeff info into class for auto equivalent form
        elif auto_equivalent_flag and not equivalent_flag:
            
            # Set coeff based on FF class
            if ff_class == 0:
                coeff = [angle_coeff.k2, angle_coeff.theta0]
            elif ff_class == 1 or ff_class == 'd':
                coeff = [angle_coeff.k2, angle_coeff.theta0]
            elif ff_class == 2:
                coeff = [angle_coeff.theta0, angle_coeff.k2, 0.0, 0.0]
            
            t = Type()
            t.type = order
            t.match = match
            t.equivalent = equiv
            t.coeffs = coeff
            t.comments = 'auto equivalent type used'
            angle_coeffs[number] = t
            angle_map[number] = order
            
        # Save angle_coeff info into class for assumed auto fill option
        elif use_assumed_auto_fill and assumed_equivalent_flag or assumed_auto_equivalent_flag:
            if not skip_printouts: log.warn('WARNING unable to find angle info for Angle Coeff {} {} {} {}, but assumed auto fill option filled in with assumed type'.format(number, type1, type2, type3))
        
            # If quartic coeffs where used (set order as type1, type2, type3)
            if assumed_equivalent_flag:
                
                # Set coeff based on FF class
                if ff_class == 0:
                    coeff = [angle_coeff.k2, angle_coeff.theta0]
                elif ff_class == 1:
                    coeff = [angle_coeff.k2, angle_coeff.theta0]
                elif ff_class == 2:
                    coeff = [angle_coeff.theta0, angle_coeff.k2, angle_coeff.k3, angle_coeff.k4]
                
                t = Type()
                t.type = order
                t.match = match
                t.equivalent = equiv
                t.coeffs = coeff
                t.comments = '{:^10} {:2} {:2} {:2}'.format('ASSUMED equivalent type used', assumed1, assumed2, assumed3) 
                angle_coeffs[number] = t
                angle_map[number] = order
                
            # elif quadratic coeffs where used (set order as type1, type2, type3)
            elif assumed_auto_equivalent_flag:
                
                # Set coeff based on FF class
                if ff_class == 0:
                    coeff = [angle_coeff.k2, angle_coeff.theta0]
                elif ff_class == 1 or ff_class == 'd':
                    coeff = [angle_coeff.k2, angle_coeff.theta0]
                elif ff_class == 2:
                    coeff = [angle_coeff.theta0, angle_coeff.k2, 0.0, 0.0]
                
                t = Type()
                t.type = order
                t.match = match
                t.equivalent = equiv
                t.coeffs = coeff
                t.comments = '{:^10} {:2} {:2} {:2}'.format('ASSUMED auto equivalent type used', assumed1, assumed2, assumed3)
                angle_coeffs[number] = t 
                angle_map[number] = order
            
        # Save angle_coeff info into class as zeros if all attempts fail
        else:
            # If this block runs due to assumed class being used, but missing parameters in file let user know
            if add_onto_assumed_file_flag:
                if not skip_printouts: log.warn('WARNING unable to find angle info for Angle Coeff {} {} {} {}, assumed_auto_fill was used but currently {} {} {} coeff is not in the file, please add to file'.format(number, type1, type2, type3, element1, element2, element3))
            else:
                if not skip_printouts: log.warn('WARNING unable to find angle info for Angle Coeff {} {} {} {}'.format(number, type1, type2, type3))
                
            # Set coeff based on FF class
            if ff_class == 0:
                coeff = 2*[0.0]
            elif ff_class == 1 or ff_class == 'd':
                coeff = 2*[0.0]
            elif ff_class == 2:
                coeff = 4*[0.0]
            
            # set order as type1, type2, type3
            order = (type1, type2, type3)
            t = Type()
            t.type = order
            t.match = ['N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A']
            t.coeffs = coeff
            t.comments = 'UNABLE to find coeff parameters'
            angle_coeffs[number] = t
            angle_map[number] = order         
    return angle_coeffs, angle_map, angle_comment


###############################
# Function for finding angles #
###############################
def find_angles(nta, BADI, angle_map, sort_remap_atomids, ff_class, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html
    """
    
    # Dictionaries to add info to
    angles = {}   # {angle number : angle object}
    
    # Find angles info
    for n, (id1, id2, id3) in enumerate(BADI.angles):
        type1 = nta[id1]; type2 = nta[id2]; type3 = nta[id3];
        
        # Try to find type in forward order 1st
        if (type1, type2, type3) in BADI.angle_types_dict:
            number = BADI.angle_types_dict[(type1, type2, type3)]
        elif (type3, type2, type1) in BADI.angle_types_dict:
            number = BADI.angle_types_dict[(type3, type2, type1)]
        elif ('X', type2, 'X') in BADI.angle_types_dict and ff_class == 'd':
            number = BADI.angle_types_dict[('X', type2, 'X')]
        else:
            log.error('ERROR Finding Angle Type Number ID based on atom-id ordering')
        
        # Find angle ordering and set angle_type and atomids
        if sort_remap_atomids:
            ordering = angle_map[number]
            if ordering == (type1, type2, type3):
                angle_type = (type1, type2, type3)
                atomids = [id1, id2, id3]
            elif ordering == (type3, type2, type1):
                angle_type = (type3, type2, type1)
                atomids = [id3, id2, id1]
            elif ordering == ('X', type2, 'X') and ff_class == 'd':
                angle_type = ('X', type2, 'X')
                atomids = [id1, id2, id3]     
            else:
                log.error('ERROR Finding Angle Type for angles ordering')
        # if not sort_remap_atomids use orginal ordering
        else:
            angle_type = (type1, type2, type3)
            atomids = [id1, id2, id3]
            
        # Save angle info into class
        a = Angle()
        a.type = number
        a.symbol = angle_type
        a.atomids = atomids
        angles[n+1] = a
    return angles


#######################################
# Function for finding dihedral types #
#######################################
def find_dihedral_parameters(frc, BADI, use_auto_equivalence, use_assumed_auto_fill, aafc, ff_class, skip_printouts, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html

    Dihedrals
    The style of dihedral potential is specified in the input command file.
    
    (1) harmonic
    
      E = K [1 + d * cos (n * phi) ]
    
      phi = radians (computed by LAMMPS)
    
      coeff1 = K (energy)
      coeff2 = d (always +1 or -1) (1=0 degrees and -1=180 degrees)
      coeff3 = n (1,2,3,4,6)
    
      Cautions when comparing to other force fields:
    
      some force fields reverse the sign convention on d so that
        E = K [1 - d * cos(n*phi)]
      some force fields divide/multiply K by the number of multiple
        torsions that contain the j-k bond in an i-j-k-l torsion
      some force fields let n be positive or negative which 
        corresponds to d = 1,-1
      in the LAMMPS force field, the trans position = 180 degrees, while
        in some force fields trans = 0 degrees
     
      3 coeffs are listed in data file
    (2) class2
    
      E = SUM(n=1,3) { K_n [ 1 - cos( n*Phi - Phi0_n ) ] }
    
      phi = radians (computed by LAMMPS)
    
      coeff1 = K_1 (energy)
      coeff2 = Phi0_1 (degrees) (converted to radians within LAMMPS)
      coeff3 = K_2 (energy)
      coeff4 = Phi0_2 (degrees) (converted to radians within LAMMPS)
      coeff5 = K_3 (energy)
      coeff6 = Phi0_3 (degrees) (converted to radians within LAMMPS)
    
      6 coeffs are listed in data file
    """
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # Dihedral coeffs dictionary
    dihedral_coeffs = {} # {dihedral type : list of coeffs}
    dihedral_map = {} # { Ordered tuple(type1, type2, type3, type4) : numeric id }
 
    # Set comments, equivs/auto coeff dicts and equivs/auto mapping dicts
    if ff_class == 0:
        dihedral_comment = 'opls' # opls Header comment
        equivalents = frc.equivalences # class1 equivalents dict
        auto_equivs = frc.auto_equivalences # class1 auto-equivalents dict
        equivalent_coeffs = frc.torsion_1_opls # class1 equivalent coeffs
        auto_equiv_coeffs = frc.torsion_1_opls # class1 auto-equivalent coeffs
    elif ff_class == 1 or ff_class == 'd':
        dihedral_comment = 'harmonic' # class1 Header comment
        equivalents = frc.equivalences # class1 equivalents dict
        auto_equivs = frc.auto_equivalences # class1 auto-equivalents dict
        equivalent_coeffs = frc.torsion_1 # class1 equivalent coeffs
        auto_equiv_coeffs = frc.torsion_1_auto # class1 auto-equivalent coeffs
    elif ff_class == 2:
        dihedral_comment = 'class2'  # class2 Header comment
        equivalents = frc.equivalences # class2 equivalents dict
        auto_equivs = frc.auto_equivalences # class2 auto-equivalents dict
        equivalent_coeffs = frc.torsion_3 # class2 equivalent coeffs
        auto_equiv_coeffs = frc.torsion_1 # class2 auto-equivalent coeffs
        
    if ff_class in ['s1', 's2']:
        dihedral_comment = 'Skeleton: {}'.format(ff_class) # skeleton Header comment

    ########################
    # Find dihedral coeffs #
    ########################
    for i in BADI.dihedral_types_lst:
        
        # Flag things
        equivalent_flag = False
        auto_equivalent_flag = False
        add_onto_assumed_file_flag = False
        assumed_equivalent_flag = False
        assumed_auto_equivalent_flag = False

        
        # bonded atom id types
        type1 = i[0]; type2 = i[1]; type3 = i[2]; type4 = i[3];
        
        # Set order/equiv/match and update later if found
        order = (type1, type2, type3, type4)
        equiv = (type1, type2, type3, type4)
        match = (type1, type2, type3, type4)
        
        # dihedral number
        number = BADI.dihedral_types_dict[i]
            
  		# if ff_class 's1' or 's2' log generic and continue to next iterations
        if ff_class in ['s1', 's2']:
            # set order as type1, type2, type3, type4
            order = (type1, type2, type3, type4)
            t = Type()
            t.type = order
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = []
            t.comments = 'Skeleton option used'
            dihedral_coeffs[number] = t
            dihedral_map[number] = order
            continue
        
        # Try matching equivalent dihedrals in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
        if ff_functions.match_4_body(type1, type2, type3, type4, log, equivalent_coeffs, equivalences=equivalents, wildcard_search=True, form='equiv')[0]:
            boolean, match, order, equiv = ff_functions.match_4_body(type1, type2, type3, type4, log, equivalent_coeffs, equivalences=equivalents, wildcard_search=True, form='equiv')
            dihedral_coeff = equivalent_coeffs[match]
            equivalent_flag = True
            
        # Try matching auto equivalent dihedrals in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences, 3) with equivalences and wild cards
        elif ff_functions.match_4_body(type1, type2, type3, type4, log, auto_equiv_coeffs, equivalences=auto_equivs, wildcard_search=True, form='auto-equiv')[0] and use_auto_equivalence:
            boolean, match, order, equiv = ff_functions.match_4_body(type1, type2, type3, type4, log, auto_equiv_coeffs, equivalences=auto_equivs, wildcard_search=True, form='auto-equiv')
            dihedral_coeff = auto_equiv_coeffs[match]
            auto_equivalent_flag = True

        # Try matching with assumed types if user wants and equivalent and auto equivalent fails
        elif use_assumed_auto_fill:
            # Find basic elements in coeff and see if they exists in the assumed auto fill types
            if type1 in frc.atom_types and type2 in frc.atom_types and type3 in frc.atom_types and type4 in frc.atom_types:
                element1 = frc.atom_types[type1].element
                element2 = frc.atom_types[type2].element
                element3 = frc.atom_types[type3].element
                element4 = frc.atom_types[type4].element
                
                # Find assumed types
                if (element1, element2, element3, element4) in aafc.dihedral_coeffs:
                    assumed1, assumed2, assumed3, assumed4 = aafc.dihedral_coeffs[(element1, element2, element3, element4)]
                elif (element4, element3, element2, element1) in aafc.dihedral_coeffs:
                    assumed4, assumed3, assumed2, assumed1  = aafc.dihedral_coeffs[(element4, element3, element2, element1)]
                else:
                    add_onto_assumed_file_flag = True
                    
                # Try using assumed1, assumed2, and assumed3 if they were in assumed types (set equivalences to False since they are already in equivalent form)
                if not add_onto_assumed_file_flag:
                    
                    # Function to remap assumed types back onto real types
                    def assumed_dihedraltypes_order_remap2_realtypes_order(assumed_order, type1, type2, type3, type4, log):
                        if assumed_order == (assumed1, assumed2, assumed3, assumed4):
                            order = (type1, type2, type3, type4)
                        elif assumed_order == (assumed4, assumed3, assumed2, assumed1):
                            order = (type4, type3, type2, type1)
                        else:
                            order = (type1, type2, type3, type4)
                            if not skip_printouts: log.warn('WARNING: assumed auto fill option used for Dihedral Coeff {} {} {} {} {}, but ordering of real types COULD NOT BE REMAPPED from assumed types. DEFAULT ORDERING USED.'.format(number, type1, type2, type3, type4))
                        return order
                    
                    # Try matching assumed type with equivalent dihedrals
                    if ff_functions.match_4_body(assumed1, assumed2, assumed3, assumed4, log, equivalent_coeffs, equivalences=False, wildcard_search=True, form='equiv')[0]:
                        boolean, match, aa_order, equiv = ff_functions.match_4_body(assumed1, assumed2, assumed3, assumed4, log, equivalent_coeffs, equivalences=False, wildcard_search=True, form='equiv')
                        dihedral_coeff = equivalent_coeffs[match]
                        assumed_equivalent_flag = True
                        
                        # Find correct order of real types based on assumed type ordering (aa_order)
                        order = assumed_dihedraltypes_order_remap2_realtypes_order(aa_order, type1, type2, type3, type4, log)
                        
                    # Try matching assumed type with auto equivalent dihedrals (set equivalences to False since they are already in equivalent form)
                    elif ff_functions.match_4_body(assumed1, assumed2, assumed3, assumed4, log, auto_equiv_coeffs, equivalences=False, wildcard_search=True, form='auto-equiv')[0]:
                        boolean, match, aa_order, equiv = ff_functions.match_4_body(assumed1, assumed2, assumed3, assumed4, log, auto_equiv_coeffs, equivalences=False, wildcard_search=True, form='auto-equiv')
                        dihedral_coeff = auto_equiv_coeffs[match]
                        assumed_auto_equivalent_flag = True
                    
                        # Find correct order of real types based on assumed type ordering (aa_order)
                        order = assumed_dihedraltypes_order_remap2_realtypes_order(aa_order, type1, type2, type3, type4, log)

                    
        # Save dihedral_coeff info into class for equivalent form
        if equivalent_flag and not auto_equivalent_flag:
            
            # Set coeff based on FF class
            if ff_class == 0:
                coeff = [dihedral_coeff.k1, dihedral_coeff.k2, dihedral_coeff.k3, dihedral_coeff.k4]
            elif ff_class == 1 or ff_class == 'd':
                coeff = ff_functions.dihedral_autoequiv(dihedral_coeff)
            elif ff_class == 2:
                coeff = [dihedral_coeff.v1, dihedral_coeff.phi1, dihedral_coeff.v2, dihedral_coeff.phi2, dihedral_coeff.v3, dihedral_coeff.phi3]
            
            t = Type()
            t.type = order
            t.match = match
            t.equivalent = equiv
            t.coeffs = coeff
            t.comments = 'equivalent type used'
            dihedral_coeffs[number] = t
            dihedral_map[number] = order
            
        # Save dihedral_coeff info into class for auto equivalent form
        elif auto_equivalent_flag and not equivalent_flag:
            
            # Set coeff based on FF class
            if ff_class == 0:
                coeff = [dihedral_coeff.k1, dihedral_coeff.k2, dihedral_coeff.k3, dihedral_coeff.k4]
            elif ff_class == 1 or ff_class == 'd':
                coeff = ff_functions.dihedral_autoequiv(dihedral_coeff)
            elif ff_class == 2:
                coeff = ff_functions.dihedral_autoequiv2equiv(dihedral_coeff)
            
            t = Type()
            t.type = order
            t.match = match
            t.equivalent = equiv
            t.coeffs = coeff
            t.comments = 'auto equivalent type used'
            dihedral_coeffs[number] = t
            dihedral_map[number] = order
            
        # Save dihedral_coeff info into class for assumed auto fill option
        elif use_assumed_auto_fill and assumed_equivalent_flag or assumed_auto_equivalent_flag:
            if not skip_printouts: log.warn('WARNING unable to find dihedral info for Dihedral Coeff {} {} {} {} {}, but assumed auto fill option filled in with assumed type'.format(number, type1, type2, type3, type4))
        
            # If equivalent coeffs where used (set order as type1, type2, type3, type4)
            if assumed_equivalent_flag:
                
                # Set coeff based on FF class
                if ff_class == 0:
                    coeff = [dihedral_coeff.k1, dihedral_coeff.k2, dihedral_coeff.k3, dihedral_coeff.k4]
                elif ff_class == 1 or ff_class == 'd':
                    coeff = ff_functions.dihedral_autoequiv(dihedral_coeff)
                elif ff_class == 2:
                    coeff = [dihedral_coeff.v1, dihedral_coeff.phi1, dihedral_coeff.v2, dihedral_coeff.phi2, dihedral_coeff.v3, dihedral_coeff.phi3]
                
                t = Type()
                t.type = order
                t.match = match
                t.equivalent = equiv
                t.coeffs = coeff
                t.comments = '{:^10} {:2} {:2} {:2} {:2}'.format('ASSUMED equivalent type used', assumed1, assumed2, assumed3, assumed4) 
                dihedral_coeffs[number] = t
                dihedral_map[number] = order
                
            # elif auto equivalent coeffs where used (set order as type1, type2, type3, type4)
            elif assumed_auto_equivalent_flag:
                
                # Set coeff based on FF class
                if ff_class == 0:
                    coeff = [dihedral_coeff.k1, dihedral_coeff.k2, dihedral_coeff.k3, dihedral_coeff.k4]
                elif ff_class == 1:
                    coeff = ff_functions.dihedral_autoequiv(dihedral_coeff)
                elif ff_class == 2:
                    coeff = ff_functions.dihedral_autoequiv2equiv(dihedral_coeff)
                
                t = Type()
                t.type = order
                t.match = match
                t.equivalent = equiv
                t.coeffs = coeff
                t.comments = '{:^10} {:2} {:2} {:2} {:2}'.format('ASSUMED auto equivalent type used', assumed1, assumed2, assumed3, assumed4)
                dihedral_coeffs[number] = t 
                dihedral_map[number] = order
            
        # Save dihedral_coeff info into class as zeros if all attempts fail
        else:
            # If this block runs due to assumed class being used, but missing parameters in file let user know
            if add_onto_assumed_file_flag:
                if not skip_printouts: log.warn('WARNING unable to find dihedral info for Dihedral Coeff {} {} {} {} {}, assumed_auto_fill was used but currently {} {} {} {} coeff is not in the file, please add to file'.format(number, type1, type2, type3, type4, element1, element2, element3, element4))
            else:
                if not skip_printouts: log.warn('WARNING unable to find dihedral info for Dihedral Coeff {} {} {} {} {}'.format(number, type1, type2, type3, type4))
                
            # Set coeff based on FF class
            if ff_class == 0:
                coeff = 4*[0.0]
            elif ff_class == 1 or ff_class == 'd':
                coeff = [0.0, 1, 0] # int vs float will be used when writing coeffs
            elif ff_class == 2:
                coeff = 6*[0.0]
            
            # set order as type1, type2, type3, type4
            order = (type1, type2, type3, type4)
            t = Type()
            t.type = order
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = coeff
            t.comments = 'UNABLE to find coeff parameters'
            dihedral_coeffs[number] = t
            dihedral_map[number] = order 
    return dihedral_coeffs, dihedral_map, dihedral_comment


##################################
# Function for finding Dihedrals #
##################################
def find_dihedrals(nta, BADI, dihedral_map, sort_remap_atomids, ff_class, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html
    """
    
    # Dictionaries to add info to
    dihedrals = {}   # {dihedral number : dihedral object}
    
    # Find angles info
    for n, (id1, id2, id3, id4) in enumerate(BADI.dihedrals):
        type1 = nta[id1]; type2 = nta[id2]; type3 = nta[id3]; type4 = nta[id4];
        
        # Try to find type in forward order 1st
        if (type1, type2, type3, type4) in BADI.dihedral_types_dict:
            number = BADI.dihedral_types_dict[(type1, type2, type3, type4)]
        elif (type4, type3, type2, type1) in BADI.dihedral_types_dict:
            number = BADI.dihedral_types_dict[(type4, type3, type2, type1)]
        elif ('X', type2, type3, 'X') in BADI.dihedral_types_dict and ff_class == 'd':
            number = BADI.dihedral_types_dict[('X', type2, type3, 'X')]
        elif ('X', type3, type2, 'X') in BADI.dihedral_types_dict and ff_class == 'd':
            number = BADI.dihedral_types_dict[('X', type3, type2, 'X')]
        else:
            log.error('ERROR Finding Dihedral Type Number ID based on atom-id ordering')
        
        # Find dihedral ordering and set dihedral_type and atomids
        if sort_remap_atomids:
            ordering = dihedral_map[number]
            if ordering == (type1, type2, type3, type4):
                dihedral_type = (type1, type2, type3, type4)
                atomids = [id1, id2, id3, id4]
            elif ordering == (type4, type3, type2, type1):
                dihedral_type = (type4, type3, type2, type1)
                atomids = [id4, id3, id2, id1]
            elif ordering == ('X', type2, type3, 'X') and ff_class == 'd':
                dihedral_type = ('X', type2, type3, 'X')
                atomids = [id1, id2, id3, id4]
            elif ordering == ('X', type3, type2, 'X') and ff_class == 'd':
                dihedral_type = ('X', type3, type2, 'X')
                atomids = [id4, id3, id2, id1]
            else:
                log.error('ERROR Finding Dihedral Type for dihedrals ordering')   
        # if not sort_remap_atomids use orginal ordering
        else:
            dihedral_type = (type1, type2, type3, type4)
            atomids = [id1, id2, id3, id4]
        
        # Save dihedral info into class
        d = Dihedral()
        d.type = number
        d.symbol = dihedral_type
        d.atomids = atomids
        dihedrals[n+1] = d
    return dihedrals


#######################################
# Function for finding improper types #
#######################################
def find_improper_parameters(frc, BADI, use_auto_equivalence, use_assumed_auto_fill, aafc, ff_class, skip_printouts, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html

    Impropers
    The style of improper potential is specified in the input command file.
    
    (1) harmonic
    
      E = K (chi - chi0)^2
    
      chi = radians (computed by LAMMPS)
    
      coeff1 = K (energy/radian^2) (the usual 1/2 is included in the K)
      coeff2 = chi0 (degrees) (converted to radians within LAMMPS)
    
      in data file, listing of 4 atoms requires atom-1 as central atom
      some force fields (AMBER,Discover) have atom-2 as central atom - it is really
        an out-of-plane torsion, may need to treat as dihedral in LAMMPS
    
      2 coeffs are listed in data file
    (2) class2
    
      same formula, coeffs, and meaning as "harmonic" except that LAMMPS
        averages all 3 angle-contributions to chi
      in class II this is called a Wilson out-of-plane interaction
    
      2 coeffs are listed in data file
      
      Josh's comments:
          - atom type 2 is central atom and outer atoms are 1,3,4 [IE I='h', J='c', L='n',
            and K='o';   outer_atoms = sorted(I,L,K) = ('h', 'n', 'o') and the improper then
            is ('h', 'c', 'n', 'o') or(outer_atoms[0], J, outer_atoms[1], outer_atoms[2]) ]
          - this means we need to check the following types:
              - 1-2 as constant and look for 3-4 and 4-3
              - 3-2 as constant and look for 1-4 and 4-1
              - 4-2 as constant and look for 1-3 and 3-1
    """
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # improper coeffs dictionary
    improper_coeffs = {} # {improper type : list of coeffs}   
    improper_map = {} # { Ordered tuple(type1, type2, type3, type4) : numeric id }
    
    # Set comments, equivs/auto coeff dicts and equivs/auto mapping dicts
    if ff_class == 0:
        improper_comment = 'cvff' # class1 Header comment
        equivalents = frc.equivalences # class1 equivalents dict
        auto_equivs = frc.auto_equivalences # class1 auto-equivalents dict
        equivalent_coeffs = frc.out_of_plane # class1 equivalent coeffs (both equiv and auto-equiv types will be in this dict)
        auto_equiv_coeffs = frc.out_of_plane_auto # class1 equivalent coeffs (both equiv and auto-equiv types will be in this dict)
    elif ff_class == 1:
        improper_comment = 'cvff' # class1 Header comment
        equivalents = frc.equivalences # class1 equivalents dict
        auto_equivs = frc.auto_equivalences # class1 auto-equivalents dict
        equivalent_coeffs = frc.out_of_plane # class1 equivalent coeffs (both equiv and auto-equiv types will be in this dict)
        auto_equiv_coeffs = frc.out_of_plane_auto # class1 equivalent coeffs (both equiv and auto-equiv types will be in this dict)
    elif ff_class == 2:
        improper_comment = 'class2'  # class2 Header comment
        equivalents = frc.equivalences # class2 equivalents dict
        auto_equivs = frc.auto_equivalences # class2 auto-equivalents dict
        equivalent_coeffs = frc.wilson_out_of_plane # class2 equivalent coeffs
        auto_equiv_coeffs = frc.wilson_out_of_plane_auto # class2 auto-equivalent coeffs
    if ff_class == 'd':
        improper_comment = 'umbrella' # classd Header comment
        equivalents = frc.equivalences # classd equivalents dict
        auto_equivs = {} # classd auto-equivalents dict
        equivalent_coeffs = frc.out_of_plane_DREIDING # classd equivalent coeffs (both equiv and auto-equiv types will be in this dict)
        auto_equiv_coeffs = {} # classd equivalent coeffs (both equiv and auto-equiv types will be in this dict)

    if ff_class in ['s1', 's2']:
        improper_comment = 'Skeleton: {}'.format(ff_class) # skeleton Header comment
    
    ########################
    # Find improper coeffs #
    ########################
    # flagged angleangles to skip finding improper data for (will be implemented at the end of loop)
    flagged_angleangles = BADI.flagged_angleangles    
    for n, i in enumerate(BADI.improper_types_lst):
        
        # Flag things
        equivalent_flag = False
        auto_equivalent_flag = False
        add_onto_assumed_file_flag = False
        assumed_equivalent_flag = False
        assumed_auto_equivalent_flag = False
        skip_flag = False

        # bonded atom id types
        type1 = i[0]; type2 = i[1]; type3 = i[2]; type4 = i[3];
        
        # Set order/equiv/match and update later if found
        order = (type1, type2, type3, type4)
        equiv = (type1, type2, type3, type4)
        match = (type1, type2, type3, type4)

        # Find improper number based on n, either from improper_types_dict or angleangle_types_dict
        if n < len(BADI.improper_types_dict):
            number = BADI.improper_types_dict[i]
        elif n >= len(BADI.improper_types_dict):
            number = BADI.angleangle_types_dict[i]
        else:
            log.error('ERROR finding which dictionary to find number ID for oop or angleangle')
            
            
		# if ff_class 's1' or 's2' log generic and continue to next iterations
        if ff_class in ['s1', 's2']:
            # find nb comment and tmp_comment
            if number in flagged_angleangles:
                nb_comment = 'nb!=3'
                tmp_comment = 'skipped over b/c number of bonded atoms to central atom != 3 (angleangle set)'
            else:
                nb_comment = 'nb==3'
                tmp_comment = 'Skeleton option used'
            
            # set order as type1, type2, type3, type4
            order = (type1, type2, type3, type4)
            t = Type()
            t.type = order
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = []
            t.comments = tmp_comment
            t.nb = nb_comment
            improper_coeffs[number] = t
            improper_map[number] = order
            continue
            
        
        # Try matching equivalent impropers in all 6 permuations. Attempts: 1) without equivalences, 2) with equivalences
        if ff_functions.match_4_body_oop(type1, type2, type3, type4, log, equivalent_coeffs, equivalences=equivalents, wildcard_search=True, form='equiv')[0]:
            boolean, match, order, equiv = ff_functions.match_4_body_oop(type1, type2, type3, type4, log, equivalent_coeffs, equivalences=equivalents, wildcard_search=True, form='equiv')
            improper_coeff = equivalent_coeffs[match]
            equivalent_flag = True
            
        # Try matching auto equivalent impropers in all 6 permuationse. Attempts: 1) without equivalences, 2) with equivalences, 3) with equivalences and wild cards
        elif ff_functions.match_4_body_oop(type1, type2, type3, type4, log, auto_equiv_coeffs, equivalences=auto_equivs, wildcard_search=True, form='auto-equiv')[0] and use_auto_equivalence:
            boolean, match, order, equiv = ff_functions.match_4_body_oop(type1, type2, type3, type4, log, auto_equiv_coeffs, equivalences=auto_equivs, wildcard_search=True, form='auto-equiv')
            improper_coeff = auto_equiv_coeffs[match]
            auto_equivalent_flag = True
            
        # Try matching with assumed types if user wants and equivalent and auto equivalent fails
        elif not equivalent_flag and not auto_equivalent_flag and use_assumed_auto_fill:
            # Find basic elements in coeff and see if they exists in the assumed auto fill types
            if type1 in frc.atom_types and type2 in frc.atom_types and type3 in frc.atom_types and type4 in frc.atom_types:
                element1 = frc.atom_types[type1].element
                element2 = frc.atom_types[type2].element
                element3 = frc.atom_types[type3].element
                element4 = frc.atom_types[type4].element
                
                # Try as 1-2 as constant and varying 3-4 ordering
                if (element1, element2, element3, element4) in aafc.improper_coeffs:
                    assumed1, assumed2, assumed3, assumed4 = aafc.improper_coeffs[(element1, element2, element3, element4)]
                elif (element1, element2, element4, element3) in aafc.improper_coeffs:
                        assumed1, assumed2, assumed4, assumed3 = aafc.improper_coeffs[(element1, element2, element4, element3)]
                    
                # Try as 3-2 as constant and varying 1-4 ordering
                elif (element3, element2, element1, element4) in aafc.improper_coeffs:
                    assumed3, assumed2, assumed1, assumed4 = aafc.improper_coeffs[(element3, element2, element1, element4)]
                elif (element3, element2, element4, element1) in aafc.improper_coeffs:
                        assumed3, assumed2, assumed4, assumed1 = aafc.improper_coeffs[(element3, element2, element4, element1)]
                
                # Try as 4-2 as constant and varying 1-3 ordering
                elif (element4, element2, element1, element3) in aafc.improper_coeffs:
                    assumed4, assumed2, assumed1, assumed3 = aafc.improper_coeffs[(element4, element2, element1, element3)]
                elif (element4, element2, element3, element1) in aafc.improper_coeffs:
                        assumed4, assumed2, assumed3, assumed1 = aafc.improper_coeffs[(element4, element2, element3, element1)]

                # Save failure as a flag for printing
                else:
                    add_onto_assumed_file_flag = True
                    
                # Try using assumed1, assumed2, and assumed3 if they were in assumed types 
                if not add_onto_assumed_file_flag:
                    
                    # Function to remap assumed types back onto real types
                    def assumed_impropertypes_order_remap2_realtypes_order(assumed_order, type1, type2, type3, type4, log):
                        # Try 1-2 as constant and varying 3-4 ordering
                        if assumed_order == (assumed1, assumed2, assumed3, assumed4):
                            order = (type1, type2, type3, type4)
                        elif assumed_order == (assumed1, assumed2, assumed4, assumed3):
                            order = (type1, type2, type4, type3)
                            
                        # Try as 3-2 as constant and varying 1-4 ordering
                        elif assumed_order == (assumed3, assumed2, assumed1, assumed4):
                            order = (type3, type2, type1, type4)
                        elif assumed_order == (assumed3, assumed2, assumed4, assumed1):
                            order = (type3, type2, type4, type1)
                            
                        # Try as 4-2 as constant and varying 1-3 ordering
                        elif assumed_order == (assumed4, assumed2, assumed1, assumed3):
                            order = (type4, type2, type1, type3)
                        elif assumed_order == (assumed4, assumed2, assumed3, assumed1):
                            order = (type4, type2, type3, type1)
                        else:
                            order = (type1, type2, type3, type4)
                            if not skip_printouts: log.warn('WARNING: assumed auto fill option used for Improper Coeff {} {} {} {} {}, but ordering of real types COULD NOT BE REMAPPED from assumed types. DEFAULT ORDERING USED.'.format(number, type1, type2, type3, type4))
                        return order
                    
                    # Try matching assumed type with equivalent impropers (set equivalences to False since they are already in equivalent form)
                    if ff_functions.match_4_body_oop(assumed1, assumed2, assumed3, assumed4, log, equivalent_coeffs, equivalences=False, wildcard_search=False, form='equiv')[0]:
                        boolean, match, aa_order, equiv = ff_functions.match_4_body_oop(assumed1, assumed2, assumed3, assumed4, log, equivalent_coeffs, equivalences=False, wildcard_search=False, form='equiv')
                        improper_coeff = equivalent_coeffs[match]
                        assumed_equivalent_flag = True
                        
                        # Find correct order of real types based on assumed type ordering (aa_order)
                        order = assumed_impropertypes_order_remap2_realtypes_order(aa_order, type1, type2, type3, type4, log)
                        
                    # Try matching assumed type with auto equivalent impropers (set equivalences to False since they are already in equivalent form)
                    elif ff_functions.match_4_body_oop(assumed1, assumed2, assumed3, assumed4, log, auto_equiv_coeffs, equivalences=False, wildcard_search=True, form='auto-equiv')[0] and use_auto_equivalence:
                        boolean, match, aa_order, equiv = ff_functions.match_4_body_oop(assumed1, assumed2, assumed3, assumed4, log, auto_equiv_coeffs, equivalences=False, wildcard_search=True, form='auto-equiv')
                        improper_coeff = auto_equiv_coeffs[match]
                        assumed_auto_equivalent_flag = True
                        
                        # Find correct order of real types based on assumed type ordering (aa_order)
                        order = assumed_impropertypes_order_remap2_realtypes_order(aa_order, type1, type2, type3, type4, log)

        ##################################################
        # only find impropers that are impropers and     #
        # not flagged as an angleangle set of atom types # 
        # SET ALL OTHER FLAGS AS FALSE                   #
        ##################################################
        if number in flagged_angleangles:
            skip_flag = True
            
        # Save improper_coeff info into class for if it was to be skipped
        if skip_flag:
            
            # Set coeff based on FF class
            if ff_class == 0:
                coeff = [0.0, 1, 0] # int vs float will be used when writing final file
            if ff_class == 1:
                coeff = [0.0, 1, 0] # int vs float will be used when writing final file
            elif ff_class == 2:
                coeff = 2*[0.0]
            elif ff_class == 'd':
                coeff = 2*[0.0]
            
            order = (type1, type2, type3, type4)
            t = Type()
            t.type = order
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = coeff
            t.comments = 'skipped over b/c number of bonded atoms to central atom != 3 (angleangle set)'
            t.nb = 'nb!=3'
            improper_coeffs[number] = t
            improper_map[number] = order
            
        # Save improper_coeff info into class for equivalent form
        elif equivalent_flag and not auto_equivalent_flag:
            
            # Set coeff based on FF class
            if ff_class == 0:
                chi0 = improper_coeff.chi0;
                if chi0 == 180: chi_int = -1;
                else: chi_int = 1
                coeff = [improper_coeff.kchi, chi_int, improper_coeff.n]
            elif ff_class == 1:
                chi0 = improper_coeff.chi0;
                if chi0 == 180: chi_int = -1;
                else: chi_int = 1
                coeff = [improper_coeff.kchi, chi_int, improper_coeff.n]
            elif ff_class == 2:
                coeff = [improper_coeff.kchi, improper_coeff.chi0]
            elif ff_class == 'd':
                coeff = [improper_coeff.kl, improper_coeff.phi0]
            
            t = Type()
            t.type = order
            t.match = match
            t.equivalent = equiv
            t.coeffs = coeff
            t.comments = 'equivalent type used'
            t.nb = 'nb==3'
            improper_coeffs[number] = t
            improper_map[number] = order
            
        # Save improper_coeff info into class for auto equivalent form
        elif auto_equivalent_flag and not equivalent_flag:
                        
            # Set coeff based on FF class
            if ff_class == 0:
                chi0 = improper_coeff.chi0;
                if chi0 == 180: chi_int = -1;
                else: chi_int = 1
                coeff = [improper_coeff.kchi, chi_int, improper_coeff.n]
            elif ff_class == 1:
                chi0 = improper_coeff.chi0;
                if chi0 == 180: chi_int = -1;
                else: chi_int = 1
                coeff = [improper_coeff.kchi, chi_int, improper_coeff.n]
            elif ff_class == 2:
                coeff = [improper_coeff.kchi, improper_coeff.chi0]
            elif ff_class == 'd':
                coeff = [improper_coeff.kl, improper_coeff.phi0]
                
            t = Type()
            t.type = order
            t.match = match
            t.equivalent = equiv
            t.coeffs = coeff
            t.comments = 'auto equivalent type used'
            t.nb = 'nb==3'
            improper_coeffs[number] = t
            improper_map[number] = order
            
        # Save improper_coeff info into class for assumed auto fill option
        elif use_assumed_auto_fill and assumed_equivalent_flag or assumed_auto_equivalent_flag:
            if not skip_printouts: log.warn('WARNING unable to find improper info for Improper Coeff {} {} {} {} {}, but assumed auto fill option filled in with assumed type'.format(number, type1, type2, type3, type4))
        
            # If auto equivalent coeffs where used (set order as type1, type2, type3, type4)
            if assumed_equivalent_flag:
                
                # Set coeff based on FF class
                if ff_class == 0:
                    chi0 = improper_coeff.chi0;
                    if chi0 == 180: chi_int = -1;
                    else: chi_int = 1
                    coeff = [improper_coeff.kchi, chi_int, improper_coeff.n]
                elif ff_class == 1:
                    chi0 = improper_coeff.chi0;
                    if chi0 == 180: chi_int = -1;
                    else: chi_int = 1
                    coeff = [improper_coeff.kchi, chi_int, improper_coeff.n]
                elif ff_class == 2:
                    coeff = [improper_coeff.kchi, improper_coeff.chi0]
                elif ff_class == 'd':
                    coeff = [improper_coeff.kl, improper_coeff.phi0]
                
                t = Type()
                t.type = order
                t.match = match
                t.equivalent = equiv
                t.coeffs = coeff
                t.comments = '{:^10} {:2} {:2} {:2} {:2}'.format('ASSUMED equivalent type used', assumed1, assumed2, assumed3, assumed4) 
                t.nb = 'nb==3'
                improper_coeffs[number] = t
                improper_map[number] = order
                
            # elif auto equivalent coeffs where used (set order as type1, type2, type3, type4)
            elif assumed_auto_equivalent_flag:
                
                # Set coeff based on FF class
                if ff_class == 0:
                    chi0 = improper_coeff.chi0;
                    if chi0 == 180: chi_int = -1;
                    else: chi_int = 1
                    coeff = [improper_coeff.kchi, chi_int, improper_coeff.n]
                elif ff_class == 1:
                    chi0 = improper_coeff.chi0;
                    if chi0 == 180: chi_int = -1;
                    else: chi_int = 1
                    coeff = [improper_coeff.kchi, chi_int, improper_coeff.n]
                elif ff_class == 2:
                    coeff = [improper_coeff.kchi, improper_coeff.chi0]
                elif ff_class == 'd':
                    coeff = [improper_coeff.kl, improper_coeff.phi0]
                    
                t = Type()
                t.type = order
                t.match = match
                t.equivalent = equiv
                t.coeffs = coeff
                t.comments = '{:^10} {:2} {:2} {:2} {:2}'.format('ASSUMED auto equivalent type used', assumed1, assumed2, assumed3, assumed4)
                t.nb = 'nb==3'
                improper_coeffs[number] = t 
                improper_map[number] = order
                
        # Save improper_coeff info into class as zeros if all attempts fail
        else:
            # If this block runs due to assumed class being used, but missing parameters in file let user know
            if add_onto_assumed_file_flag:
                if not skip_printouts: log.warn('WARNING unable to find improper info for Improper Coeff {} {} {} {} {}, assumed_auto_fill was used but currently {} {} {} {} coeff is not in the file, please add to file'.format(number, type1, type2, type3, type4, element1, element2, element3, element4))
            else:
                if not skip_printouts: log.warn('WARNING unable to find improper info for Improper Coeff {} {} {} {} {}'.format(number, type1, type2, type3, type4))

            # Set coeff based on FF class
            if ff_class == 0:
                coeff = [0.0, 1, 0] # int vs float will be used when writing final file
            if ff_class == 1:
                coeff = [0.0, 1, 0] # int vs float will be used when writing final file
            elif ff_class == 2:
                coeff = 2*[0.0]
            elif ff_class == 'd':
                coeff = 2*[0.0]
            
            # set order as type1, type2, type3, type4
            order = (type1, type2, type3, type4)
            t = Type()
            t.type = order
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = coeff
            t.comments = 'UNABLE to find coeff parameters'
            t.nb = 'nb==3'
            improper_coeffs[number] = t
            improper_map[number] = order 
    return improper_coeffs, improper_map, improper_comment


##################################
# Function for finding impropers #
##################################
def find_impropers(nta, BADI, improper_map, m, sort_remap_atomids, ff_class, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html
    """
    
    # Dictionaries to add info to
    impropers = {}   # {improper number : improper object}
    
    # Find angles info
    for n, (id1, id2, id3, id4) in enumerate(BADI.impropers):
        type1 = nta[id1]; type2 = nta[id2]; type3 = nta[id3]; type4 = nta[id4];
        
        # Find number of bonded atoms (if this != 3 it is not an improper set, but an angleangle set - only log for class 2)
        nb = len(BADI.graph[id2])
        
        # Set dict2search based on nb number of bonded atoms "nb" (will allow for duplicate types between OOP and angleangle sets)
        if nb == 3:
            dict2search = BADI.improper_types_dict
            nb_type = 'nb==3'
        elif nb > 3:
            dict2search = BADI.angleangle_types_dict
            nb_type = 'nb!=3'
        else:
            log.error('ERROR finding which dict2search for oop or angleangle number ID')

        # Try as 1-2 as constant and varying 3-4 ordering
        if (type1, type2, type3, type4) in dict2search:
            number = dict2search[(type1, type2, type3, type4)]
        elif (type1, type2, type4, type3) in dict2search:
            number = dict2search[(type1, type2, type4, type3)]
        # Try as 3-2 as constant and varying 1-4 ordering
        elif (type3, type2, type1, type4) in dict2search:
            number = dict2search[(type3, type2, type1, type4)]
        elif (type3, type2, type4, type1) in dict2search:
            number = dict2search[(type3, type2, type4, type1)]
        # Try as 4-2 as constant and varying 1-3 ordering
        elif (type4, type2, type1, type3) in dict2search:
            number = dict2search[(type4, type2, type1, type3)]
        elif (type4, type2, type3, type1) in dict2search:
            number = dict2search[(type4, type2, type3, type1)]  
        elif ('X', type2, 'X', 'X') in dict2search and ff_class == 'd':
            number = dict2search[('X', type2, 'X', 'X')]
        else:
            log.error('ERROR Finding Improper Type for improper ordering')
        
        
        # Find improper ordering and set improper_type and atomids
        if sort_remap_atomids:
            ordering = improper_map[number]
            # Try as 1-2 as constant and varying 3-4 ordering
            if ordering == (type1, type2, type3, type4):
                improper_type = (type1, type2, type3, type4)
                atomids = [id1, id2, id3, id4]
            elif ordering == (type1, type2, type4, type3):
                improper_type = (type1, type2, type4, type3)
                atomids = [id1, id2, id4, id3]
            # Try as 3-2 as constant and varying 1-4 ordering
            elif ordering == (type3, type2, type1, type4):
                improper_type = (type3, type2, type1, type4)
                atomids = [id3, id2, id1, id4]
            elif ordering == (type3, type2, type4, type1):
                improper_type = (type3, type2, type4, type1)
                atomids = [id3, id2, id4, id1]
            # Try as 4-2 as constant and varying 1-3 ordering
            elif ordering == (type4, type2, type1, type3):
                improper_type = (type4, type2, type1, type3)
                atomids = [id4, id2, id1, id3]
            elif ordering == (type4, type2, type3, type1):
                improper_type = (type4, type2, type3, type1)
                atomids = [id4, id2, id3, id1]
            elif ordering == ('X', type2, 'X', 'X') and ff_class == 'd':
                improper_type = ('X', type2, 'X', 'X')
                atomids = [id1, id2, id3, id4]
            else:
                log.error('ERROR Finding Improper Type for improper ordering')  
        # if not sort_remap_atomids use orginal ordering
        else:
            improper_type = (type1, type2, type3, type4)
            atomids = [id1, id2, id3, id4]
        
        # Save improper info into class
        i = Improper()
        i.type = number
        i.nb = nb_type
        i.symbol = improper_type
        i.atomids = atomids
        impropers[n+1] = i
    return impropers


#######################################
# Function for finding bondbond types #
####################################### 
def find_bondbond_parameters(frc, BADI, angle_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log):    
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html
    Bond-Bond (computed within class II angles)

      E = K (r - r0) * (r' - r0')

      r,r' = distance (computed by LAMMPS)

      coeff1 = K (energy/distance^2) from bond-bond section in .frc file
      coeff2 = r0 (distance)         from quadratic_bonds or quartic_bonds section in .frc for bond i-j (1-2) in angle i-j-k (1-2-3)
      coeff3 = r0' (distance)        from quadratic_bonds or quartic_bonds section in .frc for bond j-k (2-3) in angle i-j-k (1-2-3)

      3 coeffs are input in data file
      
      Ordering for datafile to write: [coeff1,  coeff2,  coeff3]
      
      Each of the cross terms are searched separately even though they share a given angle/bond type. This allows 
      parameters to be in different order in the forcefield for each cross term or maybe not even there. 
    """
    
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # bondbond coeffs dictionary
    bondbond_coeffs = {} # {bondbond type : list of coeffs}
    
    ########################
    # Find bondbond coeffs #
    ########################
    for i in BADI.angle_types_lst:
        
        # Flag things
        bondbond_flag = False
        
        # bonded atom id types
        type1 = i[0]; type2 = i[1]; type3 = i[2];
        
        # bondbond number
        number = BADI.angle_types_dict[i]
        
        # find angles_order from angles_map if atom ids were resorted
        if sort_remap_atomids:
            type1, type2, type3 = angle_map[number]
                
		# if ff_class 's2' log generic and continue to next iterations
        if ff_class == 's2':
            # Save bondbond_coeff info into class
            t = Type()
            t.type = (type1, type2, type3)
            t.match = ['N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A']
            t.coeffs = []
            t.comments = 'Skeleton option used'
            bondbond_coeffs[number] = t
            continue    
        
        # Try matching bondbond in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
        if ff_functions.match_3_body(type1, type2, type3, log, frc.bondbond, equivalences=frc.equivalences, wildcard_search=False, form='equiv')[0]:
            boolean, match, order, equiv = ff_functions.match_3_body(type1, type2, type3, log, frc.bondbond, equivalences=frc.equivalences, wildcard_search=False, form='equiv')
            bondbond_coeff = frc.bondbond[match]
            bondbond_flag = True
            
        # Try finding r0_12 and r0_23
        r0_12, match_12, equiv_12 = ff_functions.get_crossterm_r0(type1, type2, log, use_auto_equivalence, frc)
        r0_23, match_23, equiv_23 = ff_functions.get_crossterm_r0(type2, type3, log, use_auto_equivalence, frc)
                            
        # Set comments for r0's 12 and 23
        c_12 = 'bond 12:   {:^4} {:^4} equivs:  {:^5} {:^5} match:  {:^5} {:^5}'.format(type1, type2, equiv_12[0], equiv_12[1], match_12[0], match_12[1])
        c_23 = 'bond 23:   {:^4} {:^4} equivs:  {:^5} {:^5} match:  {:^5} {:^5}'.format(type2, type3, equiv_23[0], equiv_23[1], match_23[0], match_23[1])
            
        # Set bondbond coeff if found else tell user that if failed
        if bondbond_flag:
            kb_bp = bondbond_coeff.kb_bp
            comment = 'Bondbond type used'; match = match; equiv = equiv;
            map_order = (type1, type2, type3) # Map order to keep standard coeff and cross-term coeff comments identical
        else:
            kb_bp = 0.0; map_order = (type1, type2, type3) # Map order to keep standard coeff and cross-term coeff comments identical
            comment = 'UNABLE to find coeff parameters'; match = ['N/A', 'N/A', 'N/A']; equiv = ['N/A', 'N/A', 'N/A'];
            if not skip_printouts: log.warn('WARNING unable to find bondbond info for Bondbond Coeff {} {} {} {}'.format(number, type1, type2, type3))
            
        # Save bondbond_coeff info into class
        t = Type()
        t.type = map_order
        t.match = match
        t.equivalent = equiv
        t.coeffs = [kb_bp, r0_12, r0_23]
        t.comments = '{:<40} {:<75} {:<75}'.format(comment, c_12, c_23)
        bondbond_coeffs[number] = t  
    return bondbond_coeffs        


########################################
# Function for finding bondangle types #
######################################## 
def find_bondangle_parameters(frc, BADI, angle_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log):    
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html

    Bond-Angle (computed within class II angles for each of 2 bonds)
    
      E = K_n (r - r0_n) * (theta - theta0)
    
      r = distance (computed by LAMMPS)
      theta = radians (computed by LAMMPS)
    
      coeff1 = K_1 (energy/distance-radians) from bond-angle section in .frc file
      coeff2 = K_2 (energy/distance-radians) from bond-angle section in .frc file
      coeff3 = r0_1 (distance)               from quadratic_bonds or quartic_bonds section in .frc for bond i-j (1-2) in angle i-j-k (1-2-3)
      coeff4 = r0_2 (distance)               from quadratic_bonds or quartic_bonds section in .frc for bond i-j (2-3) in angle i-j-k (1-2-3)
    
      Note: theta0 is known from angle coeffs so don't need it specified here
    
      4 coeffs are listed in data file
      
      Ordering for datafile to write: [coeff1,  coeff2,  coeff3,  coeff4]
      
      Each of the cross terms are searched separately even though they share a given angle/bond type. This allows 
      parameters to bevin different order in the forcefield for each cross term or maybe not even there. 
    """
    
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # bondangle coeffs dictionary
    bondangle_coeffs = {} # {bondangle type : list of coeffs}
    
    #########################
    # Find bondangle coeffs #
    #########################
    for i in BADI.angle_types_lst:
        
        # Flag things
        bondangle_flag = False
        
        # bonded atom id types
        type1 = i[0]; type2 = i[1]; type3 = i[2];
        
        # bondangle number
        number = BADI.angle_types_dict[i]
        
        # find angles_order from angles_map if atom ids were resorted
        if sort_remap_atomids:
            type1, type2, type3 = angle_map[number]
                
		# if ff_class 's2' log generic and continue to next iterations
        if ff_class == 's2':
            # Save bondangle_coeff info into class
            t = Type()
            t.type = (type1, type2, type3)
            t.match = ['N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A']
            t.coeffs = []
            t.comments = 'Skeleton option used'
            bondangle_coeffs[number] = t
            continue    
        
        # Try matching bondangle in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
        if ff_functions.match_3_body(type1, type2, type3, log, frc.bondangle, equivalences=frc.equivalences, wildcard_search=False, form='equiv')[0]:
            boolean, match, order, equiv = ff_functions.match_3_body(type1, type2, type3, log, frc.bondangle, equivalences=frc.equivalences, wildcard_search=False, form='equiv')
            bondangle_coeff = frc.bondangle[match]
            bondangle_flag = True
            
        # Try finding r0_12 and r0_23
        r0_12, match_12, equiv_12 = ff_functions.get_crossterm_r0(type1, type2, log, use_auto_equivalence, frc)
        r0_23, match_23, equiv_23 = ff_functions.get_crossterm_r0(type2, type3, log, use_auto_equivalence, frc)
        
        # Set comments for r0's 12 and 23
        c_12 = 'bond 12:   {:^4} {:^4} equivs:  {:^5} {:^5} match:  {:^5} {:^5}'.format(type1, type2, equiv_12[0], equiv_12[1], match_12[0], match_12[1])
        c_23 = 'bond 23:   {:^4} {:^4} equivs:  {:^5} {:^5} match:  {:^5} {:^5}'.format(type2, type3, equiv_23[0], equiv_23[1], match_23[0], match_23[1])

        # Set bondangle coeff if found else tell user that if failed
        if bondangle_flag:
            
            # find if found in backwards order or forwards and set coeffs accordingly. else set as zeros and warn
            if order == (type1, type2, type3):
                k1 = bondangle_coeff.kb_theta
                k2 = bondangle_coeff.kbp_theta
                equiv = equiv
                comment = 'Bondangle type used'
                map_order = (type1, type2, type3) # Map order to keep standard coeff and cross-term coeff comments identical
            elif order == (type3, type2, type1):
                k2 = bondangle_coeff.kb_theta
                k1 = bondangle_coeff.kbp_theta
                equiv = tuple(reversed(equiv))
                comment = 'Bondangle type used'
                map_order = (type1, type2, type3) # Map order to keep standard coeff and cross-term coeff comments identical
            else:
                k1 = 0.0; k2 = 0.0; map_order = (type1, type2, type3) # Map order to keep standard coeff and cross-term coeff comments identical
                comment = 'UNABLE to find coeff parameters'; match = ['N/A', 'N/A', 'N/A']; equiv = ['N/A', 'N/A', 'N/A'];
                if not skip_printouts: log.warn(f'WARNING: Unable to determine bondangle {number} ordering for k1 and k2. Setting as zeros!')
                
        else:
            k1 = 0.0;  k2 = 0.0; map_order = (type1, type2, type3) # Map order to keep standard coeff and cross-term coeff comments identical
            comment = 'UNABLE to find coeff parameters'; match = ['N/A', 'N/A', 'N/A']; equiv = ['N/A', 'N/A', 'N/A'];
            if not skip_printouts: log.warn('WARNING unable to find bondangle info for Bondangle Coeff {} {} {} {}'.format(number, type1, type2, type3))

            
        # Save bondangle_coeff info into class
        t = Type()
        t.type = map_order
        t.match = match
        t.equivalent = equiv
        t.coeffs = [k1, k2, r0_12, r0_23]
        t.comments = '{:<40} {:<75} {:<75}'.format(comment, c_12, c_23)
        bondangle_coeffs[number] = t
    return bondangle_coeffs


################################################
# Function for finding AngleAngleTorsion types #
################################################
def find_angleangletorsion_parameters(frc, BADI, dihedral_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html

    Angle-Angle-Torsion (computed within class II dihedral)
    
    E = K (theta - theta0) * (theta' - theta0') * (phi - phi0)
  
    theta,theta' = radians (computed by LAMMPS)
    phi = radians (computed by LAMMPS)
  
    coeff1 = K (energy/radians^3)                                       from angle-angle-torsion section in .frc file for dihedral i-j-k-l (try to find with this ordering i-j-k-l in forward and reverse - 2 permutations)
    coeff2 = theta0 (degrees) (converted to radians within LAMMPS)      from quadratic_angles or quaratic_angles section in .frc for angle i-j-k (1-2-3) in Improper i-j-k-l (1-2-3-4)
    coeff3 = theta0' (degrees) (converted to radians within LAMMPS)     from quadratic_angles or quaratic_angles section in .frc for angle j-k-l (2-3-4) in Improper i-j-k-l (1-2-3-4)
  
    Note: phi0 is known from dihedral coeffs so don't need it specified here
  
    3 coeffs are listed in data file
    
    Ordering for datafile to write: [coeff1,  coeff2,  coeff3]
    
    Each of the cross terms are searched separately even though they share a given dihedral/angle type. This allows 
    parameters to be in different order in the forcefield for each cross term or maybe not even there.                                               

    """
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # angleangletorsion coeffs dictionary
    angleangletorsion_coeffs = {} # {angleangle_torsion type : list of coeffs}
    
    #################################
    # Find angleangletorsion coeffs #
    #################################
    for i in BADI.dihedral_types_lst:
        
        # Set flags
        angleangletorsion_flag = False
        
        # bonded atom id types
        type1 = i[0]; type2 = i[1]; type3 = i[2]; type4 = i[3];
        
        # angleangletorsion number
        number = BADI.dihedral_types_dict[i]
        
        # find dihedral order from dihedral_map if atom ids were resorted
        if sort_remap_atomids:
            type1, type2, type3, type4 = dihedral_map[number]
                
		# if ff_class 's2' log generic and continue to next iterations
        if ff_class == 's2':
            # save angleangletorsion_coeff into class
            t = Type()
            t.type = (type1, type2, type3, type4)
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = []
            t.comments = 'Skeleton option used'
            angleangletorsion_coeffs[number] = t
            continue   
        
        # Try matching frc.angleangletorsion in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
        if ff_functions.match_4_body(type1, type2, type3, type4, log, frc.angleangletorsion, equivalences=frc.equivalences, wildcard_search=False, form='equiv')[0]:
            boolean, match, order, equiv = ff_functions.match_4_body(type1, type2, type3, type4, log, frc.angleangletorsion, equivalences=frc.equivalences, wildcard_search=False, form='equiv')
            angleangletorsion_coeff = frc.angleangletorsion[match]
            angleangletorsion_flag = True
            
        # Set angleangletorsion_coeff if found else tell user that if failed
        if angleangletorsion_flag:
            k_ang_ang_tor = angleangletorsion_coeff.k_ang_ang_tor
            comment = 'Angleangletorsion type used'; match = match; equiv = equiv;
            map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
        else:
            k_ang_ang_tor = 0.0; map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
            comment = 'UNABLE to find coeff parameters'; match = ['N/A', 'N/A', 'N/A', 'N/A']; equiv = ['N/A', 'N/A', 'N/A', 'N/A'];
            if not skip_printouts: log.warn('WARNING unable to find angleangletorsion info for Angleangletorsion Coeff {} {} {} {} {}'.format(number, type1, type2, type3, type4))
            
        # Try finding theta0_123 and theta0_234
        theta0_123, match_123, equiv_123 = ff_functions.get_crossterm_theta0(type1, type2, type3, log, use_auto_equivalence, frc)
        theta0_234, match_234, equiv_234 = ff_functions.get_crossterm_theta0(type2, type3, type4, log, use_auto_equivalence, frc)

        # theta0_123 and theta0_234 comments
        c_123 = 'angle 123:   {:^6} {:^6} {:^6} equivs:  {:^6} {:^6} {:^6} match:  {:^6} {:^6} {:^6}'.format(type1, type2, type3, equiv_123[0], equiv_123[1], equiv_123[2], match_123[0], match_123[1], match_123[2])
        c_234 = 'angle 234:   {:^6} {:^6} {:^6} equivs:  {:^6} {:^6} {:^6} match:  {:^6} {:^6} {:^6}'.format(type2, type3, type4, equiv_234[0], equiv_234[1], equiv_234[2], match_234[0], match_234[1], match_234[2])
    
        # Build final angleangletorsion_coeff
        coeff_123 = [k_ang_ang_tor, theta0_123, theta0_234] 
        
        # save angleangletorsion_coeff into class
        t = Type()
        t.type = map_order
        t.match = match
        t.equivalent = equiv
        t.coeffs = coeff_123
        t.comments = '{:<50} {:<100} {:<100}'.format(comment, c_123, c_234)
        angleangletorsion_coeffs[number] = t
    return angleangletorsion_coeffs


#############################################
# Function for finding EndBondTorsion types #
#############################################
def find_endbondtorsion_parameters(frc, BADI, dihedral_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log): 
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html

    End-Bond-Torsion (computed within class II dihedral for each of 2 bonds)
    
      E = (r - r0_n) * [ F1_n*cos(phi) + F2_n*cos(2*phi) + F3_n*cos(3*phi) ]
    
      r = distance (computed by LAMMPS)
      phi = radians (computed by LAMMPS)
    
      coeff1 = F1_1 (energy/distance) from end-bond-torsion section in .frc file for dihedral i-j-k-l (try to find with this ordering i-j-k-l in forward and reverse - 2 permutations)
      coeff2 = F2_1 (energy/distance) from end-bond-torsion section in .frc file for dihedral i-j-k-l (try to find with this ordering i-j-k-l in forward and reverse - 2 permutations)
      coeff3 = F3_1 (energy/distance) from end-bond-torsion section in .frc file for dihedral i-j-k-l (try to find with this ordering i-j-k-l in forward and reverse - 2 permutations)
      coeff4 = F1_2 (energy/distance) from end-bond-torsion section in .frc file for dihedral i-j-k-l (try to find with this ordering i-j-k-l in forward and reverse - 2 permutations)
      coeff5 = F2_3 (energy/distance) from end-bond-torsion section in .frc file for dihedral i-j-k-l (try to find with this ordering i-j-k-l in forward and reverse - 2 permutations)
      coeff6 = F3_3 (energy/distance) from end-bond-torsion section in .frc file for dihedral i-j-k-l (try to find with this ordering i-j-k-l in forward and reverse - 2 permutations)
      coeff7 = r0_1 (distance)        from quadratic_bonds or quaratic_bonds section in .frc for bond i-j (1-2) in dihedral i-j-k-l (1-2-3-4) (try to find with this ordering i-j in forward and reverse - 2 permutations)
      coeff8 = r0_2 (distance)        from quadratic_bonds or quaratic_bonds section in .frc for bond k-l (3-4) in dihedral i-j-k-l (1-2-3-4) (try to find with this ordering k-l in forward and reverse - 2 permutations)
    
      8 coeffs are listed in data file
    
      Ordering for datafile to write: [coeff1,  coeff2,  coeff3,  coeff4,  coeff5,  coeff6,  coeff7,  coeff8]
    
    Each of the cross terms are searched separately even though they share a given dihedral/angle type. This allows 
    parameters to be in different order in the forcefield for each cross term or maybe not even there.                                               

    """
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # endbondtorsion coeffs dictionary
    endbondtorsion_coeffs = {} # {endbondtorsion type : list of coeffs}
    
    ##############################
    # Find endbondtorsion coeffs #
    ##############################
    for i in BADI.dihedral_types_lst:
        
        # Set flags
        endbondtorsion_flag = False
        
        # bonded atom id types
        type1 = i[0]; type2 = i[1]; type3 = i[2]; type4 = i[3];
        
        # endbondtorsion number
        number = BADI.dihedral_types_dict[i]
                   
        # find dihedral order from dihedral_map if atom ids were resorted
        if sort_remap_atomids:
            type1, type2, type3, type4 = dihedral_map[number]
                
		# if ff_class 's2' log generic and continue to next iterations
        if ff_class == 's2':
            # save endbondtorsion_coeff into class
            t = Type()
            t.type = (type1, type2, type3, type4)
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = []
            t.comments = 'Skeleton option used'
            endbondtorsion_coeffs[number] = t 
            continue   
        
        # Try matching frc.endbondtorsion in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
        if ff_functions.match_4_body(type1, type2, type3, type4, log, frc.endbondtorsion, equivalences=frc.equivalences, wildcard_search=False, form='equiv')[0]:
            boolean, match, order, equiv = ff_functions.match_4_body(type1, type2, type3, type4, log, frc.endbondtorsion, equivalences=frc.equivalences, wildcard_search=False, form='equiv')
            ebt = frc.endbondtorsion[match]
            endbondtorsion_flag = True
            
        # Set endbondtorsion_coeff if found else tell user that if failed
        if endbondtorsion_flag:
            # find if found in backwards order or forwards and set coeffs accordingly. else set as zeros and warn
            if order == (type1, type2, type3, type4):
                f1 = ebt.l_f1
                f2 = ebt.l_f2
                f3 = ebt.l_f3
                f4 = ebt.r_f1
                f5 = ebt.r_f2
                f6 = ebt.r_f3
                equiv = equiv
                comment = 'Endbondtorsion type used';
                map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
            elif order ==  (type4, type3, type2, type1):
                f4 = ebt.l_f1
                f5 = ebt.l_f2
                f6 = ebt.l_f3
                f1 = ebt.r_f1
                f2 = ebt.r_f2
                f3 = ebt.r_f3
                equiv = tuple(reversed(equiv))
                comment = 'Endbondtorsion type used';
                map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
            else:
                f1 = 0.0; f2 = 0.0; f3 = 0.0; f4 = 0.0; f5 = 0.0; f6 = 0.0; map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
                comment = 'UNABLE to find coeff parameters'; match = ['N/A', 'N/A', 'N/A', 'N/A']; equiv = ['N/A', 'N/A', 'N/A', 'N/A'];
                if not skip_printouts: log.warn(f'WARNING: Unable to determine endbondtorsion {number} ordering for f1, f2, f3, f4, f5, f6. Setting as zeros!')      
        else:
            f1 = 0.0; f2 = 0.0; f3 = 0.0; f4 = 0.0; f5 = 0.0; f6 = 0.0; map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
            comment = 'UNABLE to find coeff parameters'; match = ['N/A', 'N/A', 'N/A', 'N/A']; equiv = ['N/A', 'N/A', 'N/A', 'N/A'];
            if not skip_printouts: log.warn('WARNING unable to find endbondtorsion info for Endbondtorsion Coeff {} {} {} {} {}'.format(number, type1, type2, type3, type4))
            
        # Try finding r0_12 and r0_34
        r0_12, match_12, equiv_12 = ff_functions.get_crossterm_r0(type1, type2, log, use_auto_equivalence, frc)
        r0_34, match_34, equiv_34 = ff_functions.get_crossterm_r0(type3, type4, log, use_auto_equivalence, frc)
    
        # Set comments for r0's 12 and 34
        c_12 = 'bond 12:   {:^4} {:^4} equivs:  {:^5} {:^5} match:  {:^5} {:^5}'.format(type1, type2, equiv_12[0], equiv_12[1], match_12[0], match_12[1])
        c_34 = 'bond 34:   {:^4} {:^4} equivs:  {:^5} {:^5} match:  {:^5} {:^5}'.format(type3, type4, equiv_34[0], equiv_34[1], match_34[0], match_34[1])
        
        # Build final endbondtorsion_coeff
        coeff_12345678 = [f1, f2, f3, f4, f5, f6, r0_12, r0_34]
        
        # save endbondtorsion_coeff into class
        t = Type()
        t.type = map_order
        t.match = match
        t.equivalent = equiv
        t.coeffs = coeff_12345678
        t.comments = '{:<50} {:<75} {:<75}'.format(comment, c_12, c_34)
        endbondtorsion_coeffs[number] = t    
    return endbondtorsion_coeffs


################################################
# Function for finding MiddleBondTorsion types #
################################################
def find_middlebondtorsion_parameters(frc, BADI, dihedral_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html
    
    Middle-Bond-Torsion (computed within class II dihedral)

      E = (r - r0) * [ F1*cos(phi) + F2*cos(2*phi) + F3*cos(3*phi) ]
    
      r = distance (computed by LAMMPS)
      phi = radians (computed by LAMMPS)
    
      coeff1 = F1 (energy/distance) from middle-bond-torsion section in .frc file for dihedral i-j-k-l (try to find with this ordering i-j-k-l in forward and reverse - 2 permutations)
      coeff2 = F2 (energy/distance) from middle-bond-torsion section in .frc file for dihedral i-j-k-l (try to find with this ordering i-j-k-l in forward and reverse - 2 permutations)
      coeff3 = F3 (energy/distance) from middle-bond-torsion section in .frc file for dihedral i-j-k-l (try to find with this ordering i-j-k-l in forward and reverse - 2 permutations)
      coeff4 = r0 (distance)        from quadratic_bonds or quaratic_bonds section in .frc for bond j-k (1-2) in dihedral i-j-k-l (1-2-3-4) (try to find with this ordering j-k in forward and reverse - 2 permutations)
    
      4 coeffs are listed in data file
    
      Ordering for datafile to write: [coeff1,  coeff2,  coeff3,  coeff4]
    
    Each of the cross terms are searched separately even though they share a given dihedral/angle type. This allows 
    parameters to be in different order in the forcefield for each cross term or maybe not even there.                                               

    """
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # middlebondtorsion coeffs dictionary
    middlebondtorsion_coeffs = {} # {middlebondtorsion type : list of coeffs}
    
    #################################
    # Find middlebondtorsion coeffs #
    #################################
    for i in BADI.dihedral_types_lst:
        
        # Set flags
        middlebondtorsion_flag = False
        
        # bonded atom id types
        type1 = i[0]; type2 = i[1]; type3 = i[2]; type4 = i[3];
        
        # middlebondtorsion number
        number = BADI.dihedral_types_dict[i]
            
        # find dihedral order from dihedral_map if atom ids were resorted
        if sort_remap_atomids:
            type1, type2, type3, type4 = dihedral_map[number]
                
		# if ff_class 's2' log generic and continue to next iterations
        if ff_class == 's2':
            # save middlebondtorsion_coeff into class
            t = Type()
            t.type = (type1, type2, type3, type4)
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = []
            t.comments = 'Skeleton option used'
            middlebondtorsion_coeffs[number] = t 
            continue   
        
        # Try matching frc.middlebondtorsion in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
        if ff_functions.match_4_body(type1, type2, type3, type4, log, frc.middlebondtorsion, equivalences=frc.equivalences, wildcard_search=False, form='equiv')[0]:
            boolean, match, order, equiv = ff_functions.match_4_body(type1, type2, type3, type4, log, frc.middlebondtorsion, equivalences=frc.equivalences, wildcard_search=False, form='equiv')
            mbt = frc.middlebondtorsion[match]
            middlebondtorsion_flag = True
            
        # Set endbondtorsion_coeff if found else tell user that if failed
        if middlebondtorsion_flag:
            # find if found in backwards order or forwards and set coeffs accordingly. else set as zeros and warn
            if order == (type1, type2, type3, type4):
                f1 = mbt.f1
                f2 = mbt.f2
                f3 = mbt.f3
                equiv = equiv
                comment = 'Middlebondtorsion type used';
                map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
            elif order == (type4, type3, type2, type1):
                f1 = mbt.f1
                f2 = mbt.f2
                f3 = mbt.f3
                equiv = equiv
                comment = 'Middlebondtorsion type used';
                map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
            else:
                f1 = 0.0; f2 = 0.0; f3 = 0.0; map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
                comment = 'UNABLE to find coeff parameters'; match = ['N/A', 'N/A', 'N/A', 'N/A']; equiv = ['N/A', 'N/A', 'N/A', 'N/A'];
                if not skip_printouts: log.warn(f'WARNING: Unable to determine middlebondtorsion {number} ordering for f1, f2, f3. Setting as zeros!')
        else:
            f1 = 0.0; f2 = 0.0; f3 = 0.0; map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
            comment = 'UNABLE to find coeff parameters'; match = ['N/A', 'N/A', 'N/A', 'N/A']; equiv = ['N/A', 'N/A', 'N/A', 'N/A'];
            if not skip_printouts: log.warn('WARNING unable to find middlebondtorsion info for Middlebondtorsion Coeff {} {} {} {} {}'.format(number, type1, type2, type3, type4))
            
        # Try finding r0_23 and set comments
        r0_23, match_23, equiv_23 = ff_functions.get_crossterm_r0(type2, type3, log, use_auto_equivalence, frc)
        c_23 = 'bond 23:   {:^4} {:^4} equivs:  {:^5} {:^5} match:  {:^5} {:^5}'.format(type2, type3, equiv_23[0], equiv_23[1], match_23[0], match_23[1])
        
        # Build final middlebondtorsion_coeff
        coeff_1234 = [f1, f2, f3, r0_23]
        
        # save middlebondtorsion_coeff into class
        t = Type()
        t.type = map_order
        t.match = match
        t.equivalent = equiv
        t.coeffs = coeff_1234
        t.comments = '{:<50} {:<75}'.format(comment, c_23)
        middlebondtorsion_coeffs[number] = t        
    return middlebondtorsion_coeffs


#########################################
# Function for finding BondBond13 types #
#########################################
def find_bondbond13_parameters(frc, BADI, dihedral_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html
    
    Bond-Bond-13-Torsion (computed within class II dihedral)
    
      (undocumented)
      
     Josh's notes from PCFF/IFF .frc file
     #bond-bond_1_3        cff91
    
     E = K(b,b') * (R - R0) * (R' - R0')
      
    coeff1 = K (energy/distance^2) from bond-bond-1-3 section in .frc file
    coeff2 = r0 (distance)        from quadratic_bonds or quaratic_bonds section in .frc for bond i-j (1-2) in angle i-j-k-l (1-2-3-4) (try to find with this ordering i-j in forward and reverse - 2 permutations)
    coeff3 = r0 (distance)        from quadratic_bonds or quaratic_bonds section in .frc for bond k-l (3-4) in angle i-j-k-l (1-2-3-4) (try to find with this ordering k-l in forward and reverse - 2 permutations)
  
    3 coeffs are listed in data file
  
    Ordering for datafile to write: [coeff1,  coeff2,  coeff3]
    
        
    Each of the cross terms are searched separately even though they share a given dihedral/angle type. This allows 
    parameters to be in different order in the forcefield for each cross term or maybe not even there.                                               

    """
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # middlebondtorsion coeffs dictionary
    bondbond13_coeffs = {} # {bondbond13 type : list of coeffs}
    
    ##########################
    # Find bondbond13 coeffs #
    ##########################
    for i in BADI.dihedral_types_lst:
        
        # Set flags
        bondbond13_flag = False
        
        # bonded atom id types
        type1 = i[0]; type2 = i[1]; type3 = i[2]; type4 = i[3];
        
        # bondbond13 number
        number = BADI.dihedral_types_dict[i]
        
        # find dihedral order from dihedral_map if atom ids were resorted
        if sort_remap_atomids:
            type1, type2, type3, type4 = dihedral_map[number]
                
		# if ff_class 's2' log generic and continue to next iterations
        if ff_class == 's2':
            # save bondbond13_coeff into class
            t = Type()
            t.type = (type1, type2, type3, type4)
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = []
            t.comments = 'Skeleton option used'
            bondbond13_coeffs[number] = t 
            continue   
        
        # Try matching frc.bondbond13 in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
        if ff_functions.match_4_body(type1, type2, type3, type4, log, frc.bondbond13, equivalences=frc.equivalences, wildcard_search=False, form='equiv')[0]:
            boolean, match, order, equiv = ff_functions.match_4_body(type1, type2, type3, type4, log, frc.bondbond13, equivalences=frc.equivalences, wildcard_search=False, form='equiv')
            bondbond13_coeff = frc.bondbond13[match]
            bondbond13_flag = True
            
        # Set bondbond13_coeff if found else tell user that if failed
        if bondbond13_flag:
            bondbond13 = bondbond13_coeff.kb_bp
            comment = 'Bondbond13 type used'; match = match; equiv = equiv;
            map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
        else:
            bondbond13 = 0.0; map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
            comment = 'UNABLE to find coeff parameters'; match = ['N/A', 'N/A', 'N/A', 'N/A']; equiv = ['N/A', 'N/A', 'N/A', 'N/A'];
            if not skip_printouts: log.warn('WARNING unable to find bondbond13 info for Bondbond13 Coeff {} {} {} {} {}'.format(number, type1, type2, type3, type4))
            
        # Try finding r0_12 and r0_34
        r0_12, match_12, equiv_12 = ff_functions.get_crossterm_r0(type1, type2, log, use_auto_equivalence, frc)
        r0_34, match_34, equiv_34 = ff_functions.get_crossterm_r0(type3, type4, log, use_auto_equivalence, frc)
    
        # Set comments for r0's 12 and 34
        c_12 = 'bond 12:   {:^4} {:^4} equivs:  {:^5} {:^5} match:  {:^5} {:^5}'.format(type1, type2, equiv_12[0], equiv_12[1], match_12[0], match_12[1])
        c_34 = 'bond 34:   {:^4} {:^4} equivs:  {:^5} {:^5} match:  {:^5} {:^5}'.format(type3, type4, equiv_34[0], equiv_34[1], match_34[0], match_34[1])
        
        # Build final middlebondtorsion_coeff
        coeff_123 = [bondbond13, r0_12, r0_34]
        
        # save bondbond13_coeff into class
        t = Type()
        t.type = map_order
        t.match = match
        t.equivalent = equiv
        t.coeffs = coeff_123
        t.comments = '{:<50} {:<75} {:<75}'.format(comment, c_12, c_34)
        bondbond13_coeffs[number] = t  
    return bondbond13_coeffs


###########################################
# Function for finding AngleTorsion types #
###########################################
def find_angletorsion_parameters(frc, BADI, dihedral_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html
    
    Angle-Torsion (computed within class II dihedral for each of 2 angles)
    
      E = (theta - theta0) * [ F1_n*cos(phi) + F2_n*cos(2*phi) + F3_n*cos(3*phi) ]
    
      theta = radians (computed by LAMMPS)
      phi = radians (computed by LAMMPS)
    
      coeff1 = F1_1 (energy/radians)
      coeff2 = F2_1 (energy/radians)
      coeff3 = F3_1 (energy/radians)
      coeff4 = F1_2 (energy/radians)
      coeff5 = F2_3 (energy/radians)
      coeff6 = F3_3 (energy/radians)
      coeff7 = theta0_1 (degrees) (converted to radians within LAMMPS) from quadratic_angles or quaratic_angles section in .frc for angle i-j-k (1-2-3) in Improper i-j-k-l (1-2-3-4)
      coeff8 = theta0_2 (degrees) (converted to radians within LAMMPS) from quadratic_angles or quaratic_angles section in .frc for angle j-k-l (2-3-4) in Improper i-j-k-l (1-2-3-4)
    
      8 coeffs are listed in data file
  
    Ordering for datafile to write: [coeff1,  coeff2,  coeff3,  coeff4,  coeff5,  coeff6,  coeff7,  coeff8]
            
    Each of the cross terms are searched separately even though they share a given dihedral/angle type. This allows 
    parameters to be in different order in the forcefield for each cross term or maybe not even there.                                               

    """
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # angletorsion coeffs dictionary
    angletorsion_coeffs = {} # {angletorsion type : list of coeffs}
    
    ############################
    # Find angletorsion coeffs #
    ############################
    for i in BADI.dihedral_types_lst:
        
        # Set flags
        angletorsion_flag = False
        
        # bonded atom id types
        type1 = i[0]; type2 = i[1]; type3 = i[2]; type4 = i[3];
        
        # angletorsion number
        number = BADI.dihedral_types_dict[i]
        
        # find dihedral order from dihedral_map if atom ids were resorted
        if sort_remap_atomids:
            type1, type2, type3, type4 = dihedral_map[number]
                
		# if ff_class 's2' log generic and continue to next iterations
        if ff_class == 's2':
            # save angletorsion_coeffs_coeff into class
            t = Type()
            t.type = (type1, type2, type3, type4)
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = []
            t.comments = 'Skeleton option used'
            angletorsion_coeffs[number] = t
            continue   
        
        # Try matching frc.angletorsion in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
        if ff_functions.match_4_body(type1, type2, type3, type4, log, frc.angletorsion, equivalences=frc.equivalences, wildcard_search=False, form='equiv')[0]:
            boolean, match, order, equiv = ff_functions.match_4_body(type1, type2, type3, type4, log, frc.angletorsion, equivalences=frc.equivalences, wildcard_search=False, form='equiv')
            at = frc.angletorsion[match]
            angletorsion_flag = True
            
        # Set angletorsion_coeff if found else tell user that if failed
        if angletorsion_flag:
            # find if found in backwards order or forwards and set coeffs accordingly. else set as zeros and warn
            if order == (type1, type2, type3, type4):
                f1 = at.l_f1
                f2 = at.l_f2
                f3 = at.l_f3
                f4 = at.r_f1
                f5 = at.r_f2
                f6 = at.r_f3
                equiv = equiv
                comment = 'Angletorsion type used';
                map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
            elif order == (type4, type3, type2, type1):
                f4 = at.l_f1
                f5 = at.l_f2
                f6 = at.l_f3
                f1 = at.r_f1
                f2 = at.r_f2
                f3 = at.r_f3
                equiv = tuple(reversed(equiv))
                comment = 'Angletorsion type used';
                map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
            else:
                f1 = 0.0; f2 = 0.0; f3 = 0.0; f4 = 0.0; f5 = 0.0; f6 = 0.0; map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
                log.warn(f'WARNING: Unable to determine angletorsion {number} ordering for f1, f2, f3, f4, f5, f6. Setting as zeros!')
                comment = 'UNABLE to find coeff parameters'; match = ['N/A', 'N/A', 'N/A', 'N/A']; equiv = ['N/A', 'N/A', 'N/A', 'N/A'];

        else:
            f1 = 0.0; f2 = 0.0; f3 = 0.0; f4 = 0.0; f5 = 0.0; f6 = 0.0; map_order = (type1, type2, type3, type4) # Map order to keep standard coeff and cross-term coeff comments identical
            comment = 'UNABLE to find coeff parameters'; match = ['N/A', 'N/A', 'N/A', 'N/A']; equiv = ['N/A', 'N/A', 'N/A', 'N/A'];
            if not skip_printouts: log.warn('WARNING unable to find angletorsion info for Angletorsion Coeff {} {} {} {} {}'.format(number, type1, type2, type3, type4))
            
        # Try finding theta0_123 and theta0_234
        theta0_123, match_123, equiv_123 = ff_functions.get_crossterm_theta0(type1, type2, type3, log, use_auto_equivalence, frc)
        theta0_234, match_234, equiv_234 = ff_functions.get_crossterm_theta0(type2, type3, type4, log, use_auto_equivalence, frc)

        # theta0_123 and theta0_234 comments
        c_123 = 'angle 123:   {:^6} {:^6} {:^6} equivs:  {:^6} {:^6} {:^6} match:  {:^6} {:^6} {:^6}'.format(type1, type2, type3, equiv_123[0], equiv_123[1], equiv_123[2], match_123[0], match_123[1], match_123[2])
        c_234 = 'angle 234:   {:^6} {:^6} {:^6} equivs:  {:^6} {:^6} {:^6} match:  {:^6} {:^6} {:^6}'.format(type2, type3, type4, equiv_234[0], equiv_234[1], equiv_234[2], match_234[0], match_234[1], match_234[2])
    
        # Build final angletorsion_coeff
        coeff_12345678 = [f1, f2, f3, f4, f5, f6, theta0_123, theta0_234]
        
        # save angletorsion_coeffs_coeff into class
        t = Type()
        t.type = map_order
        t.match = match
        t.equivalent = equiv
        t.coeffs = coeff_12345678
        t.comments = '{:<50} {:<100} {:<100}'.format(comment, c_123, c_234)
        angletorsion_coeffs[number] = t
    return angletorsion_coeffs


#########################################
# Function for finding angleangle types #
#########################################
def find_angleangle_parameters(frc, BADI, improper_map, use_auto_equivalence, sort_remap_atomids, ff_class, skip_printouts, log):
    """
    https://docs.lammps.org/99/force_fields.html   or    https://docs.lammps.org/Manual.html

    Angle-Angle (computed within class II improper for each of 3 pairs of angles)
    
      E = K_n (theta - theta0_n) * (theta' - theta0_n')
    
      theta,theta' = radians (computed by LAMMPS)
    
      coeff1 = K_1 (energy/radians^2)                                   from angle-angle section in .frc file for Improper
      coeff2 = K_2 (energy/radians^2)                                   from angle-angle section in .frc file for Improper
      coeff3 = K_3 (energy/radians^2)                                   from angle-angle section in .frc file for Improper
      coeff4 = theta0_1 (degrees) (converted to radians within LAMMPS)  from quadratic_angles or quartic_angles section in .frc for angle i-j-k (1-2-3) in Improper i-j-k-l (1-2-3-4)
      coeff5 = theta0_2 (degrees) (converted to radians within LAMMPS)  from quadratic_angles or quartic_angles section in .frc for angle i-j-l (1-2-4) in Improper i-j-k-l (1-2-3-4)
      coeff6 = theta0_3 (degrees) (converted to radians within LAMMPS)  from quadratic_angles or quartic_angles section in .frc for angle k-j-l (3-2-4) in Improper i-j-k-l (1-2-3-4)
    
      6 coeffs are listed in data file
      
      Ordering for datafile to write: [coeff1,  coeff2,  coeff3,  coeff4,  coeff5,  coeff6]
      
      This is perhaps one of the more complex coeffs to find all data for since each K_n term needs to be searched in
      many different permutations, so it was decided to code functions inside the ff_functions.get_angleangle_data function
      to assist and help with code readability and easy of accessing each K_n parameter. 
      
      Each of the cross terms are searched separately even though they share a given improper/angle type. This allows 
      parameters to be in different order in the forcefield for each cross term or maybe not even there.                                               
      
      The Theta0_n will be found, by 1st attempting to use quartic angles and if it can't be found it will try to be 
      found using quadratic angles. The only reason for doing this is just to try to find as many parameters as possible
      even if the K_n value is zero anyways. This would help if creating assumed_auto_fill features in the future.
    """
    
    # Add printing buffer
    if not skip_printouts: log.out('')

    # angleangle coeffs dictionary
    angleangle_coeffs = {} # {angleangle type : list of coeffs}
    
    ##########################
    # Find angleangle coeffs #
    ##########################
    # flagged angleangles to find angleangle data for (will be implemented at the end of loop)
    flagged_angleangles = BADI.flagged_angleangles    
    for n, i in enumerate(BADI.improper_types_lst):
        
        # Flag things
        skip_flag = False
        
        # bonded atom id types
        type1 = i[0]; type2 = i[1]; type3 = i[2]; type4 = i[3];
        
        # Find improper number based on n, either from improper_types_dict or angleangle_types_dict
        if n < len(BADI.improper_types_dict):
            number = BADI.improper_types_dict[i]
        elif n >= len(BADI.improper_types_dict):
            number = BADI.angleangle_types_dict[i]
        else:
            log.error('ERROR finding which dictionary to find number ID for oop or angleangle')       
        
        # find improper order from improper_map if atom ids were resorted
        if sort_remap_atomids:
            type1, type2, type3, type4 = improper_map[number]
                
		# if ff_class 's2' log generic and continue to next iterations
        if ff_class == 's2':
            # find nb comment and tmp_comment
            if number in flagged_angleangles:
                nb_comment = 'nb!=3'
                tmp_comment = 'Skeleton option used'
            else:
                nb_comment = 'nb==3'
                tmp_comment = 'skipped over b/c number of bonded atoms to central atom == 3 (improper set)'
            
            # save angleangle_coeffs into class
            t = Type()
            t.type = (type1, type2, type3, type4)
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = []
            t.nb = nb_comment
            t.comments = tmp_comment
            angleangle_coeffs[number] = t
            continue  
            
        # Try matching angleangle in all 6 permuations. Attempts: 1) without equivalences, 2) with equivalences
        coeff_123, angleangle_flag, equivalent = ff_functions.get_angleangle_data(type1, type2, type3, type4, log, frc)
        
        # Try finding theta0_123, theta0_124, and theta0_324
        theta0_123, match_123, equiv_123 = ff_functions.get_crossterm_theta0(type1, type2, type3, log, use_auto_equivalence, frc)
        theta0_124, match_124, equiv_124 = ff_functions.get_crossterm_theta0(type1, type2, type4, log, use_auto_equivalence, frc)
        theta0_324, match_324, equiv_324 = ff_functions.get_crossterm_theta0(type3, type2, type4, log, use_auto_equivalence, frc)

        # theta0_123, theta0_124, and theta0_324
        c_123 = 'angle 123:   {:^6} {:^6} {:^6} equivs:  {:^6} {:^6} {:^6} match:  {:^6} {:^6} {:^6}'.format(type1, type2, type3, equiv_123[0], equiv_123[1], equiv_123[2], match_123[0], match_123[1], match_123[2])
        c_124 = 'angle 124:   {:^6} {:^6} {:^6} equivs:  {:^6} {:^6} {:^6} match:  {:^6} {:^6} {:^6}'.format(type1, type2, type4, equiv_124[0], equiv_124[1], equiv_124[2], match_124[0], match_124[1], match_124[2])
        c_324 = 'angle 324:   {:^6} {:^6} {:^6} equivs:  {:^6} {:^6} {:^6} match:  {:^6} {:^6} {:^6}'.format(type3, type2, type4, equiv_324[0], equiv_324[1], equiv_324[2], match_324[0], match_324[1], match_324[2])

        # Build final angleangle_coeff
        coeff_123456 = [coeff_123[0], coeff_123[1], coeff_123[2], theta0_123, theta0_124, theta0_324] 
        
        # Map order to keep standard coeff and cross-term coeff comments identical
        map_order = (type1, type2, type3, type4)
        
        ###########################################################################
        # Only find angleangle types that are not in flagged_angleangle types     #
        ###########################################################################
        if number not in flagged_angleangles: skip_flag = True
            
        # IF skip flag save angleangle_coeff info into class for if it was to be skipped
        if skip_flag:
            t = Type()
            t.type = map_order
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = [0.0, 0.0, 0.0, theta0_123, theta0_324, theta0_124]
            t.comments = '{:<75} {:<100} {:<100} {:<100}'.format('skipped over b/c number of bonded atoms to central atom == 3 (improper set)', c_123, c_124, c_324)
            t.nb = 'nb==3'
            angleangle_coeffs[number] = t
            
        # If angleangle_flag save angleangle coeff into class
        elif angleangle_flag:
            t = Type()
            t.type = map_order
            t.match = equivalent
            t.equivalent = equivalent
            t.coeffs = coeff_123456
            t.comments = '{:<75} {:<100} {:<100} {:<100}'.format('Angleangle type used', c_123, c_124, c_324)
            t.nb = 'nb!=3'
            angleangle_coeffs[number] = t
            
        # else set as zeros and tell user not found
        else:
            if not skip_printouts: log.warn('WARNING unable to find angleangle info for Angleangle Coeff {} {} {} {} {}'.format(number, type1, type2, type3, type4))
            t = Type()
            t.type = map_order
            t.match = ['N/A', 'N/A', 'N/A', 'N/A']
            t.equivalent = ['N/A', 'N/A', 'N/A', 'N/A']
            t.coeffs = [0.0, 0.0, 0.0, theta0_123, theta0_324, theta0_124]
            t.comments = '{:<75} {:<100} {:<100} {:<100}'.format('UNABLE to find coeff parameters', c_123, c_124, c_324)
            t.nb = 'nb!=3'
            angleangle_coeffs[number] = t
    return angleangle_coeffs