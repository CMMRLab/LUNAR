# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
October 23rd, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.sheet_builder.add_functional_groups as add_functional_groups
import src.sheet_builder.misc_functions as misc_functions
import src.sheet_builder.add_pi_electrons as ape
import src.sheet_builder.atoms2lmp as atoms2lmp
import src.mol2SYBYL2lmp as mol2SYBYL2lmp
import src.write_lmp as write_lmp
import src.read_lmp as read_lmp
import src.mol2lmp as mol2lmp
import src.pdb2lmp as pdb2lmp
import math
import copy
import os



#################################################################################
# Function to parse out percentage BondingType, percentage, direction and molID #
#################################################################################
def parse_bonding_type(string, log):
    # Set defaults
    bonding_type = ''; bounded1 = ''; bounded2 = ''
    
    # Start parsing
    starting_char = '<'
    starting_count = 0
    ending_char = '>'
    ending_count = 0
    for i in string:
        if i == ' ': continue
        if i == starting_char:
            starting_count += 1
            continue
        if i == ending_char: 
            ending_count += 1
            continue
        if starting_count == 0 and ending_count == 0:
            bonding_type = bonding_type + i
        if starting_count == 1 and ending_count <= 1:
            bounded1 = bounded1 + i
        if starting_count == 2 and ending_count < 2:
            bounded2 = bounded2 + i
            
    # Get strings from bounded1
    bonding_type = bonding_type.strip()
    bounded1 = bounded1.split(',')
    if len(bounded1) == 3:
        percentage = bounded1[0]
        direction = bounded1[1]
        try: molID = int(bounded1[2])
        except: molID = bounded1[2]
    else:
        percentage = ''
        direction = '+'
        molID = '*'

    # Get strings, floats, or ints from bounded2          
    atomIDs = []
    for i in bounded2.split(','):
        atomID = i.strip()
        if atomID == '': continue
        try: 
            atomID = float(atomID)
            if atomID.is_integer():
                atomID = int(atomID)
        except: pass
        atomIDs.append(atomID)

    # Process generated strings
    try: percentage = float(percentage.strip())
    except: log.error(f'ERROR functional_atoms {string} is not in the format "BondingType<MaxPercent,Direction,MolID><ID1, ID2, IDN>". As MaxPercent is not a float or int value.')
    return bonding_type, percentage, direction, molID, atomIDs


#################################################
# Function to find a cylinder that encapsulates #
# the orientated molecule in the Z-direction.   #
#################################################
def find_encapsuling_cylinder(m, orientation_vector, log, write_data=False):
    # Rotate molecule to point in the Z-direction
    m.atoms = misc_functions.rotate_molecule_from_a_to_b(m.atoms, orientation_vector, (0, 0, -1))
    
    # Find all new positions of rotated molecule
    x = []; y = []; z = []
    for i in m.atoms:
        atom = m.atoms[i]
        x.append(atom.x)
        y.append(atom.y)
        z.append(atom.z)
    
    # Find new box dimensions
    xlo = min(x); xhi = max(x)
    ylo = min(y); yhi = max(y)
    zlo = min(z); zhi = max(z)
    box = {'xlo':xlo, 'xhi': xhi,
           'ylo':ylo, 'yhi': yhi,
           'zlo':zlo, 'zhi': zhi}
    
    # Find box that encapsults the atoms, then derive
    # the cylinder parameters from the box dimensions
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']
    
    # Define cylinder parameters from box dimensions
    h = lz
    d = (lx**2 + ly**2)**0.5
    
    # Write an example datafile if desired
    if write_data:
        # Create m-object like read_lmp would generate
        atoms = m.atoms
        bonds = [m.bonds[i].atomids for i in m.bonds]
        masses = {}
        header = 'sheet_builder.find_encapsuling_cylinder()'
        basename = m.filename[:m.filename.rfind('.')] + '_find_encapsuling_cylinder'
        mm = atoms2lmp.Molecule_File(basename, atoms, bonds, box, masses, header)
        
        # Write LAMMPS datafile
        atom_style = 'full'
        include_type_labels = False
        write_lmp.file(mm, basename+'.data', header, atom_style, include_type_labels, log)
    return d, h


#####################################
# Function to add terminating atoms #
#####################################
def add(atoms, bonds, box, run_mode, grafting_files, seed, functional_atoms, boundary, bond_length, minimum_distance, molID_attributes, log):
    
    # We will use add_function_groups to find where to add in the files. So we 
    # will need to generate a valid functional_atoms string, where the atomIDs
    # are removed. The Type1 will be set as the filename for easy passing of
    # information in this script
    dummy_functional_atoms = []
    grafting = {} # {'graftID:N' : (atomIDs, orientation_vector, mol-object) }
    cylinders = {} # {'graftID:N' : (d, h)}
    for graftID, atom in enumerate(grafting_files.split(';')):
        atom = atom.strip()
        types = atom.split('|')
        if len(types) >= 2:
            atomtype = 'graftID:{}'.format(graftID)
            filename = types[1]
            bonding_type, percentage, direction, molID, atomIDs = parse_bonding_type(types[0], log)
            if len(atomIDs) == 1:
                tmp_functional_atoms = '{}<{},{},{}>|{}'.format(bonding_type, percentage, direction, molID, atomtype)
            elif len(atomIDs) == 2:
                tmp_functional_atoms = '{}<{},{},{}>|{}|'.format(bonding_type, percentage, direction, molID, atomtype)
            else:
                log.error(f'  ERROR grafting atoms string: {atom} has {len(atomIDs)} defined. Currently only 1 atomID or 2 atomIDs can be defined')
        
            # Start reading in files
            if filename.endswith('mol') or filename.endswith('sdf'):
                if os.path.isfile(filename):
                    m = mol2lmp.Molecule_File(filename)
                    log.out(f'  Read in {filename} chemdraw .mol or .sdf file')
                    
                    # We will "monkey paste" a per-atom attribute of .atomtype from the
                    # element attribute. The rest of this code will then use the .atomtype
                    # attribute to represent the added "atom types" from this file.
                    for i in m.atoms:
                        m.atoms[i].atomtype = m.atoms[i].element
                else: log.error(f'  ERROR .mol or .sdf file: {filename} does not exist')
                
            # Read .mol2 file (VMD MDL MOL2 file)
            elif filename.endswith('mol2'):
                if os.path.isfile(filename):
                    m = mol2SYBYL2lmp.Molecule_File(filename)
                    log.out(f'  Read in {filename} SYBYL MOL2 file')
                    
                    # We will "monkey paste" a per-atom attribute of .atomtype from the
                    # element attribute. The rest of this code will then use the .atomtype
                    # attribute to represent the added "atom types" from this file.
                    for i in m.atoms:
                        m.atoms[i].atomtype = m.atoms[i].element
                else: log.error(f'  ERROR .mol2 file: {filename} does not exist')
                
            # Read .pdb file
            elif filename.endswith('pdb'):
                if os.path.isfile(filename):
                    m = pdb2lmp.Molecule_File(filename)
                    log.out(f'  Read in {filename} pdb file')
                    
                    # We will "monkey paste" a per-atom attribute of .atomtype from the
                    # element attribute. The rest of this code will then use the .atomtype
                    # attribute to represent the added "atom types" from this file.
                    for i in m.atoms:
                        m.atoms[i].atomtype = m.atoms[i].atom_name
                else: log.error(f'  ERROR .pdb file: {filename} does not exist')
                
            elif filename.endswith('data') or filename.endswith('dat') or filename.endswith('data.gz') or filename.endswith('dat.gz'):
                if os.path.isfile(filename):
                    m = read_lmp.Molecule_File(filename, method='forward', sections = ['Atoms', 'Bonds'])
                    log.out(f'Read in {filename} LAMMPS datafile')

                    # We will "monkey paste" a per-atom attribute of .atomtype from the
                    # type label (if present) or comment from masses section (if a type
                    # label does not exist). The rest of this code will then use the
                    # .atomtype attribute to represent the added "atom types" from this
                    # file.
                    methods = set()
                    for i in m.atoms:
                        atom = m.atoms[i]
                        if atom.type in m.atom_type_labels_reverse:
                            monkey_type = m.atom_type_labels_reverse[atom.type]
                            methods.add('type labels')
                        elif atom.type in m.masses: 
                            monkey_type = m.masses[atom.type].type
                            if monkey_type == 'N/A': monkey_type = str(atom.type)
                            methods.add('masses comment')
                        else:
                            monkey_type = str(atom.type)
                            methods.add('numeric LAMMPS atomTypeID')
                        m.atoms[i].atomtype = monkey_type
                    log.out(f'  The atom type was set based on the following method(s): "{" ".join(methods)}"')
                else: log.error(f'ERROR lammps datafile: {filename} does not exist')
                
            else: log.error(f'  ERROR filename {filename} has an unsupported file extension')
            
            # We will need to define the postion vector(s) from the atomIDs, which
            # will be done by finding another atomID in the molecule that is the 
            # furtherest from the current position (either an atom position or centered
            # between two atom positions)
            if len(atomIDs) == 1:
                atom = m.atoms[atomIDs[0]]
                current_position = (atom.x, atom.y, atom.z)
            if len(atomIDs) == 2:
                atom1 = m.atoms[atomIDs[0]]
                atom2 = m.atoms[atomIDs[1]]
                current_position = ( (atom1.x + atom2.x)/2, (atom1.y + atom2.y)/2, (atom1.z + atom2.z)/2 )
                
            # Compute unit position vector
            orientation_vector = misc_functions.find_average_orientation_from_a_point(m, current_position)
            # furtherest_atom = m.atoms[misc_functions.furtherest_away_atomID(m, current_position)]
            # reference_postion = (furtherest_atom.x, furtherest_atom.y, furtherest_atom.z)
            # orientation_vector = misc_functions.compute_unit_position_vector(current_position, reference_postion)
            
            # We should shift the entire molecule to the "current_position", so
            # that all future rotations will happen "around" that point in space
            shiftx = 0 - current_position[0]
            shifty = 0 - current_position[1]
            shiftz = 0 - current_position[2]
            m.atoms = misc_functions.shift_molecule(m.atoms, shiftx, shifty, shiftz)
            
            # Find the cylinder parameters that encapsolute the orinetated fragment
            d, h = find_encapsuling_cylinder(copy.deepcopy(m), orientation_vector, log, write_data=False)
            
            # Finally log some information for the next step
            dummy_functional_atoms.append(tmp_functional_atoms)
            grafting[atomtype] = (atomIDs, orientation_vector, m)
            cylinders[atomtype] = (d, h)
    
    
    # We will add in the functional_atoms atoms here as well, so users
    # can functionalize and graft sheets or tubes all at the same time
    for i in functional_atoms.split(';'):
        dummy_functional_atoms.append(i.strip())
        
        
    # Derived a new minimum_distance if user wants
    log.out('')
    if 'cylinder' in str(minimum_distance):
        # Find larget diameter
        key, value = max(cylinders.items(), key=lambda kv:kv[1][0])
        largest_diameter = value[0]
        try: scale_factor = float(minimum_distance.split(':')[-1])
        except: scale_factor = 1.0
        log.out('  Using cylinder minimum_distance option with the following parameters:')
        log.out(f'    largest_diameter = {largest_diameter}')
        log.out(f'    scale_factor     = {scale_factor}')
        log.out( '    minimum_distance = scale_factor*largest_diameter')
        old_minimum_distance = minimum_distance
        minimum_distance = scale_factor*largest_diameter
        log.out(f'  Updated minimum_distance from "{old_minimum_distance}" to {minimum_distance}')
        
        
    # Compute the theortical maximum percent that can be functionalized with the minimum distance constraint
    log.out('')
    log.out('  Theoretical maximum grafting percent per molID based on hexagonal packing & minimum_distance')
    log.out('    Nlength = floor(length/minimum_distance)')
    log.out('    Nwidth  = floor(width/minimum_distance)')
    log.out('    Ntotal  = (2*Nlength*Nwidth) - Nlength')
    log.out('    Max (%) = 100*(Ntotal/Natoms)')
    log.out('   -----------------------------------------------------------------------------------')
    log.out('    {:^6} {:^10} {:^10} {:^10} {:^10} {:^10} {:^10} {:^10}'.format('molID', 'length', 'width', 'Nlength', 'Nwidth', 'Ntotal', 'Natoms', 'Max (%)'))
    log.out('   -----------------------------------------------------------------------------------')  
    for molID in molID_attributes:
        length, width, natoms = molID_attributes[molID]
        nlength = math.floor(length/minimum_distance)
        nwidth = math.floor(width/minimum_distance)

        # Compute the number of circles that can fit in a NxM area of hexagons
        # https://medium.com/science-spectrum/the-art-of-circle-packing-3b025076cddd
        n_rows = nlength
        n_columns = nwidth
        ntotal = (2*n_rows*n_columns) - n_rows
        max_percent = '{:.2f}'.format(100*(ntotal/natoms))
        log.out('    {:^6} {:^10.4f} {:^10.4f} {:^10} {:^10} {:^10} {:^10} {:^10}'.format(molID, length, width, nlength, nwidth, ntotal, natoms, max_percent))
    log.out('')
        
    
    # Use add_function_groups to find where to add in the files
    dummy_functional_atoms = '; '.join(dummy_functional_atoms)
    atoms, bonds, box = add_functional_groups.add(atoms, bonds, box, run_mode, dummy_functional_atoms, seed, boundary, bond_length, minimum_distance, log)
    
    # Generate graph
    graph = {i:[] for i in atoms}
    for id1, id2 in bonds:  
        graph[id1].append(id2)
        graph[id2].append(id1)
    
    # Loop through atoms and find atoms of graftID:N, which are functional groups
    # that will be termed as "stud_atoms" as that will be the "stud" we use to 
    # attach the grafting atoms to
    atoms_count = max([max(atoms), len(atoms)]) # Use max of max or len to allow for non-contiguous atomIDs
    add_atoms = {} # { atomID : atom-object }
    add_bonds = [] # [(id1, id2), n-new-bonds]
    del_atoms = [] # [id1, id2, n-del-atoms]
    del_bonds = [] # [(id1, id2), n-new-bonds]
    deltas = {'x':set([0]), 'y':set([0]), 'z':set([0])}
    zero_tolerance = 0.1 # Will be used to determine if unit vector of "stud" is near zero for box modifications
    for i in atoms:
        stud_atom = atoms[i]

        if stud_atom.atomtype in grafting:
            # Determine the "Studded functional group" unit position vector
            if not graph[i]: continue
            neigh_positions = ([], [], []) # ([x1, x2, ...], [y1, y2, ...], [z1, z2, ...])
            for neighID in graph[i]:
                sheet_atom = atoms[neighID]
                neigh_positions[0].append(sheet_atom.x)
                neigh_positions[1].append(sheet_atom.y)
                neigh_positions[2].append(sheet_atom.z)
            current_position = tuple([sum(cartesian)/len(cartesian) for cartesian in neigh_positions])
            reference_postion = (stud_atom.x, stud_atom.y, stud_atom.z)            
            stud_vector = misc_functions.compute_unit_position_vector(current_position, reference_postion)
            
            # Access the grafting vector molecule data
            atomIDs, orientation_vector, m = grafting[stud_atom.atomtype]
            
            # Generate a copy of grafting molecules and bonds, to not disturb
            # the orginals and then the rotate grafting molecule atoms
            grafting_atoms = misc_functions.rotate_molecule_from_a_to_b(copy.deepcopy(m.atoms), orientation_vector, stud_vector)
            grafting_bonds = copy.deepcopy(m.bonds)
            
            # If there are two atomIDs and the stud atom has two atomIDs we need to
            # apply one more rotation to force the "unit tangent vectors" to align
            if len(graph[i]) == 2 and len(atomIDs) == 2:
                s1, s2 = [atoms[k] for k in graph[i]]
                g2, g1 = [grafting_atoms[k] for k in atomIDs] # switch g1 and g2, so vector points in opposite direction (helps with setting correct bonds)
                stud_tangent_vector = misc_functions.compute_unit_position_vector((s1.x, s1.y, s1.z), (s2.x, s2.y, s2.z))
                grafting_tangent_vector = misc_functions.compute_unit_position_vector((g1.x, g1.y, g1.z), (g2.x, g2.y, g2.z))
                grafting_atoms = misc_functions.rotate_molecule_from_a_to_b(grafting_atoms, grafting_tangent_vector, stud_tangent_vector)
                
            
            # Start adding the grafting atoms to the molecular system
            atomID_map = {} # {graft_atomID : added_atomID }
            for j in grafting_atoms:
                grafting_atom = grafting_atoms[j]
                
                # Move atom accordingly
                tx = stud_atom.x - grafting_atom.x 
                ty = stud_atom.y - grafting_atom.y 
                tz = stud_atom.z - grafting_atom.z
                
                # create atom and bonds
                atoms_count += 1
                typeint = 1
                add_atoms[atoms_count] = ape.create_atoms(stud_atom.molid, typeint, grafting_atom.atomtype, tx, ty, tz, stud_atom.ix, stud_atom.iy, stud_atom.iz)
                atomID_map[j] = atoms_count
                
                # log deltas to adjust simulation cell after
                if abs(stud_vector[0]) > zero_tolerance: deltas['x'].add(tx - stud_atom.x)
                if abs(stud_vector[1]) > zero_tolerance: deltas['y'].add(ty - stud_atom.y)
                if abs(stud_vector[2]) > zero_tolerance: deltas['z'].add(tz - stud_atom.z)
                
            # Add bonds that are defined in the grafting molecule to the system
            for k in grafting_bonds:
                id1, id2 = grafting_bonds[k].atomids
                add_bonds.append((atomID_map[id1], atomID_map[id2]))
                
            # Add the bonds that will bond the molecule to the system
            for neighID, id2 in zip(graph[i], atomIDs):
                add_bonds.append((neighID, atomID_map[id2]))  
                
            # Remove the "Stud atom" and the "Stud bond" from the system
            del_atoms.append(i)
            for neighID in graph[i]:
                del_bonds.append((i, neighID))
                del_bonds.append((neighID, i))
                
    # Go through and add all new atoms and bonds
    for i in add_atoms:
        atoms[i] = add_atoms[i]
    for bond in add_bonds:
        bonds.append(bond)
    
    # Go through and delete the "Stud atoms and bonds"
    for i in del_atoms:
        try: del atoms[i]   
        except: pass
    for bond in del_bonds:
        try: bonds.remove(bond)  
        except: pass   
    
    # Adjust box based on changes in relative atom positions
    scale = 1.0
    try: max_dx = scale*max([abs(i) for i in deltas['x']])
    except: max_dx = 0
    try: max_dy = scale*max([abs(i) for i in deltas['y']])
    except: max_dy = 0
    try: max_dz = scale*max([abs(i) for i in deltas['z']])
    except: max_dz = 0
    box['xlo'] -= max_dx
    box['xhi'] += max_dx
    box['ylo'] -= max_dy
    box['yhi'] += max_dy
    box['zlo'] -= max_dz
    box['zhi'] += max_dz
    
    # Go through and wrap atoms
    atoms = misc_functions.wrap_atoms(atoms, box)
    
    return atoms, bonds, box