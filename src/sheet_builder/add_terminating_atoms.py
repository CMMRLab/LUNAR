# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
September 21st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
###################################
# Function to create atoms object #
###################################
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz .comment .atomtype
def create_atoms(molid, typeint, atomtype, x, y, z, ix, iy, iz):
    a = Atom()
    a.molid = molid
    a.type = typeint
    a.atomtype = atomtype
    a.comment = atomtype
    a.charge = 0.0
    a.x = x
    a.y = y
    a.z = z
    a.ix = ix
    a.iy = iy
    a.iz = iz
    return a


#####################################
# Function to add terminating atoms #
#####################################
def add(atoms, bonds, box, terminating_atoms, log):
    
    # Determine atoms to add
    terminating = terminating_atoms.split(';')
    terminators = {} # {BondingType : [list of terminating atoms]}
    for atom in terminating:
        types = atom.split('|')
        if len(types) >= 2:
            bonding_type = types[0].strip()
            adding_types = [i.strip() for n, i in enumerate(types) if n > 0]
            if bonding_type in terminators:
                log.warn(f'  WARNING BondingType "{bonding_type}" defined twice. Will use last defined BondingType in terminating_atoms.')
            terminators[bonding_type] = adding_types
            
    # Determine newtype integer and log found types
    current_types = sorted(list({atoms[i].type for i in atoms}))
    type_offset = max(current_types)
    types = {} # {(BondingType, TerminatingType) : Type Intger}
    for bonding_type in terminators:
        adding_types = terminators[bonding_type]
        log.out(f'  Will attempt to add {"-".join(adding_types)} to "{bonding_type}" types.')
        for terminating_type in adding_types:
            type_offset += 1
            types[(bonding_type, terminating_type)] = type_offset         
    
    # Generate graph
    graph = {i:[] for i in atoms}
    for id1, id2 in bonds:  
        graph[id1].append(id2)
        graph[id2].append(id1)
        
    # Start terminating atoms (assume non-periodic as that is when terminating atoms is required)
    deltas = {'x':set([0]), 'y':set([0]), 'z':set([0])}
    atoms_count = len(atoms); bonds_count = len(bonds); nbonds = len(bonds)
    terminator_count = {i:0 for i in terminators}
    bond_length = 1.0 # set all bond lengths as 1.0 for the time being ...
    for id1 in graph:
        nb1 = len(graph[id1])
        atom1 = atoms[id1]
        atom_type1 = atom1.atomtype
        x1 = atom1.x; y1 = atom1.y; z1 = atom1.z;
        
        # Start terminating atoms if bonding condition and terminators are described for the atom
        if nb1 < 3 and atom_type1 in terminators:
            
            # Compute normalized position vector
            position_vector = [0, 0, 0]
            for id2 in graph[id1]:
                atom2 = atoms[id2]
                x2 = atom2.x; y2 = atom2.y; z2 = atom2.z;
                position_vector[0] += x1 - x2
                position_vector[1] += y1 - y2
                position_vector[2] += z1 - z2
            vector_length = sum([i*i for i in position_vector])**0.5
            position_vector = [i/vector_length for i in position_vector]
            
            # Start adding terminators
            tmpid = id1
            molid = atom1.molid
            adding_types = terminators[atom_type1]
            for n, atom_type2 in enumerate(adding_types, 1):
                atoms_count += 1
                tx = x1 + bond_length*(n*position_vector[0])
                ty = y1 + bond_length*(n*position_vector[1])
                tz = z1 + bond_length*(n*position_vector[2])
                typeint = types[(atom_type1, atom_type2)]
                atoms[atoms_count] = create_atoms(molid, typeint, atom_type2, tx, ty, tz, atom1.ix, atom1.iy, atom1.iz)
                bonds_count += 1
                bonds.append(tuple(sorted([tmpid, atoms_count])))
                tmpid = atoms_count
                
                # log deltas to adjust simulation cell after
                deltas['x'].add(tx - x1)
                deltas['y'].add(ty - y1)
                deltas['z'].add(tz - z1)
                
            # Tally terminator_count
            terminator_count[atom_type1] += 1
            
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
                
    # Print number of created atoms and bonds
    log.out('  Added the following:')
    for i in terminator_count:
        log.out('  {:>6} atoms to "{}"'.format(terminator_count[i], i))
    log.out('  {:>6} bonds'.format(bonds_count - nbonds))
    
    # Print recommend startup method
    log.out('')
    log.out('  The following startup method is recommended as the addition of terminating atoms uses simple vector math')
    log.out('  and not chemical intuitions based on element types to place the termaniting atoms.')
    log.out('    Step 1:')
    log.out('      Minizime the simulation in LAMMPS with a command like "minimize     1.0e-4 1.0e-6 1000 100000"')
    log.out('')
    log.out('    Step 2:')
    log.out('      Run an small nve/limit simulation to allow the force field to slowly adjust the added atoms')
    log.out('      dynamically. This step may not be required if the minimization was able to move the added atoms')
    log.out('      close to the energy minima of the force field being used to describe the atomic interations.')
    log.out('      Example portion of LAMMPS script to run an nve/limit simulaiton:')
    log.out('        timestep       0.5 # may need to be changed to 0.1')
    log.out('        fix            1 all nve/limit 0.1')
    log.out('        run            20000 # (10000/0.5 = 20000) = 10ps with 0.5 dt')
    log.out('        write_data     added_atoms_initialization.data')
    log.out('        unfix          1')
    return atoms, bonds, box