# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
September 21st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.sheet_builder.misc_functions as misc_functions
import src.sheet_builder.add_pi_electrons as ape
import random
import math


####################################################
# Function to parse out percentage and BondingType #
####################################################
def parse_bonding_type(string, log):
    # Set defaults
    bonding_type = ''; percentage = ''
    
    # Start parsing
    starting_char = '<'
    starting_flag = False
    ending_char = '>'
    ending_flag = False
    for i in string:
        if i == ' ': continue
        if i == starting_char:
            starting_flag = True
            continue
        if i == ending_char: 
            ending_flag = True
            continue
        if not starting_flag and not ending_flag:
            bonding_type = bonding_type + i
        if starting_flag and not ending_flag:
            percentage = percentage + i
    
    # Process generated strings
    try: percentage = float(percentage.strip())
    except: log.error(f'ERROR functional_atoms {string} is not in the format "BondingType<MaxPercent>". As MaxPercent is not a float or int value.')
    bonding_type = bonding_type.strip()
    return bonding_type, percentage

#####################################
# Function to add terminating atoms #
#####################################
def add(atoms, bonds, box, run_mode, functional_atoms, functional_seed, log):
    if functional_seed > 0: random.seed(functional_seed)
    
    # Define directions to place functional groups. The following options are available for each type:
    #    sheets_direction
    #      - 'positive' which will point the functional groups in the "positive" direction
    #      - 'negative' which will point the functional groups in the "negative" direction
    #    tubes_direction
    #      - 'inward' which will point the functional groups in the "inward" direction 
    #      - 'outward' which will point the functional groups in the "outward" direction 
    sheets_direction = 'positive'
    tubes_direction = 'outward'
    
    # Find simulation cell size
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']
    cx = (box['xhi'] + box['xlo'])/2
    cy = (box['yhi'] + box['ylo'])/2
    cz = (box['zhi'] + box['zlo'])/2
    
    # set max_x, max_y, max_z w/ minimum image convention
    max_x = lx/2; max_y = ly/2; max_z = lz/2;
    
    # Determine atoms to add
    percents = {} # {BondingType : MaxPercent}
    ringed = {} # {BondingType : True or False}
    groups = {} # {BondingType : [list of group atoms]}
    for atom in functional_atoms.split(';'):
        types = atom.split('|')
        if len(types) >= 2:
            bonding_type, percentage = parse_bonding_type(types[0], log)
            adding_types = [i.strip() for n, i in enumerate(types) if n > 0 and len(i) > 0]
            if bonding_type in groups:
                log.warn(f'  WARNING BondingType "{bonding_type}" defined twice. Will use last defined BondingType in functional_atoms.')
            groups[bonding_type] = adding_types
            percents[bonding_type] = percentage
            
            # Determine if the functional group is meant to create a 3-member ring like an epoxide group
            tmp = atom.strip()
            if len(adding_types) and tmp.endswith('|'):
                ringed[bonding_type] = True
            else: ringed[bonding_type] = False            
            
    # Determine newtype integer and log found types
    current_types = sorted(list({atoms[i].type for i in atoms}))
    type_offset = max(current_types)
    types = {} # {(BondingType, TerminatingType) : Type Intger}
    for bonding_type in groups:
        percent = percents[bonding_type]
        adding_types = groups[bonding_type]
        log.out(f'  Will attempt to add {"-".join(adding_types)} to {percent}% of "{bonding_type}" types.')
        for terminating_type in adding_types:
            type_offset += 1
            types[(bonding_type, terminating_type)] = type_offset         
    
    # Generate graph
    graph = {i:[] for i in atoms}
    for id1, id2 in bonds:  
        graph[id1].append(id2)
        graph[id2].append(id1)
        
    # Find which atoms can be functionalized based on groups definition
    functionalizable = {i:[] for i in groups} # {BondingType: [list of possible atoms]}
    for i in atoms:
        nb = len(graph[i])
        atom_type = atoms[i].atomtype
        if atom_type in functionalizable and nb in [2, 3]:
            functionalizable[atom_type].append(i)
            
    # Find quantities of atoms to add
    available = {} # {BondingType : number of available atoms to add functional groups to}
    quantities = {} # {BondingType : number of functional groups to add to the atoms}
    for bonding_type in percents:
        percent = percents[bonding_type]
        navailable = len(functionalizable[bonding_type])
        available[bonding_type] = navailable
        if percent <= 0:
            quantities[bonding_type] = 0
        elif percent >= 100:
            quantities[bonding_type] = navailable
        else:
            quantities[bonding_type] = math.floor( (percent/100)*navailable )
            
    # Go through and start finding which atoms will be functionalized
    atoms2functionalize = {i:[] for i in quantities} # {BondingType : [list of atoms to functionalize]}
    rings2functionalize = {} # {atomID : neighboringID to create ring }
    for bonding_type in quantities:
        qty = quantities[bonding_type]
        ids = functionalizable[bonding_type]#.copy()
        if ringed[bonding_type]:
            increment = 2
        else: increment = 1
        for n in range(0, qty, increment):
            random_index = random.randint(0, len(ids)-1)
            atomid = ids[random_index]
            atoms2functionalize[bonding_type].append(atomid)
            del ids[random_index]
            
            # Find neighboring atom to use to create ring
            if ringed[bonding_type]:
                bonded = [i for i in graph[atomid] if atoms[i].atomtype in groups and i in ids]
                if bonded:
                    random_index = random.randint(0, len(bonded)-1)
                    neighid = bonded[random_index]
                    rings2functionalize[atomid] = neighid
                    
                    # Remove neighid from current ids if it is in ids
                    if neighid in ids:
                        ids.remove(neighid)
                    
                    # Go through and remove neighid from functionalizable
                    for i in functionalizable:
                        if neighid in functionalizable[i]:
                            functionalizable[i].remove(neighid)            
        
    # Go through and start functionalizing atoms
    atoms_count = len(atoms); bonds_count = len(bonds); nbonds = len(bonds)
    deltas = {'x':set([0]), 'y':set([0]), 'z':set([0])}
    group_count = {i:0 for i in groups}
    bond_length = 1.0 # set all bond lengths as 1.0 for the time being ...
    for bonding_type in atoms2functionalize:
        for id1 in atoms2functionalize[bonding_type]:
            atom1 = atoms[id1]
            atom_type1 = atom1.atomtype
            x1 = atom1.x; y1 = atom1.y; z1 = atom1.z;
            
            # Loop through 1st-neighs to get xyz data (remove PBCs as well)
            xyz = []
            for id2 in graph[id1]:
                atom2 = atoms[id2]
                x2 = atom2.x; y2 = atom2.y; z2 = atom2.z;
                
                # Shift x, y and z-direction if periodic
                x11, x22 = ape.shift_pbc_atom(x1, x2, lx, max_x)
                y11, y22 = ape.shift_pbc_atom(y1, y2, ly, max_y)
                z11, z22 = ape.shift_pbc_atom(z1, z2, lz, max_z)
                
                # Log I22-postions
                xyz.append([x22, y22, z22])
                
            # If len(xyz) <= 2 add I11-postion to xyz lst (to add more atoms to fit plane to)
            if len(xyz) <= 2: xyz.append([x11, y11, z11])
            
            # Fitting plane to 1st neighbors non-periodic positions
            c, normal = ape.fitplane(xyz)
            
            # Start adding terminators
            tmpid = id1
            molid = atom1.molid
            adding_types = groups[atom_type1]
            for n, atom_type2 in enumerate(adding_types, 1):
                # Point functional groups accordingly based on inputs and direction of normal vector
                tx_pos = x1 + bond_length*(n*normal[0])
                ty_pos = y1 + bond_length*(n*normal[1])
                tz_pos = z1 + bond_length*(n*normal[2])
                tx_neg = x1 - bond_length*(n*normal[0])
                ty_neg = y1 - bond_length*(n*normal[1])
                tz_neg = z1 - bond_length*(n*normal[2])
                if run_mode in ['symmetric-tube', 'chiral-tube']:
                    distance_pos = misc_functions.compute_distance(cx, cy, cz, tx_pos, ty_pos, tz_pos)
                    distance_neg = misc_functions.compute_distance(cx, cy, cz, tx_neg, ty_neg, tz_neg)
                    if distance_pos > distance_neg:
                        if tubes_direction == 'outward': tx = tx_pos; ty = ty_pos; tz = tz_pos
                        elif tubes_direction == 'inward': tx = tx_neg; ty = ty_neg; tz = tz_neg
                        else: raise Exception(f'ERROR tubes_direction {tubes_direction} not supported. Supported directions are "inward" or "outward"')
                    else:
                        if tubes_direction == 'outward': tx = tx_neg; ty = ty_neg; tz = tz_neg
                        elif tubes_direction == 'inward': tx = tx_pos; ty = ty_pos; tz = tz_pos
                        else: raise Exception(f'ERROR tubes_direction {tubes_direction} not supported. Supported directions are "inward" or "outward"')
                elif run_mode == 'sheet':
                    if sheets_direction == 'positive': tx = tx_pos; ty = ty_pos; tz = tz_pos
                    elif sheets_direction == 'negative': tx = tx_neg; ty = ty_neg; tz = tz_neg
                    else: raise Exception(f'ERROR sheets_direction {sheets_direction} not supported. Supported directions are "positive" or "negative"')
                else:
                    tx = tx_pos; ty = ty_pos; tz = tz_pos
                
                # Special "ringed functionalization" to create functional groups like epoxide rings in a graphene sheet
                if id1 in rings2functionalize:
                    id3 = rings2functionalize[id1]
                    atom3 = atoms[id3]
                    x3 = atom3.x; y3 = atom3.y; z3 = atom3.z;
                    
                    # Move atom accordingly
                    tx += (x3 - x1)/2
                    ty += (y3 - y1)/2
                    tz += (z3 - z1)/2
                    
                    # create atom and bonds
                    atoms_count += 1
                    typeint = types[(atom_type1, atom_type2)]
                    atoms[atoms_count] = ape.create_atoms(molid, typeint, atom_type2, tx, ty, tz, atom1.ix, atom1.iy, atom1.iz)
                    bonds_count += 2
                    bonds.append(tuple(sorted([tmpid, atoms_count])))
                    bonds.append(tuple(sorted([id3, atoms_count])))
                    tmpid = atoms_count
                    
                    # set natoms increment for tallying the number of functionalized atoms
                    natoms_increment = 2
                    
                # Normal "linear functionalization"
                else:
                    atoms_count += 1
                    typeint = types[(atom_type1, atom_type2)]
                    atoms[atoms_count] = ape.create_atoms(molid, typeint, atom_type2, tx, ty, tz, atom1.ix, atom1.iy, atom1.iz)
                    bonds_count += 1
                    bonds.append(tuple(sorted([tmpid, atoms_count])))
                    tmpid = atoms_count
                    
                    # set natoms increment for tallying the number of functionalized atoms
                    natoms_increment = 1
                
                    # log deltas to adjust simulation cell after
                    deltas['x'].add(tx - x1)
                    deltas['y'].add(ty - y1)
                    deltas['z'].add(tz - z1)
                
            # Tally number of functionalized atoms
            group_count[atom_type1] += natoms_increment
            
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
    for i in group_count:
        added = group_count[i]
        total = available[i]
        if total > 0: percent = 100*(added/total)
        else: percent = 0
        log.out('  {:>6} atoms to "{}" ({}/{} = {:.4f}%)'.format(added, i, added, total, percent))
    log.out('  {:>6} bonds'.format(bonds_count - nbonds))
    
    # Print recommend startup method
    log.out('')
    log.out('  The following startup method is recommended as the addition of funcational atoms uses simple vector math')
    log.out('  and not chemical intuitions based on element types to place the functionalizing atoms.')
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