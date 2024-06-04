# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
May 18th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.sheet_builder.misc_functions as misc_functions
import src.sheet_builder.build_sheets as build_sheets
import math




##################################################################
# Function to generate multiwall nanotubes in armchair or zigzag #
##################################################################
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz .comment .atomtype
def generate_MWCNT(length, diameter, r0, types, edgetype, layer_spacing, ntubes, axis, periodic_bonds, log):
    log.out('\n\nGenerating symmetric nanotube(s) ...')
    
    # Determine diameters
    diameters = [diameter + 2*n*layer_spacing for n in range(ntubes)]
    atoms = {} # { atomID : Atoms Object }
    ID = 0
    for n, diameter in enumerate(diameters, 1):
        # CNT will be generated in Z-dir, where the starting sheet will be generated in
        # the 'xz' plane. Thus lx controls circumference and ly controls nano tube length
        plane = 'xz'
        lx = diameter*math.pi # lx in nanotube
        ly = length # lz in nanotube
        
        # Build sheet (armchair vs zigzag are flipped for sheet generation as sheet gets rotated)
        if edgetype == 'armchair': sheet_edgetype = 'zigzag'
        if edgetype == 'zigzag': sheet_edgetype = 'armchair'
        layer_spacing = 3.354
        nlayers = 1
        stacking = 'AA'
        periodic_bonds_sheets = True
        pflag = False
        sheet_atoms, sheet_box = build_sheets.generate(lx, ly, r0, sheet_edgetype, types, layer_spacing, nlayers, stacking, plane, periodic_bonds_sheets, pflag, log)
        circumference = sheet_box['xhi'] - sheet_box['xlo']
        length = sheet_box['zhi'] - sheet_box['zlo']
        log.out('  edge={} with a diameter={:.4f} and a length={:.4f}'.format(edgetype, diameter, length))
        
        # Start wrapping atoms (Z-coords will stay untouched)
        radius = circumference/(2*math.pi)
        for i in sheet_atoms:
            ID += 1
            atom = sheet_atoms[i]
            arc = atom.x
            theta = arc/radius
            newx = radius*math.cos(theta)
            newy = radius*math.sin(theta)
            newz = atom.z
            a = Atom()
            a.molid = n
            a.type = atom.type
            a.atomtype = atom.atomtype
            a.comment = atom.atomtype
            a.charge = 0.0
            a.x = newx
            a.y = newy
            a.z = newz
            a.ix = 0
            a.iy = 0
            a.iz = 0
            atoms[ID] = a

    # Create box
    true_length = sheet_box['zhi'] - sheet_box['zlo']
    true_diameter = circumference/math.pi
    radial_box_edge = true_diameter/2 + 2*r0 # +- 2*r0 buffer on radial dimension of box
    box = {'xlo': -radial_box_edge,
           'xhi':  radial_box_edge,
           'ylo': -radial_box_edge,
           'yhi':  radial_box_edge,
           'zlo': -true_length/2,
           'zhi':  true_length/2}
    
    # If axis is 'x' or 'y' rotate atoms
    if axis in ['x', 'y']:
        phi = 0; theta = 0; psi = 0;
        if axis == 'x': theta = 90
        if axis == 'y': phi = 90
        atoms = misc_functions.rotate_molecule(atoms, phi, theta, psi)
        
        # Reset box based on rotation
        if axis == 'x': # z->x
            box = {'xlo': -true_length/2,
                   'xhi':  true_length/2,
                   'ylo': -radial_box_edge,
                   'yhi':  radial_box_edge,
                   'zlo': -radial_box_edge,
                   'zhi':  radial_box_edge}
        if axis == 'y': #z->y
            box = {'xlo': -radial_box_edge,
                   'xhi':  radial_box_edge,
                   'ylo': -true_length/2,
                   'yhi':  true_length/2,
                   'zlo': -radial_box_edge,
                   'zhi':  radial_box_edge}
    
    # if system is not period increase box size
    if not periodic_bonds:
        increase = 2*r0
        if axis == 'x':
            box['xlo'] -= increase
            box['xhi'] += increase
        if axis == 'y':
            box['ylo'] -= increase
            box['yhi'] += increase
        if axis == 'z':
            box['zlo'] -= increase
            box['zhi'] += increase
    return atoms, box




#####################################################################
# Function to generate single wall nano tubes that have chirallity  #
#####################################################################
def generate_chiral(n, m, desried_length, r0, types, axis, periodic_bonds, log):  
    log.out('\n\nGenerating a chiral nanotube ...')
    
    # Find circumference and diameter
    circumference = math.sqrt(3)*r0*math.sqrt(n**2 + m**2 + n*m) # will be X-dir
    diameter = (math.sqrt(3)*r0/math.pi)*math.sqrt(n**2 + m**2 + n*m)
    radius = circumference/(2*math.pi)
    
    # Find tube length
    gcd = math.gcd(m, n)
    if (n - m)%3 != 0:
        dr = gcd
    else:
        dr = 3*gcd
    length = 3*r0*math.sqrt(n**2 + m**2 + n*m)/dr # will be Z-dir
    
    # Find alpha and log
    alpha = math.degrees(math.atan( (math.sqrt(3)*m)/(m + 2*n) ))
    log.out('  n={} and m={} generates a tube with a diameter={:.4f}, a length={:.4f}, and an alpha={:.4f}'.format(n, m, diameter, length, alpha))
    
    # Generate graphene like lattice in xz plane (will be geometrically centered about (0, 0, 0))
    max_length = 3*max([circumference, length]) # 3x larger on each side (centered about (0, 0, 0))
    layer_spacing = 3.354
    nlayers = 1
    plane = 'xz'
    stacking = 'AA'
    sheet_edgetype = 'armchair'
    periodic_bonds_sheets = True
    pflag = False
    sheet_atoms, sheet_box = build_sheets.generate(max_length, max_length, r0, sheet_edgetype, types, layer_spacing, nlayers, stacking, plane, periodic_bonds_sheets, pflag, log)
    
    # Rotate the sheet atoms about Y-axis
    phi = 0; theta = -alpha; psi = 0;
    rotated_atoms = misc_functions.rotate_molecule(sheet_atoms, phi, theta, psi)
    
    # Find lower-left atom to "anchor" all other atoms for cutting out a piece of the
    # rotated atoms sheet. The atom should type 2 or 3 (horizontal armchair atoms)
    anchoring_atoms = {} # { atomID : dist2center }
    for i in rotated_atoms:
        atom = rotated_atoms[i]
        if atom.type in [2, 3] and atom.x < -length/2 and atom.z < -circumference/2: 
            anchoring_atoms[i] = misc_functions.compute_distance(atom.x, atom.y, atom.z, 0, 0, 0)
    anchoringID = min(anchoring_atoms, key=anchoring_atoms.get)
    atom = rotated_atoms[anchoringID]
    
    # Find the atoms to "cut out" of the rotated atoms
    cutout_atoms = {} # { atomID : Atoms Object }
    atom1 = rotated_atoms[anchoringID]
    x1 = atom1.x; z1 = atom1.z; ID = 0
    tolerance = r0/6 # += on "cut out" atoms. Can be fairly large since duplicate atomic positions will be checked while "rolling" the tube
    min_x = x1 - tolerance
    min_z = z1 - tolerance
    max_x = x1 + circumference + tolerance
    max_z = z1 + length + tolerance
    zpositions = set()
    for i in rotated_atoms:
        atom2 = rotated_atoms[i]
        x2 = atom2.x; z2 = atom2.z;
        if min_x <= x2 < max_x and min_z <= z2 < max_z:
            ID += 1
            cutout_atoms[ID] = atom2
            zpositions.add(z2)
            
    # Atom positions will be checked if they already exists so setup a domain decomposition to speed-up calculations
    lz = max(zpositions) - min(zpositions); domain_size = 2*r0
    nzz = math.ceil(lz/domain_size)
    if nzz == 0: nzz = 1
    dz = lz/nzz; ID = 0
    zadd = dz/2 +  min(zpositions)
    domain = {} # { domainID : (zlo, zhi) }
    for nz in range(nzz):
        zc = nz*dz + zadd
        zlo_sub = zc - dz/2
        zhi_sub = zc + dz/2
        domain[nz+1] = (zlo_sub, zhi_sub)
        
    # Generate domain_graph
    min_domainID = min(domain); max_domainID = max(domain);
    domain_graph = {} # { domainID : [nearest domainIDs]}
    domain_periodicity = {} # { domainID: near edge flag}
    for domainID in domain:
        domains = [domainID]
        if domainID == min_domainID:
            if min_domainID+1 in domain: domains.append(min_domainID+1)
            if max_domainID in domain: domains.append(max_domainID)
            near_edge = True
        elif domainID == max_domainID:
            if max_domainID-1 in domain: domains.append(max_domainID-1)
            if min_domainID in domain: domains.append(min_domainID)
            near_edge = True
        else:
            if domainID-1 in domain: domains.append(domainID-1)
            if domainID+1 in domain: domains.append(domainID+1)
            near_edge = False
        domain_graph[domainID] = domains 
        domain_periodicity[domainID] = near_edge
        
    # Start wrapping atoms (Z-coords will stay untouched)
    atoms = {} # { atomID : Atoms Object }
    spans = {'x':[], 'y':[], 'z':[]}; guessID = 1; ID = 0
    positions = {i:[] for i in domain} # { domainID : [(x, y, z), ...] }
    for i in cutout_atoms:
        # roll the tube
        atom = cutout_atoms[i]
        arc = atom.x
        theta = arc/radius
        newx = radius*math.cos(theta)
        newy = radius*math.sin(theta)
        newz = atom.z
        spans['x'].append(newx)
        spans['y'].append(newy)
        spans['z'].append(newz)
        position = (newx, newy, newz)
        
        # Check to if atom location already exists
        domainID, guessID = get_domainID(domain, guessID, newz)
        exists = check_if_position_exits(position, domainID, positions, tolerance, length, domain_graph, domain_periodicity)
        
        # if atom position doesnt exist, save atom position
        if not exists:
            positions[domainID].append(position)
            ID += 1
            a = Atom()
            a.molid = n
            a.type = atom.type
            a.atomtype = atom.atomtype
            a.comment = atom.atomtype
            a.charge = 0.0
            a.x = newx
            a.y = newy
            a.z = newz
            a.ix = 0
            a.iy = 0
            a.iz = 0
            atoms[ID] = a

    # Center atoms
    xavg = sum(spans['x'])/len(spans['x'])
    yavg = sum(spans['y'])/len(spans['y'])
    zavg = sum(spans['z'])/len(spans['z'])
    zpositions = []
    for i in atoms:
        atom = atoms[i]
        atom.x -= xavg
        atom.y -= yavg
        atom.z -= zavg
        zpositions.append(atom.z)
    
    # Generate box
    radial_box_edge = radius + 2*r0 # 2*r0 buffer due to non-periodic bonding dimension
    zlo = min(zpositions); zhi = max(zpositions)
    lz = zhi - zlo
    if length - lz > 0:
        buffer = (length - lz)/2
    else: buffer = r0/8 # r0/8 buffer due to floating point errors
    axial_box_lo = zlo - buffer
    axial_box_hi = zhi + buffer
    box = {'xlo': -radial_box_edge,
           'xhi':  radial_box_edge,
           'ylo': -radial_box_edge,
           'yhi':  radial_box_edge,
           'zlo':  axial_box_lo,
           'zhi':  axial_box_hi}
    atoms = misc_functions.wrap_atoms(atoms, box) 
    
    # If axis is 'x' or 'y' rotate atoms
    if axis in ['x', 'y']:
        phi = 0; theta = 0; psi = 0;
        if axis == 'x': theta = 90
        if axis == 'y': phi = 90
        atoms = misc_functions.rotate_molecule(atoms, phi, theta, psi)
        
        # Reset box based on rotation
        if axis == 'x': # z->x
            box = {'xlo':  axial_box_lo,
                   'xhi':  axial_box_hi,
                   'ylo': -radial_box_edge,
                   'yhi':  radial_box_edge,
                   'zlo': -radial_box_edge,
                   'zhi':  radial_box_edge}
        if axis == 'y': #z->y
            box = {'xlo': -radial_box_edge,
                   'xhi':  radial_box_edge,
                   'ylo':  axial_box_lo,
                   'yhi':  axial_box_hi,
                   'zlo': -radial_box_edge,
                   'zhi':  radial_box_edge}
        
    # Adjust height to user perference as best as possible
    nunits = int(round(desried_length/length))
    if nunits > 1:
        # Replicate system in X- or Y- or Z-direction by nunits
        log.out('  Replicating tube by {} units. Creating a final tube length of {:.4f}'.format(nunits, nunits*length))
        nx = 1; ny = 1; nz = 1; center = True
        if axis == 'x': nx = nunits
        if axis == 'y': ny = nunits
        if axis == 'z': nz = nunits
        atoms, box = misc_functions.replicate_atoms(atoms, box, nx, ny, nz, center, increment_molid=False)
        
    # if system is not periodic wrap any atoms with a single non-periodic bond to the other side of the box and increase box size
    if not periodic_bonds:
        atoms = wrap_single_bond_atoms_on_chiral_tube(atoms, box, r0, axis, log)
        increase = 2*r0
        if axis == 'x':
            box['xlo'] -= increase
            box['xhi'] += increase
        if axis == 'y':
            box['ylo'] -= increase
            box['yhi'] += increase
        if axis == 'z':
            box['zlo'] -= increase
            box['zhi'] += increase
    return atoms, box


########################################################
# Function to get position domainID and update guessID #
########################################################
def get_domainID(domain, guessID, newz):
    zlo_sub, zhi_sub = domain[guessID]
    if zlo_sub <= newz < zhi_sub:
        domainID = guessID
    else:
        domainID = guessID
        for j in domain:
            zlo_sub, zhi_sub = domain[j]
            if zlo_sub <= newz <= zhi_sub:
                domainID = j; guessID = j; break
    return domainID, guessID


################################################
# Function to check if position already exists #
################################################
def check_if_position_exits(position, domainID, positions, tolerance, length, domain_graph, domain_periodicity):
    # Setup z-image flags (tube axis is aligned in z, thus is the only periodic face)
    if domain_periodicity[domainID]:
        zimages = [(0, 0, 0), (0, 0, -1), (0, 0, 1)]
    else: zimages = [(0, 0, 0)]
        
    # Compute distances
    distances = []
    for ix, iy, iz in zimages:
        x1 = position[0] + ix*0
        y1 = position[1] + iy*0
        z1 = position[2] + iz*length
        for ID in domain_graph[domainID]:
            for x2, y2, z2 in positions[ID]:
                if abs(x1 - x2) > tolerance: continue
                elif abs(y1 - y2) > tolerance: continue
                elif abs(z1 - z2) > tolerance: continue
                distances.append(misc_functions.compute_distance(x1, y1, z1, x2, y2, z2))
                
    # check if atom is close to another to see if it already exists
    exists = False
    if distances:
        if min(distances) <= tolerance:
            exists = True
    return exists

##########################################################
# Function to wrap upper atoms in a chiral nanotube with #
# a single bond down to the other size to create the a   #
# "stair-step" ring appearence when boundary is 'f f f'  #
##########################################################
def wrap_single_bond_atoms_on_chiral_tube(atoms, box, r0, axis, log):
    # First find bonds for non-periodic system
    pflag = True # do not print bonding results
    max_bonds_per_atom = 3 # only use the first 3 closest atoms to generate bonds
    tolerance = r0/2 # max_distance = r0*tolerance, which tolerance can be large since code will use 3 closets atoms to create bonds
    domain_size = 5*r0 # generate a moderate domain size for domain decomposition (best of all worlds for performance)
    boundary = 'f f f'
    bonds = misc_functions.find_bonds(atoms, box, boundary, r0, tolerance, max_bonds_per_atom, domain_size, pflag, log)
    
    # Generate graph
    graph = {i:[] for i in atoms}
    for id1, id2 in bonds:     
        graph[id1].append(id2)
        graph[id2].append(id1)
        
    # Find atoms that are "near the top of the tube and have a single bond
    atoms2wrap = []
    for i in atoms:
        atom = atoms[i]
        if len(graph[i]) == 1:
            if axis == 'x' and atom.x >= box['xhi'] - 3*r0:
                atoms2wrap.append(i)
            if axis == 'y' and atom.y >= box['yhi'] - 3*r0:
                atoms2wrap.append(i)
            if axis == 'z' and atom.z >= box['zhi'] - 3*r0:
                atoms2wrap.append(i)
    
    # wrap atoms
    if axis == 'x': length = box['xhi'] - box['xlo']
    if axis == 'y': length = box['yhi'] - box['ylo']
    if axis == 'z': length = box['zhi'] - box['zlo']
    for i in atoms2wrap:
        if axis == 'x': atoms[i].x -= length
        if axis == 'y': atoms[i].y -= length
        if axis == 'z': atoms[i].z -= length
    return atoms