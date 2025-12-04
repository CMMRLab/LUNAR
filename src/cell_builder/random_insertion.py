# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
December 4, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.cell_builder.misc_functions as misc_functions
from itertools import product
import math
import time


################################################
# Generate images via minimum image convention #
################################################
def generate_iflags(boundary, log):
    # Boundary conditions
    pflags = boundary.split() # split boundary
    count = pflags.count('f') + pflags.count('p')
    if len(pflags) != 3 and count != 3:
        log.error('ERROR boundary does not contain 3-dimensions or boundary flags are not f or p')
    
    # Image information
    nimages = 1 # only use minimum image convention
    images = set([]) # set to hold unique images
    # Loop through nimages to generate image flags
    for ix in range(-nimages, nimages+1):
        for iy in range(-nimages, nimages+1):
            for iz in range(-nimages, nimages+1):
                
                # Update image based on boundary conditions
                if ix != 0 and pflags[0] == 'f': ix = 0
                if iy != 0 and pflags[1] == 'f': iy = 0
                if iz != 0 and pflags[2] == 'f': iz = 0
            
                # Log image
                images.add( (ix, iy, iz) )
                
    # Sort images in ascending order of magnitude to reduce run time
    images = sorted(images, key=lambda x: sum([abs(i) for i in x]) )
    return images, pflags


###########################################################
# Function to find N-number of possible periodic postions #
###########################################################
def find_periodic_postions(scaled_images, x1, y1, z1, cx, cy, cz, Npos=14):
    postions_distance = {} # {distance from center : (pbc-x, pbc-y, pbc-z) }
    for ixlx, iyly, izlz in scaled_images:
        x1i = x1+ixlx; y1i = y1+iyly; z1i = z1+izlz;
        dist_from_center = misc_functions.compute_distance(x1i, y1i, z1i, cx, cy, cz)
        postions_distance[dist_from_center] = (x1i, y1i, z1i)
    postions_distance = dict(sorted(postions_distance.items(), key=lambda x:abs(x[0]) )) # [0=keys;1=values]
    positions = []
    for N, dist in enumerate(postions_distance):
        if N < Npos: positions.append(postions_distance[dist])
        else: break
    return positions


##################################################################################
# Function to check if the atom is near a box edge. Then build near edge dict to #
# use as a look up table to aviod multiple re-calculations as the loops progress #
##################################################################################
def check_near_edge(x, y, z, max_from_edge, xlo, xhi, ylo, yhi, zlo, zhi):
    if abs(x - xlo) < max_from_edge or abs(x + xlo) < max_from_edge: pbcflag = True
    elif abs(x - xhi) < max_from_edge or abs(x + xhi) < max_from_edge: pbcflag = True        
    elif abs(y - ylo) < max_from_edge or abs(y + ylo) < max_from_edge: pbcflag = True
    elif abs(y - yhi) < max_from_edge or abs(y + yhi) < max_from_edge: pbcflag = True
    elif abs(z - zlo) < max_from_edge or abs(z + zlo) < max_from_edge: pbcflag = True
    elif abs(z - zhi) < max_from_edge or abs(z + zhi) < max_from_edge: pbcflag = True
    else: pbcflag = False
    return pbcflag


##################################################
# Function to generate domain from system "self" #
##################################################
def generate_domain(sys, domain_size, scaled_images, pflag, log):
    #-----------------#
    # Generate domain #
    #-----------------#
    nxx = math.ceil(sys.lx/domain_size)
    nyy = math.ceil(sys.ly/domain_size)
    nzz = math.ceil(sys.lz/domain_size)
    if nxx == 0: nxx = 1
    if nyy == 0: nyy = 1
    if nzz == 0: nzz = 1    
    dxx = sys.lx/nxx
    dyy = sys.ly/nyy
    dzz = sys.lz/nzz
    if pflag: log.out('Using {} x {} x {} sub domains of size {:.2f} x {:.2f} x {:.2f} to perform domain decomposition'.format(nxx, nyy, nzz, dxx, dyy, dzz))
    if pflag: log.out('Finding domain to use cell linked list algorithm for interatomic distance calculations ...')
   
    # Domain-coordinates are indexed starting from
    # 1-to-nii in the xx, yy, and zz-directions
    start_time = time.time()
    coords = list(product(range(1, nxx+1),
                          range(1, nyy+1),
                          range(1, nzz+1)))

    # Generate a map from domainID to Domain-coordinates
    indexes_forward = {} #{ domainID : (nx, ny, nz) }
    indexes_reverse = {} #{ (nx, ny, nz) : domainID }
    for ID, xyz in enumerate(coords, start=1):
        xyz = tuple(xyz)
        indexes_forward[ID] = xyz
        indexes_reverse[xyz] = ID
    ndomain = len(coords)
    execution_time = (time.time() - start_time)
    if pflag: log.out('    Time in seconds to generate domain: ' + str(execution_time))
                
    #----------------#
    # Generate graph #
    #----------------#
    if pflag: log.out('Finding cell linked graph for interatomic distance calculations ...')
    start_time = time.time()
    
    # self+13 neighboring domains. down, east, south priority ordering.
    neighbor_shifts = [(0, 0, 0),    (0, 0, -1), (0, -1, 0),  (0, -1, -1), (-1, 0, 0), (-1, 0, -1), (-1, -1, 0),
                       (-1, -1, -1), (-1, 1, 0), (-1, 1, -1), (0, 1, -1),  (1, 1, -1), (1, 0, -1),  (1, -1, -1)]
    neighbor_shifts = neighbor_shifts[1:] # remove zero shift as we dont need it (part of this comes from Tristan)
    
    domain_graph = {ID:set() for ID in indexes_forward} # { ID : set(bonded IDs) }
    domain_graph[0] = {ID for ID in indexes_forward}
    progress_increment = 10; count = 0;
    for id1 in indexes_forward:
        nx, ny, nz = indexes_forward[id1]
        # Optional printing of progress
        if pflag:
            count += 1
            if 100*count/ndomain % progress_increment == 0:
                log.out('    progress: {} %'.format(int(100*count/ndomain)))
        
        # Start shifting around Domain-coordinates
        for ix, iy, iz in neighbor_shifts:
            mx = nx + ix
            my = ny + iy
            mz = nz + iz
            if mx < 1:
                mx = nxx
            elif mx > nxx:
                mx = 1
                
            if my < 1:
                my = nyy
            elif my > nyy:
                my = 1
                
            if mz < 1:
                mz = nzz
            elif mz > nzz:
                mz = 1
                
            id2 = indexes_reverse[(mx, my, mz)]
            domain_graph[id1].add(id2)
 
    execution_time = (time.time() - start_time)
    if pflag: log.out('    Time in seconds to generate domain connectivity: ' + str(execution_time))
                            
    # setup dict to hold info needed to assign atoms to domain
    atoms2domain = {'indexes_forward':indexes_forward,
                    'indexes_reverse':indexes_reverse,
                    'deltas': [dxx, dyy, dzz],
                    'xlo': sys.xlo,
                    'ylo': sys.ylo,
                    'zlo': sys.zlo}
    return domain_graph, atoms2domain


#####################################
# Function to assing atom to domain #
#####################################
def assign_atom_a_domainID(x, y, z, atoms2domain):
    try:
        xlo = atoms2domain['xlo']
        ylo = atoms2domain['ylo']
        zlo = atoms2domain['zlo']
        dxx, dyy, dzz = atoms2domain['deltas']
        indexes_reverse = atoms2domain['indexes_reverse']
        nx = math.ceil( (x - xlo)/dxx )
        ny = math.ceil( (y - ylo)/dyy )
        nz = math.ceil( (z - zlo)/dzz )
        domainID = indexes_reverse[(nx, ny, nz)]
    except: domainID = 0
    return domainID


#############################
# Function to mix LJ sigmas #
#############################
def mix_LJ_sigmas(sigma1, sigma2, mixing_rule, tolerance):
    # 'geometric' or 'arithmetic' or 'sixthpower'
    mixed_atomsize = tolerance
    if 'arithmetic' in mixing_rule:
        mixed_atomsize = 0.5*(sigma1 + sigma2)
        # typically used with 12-6 LJ, which has minimum energy at 2**(1/6) or 1.12 of sigma
        if 'min' in mixing_rule:
            mixed_atomsize = (2**(1/6))*mixed_atomsize
    elif 'geometric' in mixing_rule:
        mixed_atomsize = math.sqrt(sigma1*sigma2)
        # typically used with 12-6 LJ, which has minimum energy at 2**(1/6) or 1.12 of sigma
        if 'min' in mixing_rule:
            mixed_atomsize = (2**(1/6))*mixed_atomsize
    elif 'sixthpower' in mixing_rule:
        mixed_atomsize = (0.5*(sigma1**6 + sigma2**6))**(1/6)
        # typically used with 9-6 LJ, which has minimum energy at (1.5)**(1/3) or 1.14 of sigma
        if 'min' in mixing_rule:
            mixed_atomsize = (1.5)**(1/3)*mixed_atomsize
    else:
        raise Exception(f"ERROR mixing_rule = {mixing_rule} is not supported. Currently supported rules: 'geometric' or 'arithmetic' or 'sixthpower' or 'geometric-min' or 'arithmetic-min' or 'sixthpower-min'")
    return mixed_atomsize


####################################################################
# Function to check if atom overlaps any in current system: serial #
####################################################################
def overlap_check(sys, m, linked_lst, domain_graph, xshift, yshift, zshift, phi, theta, psi, tolerance, mix_sigma, mixing_rule, boundary_conditions, scaled_images, atoms2domain):
    # Set default overlap, inside Boolean, and insert_molecule
    overlap = False; inside_box = True; insert_molecule = True
    
    # Set mixed and half atom size from tolerance and update later if mix_sigma Boolean Pair Coeffs
    mixed_atomsize = tolerance; half_atomsize = tolerance/2
    
    # Find rotational matrix
    RzRy = misc_functions.matrix_by_matrix(misc_functions.Rz(psi), misc_functions.Ry(theta))
    RzRyRx = misc_functions.matrix_by_matrix(RzRy, misc_functions.Rx(phi))

    # Check if molecule is inside the box
    positions = {} # {atomID:(x, y, z)}
    for id1 in m.atoms:
        atom1 = m.atoms[id1]
        if mix_sigma:
            pair_coeff1 = m.pair_coeffs[atom1.type].coeffs
            half_atomsize = pair_coeff1[tolerance]/2
        x1, y1, z1 = misc_functions.vector_by_matrix(RzRyRx, [atom1.x, atom1.y, atom1.z])
        x1 += xshift
        y1 += yshift
        z1 += zshift
        positions[id1] = (x1, y1, z1)
        if x1 <= sys.xlo + half_atomsize:
            if boundary_conditions[0] == 'f': insert_molecule = False
            inside_box = False
        if x1 >= sys.xhi - half_atomsize:
            if boundary_conditions[0] == 'f': insert_molecule = False
            inside_box = False
        if y1 <= sys.ylo + half_atomsize:
            if boundary_conditions[1] == 'f': insert_molecule = False
            inside_box = False
        if y1 >= sys.yhi - half_atomsize:
            if boundary_conditions[1] == 'f': insert_molecule = False
            inside_box = False
        if z1 <= sys.zlo + half_atomsize:
            if boundary_conditions[2] == 'f': insert_molecule = False
            inside_box = False
        if z1 >= sys.zhi - half_atomsize:
            if boundary_conditions[2] == 'f': insert_molecule = False
            inside_box = False
        if not insert_molecule: break
    
    # Check for overlaps
    if sys.atoms and insert_molecule:
        for id1 in m.atoms:
            if overlap: break
            atom1 = m.atoms[id1]
            if mix_sigma:
                pair_coeff1 = m.pair_coeffs[atom1.type].coeffs
                sigma1 = pair_coeff1[tolerance]
                half_atomsize = pair_coeff1[tolerance]/2
            x1, y1, z1 = positions[id1] 
            if not inside_box:
                if x1 <= sys.xlo: x1 += sys.lx
                if x1 >= sys.xhi: x1 -= sys.lx
                if y1 <= sys.ylo: y1 += sys.ly
                if y1 >= sys.yhi: y1 -= sys.ly
                if z1 <= sys.zlo: z1 += sys.lz
                if z1 >= sys.zhi: z1 -= sys.lz
            edgeflag = check_near_edge(x1, y1, z1, 1.1*half_atomsize, sys.xlo, sys.xhi, sys.ylo, sys.yhi, sys.zlo, sys.zhi)
            if edgeflag: periodic_postions = find_periodic_postions(scaled_images, x1, y1, z1, sys.cx, sys.cy, sys.cz, Npos=14)
            domainID = assign_atom_a_domainID(x1, y1, z1, atoms2domain)
            domains = list(domain_graph[domainID]) + [domainID]
            for domainID_linked in domains:
                if overlap: break
                neighboring_atoms = linked_lst[domainID_linked]
                for id2 in neighboring_atoms:
                    atom2 = sys.atoms[id2]
                    x2 = atom2.x
                    y2 = atom2.y
                    z2 = atom2.z
                    if mix_sigma:
                        pair_coeff2 = atom2.pair_coeff
                        sigma2 = pair_coeff2[tolerance]
                        mixed_atomsize = mix_LJ_sigmas(sigma1, sigma2, mixing_rule, tolerance)
                    if edgeflag:
                        for x1i, y1i, z1i in periodic_postions:
                            if abs(x1i - x2) > mixed_atomsize: continue
                            elif abs(y1i - y2) > mixed_atomsize: continue
                            elif abs(z1i - z2) > mixed_atomsize: continue
                            distance = misc_functions.compute_distance(x1i, y1i, z1i, x2, y2, z2)
                            if distance <= mixed_atomsize:
                                overlap = True
                                break
                    else:
                        if abs(x1 - x2) > mixed_atomsize: continue
                        elif abs(y1 - y2) > mixed_atomsize: continue
                        elif abs(z1 - z2) > mixed_atomsize: continue
                        distance = misc_functions.compute_distance(x1, y1, z1, x2, y2, z2)
                        if distance <= mixed_atomsize:
                            overlap = True
                            break
    return overlap, inside_box, insert_molecule


################################################################################
# Function to determine the amount of reshifting to apply to "stack molecules" #
################################################################################
def insert_reshift(sys, m, linked_lst, domain_graph, xshift, yshift, zshift, phi, theta, psi, tolerance, mix_sigma, 
                   mixing_rule, boundary_conditions, scaled_images, atoms2domain, reshift_cutoff, domain_depth):
    # Set default new shifts
    max_shiftx = 0
    max_shifty = 0
    max_shiftz = 0
    
    # Set mixed and half atom size from tolerance and update later if mix_sigma Boolean Pair Coeffs
    mixed_atomsize = tolerance; half_atomsize = tolerance/2
    
    # Find rotational matrix
    RzRy = misc_functions.matrix_by_matrix(misc_functions.Rz(psi), misc_functions.Ry(theta))
    RzRyRx = misc_functions.matrix_by_matrix(RzRy, misc_functions.Rx(phi))
    
    # Distances for shifting atoms
    if sys.atoms:
        if domain_depth <= 1:
            depth_domain_graph = domain_graph
        else:
            # Generate a domain graph that is "deep enough" for reshift_cutoff
            depth_domain_graph = {} # { ID : set(bonded IDs up to max-depth) }
            for domainID in domain_graph:
                # domainID=0 is all domains (in-case an atom is not in any domain).
                # This means we do not need to find cummulative neighbors and only
                # would slow things down.
                if domainID == 0:
                    depth_domain_graph[domainID] = domain_graph[domainID]
                else:
                    neighbors = find_cumulative_neighs(domain_graph, domainID, domain_depth)
                    bonded = set().union(*neighbors.values())
                    depth_domain_graph[domainID] = bonded
                    
        if not isinstance(reshift_cutoff, (int, float)):
            reshift_cutoff = float('inf') # reset to something extremely large
        
        # Start finding distances
        distances = {} # {(insert_atomID, system_atomID):(dx, dy, dz, distance)}
        mixed_atomsizes = {} # {(insert_atomID, system_atomID):mixed_atomsize}
        for id1 in m.atoms:
            atom1 = m.atoms[id1]
            if mix_sigma:
                pair_coeff1 = m.pair_coeffs[atom1.type].coeffs
                sigma1 = pair_coeff1[tolerance]
                half_atomsize = pair_coeff1[tolerance]/2
            x1, y1, z1 = misc_functions.vector_by_matrix(RzRyRx, [atom1.x, atom1.y, atom1.z])
            x1 += xshift
            y1 += yshift
            z1 += zshift
            
            inside_box = True
            if x1 <= sys.xlo + half_atomsize:
                inside_box = False
            if x1 >= sys.xhi - half_atomsize:
                inside_box = False
            if y1 <= sys.ylo + half_atomsize:
                inside_box = False
            if y1 >= sys.yhi - half_atomsize:
                inside_box = False
            if z1 <= sys.zlo + half_atomsize:
                inside_box = False
            if z1 >= sys.zhi - half_atomsize:
                inside_box = False
            if not inside_box:
                if x1 <= sys.xlo: x1 += sys.lx
                if x1 >= sys.xhi: x1 -= sys.lx
                if y1 <= sys.ylo: y1 += sys.ly
                if y1 >= sys.yhi: y1 -= sys.ly
                if z1 <= sys.zlo: z1 += sys.lz
                if z1 >= sys.zhi: z1 -= sys.lz
            edgeflag = check_near_edge(x1, y1, z1, 1.1*half_atomsize, sys.xlo, sys.xhi, sys.ylo, sys.yhi, sys.zlo, sys.zhi)
            if edgeflag: periodic_postions = find_periodic_postions(scaled_images, x1, y1, z1, sys.cx, sys.cy, sys.cz, Npos=14)
            domainID = assign_atom_a_domainID(x1, y1, z1, atoms2domain)
            domains = list(depth_domain_graph[domainID]) + [domainID]
            for domainID_linked in domains:
                neighboring_atoms = linked_lst[domainID_linked]
                for id2 in neighboring_atoms:
                    atom2 = sys.atoms[id2]
                    x2 = atom2.x
                    y2 = atom2.y
                    z2 = atom2.z
                    if mix_sigma:
                        pair_coeff2 = atom2.pair_coeff
                        sigma2 = pair_coeff2[tolerance]
                        mixed_atomsize = mix_LJ_sigmas(sigma1, sigma2, mixing_rule, tolerance)
                    if edgeflag:
                        periodic_distances = {} # {distance:(dx, dy, dz)}
                        for x1i, y1i, z1i in periodic_postions:
                            if abs(x1i - x2) > reshift_cutoff: continue
                            elif abs(y1i - y2) > reshift_cutoff: continue
                            elif abs(z1i - z2) > reshift_cutoff: continue
                            dx = x2 - x1i
                            dy = y2 - y1i
                            dz = z2 - z1i
                            distance = math.sqrt(dx*dx + dy*dy + dz*dz)
                            periodic_distances[distance] = (dx, dy, dz)
                            
                        # Select closest periodic image
                        if not periodic_distances: continue
                        distance = min(periodic_distances.keys())
                        if distance <= reshift_cutoff:
                            dx, dy, dz = periodic_distances[distance]
                            distances[(id1, id2)] = (dx, dy, dz, distance)
                            mixed_atomsizes[(id1, id2)] = mixed_atomsize
                    else:
                        if abs(x1 - x2) > reshift_cutoff: continue
                        elif abs(y1 - y2) > reshift_cutoff: continue
                        elif abs(z1 - z2) > reshift_cutoff: continue
                        dx = x2 - x1
                        dy = y2 - y1
                        dz = z2 - z1
                        distance = math.sqrt(dx*dx + dy*dy + dz*dz)
                        if distance <= reshift_cutoff:
                            distances[(id1, id2)] = (dx, dy, dz, distance)
                            mixed_atomsizes[(id1, id2)] = mixed_atomsize

            
        # Define the direction to shift the molecule based on the smallest dx, dy, and dz direciton
        if distances:
            # Find the two closests atoms
            closest_atoms, closests_vector = min(distances.items(), key=lambda kv: abs(kv[1][3]) )
            dx, dy, dz, distance = closests_vector # d0 = system - insert
            ux, uy, uz = (-dx/distance, -dy/distance, -dz/distance)
            
            # Check position-vector is pointing in the correct direction
            dot = ux*dx + uy*dy + uz*dz  # = |d| in your current convention
            if dot > 0:
                ux, uy, uz = -ux, -uy, -uz
            
            t_candidates = []
            max_shift = None
            for (id1, id2), (dx, dy, dz, dist) in distances.items():
                rmin = mixed_atomsizes[(id1, id2)] 
                rmin2 = rmin*rmin
            
                # Current squared distance
                d2 = dx*dx + dy*dy + dz*dz
            
                # If already too close, can't move forward at all
                if d2 <= rmin2:
                    max_shift = 0.0
                    break
            
                # Quadratic coefficients for D^2(t) = rmin^2
                # D^2(t) = |d0 + t*u|^2 = t^2 + 2 (u·d0) t + |d0|^2
                dot = ux*dx + uy*dy + uz*dz
                if dot >= 0.0:
                    # If dot >= 0, moving away or sideways; can't create a new overlap
                    continue
                
                A = 1.0
                B = 2.0 * dot
                C = d2 - rmin2
            
                disc = B*B - 4.0*A*C
                if disc < 0.0:
                    # No real intersection → never reaches exactly rmin along this direction
                    continue
            
                sqrt_disc = math.sqrt(disc)
                t1 = (-B - sqrt_disc) / (2.0*A)
                t2 = (-B + sqrt_disc) / (2.0*A)
            
                # We’re currently at t = 0 with D^2(0) > rmin^2.
                # The unsafe region (distance < rmin) lies between the two roots.
                # So the first time we become unsafe is the smallest positive root.
                if t1 > 0.0: t_candidates.append(t1)
                if t2 > 0.0: t_candidates.append(t2)
            
            if t_candidates and max_shift != 0.0:
                max_shift = min(t_candidates)
            else: max_shift = 0 
            
            # Now decompose the shift into components along x,y,z
            max_shiftx = max_shift*ux
            max_shifty = max_shift*uy
            max_shiftz = max_shift*uz
    return max_shiftx, max_shifty, max_shiftz


###########################################################################
# Function to find cumulative neighbors from a node up to a maximum depth #
###########################################################################
def find_cumulative_neighs(graph, node, max_depth):
    """
    Finds the cummulative neighbors of a graph starting at
    a given node and traversing to a maximum depth of max_depth.
    Depth is defined as follows:
        
    1. first neighbors
    2. second neighbors
    3. thrid neighbors
    """
    neighbors = {i+1 : set() for i in range(max_depth)} # { depth : {id1, id2, ...} }
    neighbors[1] = set(graph[node]) 
    visited = set(list(graph[node]) + [node])
    for depth in range(1, max_depth):
        for i in neighbors[depth]:
            for j in graph[i]:
                if j in visited: continue  
                neighbors[depth+1].add(j)
                visited.add(j)
    return neighbors
