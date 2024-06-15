# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
June 14th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.cell_builder.misc_functions as misc_functions
import math


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
def find_periodic_postions(scaled_images, x1, y1, z1, cx, cy, cz, Npos=12):
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
    # Generate domain
    nxx = math.ceil(sys.lx/domain_size)
    nyy = math.ceil(sys.ly/domain_size)
    nzz = math.ceil(sys.lz/domain_size)
    cx = (sys.xhi + sys.xlo)/2
    cy = (sys.yhi + sys.ylo)/2
    cz = (sys.zhi + sys.zlo)/2
    if nxx == 0: nxx = 1
    if nyy == 0: nyy = 1
    if nzz == 0: nzz = 1
    if pflag: log.out(f'Using {nxx} x {nyy} x {nzz} sub domains of size {domain_size} to perform domain decomposition')
    dx = sys.lx/nxx; dy = sys.ly/nyy; dz = sys.lz/nzz;
    halfdx = dx/2; halfdy = dy/2; halfdz = dz/2; ID = 0;
    xadd = halfdx + sys.xlo; yadd = halfdy + sys.ylo; zadd = halfdz + sys.zlo;
    domain = {} # { domainID : (xlo, xhi, ylo, yhi, zlo, zhi, r, edgeflag) }
    if pflag: log.out('Finding domain to use cell linked list algorithm for interatomic distance calculations ...')
    ndomain = nxx*nyy*nzz; progress_increment = 10; count = 0;
    for nx in range(nxx):
        for ny in range(nyy):
            for nz in range(nzz):
                # Optional printing of progress
                if pflag:
                    count += 1
                    if 100*count/ndomain % progress_increment == 0:
                        if pflag: log.out('    progress: {} %'.format(int(100*count/ndomain)))
                
                # Generate domain center
                ID += 1
                xc = nx*dx + xadd
                yc = ny*dy + yadd
                zc = nz*dz + zadd
                xlo_sub = xc - halfdx
                xhi_sub = xc + halfdx
                ylo_sub = yc - halfdy
                yhi_sub = yc + halfdy
                zlo_sub = zc - halfdz
                zhi_sub = zc + halfdz
                r = misc_functions.compute_distance(xc, yc, zc, cx, cy, cz)
                edgeflag = check_near_edge(xc, yc, zc, domain_size, sys.xlo, sys.xhi, sys.ylo, sys.yhi, sys.zlo, sys.zhi)
                domain[ID] = (xlo_sub, xhi_sub, ylo_sub, yhi_sub, zlo_sub, zhi_sub, xc, yc, zc, r, edgeflag)
                
    # Find domain connectivity (graph)
    domain_graph = {ID:set() for ID in domain} # { ID : set(bonded IDs) }
    domain_graph[0] = {ID for ID in domain}
    sub_domain = {i for i in domain}
    if pflag: log.out('Finding cell linked graph for interatomic distance calculations ...')
    progress_increment = 10; count = 0;
    for id1 in domain:
        # Optional printing of progress
        if pflag:
            count += 1
            if 100*count/ndomain % progress_increment == 0:
                if pflag: log.out('    progress: {} %'.format(int(100*count/ndomain)))
        d1 = domain[id1]
        x1 = d1[6]; y1 = d1[7]; z1 = d1[8]; r1 = d1[9]; edge1 = d1[10];
        min_radius = r1 - 2*domain_size
        max_radius = r1 + 2*domain_size
        sub_domain.remove(id1)
        if edge1: periodic_postions = find_periodic_postions(scaled_images, x1, y1, z1, cx, cy, cz, Npos=12)
        for id2 in sub_domain:
            if min_radius < domain[id2][9] < max_radius:
                d2 = domain[id2]
                x2 = d2[6]; y2 = d2[7]; z2 = d2[8]; edge2 = d2[10];
                if edge1 and edge2: # periodic
                    for x1i, y1i, z1i in periodic_postions:
                        if abs(x1i - x2) > domain_size: continue
                        elif abs(y1i - y2) > domain_size: continue
                        elif abs(z1i - z2) > domain_size: continue
                        domain_graph[id1].add(id2)
                        domain_graph[id2].add(id1)
                        break
                else: # non-periodic
                    if abs(x1 - x2) > domain_size: continue
                    elif abs(y1 - y2) > domain_size: continue
                    elif abs(z1 - z2) > domain_size: continue
                    domain_graph[id1].add(id2)
                    domain_graph[id2].add(id1)
    return domain, domain_graph


#####################################
# Function to assing atom to domain #
#####################################
def assign_atom_a_domainID(x, y, z, guess, domain):
    domainID = 0
    xlo_sub, xhi_sub, ylo_sub, yhi_sub, zlo_sub, zhi_sub, xc, yc, zc, r2, edgeflag = domain[guess]
    if xlo_sub <= x <= xhi_sub and ylo_sub <= y <= yhi_sub and zlo_sub <= z <= zhi_sub:
        domainID = guess
    else:    
        for j in domain:
            xlo_sub, xhi_sub, ylo_sub, yhi_sub, zlo_sub, zhi_sub, xc, yc, zc, r2, edgeflag = domain[j]
            if xlo_sub <= x <= xhi_sub and ylo_sub <= y <= yhi_sub and zlo_sub <= z <= zhi_sub: 
                domainID = j; guess = j; break
    return domainID, guess


#############################
# Function to mix LJ sigmas #
#############################
def mix_LJ_sigmas(sigma1, sigma2, mixing_rule, tolerance):
    # 'geometric' or 'arithmetic' or 'sixthpower'
    mixed_atomsize = tolerance
    if 'arithmetic' in mixing_rule:
        mixed_atomsize = 0.5*(sigma1 + sigma2)
    elif 'geometric' in mixing_rule:
        mixed_atomsize = math.sqrt(sigma1*sigma2)
    elif 'sixthpower' in mixing_rule:
        mixed_atomsize = (0.5*(sigma1**6 + sigma2**6))**(1/6)
    else:
        raise Exception(f"ERROR mixing_rule = {mixing_rule} is not supported. Currently supported rules: 'geometric' or 'arithmetic' or 'sixthpower' or 'geometric-min' or 'arithmetic-min' or 'sixthpower-min'")
    if 'min' in mixing_rule:
        # both 9-6 and 12-6 LJ have minimum energy at 2**(1/6) or 1.12 of sigma as the sixth power is the respulive portion of the LJ-potential
        mixed_atomsize = (2**(1/6))*mixed_atomsize
    return mixed_atomsize


############################################################
# Function to check if atom overlaps any in current system #
############################################################
def check_for_overlap_and_inside_box(sys, m, linked_lst, domain, domain_graph, xshift, yshift, zshift, phi, theta, psi, tolerance, mix_sigma, mixing_rule, boundary_conditions, scaled_images):
    # Set default overlap, inside Boolean, and insert_molecule
    overlap = False; inside_box = True; insert_molecule = True
    
    # Set mixed and half atom size from tolerance and update later if mix_sigma Boolean Pair Coeffs
    mixed_atomsize = tolerance; half_atomsize = tolerance/2
    
    # Find rotational matrix
    RzRy = misc_functions.matrix_by_matrix(misc_functions.Rz(psi), misc_functions.Ry(theta))
    RzRyRx = misc_functions.matrix_by_matrix(RzRy, misc_functions.Rx(phi))

    # Check if molecule is inside the box
    for id1 in m.atoms:
        atom1 = m.atoms[id1]
        if mix_sigma:
            pair_coeff1 = m.pair_coeffs[atom1.type].coeffs
            half_atomsize = pair_coeff1[tolerance]/2
        x1, y1, z1 = misc_functions.vector_by_matrix(RzRyRx, [atom1.x, atom1.y, atom1.z])
        x1 += xshift
        y1 += yshift
        z1 += zshift
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
        guess = 1 # domainID guess
        for id1 in m.atoms:
            if overlap: break
            atom1 = m.atoms[id1]
            if mix_sigma:
                pair_coeff1 = m.pair_coeffs[atom1.type].coeffs
                sigma1 = pair_coeff1[tolerance]
                half_atomsize = pair_coeff1[tolerance]/2
            x1, y1, z1 = misc_functions.vector_by_matrix(RzRyRx, [atom1.x, atom1.y, atom1.z])
            x1 += xshift
            y1 += yshift
            z1 += zshift
            if not inside_box:
                if x1 <= sys.xlo: x1 += sys.lx
                if x1 >= sys.xhi: x1 -= sys.lx
                if y1 <= sys.ylo: y1 += sys.ly
                if y1 >= sys.yhi: y1 -= sys.ly
                if z1 <= sys.zlo: z1 += sys.lz
                if z1 >= sys.zhi: z1 -= sys.lz
            edgeflag = check_near_edge(x1, y1, z1, 1.1*half_atomsize, sys.xlo, sys.xhi, sys.ylo, sys.yhi, sys.zlo, sys.zhi)
            if edgeflag: periodic_postions = find_periodic_postions(scaled_images, x1, y1, z1, sys.cx, sys.cy, sys.cz, Npos=12)
            domainID, guess = assign_atom_a_domainID(x1, y1, z1, guess, domain)
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
