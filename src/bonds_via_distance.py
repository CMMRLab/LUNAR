# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
December 12th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import time
import math


##################################
# vdw radius to search for bonds #
##################################
# Set elemental vdw radius for bond search (User can adjust this as needed)
# Ref1 = Batsanov, Stepan S. "Van der Waals radii of elements." Inorganic materials 37.9 (2001): 871-885.
# Ref2 = https://en.wikipedia.org/wiki/Van_der_Waals_radius
# The goal is to set via Ref1, but if not found in Ref1 use wikipedia Ref2
vdw_radius = {'C':  1.70, # Ref1
              'H':  1.20, # Ref1
              'O':  1.55, # Ref1
              'N':  1.60, # Ref1
              'S':  1.80, # Ref1
              'F':  1.50, # Ref1
              'Si': 2.10, # Ref1
              'Xe': 2.16, # Ref2
              'Ne': 1.54, # Ref2
              'Kr': 2.02, # Ref2
              'He': 1.40, # Ref2
              'D':  2.40, # Ref1 (2*H vdw radius)
              'Cl': 1.80, # Ref1
              'Ca': 2.40, # Ref1
              'Br': 1.90, # Ref1
              'Ar': 1.88, # Ref2
              'P':  1.90, # Ref1
              'Al': 2.10, # Ref1
              'Mg': 2.20, # Ref1
              'Li': 2.20, # Ref1
              'Fe': 2.05, # Ref1
              'Na': 2.40, # Ref1
              'K':  2.80, # Ref1
              'Cs': 3.00, # Ref1
              'Ba': 2.70, # Ref1
              'Sr': 2.55, # Ref1
              'Pb': 2.30, # Ref1
              }

# Function to get vdw radii
def get_vdw_radii(element, vdw_radius, log):
    if element in vdw_radius:
        return vdw_radius[element]
    else: 
        log.out(f'\nERROR {element} element not in vdw_radius dictionary')
        log.error(f'inside bonds_via_distance.py PLEASE ADD {element}')
        
        
################################################################
# Statistic's functions to use for analyzing bond stats. *NOTE #
# not using numpy to make this code have zero dependancies*    #
################################################################
def compute_mean(data):
  return sum(data)/len(data)
 
def compute_variance(data):
  mean = compute_mean(data)
  deviations = [(x - mean)**2 for x in data]
  variance = sum(deviations)/len(data)
  return variance
 
def compute_standard_deviation(data):
  variance = compute_variance(data)
  return math.sqrt(variance)


###################################################################################
# Function to compute distance when math.dist is not available. math.dist will be #
# quicker but for some reason not all python has math.dist .... compute_distance  #
# will be time optimized but still will be slower then math.dist.                 #
###################################################################################        
def compute_distance(x1, y1, z1, x2, y2, z2):
    dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
    return math.sqrt(dx*dx + dy*dy + dz*dz)

        
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


###########################################################
# Function to find N-number of possible periodic postions #
###########################################################
def find_periodic_postions(scaled_images, x1, y1, z1, cx, cy, cz, Npos=15):
    postions_distance = {} # {distance from center : (pbc-x, pbc-y, pbc-z) }
    for ixlx, iyly, izlz in scaled_images:
        x1i = x1+ixlx; y1i = y1+iyly; z1i = z1+izlz;
        dist_from_center = compute_distance(x1i, y1i, z1i, cx, cy, cz)
        postions_distance[dist_from_center] = (x1i, y1i, z1i)
    postions_distance = dict(sorted(postions_distance.items(), key=lambda x:abs(x[0]) )) # [0=keys;1=values]
    positions = []
    for N, dist in enumerate(postions_distance):
        if N < Npos: positions.append(postions_distance[dist])
        else: break
    return positions


############################################################################################
# Function to find N-number of closets neighboring domains to an atomic position (x, y, z) #
############################################################################################
def find_closets_domains(domain, domain_graph, domainID, x, y, z, Npos=15):
    postions_distance = {} # { domainID : distance from x,y,z }
    for ID in domain_graph[domainID]:
        xc, yc, zc = domain[ID]
        postions_distance[ID] = compute_distance(x, y, z, xc, yc, zc)
    
    # Find the nearest Npos
    postions_distance = dict(sorted(postions_distance.items(), key=lambda x:x[1] )) # [0=keys;1=values]
    domains = set()
    for N, ID in enumerate(postions_distance):
        if N < Npos: domains.add(ID)
        else: break
    
    # Add the current domain to domains
    domains.add(domainID)
    return domains


################################################
# cell linked list interatomic distance search #
################################################
def interatomic_cell_linked(m, maxdist_dict, edgeflags, scaled_images, xlo, xhi, ylo, yhi, zlo, zhi, log):
    # Set domain size based on largest pairwise distances
    domain_size = 2.1*(maxdist_dict[max(maxdist_dict, key=maxdist_dict.get)])
    
    # Find domain decomposition regions
    start_time = time.time()
    lx = xhi-xlo; ly = yhi-ylo; lz = zhi-zlo;
    cx = (xhi + xlo)/2; cy = (yhi + ylo)/2; cz = (zhi + zlo)/2;
    nxx = math.ceil(lx/domain_size)
    nyy = math.ceil(ly/domain_size)
    nzz = math.ceil(lz/domain_size)
    if nxx == 0: nxx = 1
    if nyy == 0: nyy = 1
    if nzz == 0: nzz = 1
    dx = lx/nxx; dy = ly/nyy; dz = lz/nzz; ID = 0;
    xadd = dx/2 + xlo; yadd = dy/2 + ylo; zadd = dz/2 + zlo;
    log.out('\n\nFinding domain to use cell linked list algorithm for interatomic distance calculations ...')
    domain = {} # { domainID : (xc, yc, zc) }
    indexes_forward = {} # { domainID : (nx, ny, nz) }
    indexes_reverse = {} # { (nx, ny, nz) : domainID }
    for nx in range(nxx):
        for ny in range(nyy):
            for nz in range(nzz):
                ID += 1
                xc = nx*dx + xadd
                yc = ny*dy + yadd
                zc = nz*dz + zadd
                domain[ID] = (xc, yc, zc)
                indexes_forward[ID] = (nx+1, ny+1, nz+1)
                indexes_reverse[(nx+1, ny+1, nz+1)] = ID
                
    # Find domain connectivity (graph)
    def get_neighboring_indexes(ni, nii):
        neighs = [ni]
        if ni-1 < 1: # lo-side wraps to hi-side for periodicity
            neighs.append(nii)
        else: neighs.append(ni-1)
        if ni+1 > nii: # hi-side wraps to lo-side for periodicity
            neighs.append(1)
        else: neighs.append(ni+1)
        return neighs
    domain_graph = {ID:set() for ID in domain} # { ID : set(bonded IDs) }
    domain_graph[0] = {ID for ID in domain}
    log.out('Finding cell linked graph for interatomic distance calculations ...')
    for id1 in domain:    
        nx, ny, nz = indexes_forward[id1]
        nxs = get_neighboring_indexes(nx, nxx)
        nys = get_neighboring_indexes(ny, nyy)
        nzs = get_neighboring_indexes(nz, nzz)
        for ix in nxs:
            for iy in nys:
                for iz in nzs:
                    id2 = indexes_reverse[(ix, iy, iz)]
                    if id1 == id2: continue
                    domain_graph[id1].add(id2)
                    domain_graph[id2].add(id1)
                        
    # Build linked list
    linked_lst = {i:set() for i in domain} # { domainID : atoms in domain }
    linked_lst[0] = set()
    atom_domain = {} # { atomID : domainID }
    edgeflags = {} # { atomID : edgeflag }
    log.out('  Assigning atoms to each sub domain for interatomic distance calculations ...')
    for i in m.atoms:
        atom = m.atoms[i]
        x = atom.x; y = atom.y; z = atom.z
        edgeflags[i] = check_near_edge(x, y, z, domain_size, xlo, xhi, ylo, yhi, zlo, zhi)
            
        # Assign to domain
        try:
            nx = math.ceil( (x-xlo)/dx )
            ny = math.ceil( (y-ylo)/dy )
            nz = math.ceil( (z-zlo)/dz )
            domainID = indexes_reverse[(nx, ny, nz)]
        except: domainID = 0
        atom_domain[i] = domainID
        linked_lst[domainID].add(i)
        
        # check if atom is outside of box and assign to every domain if it is outside of the box
        outside = False
        if atom.x <= xlo: outside = True
        if atom.x >= xhi: outside = True
        if atom.y <= ylo: outside = True
        if atom.y >= yhi: outside = True
        if atom.z <= zlo: outside = True
        if atom.z >= zhi: outside = True
        if outside:
            for j in linked_lst:
                linked_lst[j].add(i)
                    
    # Start finding interatomic distances
    possible_bonds = {} # { tuple(id1, id2): distance_fff or distance_ppp }
    bond_status = {'periodic': 0, 'non-periodic': 0} # to tally bond types
    log.out('Finding interatomic distances for bond creation (cell linked list) ....')
    progress_increment = 10; count = 0; natoms = len(m.atoms);
    for id1 in m.atoms:      
        atom1 = m.atoms[id1]; count += 1
        x1 = atom1.x; y1 = atom1.y; z1 = atom1.z
        if edgeflags[id1]: periodic_postions = find_periodic_postions(scaled_images, x1, y1, z1, cx, cy, cz, Npos=15)
        if 100*count/natoms % progress_increment == 0:
            log.out('    progress: {} %'.format(int(100*count/natoms)))
        atom_domains = find_closets_domains(domain, domain_graph, atom_domain[id1], x1, y1, z1, Npos=27)
        if not atom_domains: continue
        for domainID in atom_domains:
            for id2 in linked_lst[domainID]:
                atom2 = m.atoms[id2]
                if id1 == id2: continue
                # Find atom2 and max distance cutoff
                atom2 = m.atoms[id2];
                maxdist = maxdist_dict[(atom1.element, atom2.element)];
                x2 = atom2.x; y2 = atom2.y; z2 = atom2.z

                # If both atoms are near an edge try finding across the PBC  
                if edgeflags[id1] and edgeflags[id2]:   
                    for x1i, y1i, z1i in periodic_postions:
                        if abs(x1i - x2) > maxdist: continue
                        elif abs(y1i - y2) > maxdist: continue
                        elif abs(z1i - z2) > maxdist: continue
                        distance_ppp = compute_distance(x1i, y1i, z1i, x2, y2, z2)
                        if distance_ppp <= maxdist: 
                            if id1 < id2: bond = (id1, id2)
                            else: bond = (id2, id1)
                            possible_bonds[bond] = distance_ppp
                            bond_status['periodic'] += 1
                            break
                            
                else: # else compute none-periodic distance
                    if abs(x1 - x2) > maxdist: continue
                    elif abs(y1 - y2) > maxdist: continue
                    elif abs(z1 - z2) > maxdist: continue
                    distance_fff = compute_distance(x1, y1, z1, x2, y2, z2)
                    if distance_fff <= maxdist:
                        # Generate bonding pairs in ascending order
                        if id1 < id2: bond = (id1, id2)
                        else: bond = (id2, id1)
                        
                        possible_bonds[bond] = distance_fff
                        bond_status['non-periodic'] += 1
                        
    execution_time = (time.time() - start_time)
    log.out(f'    Bond generation execution time: {execution_time} (seconds)')
    return possible_bonds, bond_status 
    
    
##############################################################
# Class to generate bonds via interatomic distance searching #
##############################################################
class Stats: pass # .count .avg .min .max .std .cutoff
class generate:
    def __init__(self, m, boundary, vdw_radius_scale, maxbonded, log):
        self.bonds = [] # lst of viable bonds
        self.flagged_bonds = [] # lst of flagged bonds
        self.statistics = {} # { bondtype : Stats object}
        self.nb_stats = {} # { element : Element_nb object }
        self.images = [] # lst of images searched
        self.nb_count = {} # { element : nb_dict }
    
    
        ################################################
        # Generate images via minimum image convention #
        ################################################
        # Boundary conditions
        pflags = boundary.split() # split boundary
        count = pflags.count('f') + pflags.count('p')
        if len(pflags) != 3 and count != 3:
            log.error('ERROR boundary does not contain 3-dimensions or boundary flags are not f or p')
        if 'p' in pflags and m.xy != 0 or m.xz != 0 or m.yz != 0:
            log.out('ERROR simulaiton cell is NOT orthogonal and a periodic boundary is trying')
            log.error('to be applied. Currently PBCs are limited to orthogonal boxes only ...')
        
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
                    
        # Sort images in ascending order of magnitude to reduce run time when 
        # finding bond distances, since most bonds will be found at (0, 0, 0)
        self.images = sorted(images, key=lambda x: sum([abs(i) for i in x]) )
    
    
        ####################################
        # Find set simulataion cell bounds #
        ####################################
        xline = m.xbox_line.split(); yline = m.ybox_line.split(); zline = m.zbox_line.split();
        xlo = float(xline[0]); xhi = float(xline[1]);
        ylo = float(yline[0]); yhi = float(yline[1]);
        zlo = float(zline[0]); zhi = float(zline[1]);
        lx = xhi-xlo; ly = yhi-ylo; lz = zhi-zlo;
        
        # Build scaled images for quicker calculations later on (Li directions already multiplied)
        scaled_images = [(ix*lx, iy*ly, iz*lz) for (ix, iy, iz) in self.images]            
        
        # Determine "nearness" to edge. Only search periodic bonds if atom is this close to an edge.
        # Only search largest possible bond found in the system and nothing more to manage run time
        elements = sorted({m.atoms[i].element for i in m.atoms})
        max_from_edge = vdw_radius_scale*max([get_vdw_radii(element, vdw_radius, log) for element in elements]) 
        edgeflags = {} # { atomID : edge flag }
        for i in m.atoms:
            x = m.atoms[i].x; y = m.atoms[i].y; z = m.atoms[i].z;
            edgeflags[i] = check_near_edge(x, y, z, max_from_edge, xlo, xhi, ylo, yhi, zlo, zhi)


        ###############################################################################################
        # Generate lookup tables for distance cutoffs based on bondtypes for quicker runtime later on #
        ###############################################################################################
        maxdist_dict = {} # { tuple(element1, element2) : cutoff distance }
        for element1 in elements:
            radius1 = vdw_radius_scale*get_vdw_radii(element1, vdw_radius, log)
            for element2 in elements:
                radius2 = vdw_radius_scale*get_vdw_radii(element2, vdw_radius, log)
                
                # Compute distance and log in forward and reverse ordering
                maxdistance = (radius1 + radius2)/2
                maxdist_dict[(element1, element2)] = maxdistance
                maxdist_dict[(element2, element1)] = maxdistance
                
                
        ##########################
        # Compute system density #
        ##########################
        system_mass_sum = 0; amu2grams = 1/6.02214076e+23; density = 0;
        for i in m.atoms:
            mass = m.masses[m.atoms[i].type]
            system_mass_sum += mass.coeffs[0]
        mass = system_mass_sum*amu2grams
        angstromcubed2cmcubed = 1e-24
        volume = lx*ly*lz*angstromcubed2cmcubed
        if volume > 0: density = mass/volume
        log.out('\n---------------------------------------')
        log.out('System Unit cell, volume, mass, density')
        log.out('---------------------------------------')
        log.out('{:<10} {:^10.5f} {:<10}'.format('Lx:', lx, 'angstrom'))
        log.out('{:<10} {:^10.5f} {:<10}'.format('Ly:', ly, 'angstrom'))
        log.out('{:<10} {:^10.5f} {:<10}'.format('Lz:', lz, 'angstrom'))
        log.out('{:<10} {:^10.4E} {:<10}'.format('volume:', volume, 'cm^3'))
        log.out('{:<10} {:^10.4E} {:<10}'.format('mass:', mass, 'grams'))
        log.out('{:<10} {:^10.5f} {:<10}' .format('density:', density, 'g/cm^3'))
        
        
        #################################################
        # Find possible bonds via interatomic distances #
        #################################################
        possible_bonds, self.bond_status = interatomic_cell_linked(m, maxdist_dict, edgeflags, scaled_images, xlo, xhi, ylo, yhi, zlo, zhi, log)
                    
        
        ####################################
        # Find bonded atoms to be able to  #
        # set max nb cut-off per atom type #
        ####################################   
        # Generate graph
        graph = {i:[] for i in m.atoms}
        for id1, id2 in possible_bonds:     
            graph[id1].append(id2)
            graph[id2].append(id1)

       
        ########################################
        # Start finding and flagging bonds     #
        # that meet cut-off values set by user #
        ########################################
        bonds = set([]); flagged_bonds = set([]);
        for id1 in graph:           
            # Creating dictionary to sort in decending order
            dist_dict = {}
            for id2 in graph[id1]:
                bond = tuple( sorted([id1, id2]) ) 
                dist_dict[bond] = possible_bonds[bond]
            dist_dict = dict( sorted(dist_dict.items(), key=lambda x:x[1]) )
           
            # Find max_nb based on element
            max_nb = maxbonded[m.atoms[id1].element]
           
            # Creating bonds based on inputs values
            count_nb = 0;
            for bond in dist_dict:
                count_nb += 1
                # If number of bonds is less than specified
                if count_nb <= max_nb:
                    bonds.add(bond)
                # If bond does not meet criteria flag the bond to not create the bond
                else:
                    flagged_bonds.add(bond)
                               
        # Removing any flagged bonds if they were created      
        bonds = list(set(bonds) - set(flagged_bonds)) 
       
        # Sorting bonds and flagged bonds                
        bonds = sorted(bonds); flagged_bonds = sorted(flagged_bonds) 
        self.bonds = bonds; self.flagged_bonds = flagged_bonds
       
        
        # generate bondtype_stats_holder
        bondtype_stats_holder = {} # { tuple(element1, element2) : [lst of distances]}
        for id1, id2 in bonds:
            distance = possible_bonds[ tuple(sorted([id1, id2])) ]
            element1 = m.atoms[id1].element
            element2 = m.atoms[id2].element
            bondtype = tuple(sorted([element1, element2]))
            
            # Log bondtype
            if bondtype in bondtype_stats_holder:
                bondtype_stats_holder[bondtype].append(distance)
            else:
                bondtype_stats_holder[bondtype] = [distance]
        
        # Find bond order stats
        for i in bondtype_stats_holder:
            lst = bondtype_stats_holder[i]
    
            # Find bond name
            bond = '{}-{}'.format(i[0], i[1])
            
            # find distance cutoff used
            cutoff = round(maxdist_dict[(i[0], i[1])], 4)
            
            # Add to stats if count if larger then zero
            if len(lst) > 0:
                s = Stats()
                s.count = int(len(lst))
                s.avg = '{:.4f}'.format( compute_mean(lst) )
                s.min = '{:.4f}'.format( min(lst) )
                s.max = '{:.4f}'.format( max(lst) )
                s.std = '{:.4f}'.format( compute_standard_deviation(lst) )
                s.cutoff = cutoff
                self.statistics[bond] = s
                
        # Find number of bonded stats
        nb_dict = {i:[] for i in m.atoms}
        for id1, id2 in bonds:
            nb_dict[id1].append(id2)
            nb_dict[id2].append(id2)
            
        # Find unique elements and generate nb_count dict to add to
        elements = sorted({m.atoms[i].element for i in m.atoms})
        self.maxbond = max([len(nb_dict[i]) for i in nb_dict]) # find max bond from nb
        for i in elements:
            self.nb_count[i] = {j:0 for j in range(self.maxbond+1)}
            
        # Add to self.nb_count
        for i in nb_dict:
            element = m.atoms[i].element
            nb = len(nb_dict[i])
            self.nb_count[element][nb] += 1