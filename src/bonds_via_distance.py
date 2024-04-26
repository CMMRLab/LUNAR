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
import sys


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


###########################################
# Brute force interatomic distance search #
###########################################
def interatomic_brute_force(m, maxdist_dict, radial, edgeflags, scaled_images, log):
    # Start finding interatomic distances
    start_time = time.time()
    possible_bonds = {} # { tuple(id1, id2): distance_fff or distance_ppp }
    bond_status = {'periodic': 0, 'non-periodic': 0} # to tally bond types
    id2_atoms = {i for i in m.atoms} # atoms to loop over in inner loop (will be reduced each outer loop iteration)
    log.out('\n\nFinding interatomic distances for bond creation (brute force) ....')
    progress_increment = 10; count = 0; natoms = len(m.atoms);
    for id1 in m.atoms:      
        # Find atom1 to access x1, y1, and z1
        atom1 = m.atoms[id1]; count += 1
        x1 = atom1.x; y1 = atom1.y; z1 = atom1.z; r1 = radial[id1]
        if edgeflags[id1]: periodic_postions = [ (x1+ixlx, y1+iyly, z1+izlz) for ixlx, iyly, izlz in scaled_images ]
        id2_atoms.remove(id1)
        if 100*count/natoms % progress_increment == 0:
            log.out('    progress: {} %'.format(int(100*count/natoms)))
        for id2 in id2_atoms:
            # Find atom2 and max distance cutoff
            atom2 = m.atoms[id2];
            maxdist = maxdist_dict[(atom1.element, atom2.element)];
            min_radius = r1 - 2*maxdist
            max_radius = r1 + 2*maxdist
             
            # Find x2, y2, and z2 for interatomic searching
            x2 = atom2.x; y2 = atom2.y; z2 = atom2.z; r2 = radial[id2]
            if min_radius < r2 < max_radius:
                # If both atoms are near an edge try finding across the PBC
                if edgeflags[id1] and edgeflags[id2]:   
                    for x1i, y1i, z1i in periodic_postions:
                        if abs(x1i - x2) > maxdist: continue
                        elif abs(y1i - y2) > maxdist: continue
                        elif abs(z1i - z2) > maxdist: continue
                        else:
                            distance_ppp = compute_distance(x1i, y1i, z1i, x2, y2, z2)
                            if distance_ppp <= maxdist: 
                                # Generate bonding pairs in ascending order
                                if id1 < id2: bond = (id1, id2)
                                else: bond = (id2, id1)
                                
                                possible_bonds[bond] = distance_ppp
                                bond_status['periodic'] += 1
                                break
                            
                else: # else compute none-periodic distance
                    if abs(x1 - x2) > maxdist: continue
                    elif abs(y1 - y2) > maxdist: continue
                    elif abs(z1 - z2) > maxdist: continue
                    else: 
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

################################################
# cell linked list interatomic distance search #
################################################
def interatomic_cell_linked(m, maxdist_dict, radial, edgeflags, scaled_images, xlo, xhi, ylo, yhi, zlo, zhi, log):
    # Find domain decomposition regions
    start_time = time.time()
    lx = xhi-xlo; ly = yhi-ylo; lz = zhi-zlo;
    cx = (xhi + xlo)/2; cy = (yhi + ylo)/2; cz = (zhi + zlo)/2;
    domain_size = 6; domain = {} # { domainID : (xlo, xhi, ylo, yhi, zlo, zhi, xc, yc, zc, r, edgeflag) }
    nxx = math.ceil(lx/domain_size)
    nyy = math.ceil(ly/domain_size)
    nzz = math.ceil(lz/domain_size)
    dx = lx/nxx; dy = ly/nyy; dz = lz/nzz;
    halfdx = dx/2; halfdy = dy/2; halfdz = dz/2; ID = 0;
    xadd = halfdx + xlo; yadd = halfdy + ylo; zadd = halfdz + zlo;
    log.out('\n\nFinding domain to use cell linked list algorithm for interatomic distance calculations ...')
    for nx in range(nxx):
        for ny in range(nyy):
            for nz in range(nzz):
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
                r = compute_distance(xc, yc, zc, cx, cy, cz)
                edgeflag = check_near_edge(xc, yc, zc, domain_size, xlo_sub, xhi_sub, ylo_sub, yhi_sub, zlo, zhi)
                domain[ID] = (xlo_sub, xhi_sub, ylo_sub, yhi_sub, zlo_sub, zhi_sub, xc, yc, zc, r, edgeflag)
                
    # Find domain connectivity (graph)
    domain_graph = {ID:set() for ID in domain} # { ID : set(bonded IDs) }
    sub_domain = {i for i in domain}
    log.out('Finding cell linked graph for interatomic distance calculations ...')
    for id1 in domain:
        d1 = domain[id1]
        x1 = d1[6]; y1 = d1[7]; z1 = d1[8]; r1 = d1[9]; edge1 = d1[10];
        min_radius = r1 - 2*domain_size
        max_radius = r1 + 2*domain_size
        sub_domain.remove(id1)
        if edge1: periodic_postions = [ (x1+ixlx, y1+iyly, z1+izlz) for ixlx, iyly, izlz in scaled_images ]
        for id2 in sub_domain:
            d2 = domain[id2]
            x2 = d2[6]; y2 = d2[7]; z2 = d2[8]; r2 = d2[9]; edge2 = d2[10];
            if min_radius < r2 < max_radius:
                if edge1 and edge2: # periodic
                    for x1i, y1i, z1i in periodic_postions:
                        if abs(x1i - x2) > domain_size: continue
                        elif abs(y1i - y2) > domain_size: continue
                        elif abs(z1i - z2) > domain_size: continue
                        else:
                            domain_graph[id1].add(id2)
                            domain_graph[id2].add(id1)
                            break
                else: # non-periodic
                    if abs(x1 - x2) > domain_size: continue
                    elif abs(y1 - y2) > domain_size: continue
                    elif abs(z1 - z2) > domain_size: continue
                    else:
                        domain_graph[id1].add(id2)
                        domain_graph[id2].add(id1)
                        
    # Build linked list
    linked_lst = {i:set() for i in domain} # { domainID : atoms in domain }
    atom_domain = {} # { atomID : domainID }
    failed_domain = set()
    for i in m.atoms:
        atom = m.atoms[i]
        x = atom.x; y = atom.y; z = atom.z
        atom_domain[i] = 0
        for j in domain:
            xlo_sub, xhi_sub, ylo_sub, yhi_sub, zlo_sub, zhi_sub, xc, yc, zc, r, edgeflag = domain[j]
            if xlo_sub < x < xhi_sub and ylo_sub < y < yhi_sub and zlo_sub < z < zhi_sub: 
                linked_lst[j].add(i); atom_domain[i] = j; break;
        if atom_domain[i] == 0: failed_domain.add(i)
                    
    # Start finding interatomic distances
    possible_bonds = {} # { tuple(id1, id2): distance_fff or distance_ppp }
    bond_status = {'periodic': 0, 'non-periodic': 0} # to tally bond types
    log.out('Finding interatomic distances for bond creation (cell linked list) ....')
    progress_increment = 10; count = 0; natoms = len(m.atoms);
    checked = {i:False for i in m.atoms}
    for id1 in m.atoms:      
        # Find atom1 to access x1, y1, and z1
        atom1 = m.atoms[id1]; count += 1
        x1 = atom1.x; y1 = atom1.y; z1 = atom1.z; r1 = radial[id1]
        if edgeflags[id1]: periodic_postions = [ (x1+ixlx, y1+iyly, z1+izlz) for ixlx, iyly, izlz in scaled_images ]
        if 100*count/natoms % progress_increment == 0:
            log.out('    progress: {} %'.format(int(100*count/natoms)))
        domain1 = atom_domain[id1]
        atom_domains = [domain1] + list(domain_graph[domain1]) + list(failed_domain)
        checked[id1] = True
        if not atom_domains: continue
        for domainID in atom_domains:
            id2_atoms = linked_lst[domainID]
            for id2 in id2_atoms:
                if checked[id2]: continue
                # Find atom2 and max distance cutoff
                atom2 = m.atoms[id2];
                maxdist = maxdist_dict[(atom1.element, atom2.element)];
                min_radius = r1 - 2*maxdist
                max_radius = r1 + 2*maxdist
                 
                # Find x2, y2, and z2 for interatomic searching
                x2 = atom2.x; y2 = atom2.y; z2 = atom2.z; r2 = radial[id2]
                if min_radius < r2 < max_radius:
                    # If both atoms are near an edge try finding across the PBC
                    if edgeflags[id1] and edgeflags[id2]:   
                        for x1i, y1i, z1i in periodic_postions:
                            if abs(x1i - x2) > maxdist: continue
                            elif abs(y1i - y2) > maxdist: continue
                            elif abs(z1i - z2) > maxdist: continue
                            else:
                                distance_ppp = compute_distance(x1i, y1i, z1i, x2, y2, z2)
                                if distance_ppp <= maxdist: 
                                    # Generate bonding pairs in ascending order
                                    if id1 < id2: bond = (id1, id2)
                                    else: bond = (id2, id1)
                                    
                                    possible_bonds[bond] = distance_ppp
                                    bond_status['periodic'] += 1
                                    break
                                
                    else: # else compute none-periodic distance
                        if abs(x1 - x2) > maxdist: continue
                        elif abs(y1 - y2) > maxdist: continue
                        elif abs(z1 - z2) > maxdist: continue
                        else: 
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
        radial = {} # { atomID : radius from center of cell }
        cx = (xhi + xlo)/2; cy = (yhi + ylo)/2; cz = (zhi + zlo)/2;
        for i in m.atoms:
            x = m.atoms[i].x; y = m.atoms[i].y; z = m.atoms[i].z;
            edgeflags[i] = check_near_edge(x, y, z, max_from_edge, xlo, xhi, ylo, yhi, zlo, zhi)
            radial[i] = compute_distance(x, y, z, cx, cy, cz)


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
        if len(m.atoms) < 7500 or density < 0.25 or lx < 18 or ly < 18 or lz < 18:
            possible_bonds, self.bond_status = interatomic_brute_force(m, maxdist_dict, radial, edgeflags, scaled_images, log)
        else: # cell linked list method is really slow for low density systems and requires at least a 3*domain_size (3*6=18) in each direction to decompose the system
            try: possible_bonds, self.bond_status = interatomic_cell_linked(m, maxdist_dict, radial, edgeflags, scaled_images, xlo, xhi, ylo, yhi, zlo, zhi, log)
            except:
                log.out('Cell linked list method failed, resorting to brute force method')
                possible_bonds, self.bond_status = interatomic_brute_force(m, maxdist_dict, radial, edgeflags, scaled_images, log)


                    
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