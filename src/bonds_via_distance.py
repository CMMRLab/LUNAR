# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.4
May 4, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.periodicity as periodicity
from itertools import product
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
        return 0
        
        
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


################################################
# cell linked list interatomic distance search #
################################################
def interatomic(m, maxdist_dict, pflags, log):
    # Get atoms count and exit if needed
    natoms = len(m.atoms)
    if natoms == 0:
        log.error('\nERROR no atoms to find interatomic distances for bond calculations')
    
    # Set domain size based on largest pairwise distances
    maxcut = max(maxdist_dict.values())
    domain_size = 1.1*maxcut
    
    # Set periodic dimension flags
    xperiodic, yperiodic, zperiodic = True, True, True
    if pflags[0] == 'f': xperiodic = False
    if pflags[1] == 'f': yperiodic = False
    if pflags[2] == 'f': zperiodic = False
    
    # Get box info
    h, h_inv, boxlo, boxhi, tilts = periodicity.get_box_parameters(m)
    lx, ly, lz, yz, xz, xy = h
    nx = max(1, math.ceil(lx / domain_size))
    ny = max(1, math.ceil(ly / domain_size))
    nz = max(1, math.ceil(lz / domain_size))
    
    # Assign atoms to domain
    tol = 1e-12
    domains = {} # {(ix, iy, iz):[id1, id2, id3 ...]}
    atom_to_domain = {} # {atomID:(ix, iy, iz)}
    start_time = time.time()
    log.out('\n\nAssigning atoms to each sub-domain for interatomic distance calculations ...')
    for id1, atom1 in m.atoms.items():
        # Convert to fractional coordinates and wrap/clamp into the primary cell
        fx, fy, fz = periodicity.pos2frac(atom1.x, atom1.y, atom1.z, h_inv, boxlo)
        if xperiodic:
            fx %= 1.0
        else: fx = min(1.0 - tol, max(0.0, fx))
        
        if yperiodic:
            fy %= 1.0
        else: fy = min(1.0 - tol, max(0.0, fy))
        
        if zperiodic:
            fz %= 1.0
        else: fz = min(1.0 - tol, max(0.0, fz))
    
        # Assigned to domain
        ix = min(nx - 1, math.floor(fx * nx)) # 0 ≤ ix ≤ nx-1
        iy = min(ny - 1, math.floor(fy * ny)) # 0 ≤ iy ≤ ny-1
        iz = min(nz - 1, math.floor(fz * nz)) # 0 ≤ iz ≤ nz-1
        domains.setdefault((ix, iy, iz), []).append(id1)
        atom_to_domain[id1] = (ix, iy, iz)
        
    execution_time = (time.time() - start_time)
    log.out(f'    Domain assigmnet execution time: {execution_time} (seconds)')
        
    # Neighbor shifting that works for very skewed triclinic cells
    if yz == 0 and xz == 0 and xy == 0:
        sx, sy, sz = 1, 1, 1
    else:
        ax = lx / nx
        by = math.sqrt((xy / ny)**2 + (ly / ny)**2)
        cz = math.sqrt((xz / nz)**2 + (yz / nz)**2 + (lz / nz)**2)
        sx = math.ceil(maxcut / ax) + 2 # could be "+ 1" but "+ 2" helps for highly skewed cells
        sy = math.ceil(maxcut / by) + 2 # could be "+ 1" but "+ 2" helps for highly skewed cells
        sz = math.ceil(maxcut / cz) + 2 # could be "+ 1" but "+ 2" helps for highly skewed cells
    neighbor_shifts = list(product(range(-sx, sx + 1), range(-sy, sy + 1), range(-sz, sz + 1)))
    
    # Start finding interatomic distances
    possible_bonds = {} # { tuple(id1, id2): distance }
    progress_increment = 10; last_progress = -1; count = 0
    start_time = time.time()
    log.out('Finding interatomic distances for bond creation (cell linked list) ....')
    for id1 in m.atoms:
        count += 1
        percent = int(100 * count / natoms)
        if percent >= last_progress + progress_increment:
            log.out(f'    progress: {percent} %')
            last_progress = percent
        
        # Atom 1 info
        checked_domains = set()
        ix, iy, iz = atom_to_domain[id1]
        atom1 = m.atoms[id1]
        for dix, diy, diz in neighbor_shifts:
            jx = (ix + dix) % nx if xperiodic else ix + dix
            if not (0 <= jx < nx): continue
            
            jy = (iy + diy) % ny if yperiodic else iy + diy
            if not (0 <= jy < ny): continue
            
            jz = (iz + diz) % nz if zperiodic else iz + diz
            if not (0 <= jz < nz): continue
            
            # If nx=ny=nz=1, or for larger skewed periodic stencils, many
            # shifts may point to the same domain (check for performance)
            jdomain = (jx, jy, jz)
            if jdomain in checked_domains:
                continue
            checked_domains.add(jdomain)
        
            # Atom 2 info
            ids2 = domains.get(jdomain, [])
            for id2 in ids2:
                if id2 <= id1: continue
                atom2 = m.atoms[id2]
                
                # See if distance is within vdw-radii maxdist
                distance = periodic_distance(atom1, atom2, h, h_inv, boxlo, xperiodic, yperiodic, zperiodic)
                maxdist = maxdist_dict[(atom1.element, atom2.element)]
                if distance <= maxdist:
                    possible_bonds[(id1, id2)] = distance

    execution_time = (time.time() - start_time)
    log.out(f'    Bond generation execution time: {execution_time} (seconds)')
    return h, possible_bonds


###########################################
# Triclinic periodic distance calculation #
###########################################
def periodic_distance(atom1, atom2, h, h_inv, boxlo, xperiodic, yperiodic, zperiodic):
    fx1, fy1, fz1 = periodicity.pos2frac(atom1.x, atom1.y, atom1.z, h_inv, boxlo)
    fx2, fy2, fz2 = periodicity.pos2frac(atom2.x, atom2.y, atom2.z, h_inv, boxlo)
    dfx = fx2 - fx1
    dfy = fy2 - fy1
    dfz = fz2 - fz1

    # Apply minimum image in fractional space
    if xperiodic:
        dfx -= round(dfx)
    if yperiodic:
        dfy -= round(dfy)
    if zperiodic:
        dfz -= round(dfz)

    dx, dy, dz = periodicity.frac2pos(dfx, dfy, dfz, h, [0.0, 0.0, 0.0])

    return math.sqrt(dx*dx + dy*dy + dz*dz)
    
    
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
    
    
        # Boundary conditions
        pflags = boundary.split() # split boundary
        count = pflags.count('f') + pflags.count('p')
        if len(pflags) != 3 or count != 3:
            log.error('ERROR boundary does not contain 3-dimensions or boundary flags are not f or p')

        # Generate lookup tables for distance cutoffs based on bondtypes for quicker runtime later on
        elements = sorted({m.atoms[i].element for i in m.atoms})
        maxdist_dict = {} # { tuple(element1, element2) : cutoff distance }
        for element1 in elements:
            radius1 = vdw_radius_scale*get_vdw_radii(element1, vdw_radius, log)
            for element2 in elements:
                radius2 = vdw_radius_scale*get_vdw_radii(element2, vdw_radius, log)
                
                # Compute distance and log in forward and reverse ordering
                maxdistance = (radius1 + radius2)/2
                maxdist_dict[(element1, element2)] = maxdistance
                maxdist_dict[(element2, element1)] = maxdistance
        
        
        # Find possible bonds via interatomic distances
        h, possible_bonds = interatomic(m, maxdist_dict, pflags, log)
        
        
        # Compute system density
        lx, ly, lz, yz, xz, xy = h
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
            # Creating dictionary to sort in ascending order
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
            nb_dict[id2].append(id1)
            
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