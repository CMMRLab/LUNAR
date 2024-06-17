# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
June 4th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import math
import time


###############################
# Function to write .nta file #
###############################
def write_nta(m, style, version, charges, filename):
    with open(filename, 'w') as f:
        datafile = filename[:filename.rfind('.')]
        f.write(f'New type assignment file for {datafile}.data created with sheet_builder {version}\n')
        
        # Write atom types
        if style == 'type':
            f.write('\n\n# Assign atom types based on atomTypeID\n')
            f.write('style type\n')
            for i in m.masses:
                coeff = m.masses[i]
                f.write('{:<6} {}\n'.format(i, coeff.type))
        if style == 'id':
            f.write('\n\n# Assign atom types based on atomID\n')
            f.write('style id\n')
            for i in m.atoms:
                atom = m.atoms[i]
                f.write('{:<6} {}\n'.format(i, atom.comment))
                
        # Write charge options
        f.write('\n\n# Assign charge based on charge nta style\n')
        f.write('charge nta\n')
        for i in m.masses:
            atomtype = m.masses[i].type
            if atomtype in charges:
                charge = charges[atomtype]
            else:
                charge = 0.0
            f.write('{:<6} {:^10.6f}\n'.format(atomtype, charge))
            
        # If any types are known IFF pi-electrons add constraints to remove unneeded dihdrals
        pi_electron_types = {'cge':'cg1', 'nbe':'nbn', 'bbe':'bbn'} # {pi-electron : bonded type}
        pi_electrons = set() 
        for i in m.masses:
            atomtype = m.masses[i].type
            if atomtype in pi_electron_types: 
                pi_electrons.add(atomtype)
        if pi_electrons:
            f.write('\n\n# Remove unneeded dihedral and improper interactions\n')
            f.write('#dihedral-remove zero  # can uncomment to remove any parms and interactions if parms are all zeros\n')
            f.write('#improper-remove zero  # can uncomment to remove any parms and interactions if parms are all zeros\n')
            
            # Write specfic dihedrals to remove
            pi_electrons = sorted(pi_electrons)
            f.write('dihedral-remove nta\n')
            for i in pi_electrons:
                f.write('{:<6} {:<6} {:<6} {:<6}\n'.format(i, '*', '*', '*'))
                
            # Write specfic impropers to remove
            f.write('improper-remove nta\n')
            for i in pi_electrons:
                central = pi_electron_types[i]
                f.write('{:<6} {:<6} {:<6} {:<6} {:<6}\n'.format('*', central, '*', '*', 'nb!=3'))
                
    return


###########################################
# Function to assign charge to the system #
###########################################
def charge_system(atoms, charges, log):
    log.out('\n\nAttempting to add charge to system based on charges in charge dictionary ...')
    system_charges = []
    for i in atoms:
        atom = atoms[i]
        if atom.type in charges:
            atom.charge = charges[atom.type]
        system_charges.append(atom.charge)  
    net_system_charge = sum(system_charges)/len(system_charges) 
    log.out('  System net charge: {:.6f}'.format(net_system_charge))
    return atoms


##################################
# Function to get sign of number #
##################################
def number_sign(number):
    if number < 0:
        return -1
    else: return 1


##########################
# Function to wrap atoms #
##########################
def wrap_atoms(atoms, box):
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']
    for i in atoms:
        atom = atoms[i]
        
        # Wrap x atom
        if atom.x < box['xlo']:
            atom.x += lx
            atom.ix -= 1
        elif atom.x > box['xhi']:
            atom.x -= lx
            atom.ix += 1
        
        # Wrap y atom
        if atom.y < box['ylo']:
            atom.y += ly
            atom.iy -= 1
        elif atom.y > box['yhi']:
            atom.y -= ly
            atom.iy += 1
            
        # Wrap z atom
        if atom.z < box['zlo']:
            atom.z += lz
            atom.iz -= 1
        elif atom.z > box['zhi']:
            atom.z -= lz
            atom.iz += 1    
    return atoms


############################
# Function to unwrap atoms #
############################
def unwrap_atoms(atoms, box):
    # Find box info
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']
    cx = (box['xhi'] + box['xlo'])/2
    cy = (box['yhi'] + box['ylo'])/2
    cz = (box['zhi'] + box['zlo'])/2
    
    # set max_x, max_y, max_z w/ minimum image convention
    max_x = lx/2; max_y = ly/2; max_z = lz/2;
    for i in atoms:
        atom = atoms[i]
        
        # Find vectors in each direction
        diffx = cx - atom.x
        diffy = cy - atom.y
        diffz = cz - atom.z
        
        # Apply minimum image convention
        if abs(diffx) > max_x:
            atom.x += number_sign(diffx)*lx
        if abs(diffy) > max_y:
            atom.y += number_sign(diffy)*ly
        if abs(diffz) > max_z:
            atom.z += number_sign(diffz)*lz  
    return atoms


###############################
# Function to replicate atoms #
###############################
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz .comment .atomtype
def replicate_atoms(atoms, box, nx, ny, nz, center, increment_molid=True):    
    # unwrap atoms by anchoring
    unwrapped = unwrap_atoms(atoms, box)
    
    # Find box info
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']
    
    # build replicated atoms
    replicated = {} # { atomID : Atoms Object }
    spans = {'x':[], 'y':[], 'z':[]}
    ID = 0; molID = 0;
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                molID += 1
                for l in unwrapped:
                    atom = unwrapped[l]
                    newx = atom.x + i*lx
                    newy = atom.y + j*ly
                    newz = atom.z + k*lz
                    spans['x'].append(newx)
                    spans['y'].append(newy)
                    spans['z'].append(newz)
                    if increment_molid: molid = molID
                    else: molid = 1
                    ID += 1
                    a = Atom()
                    a.molid = molid
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
                    replicated[ID] = a
            
    # rebuild box
    new_box = {}
    new_box['xlo'] = box['xlo']
    new_box['xhi'] = box['xhi'] + (nx-1)*lx
    new_box['ylo'] = box['ylo']
    new_box['yhi'] = box['yhi'] + (ny-1)*ly 
    new_box['zlo'] = box['zlo']
    new_box['zhi'] = box['zhi'] + (nz-1)*lz
    
    # if center, center replicated atoms
    if center:
        # Shift atoms
        xavg = sum(spans['x'])/len(spans['x'])
        yavg = sum(spans['y'])/len(spans['y'])
        zavg = sum(spans['z'])/len(spans['z'])
        for i in replicated:
            atom = replicated[i]
            atom.x -= xavg
            atom.y -= yavg
            atom.z -= zavg
        
        # Shift box
        new_box['xlo'] -= xavg
        new_box['xhi'] -= xavg
        new_box['ylo'] -= yavg
        new_box['yhi'] -= yavg
        new_box['zlo'] -= zavg
        new_box['zhi'] -= zavg
    
    # rewrap atoms
    replicated = wrap_atoms(replicated, new_box)
    return replicated, new_box 


#####################
# Rotational matrix #
#####################
cos = math.cos
sin = math.sin

def Rx(phi):
    phi = math.radians(phi)
    return [[ 1, 0,        0        ],
            [ 0, cos(phi), -sin(phi)],
            [ 0, sin(phi), cos(phi)]]
  
def Ry(theta):
    theta = math.radians(theta)
    return [[ cos(theta), 0, sin(theta)],
            [ 0,          1, 0         ],
            [-sin(theta), 0, cos(theta)]]
  
def Rz(psi):
    psi = math.radians(psi)
    return [[ cos(psi),  -sin(psi), 0 ],
            [ sin(psi),  cos(psi) , 0 ],
            [ 0,         0,         1 ]]

# Function to multiply two matrices together
def matrix_by_matrix(m1, m2):
    # Intialize matrix
    ncolumns = max( [len(m1[0]), len(m2[0])] ) 
    nrows = max( [len(m1), len(m2)] )
    result = [ncolumns*[0] for row in range(nrows)]
    
    # Iterate through rows of m1
    for i in range(len(m1)):
       # Iterate through columns of m2
       for j in range(len(m2[0])):
           # Iterate through rows of m2
           for k in range(len(m2)):
               result[i][j] += m1[i][k] * m2[k][j]
    return result

# Function to multiply a vector and a matrix together
def vector_by_matrix(m1, v1):
    return [sum([x*y for x, y in zip(v1, v2)]) for v2 in m1]


#############################################################
# Function to rotate a molecule at centered about (0, 0, 0) #
#############################################################
def rotate_molecule(atoms, phi, theta, psi):
    # Find rotational matrix
    RzRy = matrix_by_matrix(Rz(psi), Ry(theta))
    RzRyRx = matrix_by_matrix(RzRy, Rx(phi))
    for i in atoms:
        atom = atoms[i]
        
        # Rotate around atoms
        x, y, z = vector_by_matrix(RzRyRx, [atom.x, atom.y, atom.z])
        
        # Save data into m.atoms
        atom.x = x
        atom.y = y
        atom.z = z
    return atoms


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


###################################################################################
# Function to compute distance when math.dist is not available. math.dist will be #
# quicker but for some reason not all python has math.dist .... compute_distance  #
# will be time optimized but still will be slower then math.dist.                 #
###################################################################################     
def compute_distance(x1, y1, z1, x2, y2, z2):
    dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
    return math.sqrt(dx*dx + dy*dy + dz*dz)


###########################################################
# Function to find N-number of possible periodic postions #
###########################################################
def find_periodic_postions(scaled_images, x1, y1, z1, cx, cy, cz, Npos=12):
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


##########################
# Function to find bonds #
##########################
def find_bonds(atoms, box, boundary, r0, tolerance, max_bonds_per_atom, domain_size, pflag, log):
    # Generate images
    start_time = time.time()
    images, pflags = generate_iflags(boundary, log)
    if pflag: log.out(f'  Boundary flags {pflags[0]} {pflags[1]} {pflags[2]} for finding bonds')
    
    # Find domain decomposition regions
    xlo = box['xlo']
    xhi = box['xhi']
    ylo = box['ylo']
    yhi = box['yhi']
    zlo = box['zlo']
    zhi = box['zhi']
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo
    cx = (xhi + xlo)/2
    cy = (yhi + ylo)/2
    cz = (zhi + zlo)/2
    nxx = math.ceil(lx/domain_size)
    nyy = math.ceil(ly/domain_size)
    nzz = math.ceil(lz/domain_size)
    if nxx == 0: nxx = 1
    if nyy == 0: nyy = 1
    if nzz == 0: nzz = 1
    dx = lx/nxx; dy = ly/nyy; dz = lz/nzz;
    halfdx = dx/2; halfdy = dy/2; halfdz = dz/2; ID = 0;
    xadd = halfdx + xlo; yadd = halfdy + ylo; zadd = halfdz + zlo;
    if pflag: log.out('Using {} x {} x {} sub domains of size {:.2f} x {:.2f} x {:.2f} to perform domain decomposition'.format(nxx, nyy, nzz, dx, dy, dz))
    ndomain = nxx*nyy*nzz; progress_increment = 10; count = 0;
    domain = {} # { domainID : (xlo, xhi, ylo, yhi, zlo, zhi, r, edgeflag) }
    indexes_forward = {} #{ domainID : (nx, ny, nz) }
    indexes_reverse = {} #{ (nx, ny, nz) : domainID }
    cubes = [dx, dy, dz]
    for nx in range(nxx):
        for ny in range(nyy):
            for nz in range(nzz):
                # # Optional printing of progress
                # if pflag:
                #     count += 1
                #     if 100*count/ndomain % progress_increment == 0:
                #         log.out('    progress: {} %'.format(int(100*count/ndomain)))
                
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
                r = compute_distance(xc, yc, zc, cx, cy, cz)
                edgeflag = check_near_edge(xc, yc, zc, domain_size, xlo, xhi, ylo, yhi, zlo, zhi)
                domain[ID] = (xlo_sub, xhi_sub, ylo_sub, yhi_sub, zlo_sub, zhi_sub, xc, yc, zc, r, edgeflag)
                indexes_forward[ID] = (nx+1, ny+1, nz+1)
                indexes_reverse[(nx+1, ny+1, nz+1)] = ID

                
    # Find domain connectivity (graph)
    def get_neighboring_indexes(ni, nii):
        neighs = [ni]
        if ni-1 < 1:
            neighs.append(nii)
        else: neighs.append(ni-1)
        if ni+1 > nii: 
            neighs.append(1)
        else: neighs.append(ni+1)
        return sorted(neighs)
    scaled_images = [(ix*lx, iy*ly, iz*lz) for (ix, iy, iz) in images]
    domain_graph = {ID:set() for ID in domain} # { ID : set(bonded IDs) }
    domain_graph[0] = {ID for ID in domain}
    if pflag: log.out('  Finding cell linked graph for interatomic distance calculations ...')
    progress_increment = 10; count = 0;
    for id1 in domain:
        # # Optional printing of progress
        # if pflag:
        #     count += 1
        #     if 100*count/ndomain % progress_increment == 0:
        #         log.out('    progress: {} %'.format(int(100*count/ndomain)))
                
        d1 = domain[id1]
        x1 = d1[6]; y1 = d1[7]; z1 = d1[8]; r1 = d1[9]; edge1 = d1[10];
        min_radius = r1 - 2*domain_size
        max_radius = r1 + 2*domain_size
        if edge1: periodic_postions = find_periodic_postions(scaled_images, x1, y1, z1, cx, cy, cz, Npos=12)
        nx, ny, nz = indexes_forward[id1]
        nxs = get_neighboring_indexes(nx, nxx)
        nys = get_neighboring_indexes(ny, nyy)
        nzs = get_neighboring_indexes(nz, nzz)
        for ix in nxs:
            for iy in nys:
                for iz in nzs:
                    id2 = indexes_reverse[(ix, iy, iz)]
                    if id1 == id2: continue
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
                            
                        
    # Build linked list
    linked_lst = {i:set() for i in domain} # { domainID : atoms in domain }
    linked_lst[0] = set()
    atom_domain = {} # { atomID : domainID }
    edgeflags = {} # { atomID : edgeflag }
    failed_domain = set(); inner_atoms = set()
    
    if pflag: log.out('  Assigning atoms to each sub domain for interatomic distance calculations ...')
    progress_increment = 10; count = 0; natoms = len(atoms)
    for i in atoms:
        inner_atoms.add(i)
        atom = atoms[i]
        x = atom.x; y = atom.y; z = atom.z
        atom_domain[i] = 0; count += 1
        edgeflags[i] = check_near_edge(x, y, z, domain_size, xlo, xhi, ylo, yhi, zlo, zhi)
        if 100*count/natoms % progress_increment == 0:
            if pflag: log.out('    progress: {} %'.format(int(100*count/natoms)))
            
        # Assign to domain
        try:
            nx = math.ceil( (x-xlo)/cubes[0] )
            ny = math.ceil( (y-ylo)/cubes[0] )
            nz = math.ceil( (z-zlo)/cubes[0] )
            domainID = indexes_reverse[(nx, ny, nz)]
        except: domainID = 0
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
    if pflag: log.out('  Finding interatomic distances for bond creation ...')
    progress_increment = 10; count = 0; natoms = len(atoms); nfailed = 0
    maxdist = r0 + tolerance
    for id1 in atoms:      
        inner_atoms.remove(id1)
        atom1 = atoms[id1]; count += 1
        x1 = atom1.x; y1 = atom1.y; z1 = atom1.z
        if edgeflags[id1]: periodic_postions = find_periodic_postions(scaled_images, x1, y1, z1, cx, cy, cz, Npos=12)
        if 100*count/natoms % progress_increment == 0:
            if pflag: log.out('    progress: {} %'.format(int(100*count/natoms)))
        domain1 = atom_domain[id1]
        atom_domains = [domain1] + list(domain_graph[domain1]) + list(failed_domain)
        if not atom_domains: continue
        for domainID in atom_domains:
            try: id2_atoms = linked_lst[domainID]
            except: 
                id2_atoms = inner_atoms
                nfailed += 1
            for id2 in id2_atoms:
                atom2 = atoms[id2]
                if id1 == id2: continue
                if atom1.molid != atom2.molid: continue
                x2 = atom2.x; y2 = atom2.y; z2 = atom2.z
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
                            break
                else: # else compute none-periodic distance
                    if abs(x1 - x2) > maxdist: continue
                    elif abs(y1 - y2) > maxdist: continue
                    elif abs(z1 - z2) > maxdist: continue
                    distance_fff = compute_distance(x1, y1, z1, x2, y2, z2)
                    if distance_fff <= maxdist:
                        if id1 < id2: bond = (id1, id2)
                        else: bond = (id2, id1) 
                        possible_bonds[bond] = distance_fff
                          
    # Generate graph to remove any bonds that create atoms having more then 3-bonds
    graph = {i:[] for i in atoms}
    for id1, id2 in possible_bonds:     
        graph[id1].append(id2)
        graph[id2].append(id1)
    
    # if more bonds then 3 and found remove longest bonds
    bonds = set([]); flagged_bonds = set([])
    for id1 in graph:           
        # Creating dictionary to sort in decending order
        dist_dict = {}
        for id2 in graph[id1]:
            bond = tuple( sorted([id1, id2]) ) 
            dist_dict[bond] = possible_bonds[bond]
        dist_dict = dict( sorted(dist_dict.items(), key=lambda x:x[1]) )    
       
        # Creating bonds based on inputs values
        count_nb = 0;
        for bond in dist_dict:
            count_nb += 1
            # If number of bonds is less than specified
            if count_nb <= max_bonds_per_atom:
                bonds.add(bond)
            # If bond does not meet criteria flag the bond to not create the bond
            else:
                flagged_bonds.add(bond)
                           
    # Removing any flagged bonds if they were created      
    bonds = list(set(bonds) - set(flagged_bonds)) 
   
    # Sorting bonds and flagged bonds                
    bonds = sorted(bonds); flagged_bonds = sorted(flagged_bonds) 
    #bonds = sorted(list(possible_bonds.keys())) # comment/uncomment to override bond reduction scheme above
                        
    execution_time = (time.time() - start_time)
    if pflag: log.debug(f'    Number of failed domain access attempts: {nfailed}')
    if pflag: log.out(f'    Bond generation execution time: {execution_time} (seconds)')
    return bonds