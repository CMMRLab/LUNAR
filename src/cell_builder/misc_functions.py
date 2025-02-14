# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.9
February 14th, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

##############################
# Import Necessary Libraries #
##############################
import math
import os


#################################################################################################
# Function to get header info from LAMMPS datafile. The following are supported header options: #
#    Regions = {'x':[0, 0], 'y':[0, 0], 'z':[0,0]}   # {'x':[lo-iflag, hi-iflag], ... }         #
#    Rotations = {'x':10, 'y':15, 'z':20}            # {'x':max-rotation, ... }                 #
#################################################################################################
def get_dict(string, log):
    # Setup the globals namespace to limit scope of what eval() can do
    allowed_builtins = ['min','max','sum','abs','len','map','range','reversed']
    copied_builtins = globals()['__builtins__'].copy()
    globals_dict = {'__builtins__': {key:copied_builtins[key] for key in allowed_builtins} }
    
    # Setup default options
    Regions = {} # {'x':[lo-iflag, hi-iflag], 'y':[lo-iflag, hi-iflag], 'z':[lo-iflag, hi-iflag]} 
    Rotations = {} # {'x':max-rotation, 'y':max-rotation, 'z':max-rotation}
    Shift = {} # {'x':x-shift, 'y':y-shift, 'z':z-shift}
    
    
    # Function to parse header
    def get_string_dict(line):
        string_dict = '{'; bracket = False; equal = False;
        count1 = 0; count2 = 0;
        for i in line:
            if i == '=': equal = True; continue
            if i == '{': bracket = True; count1 += 1; continue;
            if i == '}': bracket = False; count2 += 1; continue;
            if count1 > 1 or count2 > 1: bracket = False
            if bracket and equal: string_dict += i
        string_dict += '}'
        return string_dict
    
    # Try getting Shift
    if 'Shift' in string:
        try:
            line = string.split('Shifts') # Find Shift string 
            string_dict = get_string_dict(line[-1])
            Shift = eval(string_dict, globals_dict)
        except: log.error(f'ERROR could not resolve Shifts = {line[-1]}')
    
    # Try getting Region
    if 'Regions' in string:
        try:
            line = string.split('Regions') # Find Region string 
            string_dict = get_string_dict(line[-1])
            Regions = eval(string_dict, globals_dict)
        except: log.error(f'ERROR could not resolve Regions = {line[-1]}')
        
    # Try getting Rotations
    if 'Rotations' in string:
        try:
            line = string.split('Rotation') # Find Rotation string 
            string_dict = get_string_dict(line[-1])
            Rotations = eval(string_dict, globals_dict)
        except: log.error(f'ERROR could not resolve Rotations = {line[-1]}')
    return Regions, Rotations, Shift


########################################################
# Function to compute 3D distance (run-time optimized) #
########################################################
def compute_distance(x1, y1, z1, x2, y2, z2):
    dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
    return math.sqrt(dx*dx + dy*dy + dz*dz)


############################################
# Function to apply shift to atoms and box #
############################################
def shift_system(m, Shift, log):
    try:
        # Update simulation cell
        xline = m.xbox_line.split(); yline = m.ybox_line.split(); zline = m.zbox_line.split();
        xlo = float(xline[0]) + Shift['x']; xhi = float(xline[1]) + Shift['x'];
        ylo = float(yline[0]) + Shift['y']; yhi = float(yline[1]) + Shift['y'];
        zlo = float(zline[0]) + Shift['z']; zhi = float(zline[1]) + Shift['z'];
        m.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(xlo, xhi, 'xlo', 'xhi')
        m.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(ylo, yhi, 'ylo', 'yhi')
        m.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(zlo, zhi, 'zlo', 'zhi')
        
        # Update atom positions
        for i in m.atoms:
            atom = m.atoms[i]
            atom.x += Shift['x']
            atom.y += Shift['y']
            atom.z += Shift['z']
        log.out(f'Applied shift {Shift} to {os.path.basename(m.filename)}')
    except: log.error(f'ERROR failed to apply shifts {Shift} to {os.path.basename(m.filename)}')
    return m
        
        
#####################################################################################
# Function to recenter atoms in m-object and find max spanning dists from (0, 0, 0) #
#####################################################################################
def center_about_000(m):
    #--------------------------------#
    # Find center of atoms and shift #
    # needed to center at (0, 0, 0)  #
    #--------------------------------#
    x = []; y = []; z = [];
    for i in m.atoms:
        atom = m.atoms[i]
        x.append(atom.x)
        y.append(atom.y)
        z.append(atom.z)
    
    # Finding geometric center
    xcenter = sum(x)/len(x); ycenter = sum(y)/len(y); zcenter = sum(z)/len(z)
    
    # Finding needed shift to center around geometric center
    xshift = (0 - xcenter); yshift = (0 - ycenter); zshift = (0 - zcenter);
    
    #-----------------------------------------------------#
    # Recenter atoms and compute distances from (0, 0, 0) #
    #-----------------------------------------------------#
    distances_from_000 = []
    for i in m.atoms:
        atom = m.atoms[i]
        
        # shift atoms
        x = atom.x + xshift
        y = atom.y + yshift
        z = atom.z + zshift
        
        # Update atom
        atom.x = x; atom.y = y; atom.z = z;
        
        # compute distance from (0, 0, 0)
        distances_from_000.append(compute_distance(0, 0, 0, x, y, z))
    
    #-----------------------------------------------------------------#
    # Add maxspan attribute to m-object as 2*max(distances_froms_000) #
    #-----------------------------------------------------------------#
    m.maxspan = 2*max(distances_from_000)
    return m


############################################
# Function to unwrap atoms via image flags #
############################################
def unwrap_atoms_via_iflags(m):
    # Find box dimensions
    xline = m.xbox_line.split();
    yline = m.ybox_line.split();
    zline = m.zbox_line.split();
    
    # Lx, Ly, and Ly
    lx = float(xline[1])-float(xline[0])
    ly = float(yline[1])-float(yline[0])
    lz = float(zline[1])-float(zline[0])
    
    #-----------------------------------------------------#
    # Shift atoms by image flags and reset them to zeroes #
    #-----------------------------------------------------#
    for i in m.atoms:
        atom = m.atoms[i]
        
        # shift atoms
        x = atom.x + lx*atom.ix
        y = atom.y + ly*atom.iy
        z = atom.z + lz*atom.iz
        
        # Update atom positions and image flags
        atom.x = x; atom.y = y; atom.z = z;
        atom.ix = 0; atom.iy = 0; atom.iz = 0;
    return m


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
def rotate_molecule(m, phi, theta, psi):
    # find rotational matrix
    RzRy = matrix_by_matrix(Rz(psi), Ry(theta))
    RzRyRx = matrix_by_matrix(RzRy, Rx(phi))
    for i in m.atoms:
        atom = m.atoms[i]
        
        # Rotate around atoms
        x, y, z = vector_by_matrix(RzRyRx, [atom.x, atom.y, atom.z])
        
        # Save data into m.atoms
        atom.x = x
        atom.y = y
        atom.z = z
    return m


###############################################
# Function to find molecules and resid molids #
###############################################
# Function to find molids and analyze molecules      
def update_molids(m, log):
    # Find clusters
    log.out('\n\nFinding new molids ....')
    # Generate graph
    graph = {i:[] for i in m.atoms}
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        graph[id1].append(id2)
        graph[id2].append(id1)
    
    # Initialize clusters and molids
    checked = {ID:False for ID in m.atoms}; clusters = set([]);
    for ID in graph:
        if checked[ID]: continue
        visited=set([ID]); queue=[ID];
        while queue:
            s = queue.pop(0) 
            for neighbor in graph[s]:
                if checked[neighbor]: continue
                visited.add(neighbor)
                queue.append(neighbor)
                checked[neighbor]=True
        clusters.add( tuple(sorted(visited)) )

    # Sort clusters in a very unique fashion
    clusters = list(clusters)
    clusters = sorted(clusters, key=lambda x: x[0]) # Sort all clusters based on 1st atomID in cluster
    clusters = sorted(clusters, key=len, reverse=True) # Sort all clusters by number of atoms
    
    # Update atom molids
    for molID, cluster in enumerate(clusters, 1):
        for ID in cluster:
            m.atoms[ID].molid = molID
        
    # add clusters attribute to m
    m.clusters = clusters
    
    # Functon to compute cluster mass                                  
    def getmass(c):                                                                                                                  
        return sum([m.masses[m.atoms[i].type].coeffs[0] for i in c])
    
    # analyze clusters
    natoms = [];  cmass = [];   
    for c in clusters:
        natoms.append(len(c))
        cmass.append(getmass(c))     
                                                                                                                  
    # # summing mass and atoms                                                                                                                                                               
    # mass_total = sum(cmass)                                                                                     
    # size_total = sum(natoms)    

    # Print cluster info and roations applied
    log.out('--------------------------------------------------------------------------Cluster Analysis-------------------------------------------------------------------')
    log.out('{:^6} {:^10} {:^10} {:>10} {:>10} {:>10} {:>14} {:>14} {:>14} {:^25}'.format('molID', 'Size', 'Mass', 'X-center',
    'Y-center', 'Z-center', 'phi (\N{DEGREE SIGN})', 'theta (\N{DEGREE SIGN})', 'psi (\N{DEGREE SIGN})', 'Filename'))
    log.out('-------------------------------------------------------------------------------------------------------------------------------------------------------------')                                 
    clusterdata = {} # {clusterID : (tuple of data)}
    for n, (size, mass) in enumerate(zip(natoms, cmass)):  
        # find cluster to get cluster formula
        cluster = list(clusters[n])
        
        # find percent size and mass
        # pmass = 0; psize = 0;
        # if mass_total > 0: pmass = round(100*mass/mass_total, 2)
        # if size_total > 0: psize = round(100*size/size_total, 2)
        
        # Find x, y, and z center location from atomID at index 0 (all atoms will have same center location)
        x, y, z = m.atoms[cluster[0]].location
        
        # Find phi, theta, psi rotations from atomID at index 0 (all atoms will have same angles) and convert to degrees
        phi, theta, psi = m.atoms[cluster[0]].rotation # all clusters should have at least 1-atom, thus index 0
        
        # Find filename from atomID at index 0 (all atoms will have the filename)
        filename = os.path.basename(m.atoms[cluster[0]].filename)
        
        # Save data into a tuple to write .txt file
        clusterdata[n+1] = (size, mass, x, y, z, phi, theta, psi, filename)
        
        # print to screen
        log.out('{:^6} {:^10} {:^10.2f} {:>10.2f} {:>10.2f} {:>10.2f} {:>14.2f} {:>14.2f} {:>14.2f} {:>5} {}'.format(n+1, size, mass, x, y, z, phi, theta, psi, '', filename))
        
    # add clusterdata attribute to m
    m.clusterdata = clusterdata
    return m


##################################################
# Function to compute system mass/volume/density #
##################################################
def compute_mass_volume_density(m):
    system_mass_sum = 0; amu2grams = 1/6.02214076e+23;
    for i in m.atoms:
        atom = m.atoms[i]
        mass = m.masses[atom.type]
        system_mass_sum += mass.coeffs[0]
        
    # convert system mass in amu to grams
    mass = system_mass_sum*amu2grams
    
    # Find box dimensions to compute volume
    angstromcubed2cmcubed = 1e-24
    xline = m.xbox_line.split();
    yline = m.ybox_line.split();
    zline = m.zbox_line.split();
    lx = float(xline[1])-float(xline[0])
    ly = float(yline[1])-float(yline[0])
    lz = float(zline[1])-float(zline[0])
    volume = lx*ly*lz*angstromcubed2cmcubed
    
    # Compute density
    density = mass/volume
    return lx, ly, lz, volume, mass, density