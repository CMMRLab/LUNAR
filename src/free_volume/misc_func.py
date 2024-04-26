# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
September 28th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import math


###################################################################################
# Function to compute distance when math.dist is not available. math.dist will be #
# quicker but for some reason not all python has math.dist .... compute_distance  #
# will be time optimized but still will be slower then math.dist.                 #
###################################################################################        
def compute_distance(x1, y1, z1, x2, y2, z2):
    dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
    return math.sqrt(dx*dx + dy*dy + dz*dz)


###########################################################################
# Function to compute distance sqrd to avoid having to take a square root #
# which a very slow calculation.compute_distance_sqrd is time optimized   #
###########################################################################  
def compute_distance_sqrd(x1, y1, z1, x2, y2, z2):
    dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
    return dx*dx + dy*dy + dz*dz


##################################
# Function to get box dimensions #
##################################
def get_box_dimensions(m):
    # Find box dimensions
    xline = m.xbox_line.split();
    yline = m.ybox_line.split();
    zline = m.zbox_line.split();
    
    # Find ilo, ihi
    xlo = float(xline[0]); xhi = float(xline[1]);
    ylo = float(yline[0]); yhi = float(yline[1]);
    zlo = float(zline[0]); zhi = float(zline[1]);
    
    
    # Lx, Ly, and Ly box dimensions
    lx = float(xline[1])-float(xline[0])
    ly = float(yline[1])-float(yline[0])
    lz = float(zline[1])-float(zline[0])
    
    # Find Cx, Cy, Cz box centers
    cx = (float(xline[1])+float(xline[0]))/2
    cy = (float(yline[1])+float(yline[0]))/2
    cz = (float(zline[1])+float(zline[0]))/2
    return lx, ly, lz, cx, cy, cz, xlo, xhi, ylo, yhi, zlo, zhi


##################################################################################
# Function to check if the atom is near a box edge. Then build near edge dict to #
# use as a look up table to aviod multiple re-calculations as the loops progress #
##################################################################################
def check_near_edge(x, y, z, max_from_edge, xlo, xhi, ylo, yhi, zlo, zhi):
    if abs(x - xlo) <= max_from_edge or abs(x + xlo) <= max_from_edge: pbcflag = True
    elif abs(x - xhi) <= max_from_edge or abs(x + xhi) <= max_from_edge: pbcflag = True  
      
    elif abs(y - ylo) <= max_from_edge or abs(y + ylo) <= max_from_edge: pbcflag = True
    elif abs(y - yhi) <= max_from_edge or abs(y + yhi) <= max_from_edge: pbcflag = True

    elif abs(z - zlo) <= max_from_edge or abs(z + zlo) <= max_from_edge: pbcflag = True
    elif abs(z - zhi) <= max_from_edge or abs(z + zhi) <= max_from_edge: pbcflag = True
    else: pbcflag = False
    return pbcflag

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


################################################
# Generate images via minimum image convention #
################################################
def generate_iflags(m, boundary, log):
    # Boundary conditions
    pflags = boundary.split() # split boundary
    count = pflags.count('f') + pflags.count('p')
    if len(pflags) != 3 and count != 3:
        log.error('ERROR boundary does not contain 3-dimensions or boundary flags are not f or p')
    if m.xy != 0 or m.xz != 0 or m.yz != 0:
        log.error('ERROR simulaiton cell is NOT orthogonal.')
    
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


#################################
# Find free volume cluster span #
#################################
def find_free_volume_cluster_span(cluster, voxel_freeIDs):
    xpos = []; ypos = []; zpos = [];
    for ID in cluster:
        x, y, z, = voxel_freeIDs[ID]
        xpos.append(x)
        ypos.append(y)
        zpos.append(z)
     
    # Find spans
    xspan = 0; yspan = 0; zspan = 0;
    if xpos: xspan = abs(max(xpos) - min(xpos))
    if ypos: yspan = abs(max(ypos) - min(ypos))
    if zpos: zspan = abs(max(zpos) - min(zpos))
    return xspan, yspan, zspan