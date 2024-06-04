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
import math




#################################
# Function to generate sheet(s) #
#################################
# https://www.allmathwords.org/en/h/hexagon.html#:~:text=center%3A%20The%20point%20at%20the,hexagon%20is%20also%20its%20centroid.
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz .comment .atomtype
def generate(lx, ly, r0, edgetype, types, layer_spacing, nlayers, stacking, plane, periodic_bonds, pflag, log):
    log.out('\n\nGenerating sheet(s) ...')
    
    # Generate starting atoms base unit to
    # duplicate for zigzag edge type
    # Y
    # ^
    # |  type2___type3
    # |      /   \
    # | type1     type4
    # |
    # +-----------------> X
    if edgetype == 'zigzag':
        base_unit = [] # [ [x, y, z, atomtype, type] ]
        dx = r0*math.cos(math.radians(60))
        dy = r0*math.sin(math.radians(60))
        shiftx = 2*r0 + 2*dx
        shifty = 2*dy
        base_unit.append( [0,         0,  0, types[1], 1] )
        base_unit.append( [dx,        dy, 0, types[2], 2] )
        base_unit.append( [r0 + dx,   dy, 0, types[3], 3] )
        base_unit.append( [r0 + 2*dx, 0,  0, types[4], 4] )
        
        # shift half a bond length to account for box +- r0/2
        for i in range(len(base_unit)):
            base_unit[i][0] += r0/2
            base_unit[i][1] += r0/2
            
        # compute nx and ny number of rings based on lx and ly
        ringx_span = 2*r0 + 2*dx
        ringy_span = 2*dy
        nx = math.floor(lx/ringx_span)
        ny = math.floor(ly/ringy_span)
            
    # Generate starting atoms base unit to
    # duplicate for armchair edge type
    # Y
    # ^
    # | type4\
    # |       \type3
    # |       |
    # |       |type2
    # |       /
    # |      /
    # | type1
    # +-------------> X
    elif edgetype == 'armchair':
        base_unit = [] # [ [x, y, z, atomtype, type] ]
        dx = r0*math.cos(math.radians(30))
        dy = r0*math.sin(math.radians(30))
        shiftx = 2*dx
        shifty = 2*r0 + 2*dy
        base_unit.append( [0,       0,     0, types[1], 1] )
        base_unit.append( [dx,      dy,    0, types[2], 2] )
        base_unit.append( [dx, r0 + dy,    0, types[3], 3] )
        base_unit.append( [0,  r0 + 2*dy,  0, types[4], 4] )
        
        # shift half a bond length to account for box +- r0/2
        for i in range(len(base_unit)):
            base_unit[i][0] += r0/2
            base_unit[i][1] += r0/2
            
        # compute nx and ny number of rings based on lx and ly
        ringx_span = 2*dx
        ringy_span = 2*r0 + 2*dy
        nx = math.floor(lx/ringx_span)
        ny = math.floor(ly/ringy_span)
    
    # Else log an error
    else:
        log.error(f"ERROR edgetype {edgetype} is not supported. Supported edgetype's are 'zigzag' or 'armchair'")
    
    # Duplicate base unit to generate sheet
    atoms = {} # { atomID : Atoms Object }
    ID = 0
    if nx == 0: nx = 1
    if ny == 0: ny = 1
    spanx = nx*ringx_span
    spany = ny*ringy_span
    height = (1 + nlayers)*layer_spacing
    for n in range(1, nlayers+1):
        xoffset = 0
        yoffset = 0
        zoffset = n*layer_spacing - height/2
        if stacking in ['AB', 'ABC'] and n%2 == 0:
            if edgetype == 'zigzag': xoffset = r0
            if edgetype == 'armchair': yoffset = r0
        if stacking == 'ABC' and n%3 == 0:
            if edgetype == 'zigzag': xoffset = 2*r0
            if edgetype == 'armchair': yoffset = 2*r0
        log.out('  edge={} with an edge length={:.4f} and a perpendicular length={:.4f}'.format(edgetype, ny*ringy_span, nx*ringx_span))
        for i in range(nx):
            for j in range(ny):
                for xx, yy, zz, atomtype, typeint in base_unit:
                    newx = i*shiftx + xx - spanx/2 # center atoms in X-dir
                    newy = j*shifty + yy - spany/2 # center atoms in Y-dir
                    newz = zz + zoffset
                    ID += 1
                    a = Atom()
                    a.molid = n
                    a.type = typeint
                    a.atomtype = atomtype
                    a.comment = atomtype
                    a.charge = 0.0
                    a.x = newx + xoffset # apply here to get "undisturbed" positions for box and then wrap atoms later
                    a.y = newy + yoffset # apply here to get "undisturbed" positions for box and then wrap atoms later
                    a.z = newz
                    a.ix = 0
                    a.iy = 0
                    a.iz = 0
                    atoms[ID] = a
                
    # Create box
    true_lx = nx*ringx_span
    true_ly = ny*ringy_span
    true_lz = nlayers*layer_spacing
    box = {'xlo': -true_lx/2,
           'xhi':  true_lx/2,
           'ylo': -true_ly/2,
           'yhi':  true_ly/2,
           'zlo': -true_lz/2,
           'zhi':  true_lz/2}  

    # Wrap atoms if stacking was AB or ABC and nylayers was greater then 1
    if stacking in ['AB', 'ABC'] and nlayers > 1:
        atoms = misc_functions.wrap_atoms(atoms, box) 
        
    # If plane is 'xz' or 'yz' rotate atoms
    if plane in ['xz', 'yz']:
        phi = 0; theta = 0; psi = 0;
        if plane == 'xz': phi = 90
        if plane == 'yz': theta = 90
        atoms = misc_functions.rotate_molecule(atoms, phi, theta, psi)
        
        # Reset box based on rotation
        if plane == 'xz': # y->z            
            box = {'xlo': -true_lx/2,
                   'xhi':  true_lx/2,
                   'ylo': -true_lz/2,
                   'yhi':  true_lz/2,
                   'zlo': -true_ly/2,
                   'zhi':  true_ly/2}  
        if plane == 'yz': # x->z
            box = {'xlo': -true_lz/2,
                   'xhi':  true_lz/2,
                   'ylo': -true_ly/2,
                   'yhi':  true_ly/2,
                   'zlo': -true_lx/2,
                   'zhi':  true_lx/2}
            
    # if system is not period increase box size
    if not periodic_bonds:
        increase = 2*r0
        if plane == 'xy':
            box['xlo'] -= increase
            box['xhi'] += increase
            box['ylo'] -= increase
            box['yhi'] += increase
        if plane == 'xz':
            box['xlo'] -= increase
            box['xhi'] += increase
            box['zlo'] -= increase
            box['zhi'] += increase
        if plane == 'yz':
            box['ylo'] -= increase
            box['yhi'] += increase
            box['zlo'] -= increase
            box['zhi'] += increase
            
    return atoms, box