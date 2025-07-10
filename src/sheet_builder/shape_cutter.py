# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
June 9, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


https://xgeometry.com/formulas/hexagon
"""
##############################
# Import Necessary Libraries #
##############################
import src.sheet_builder.misc_functions as misc_functions
import math




###############################################
# Function for cutting out atoms in a lattice #
###############################################
def cut(atoms, box, bond_length, plane, sheet_edgetype, cutter, log):    
    # Finding some defaults
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']
    cx = (box['xhi'] + box['xlo'])/2
    cy = (box['yhi'] + box['ylo'])/2
    cz = (box['zhi'] + box['zlo'])/2
    tol = bond_length/6
    a = bond_length*(math.sqrt(3)/2)
    log.out('  Box dimensions (lx, Ly, Lz)    : {:10.6f}, {:10.6f}, {:10.6f}'.format(lx, ly, lz))
    log.out('  Box center (cx, cy, cz)        : {:10.6f}, {:10.6f}, {:10.6f}'.format(cx, cy, cz))

    
    # Find the centroid of all atoms
    xs, ys, zs = [], [], []  # {molid:[xs]}, {molid:[ys]}, {molid:[zs]}
    for i in atoms:
        atom = atoms[i]
        molid = atom.molid
        xs.append(atom.x)
        ys.append(atom.y)
        zs.append(atom.z)
    xc = sum(xs)/len(atoms)
    yc = sum(ys)/len(atoms)
    zc = sum(zs)/len(atoms)
    log.out('  Atoms centroid (x, y, z)       : {:10.6f}, {:10.6f}, {:10.6f}'.format(xc, yc, zc))
    
    
    # Find atoms atoms distance from centroid to determine center ring centroid
    dist2centroid = {} # {atomID:distance-to-centroid}
    for i in atoms:
        atom = atoms[i]
        x = atom.x
        y = atom.y
        z = atom.z
        dist2centroid[i] = misc_functions.compute_distance(x, y, z, xc, yc, zc)

    # The 6 nearest atoms to the centroid are the atoms of the central ring
    dist2centroid = dict(sorted(dist2centroid.items(), key=lambda x: x[1]))
    center_ring_atoms = list(dist2centroid.keys())[:6]
    rxs, rys, rzs = [], [], []
    for i in center_ring_atoms:
        #atom = atoms[i].atomtype = 'O'
        atom = atoms[i]
        rxs.append(atom.x)
        rys.append(atom.y)
        rzs.append(atom.z)
    
    # Find centroid of central ring (default to global centroid if all no atoms where found)
    if rxs and rys and rzs:
        rxc = sum(rxs)/len(rxs)
        ryc = sum(rys)/len(rys)
        rzc = sum(rzs)/len(rzs)
    else:
        rxc = xc
        ryc = yc
        rzc = zc
    log.out('  Center ring centroid (x, y, z) : {:10.6f}, {:10.6f}, {:10.6f}'.format(rxc, ryc, rzc))

    
    current_atoms = {i for i in atoms}
    delete_atoms = set()
    for string in cutter:
        settings = get_setting(string, lx, ly, lz, cx, cy, cz, bond_length, rxc, ryc, rzc, a)
        if 'shape' not in settings:
            log.warn(f'  WARNING shape was not set in cutter string "{string}". Skipping cutting operation')
            continue

        shape = settings['shape']
        del settings['shape']
        log.out(f'\n  Starting "{shape}" cutting operation:')
        if shape in ['h', 'hex', 'hexagon']:
            # Define defaults to this cutting options
            parms = {'plane': plane,  # This NEEDS to be the first key, so updates happen in the correct ordering
                     'theta': 'hexagon',
                     'len': max(lx, ly, lz),
                     'del': 'in',
                     'tol': tol,
                     'rxc': rxc,
                     'ryc': ryc,
                     'rzc': rzc,
                     'sx': 0.0,
                     'sy': 0.0,
                     'sz': 0.0,
                     'a': 2*bond_length}
            
            # Set theta's based on plane, edge type, and desired cut-out shape
            shapes = {'h':'hexagon',
                      's': 'snowflake',
                      'hex':'hexagon',
                      'snow': 'snowflake',
                      'flake': 'snowflake'}
            thetas = {('xy', 'armchair', 'hexagon')   : 0.0,
                      ('xy', 'armchair', 'snowflake') : 30.0,
                      ('xy', 'zigzag',   'hexagon')   : 30.0,
                      ('xy', 'zigzag',   'snowflake') : 0.0,
                      
                      ('xz', 'armchair', 'hexagon')   : 0.0,
                      ('xz', 'armchair', 'snowflake') : 30.0,
                      ('xz', 'zigzag',   'hexagon')   : 30.0,
                      ('xz', 'zigzag',   'snowflake') : 0.0,
                      
                      ('yz', 'armchair', 'hexagon')   : 30.0,
                      ('yz', 'armchair', 'snowflake') : 0.0,
                      ('yz', 'zigzag',   'hexagon')   : 0.0,
                      ('yz', 'zigzag',   'snowflake') : 30.0}
            
            # Set up defaults based on plane
            if parms['plane'] == 'xy':
                parms['len'] = lz                
            elif parms['plane'] == 'xz':
                parms['len'] = ly
            elif parms['plane'] == 'yz':
                parms['len'] = lx
            else:
                log.error(f'ERROR plane={parms["plane"]} is not suppported. Supported planes are: "xy" or "xz" or "yz"')
            
            # Update defaults based on user settings
            log.out('    Center ring centroid will be subtracted from positions')
            for key in parms:
                # User-defined settings
                if key in settings:
                    value = settings[key]
                    
                    # check if theta is 'hexagon' or 'snowflake' and update accordingly
                    if key == 'theta' and not isinstance(value, (int, float)):
                        if value in shapes: value = shapes[value] # shorr cut names
                        key2thetas = (parms['plane'], sheet_edgetype, value)
                        if key2thetas in thetas:
                            value = thetas[key2thetas]
                        else: value = 0    
                    log.out('    User-defined  setting : {}={}'.format(key, value))
                    parms[key] = value
                
                # Default settings
                else:
                    value = parms[key]
                    
                    # check if theta is 'hexagon' or 'snowflake' and update accordingly
                    if key == 'theta' and not isinstance(value, (int, float)):
                        if value in shapes: value = shapes[value] # shorr cut names
                        key2thetas = (parms['plane'], sheet_edgetype, value)
                        if key2thetas in thetas:
                            value = thetas[key2thetas]
                        else: value = 0     
                        parms[key] = value
                    log.out('    Using-default setting : {}={}'.format(key, value))

            
            # Find atoms inside hex and then add to atoms to delete
            inside_atoms = find_hex_atoms(atoms, parms)
            
            # Set which atoms to delete
            if parms['del'] == 'in':
                delatoms = inside_atoms
                delete_atoms.update(delatoms) 
                log.out('    Number of found atoms : {}'.format(len(delatoms)))
            elif parms['del'] == 'out':
                delatoms = current_atoms - inside_atoms
                delete_atoms.update(delatoms)
                log.out('    Number of found atoms : {}'.format(len(delatoms)))
            else:
                log.error(f'ERROR del={parms["del"]} is not supported. Only "in" or "out" is supported')
    
    # Go through and delete atoms
    for i in delete_atoms:
        del atoms[i]
    
    return atoms


###########################################################
# Function for finding atoms inside an "extruded hexagon" #
###########################################################
def find_hex_atoms(atoms, parms):
    # Determine rotation (if any to apply to atoms)
    rotation = parms['theta']
    phi, theta, psi = 0, 0, 0
    if parms['plane'] == 'xy':
        psi = rotation
    elif parms['plane'] == 'xz':
        theta = rotation
    elif parms['plane'] == 'yz':
        phi = rotation
    RzRy = misc_functions.matrix_by_matrix(misc_functions.Rz(psi), misc_functions.Ry(theta))
    RzRyRx = misc_functions.matrix_by_matrix(RzRy, misc_functions.Rx(phi))


    # Calculate hexagon vertex distance and start labeling atoms. Hexagon's 3 principal axis:
    #    r
    #   /
    #  /
    # *------ q
    #  \
    #   \
    #    s
    # which will be used as q, r, and s below
    hex_radius = parms['a'] - parms['tol']
    hex_height = parms['len']/2
    inside_atoms = set()
    x0 = parms['rxc'] + parms['sx']
    y0 = parms['ryc'] + parms['sy']
    z0 = parms['rzc'] + parms['sz']
    for i in atoms:
        atom = atoms[i]
        x = atom.x - x0
        y = atom.y - y0
        z = atom.z - z0
        
        if rotation != 0:
            x, y, z = misc_functions.vector_by_matrix(RzRyRx, [x, y, z])
        
        if parms['plane'] == 'xy':
            # Compute axial hex coordinates
            q = (x*math.sqrt(3)/3 - y/3)
            r = (y*2/3)
            s = -q - r
            
            # Keep atoms inside hex 
            if max(abs(q), abs(r), abs(s)) <= hex_radius and abs(z) <= hex_height:
                inside_atoms.add(i)
                
        elif parms['plane'] == 'xz':
            # Compute axial hex coordinates
            q = (x*math.sqrt(3)/3 - z/3)
            r = (z*2/3)
            s = -q - r
            
            # Keep atoms inside hex 
            if max(abs(q), abs(r), abs(s)) <= hex_radius and abs(y) <= hex_height:
                inside_atoms.add(i)
                
        elif parms['plane'] == 'yz':
            # Compute axial hex coordinates
            q = (y*math.sqrt(3)/3 - z/3)
            r = (z*2/3)
            s = -q - r
            
            # Keep atoms inside hex 
            if max(abs(q), abs(r), abs(s)) <= hex_radius and abs(x) <= hex_height: 
                inside_atoms.add(i)

    return inside_atoms


########################################
# Function to get settings from string #
########################################
def get_setting(string, lx, ly, lz, cx, cy, cz, r0, rxc, ryc, rzc, a):
    # Setup the globals namespace to limit scope of what eval() can do
    allowed_builtins = ['min','max','sum','abs','len','map','range','reversed']
    copied_builtins = globals()['__builtins__'].copy()
    globals_dict = {}
    globals_dict['__builtins__'] = {key:copied_builtins[key] for key in allowed_builtins}
    globals_dict['a'] = a
    globals_dict['r0'] = r0
    globals_dict['lx'] = lx
    globals_dict['ly'] = ly
    globals_dict['lz'] = lz
    globals_dict['cx'] = cx
    globals_dict['cy'] = cy
    globals_dict['cz'] = cz

    
    # Parse misc string
    setting = {} # {keyword:float or int or Boolean}
    tmp1 = string.split(';')
    for tmp2 in tmp1:
        tmp3 = tmp2.split('=')
        if len(tmp3) >= 2:
            i = tmp3[0].strip()
            try: j = eval(tmp3[1], globals_dict)
            except: j = str(tmp3[1])
            if j == 'a': j = a
            if j == 'lx': j = lx
            if j == 'ly': j = ly
            if j == 'lz': j = lz
            if j == 'cx': j = cx
            if j == 'cy': j = cy
            if j == 'cz': j = cz
            if j == 'r0': j = r0
            if j == 'rxc': j = rxc
            if j == 'ryc': j = ryc
            if j == 'rzc': j = rzc
            setting[i] = j
    return setting