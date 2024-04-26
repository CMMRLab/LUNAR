# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 16th, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

##############################
# Import Necessary Libraries #
##############################
import math



# Find hybridization by partial topological considerations and partial geometric considerations (Works for periodic systems as well)
def partially_topological_partially_geometric(mm, log):
    log.out('Finding atoms hybridization ...')
    
    # Find set simulataion cell bounds if periodicity is needed to be removed via the minimum image convention
    xline = mm.xbox_line.split(); yline = mm.ybox_line.split(); zline = mm.zbox_line.split();
    lx = float(xline[1])-float(xline[0])
    ly = float(yline[1])-float(yline[0])
    lz = float(zline[1])-float(zline[0])
    
    # Set VSEPR_approx_angles (degrees) to find the minimized difference of avg_angle to approx_angle to find hybridization state
    VSEPR_approx_angles = {'Sp1':     180,  # Will be used to set potenial Sp1 hybridized geometries
                           'Sp2':     120,  # Will be used to set potenial Sp2 hybridized geometries
                           'Sp3':     109.5 # Will be used to set potenial Sp3 hybridized geometries
                           }
    
    # Set elements that are known to be terminating (IE only 1-bonded neighbor)
    terminal = ['H', 'F', 'Cl', 'Ca', 'Br']
       

    # Find all angles and split angles into the following categories:
    # - atomIDs at the "center" of the angle (for sp1, sp2, and sp3 type of atoms)
    # - atomIDs on the "outside" of the angle (for terminating atoms like H) 
    angles = GTangles(mm)
    id2center_angles = {i:set() for i in mm.atoms} # { atomID : [list of angles, where atomID is in the center]}
    id2outer_angles =  {i:set() for i in mm.atoms} # { atomID : [list of angles, where atomID is at the outside]}
    for angle in angles:
        id1, id2, id3 = angle
        id2outer_angles[id1].add(angle)
        id2center_angles[id2].add(angle)
        id2outer_angles[id3].add(angle)

    
    
    ######################################################
    # Loop through atoms and start finding hybridization #
    ######################################################
    for i in mm.atoms:
        
        # Find general atom info
        atom = mm.atoms[i]; nb = int(atom.nb);
        ring = int(atom.ring); element = atom.element;
        atom.hybridization = 'unknown' # Initialize atom.hybridization and update if found.
        
        
        # Get local angles by searching atomID:i to be the center atom in the angle if nb >= 2: (use only 1st neighs)
        if nb >= 2:
            local_angles = id2center_angles[i]
        else: # else get local angles by searching for atomID:i as being and outside atom
            local_angles = id2outer_angles[i]
        if not local_angles: continue
        
        
        # Find avg_angle for all angles atom may belong too, based on how local angles are found above
        atom.avg_angle = compute_avg_angle_degrees(mm, lx, ly, lz, local_angles)
        
        
        # Find hybridization start with the most sure topological methods and move to geometry when needed
        # If element is known to be terminating and nb == 1 set 'Terminal' as hybridized state
        if element in terminal:
            atom.hybridization = 'Terminal'
            
        # If either 'C' or 'N' has 4 bonded neighbors it is guaranteed to be Sp3
        elif element in ['C', 'N'] and nb >= 4: 
            atom.hybridization = 'Sp3'
            
        # if ring is greater then 0 and less then or equal to 4 it is guaranteed to be Sp3
        elif ring > 0 and ring <= 4:
            atom.hybridization = 'Sp3'
        
        # if ring is greater then or equal to 5 and nb <= 3 it is guaranteed to be Sp2
        elif ring >= 5 and nb <= 3: 
            atom.hybridization = 'Sp2'
            
        # if atom is oxygen and nb == 2 and ring == 0 it is Sp3
        elif element == 'O' and nb == 2 and ring == 0:
            atom.hybridization = 'Sp3'
            
        # if atom is oxygen and nb == 1 it is Sp2
        # elif element == 'O' and nb == 1 and avg_angle >= 110 and avg_angle <= 130:
        #     atom.hybridization = 'Sp2'
            
        # if atom is oxygen and nb == 1 it is Sp1
        elif element == 'O' and nb == 1 and atom.avg_angle > 130:
            atom.hybridization = 'Sp1'
            
        # Moving to geometric methods which will only be used if the above logic did not find the hybridized state
        else:
            atom.hybridization = hybridization_via_minimized_diff_of_angles(atom.avg_angle, VSEPR_approx_angles)
            
        # Debug print
        #print(i, element, ring, nb, atom.hybridization, atom.avg_angle, local_angles)
    return mm



# Function to minimized difference between avg_angles and VSEPR_approx_angles dict
def hybridization_via_minimized_diff_of_angles(avg_angle, VSEPR_approx_angles):
    # Initialize hybridization as 'Unknown' and update when minimized angle is found
    hybridization = 'unknown'
    
    # Iteration through VSEPR_approx_angles and tally difference from avg to ideal
    diff_log = {} # { hybridization : abs(ideal_angle - avg_angle) }
    for hybrid in VSEPR_approx_angles:
        ideal_angle = VSEPR_approx_angles[hybrid]
        diff_log[hybrid] = abs(ideal_angle - avg_angle)
    
    # Update hybridization with the minimized difference of angles
    hybridization = min(diff_log, key=diff_log.get)
    return hybridization

# Function to find global angles
def GTangles(mm):
    # Set to add angles to
    angles = set()

    # Loop through bonds and graph to generate angles
    for i in mm.bonds:        
        id1, id2 = mm.bonds[i].atomids
        # Loop through id1 adjacency list
        for id3 in mm.graph[id1]:
            id123 = set([id1, id2, id3])
            if len(id123) == 3:
                # Sort angle in ascending order
                if id2 < id3:
                    angle = (id2, id1, id3) # id1 must be next to id3
                else:
                    angle = (id3, id1, id2) # id1 must be next to id3
                    
                # log sorted angle
                #if angle not in angles and tuple(reversed(angle)) not in angles:
                angles.add(angle)
                    
        # Loop through id2 adjacency list
        for id3 in mm.graph[id2]:
            id123 = set([id1, id2, id3])
            if len(id123) == 3:
                # Sort angle in ascending order
                if id1 < id3:
                    angle = (id1, id2, id3) # id2 must by next to id3 
                else:
                    angle = (id3, id2, id1) # id2 must by next to id3
                    
                # log sorted angle
                #if angle not in angles and tuple(reversed(angle)) not in angles:
                angles.add(angle)
    return angles


# Function to compute average angle of all bonded atoms
def compute_avg_angle_degrees(mm, lx, ly, lz, angles):
    # Initialize avg_angle as zero and update later on
    avg_angle = 0; angles_degrees = []; # lst to append found angles in degrees too
    
    # Loop through angles and find positon vectors and dot products
    for angle in angles:
        
        # Shift any periodic atoms and compute dot products and magnitudes
        shifted = shift_bonded_atoms(mm, angle, lx, ly, lz)
        v1, v2 = get_pos_vectors(shifted, angle)
        dotproduct = 0; mag_v1 = 0; mag_v2 = 0;
        for i, j in zip(v1, v2):
            dotproduct += i*j
            mag_v1 += i*i
            mag_v2 += j*j
        mag_v1 = math.sqrt(mag_v1)
        mag_v2 = math.sqrt(mag_v2)
        
        # magnitudes could be zero if so error will be triggered so try computing
        try:
            angle_rad = math.acos(dotproduct/(mag_v1*mag_v2))
            angle_deg = math.degrees(angle_rad)
            angles_degrees.append(angle_deg)
        except:  pass
        
    # if angles_degrees is not empty re-compute avg_angle
    if angles_degrees:
        avg_angle = sum(angles_degrees)/len(angles_degrees)
    return avg_angle


    
# Function to get position vectors from angle
def get_pos_vectors(shifted, angle):
    # Initialize v1 and v2 as zeros and update when found
    v1 = (0, 0, 0); v2 = (0, 0, 0);
    
    # Find atomIDs from angle with id2 as the common atom
    id1 = angle[0]; id2 = angle[1]; id3 = angle[2];
    x1, y1, z1 = shifted[id1]
    x2, y2, z2 = shifted[id2]
    x3, y3, z3 = shifted[id3]
    
    # Update position vectors moving from ID2 to ID1 or ID3
    v1 = (x2-x1, y2-y1, z2-z1); v2 = (x2-x3, y2-y3, z2-z3);
    return v1, v2
    
# Function to shift x1, y1, and z1 using x2, y2, and z2 as the anchoring atom
def shift_atom_via_anchor(x1, y1, z1, x2, y2, z2, lx, ly, lz):
    # if atom is near PBC shift one atom by half of the box dimension
    diffx = x2 - x1
    if diffx > 0.5*lx:
        x1 = x1 + lx
    elif diffx < -0.5*lx:
        x1 = x1 - lx 
    diffy = y2 - y1
    if diffy > 0.5*ly:
        y1 = y1 + ly 
    elif diffy < -0.5*ly:
        y1 = y1 - ly 
    diffz = z2 - z1
    if diffz > 0.5*lz:
        z1 = z1 + lz 
    elif diffz < -0.5*lz:
        z1 = z1 - lz 
    return (x1, y1, z1)

# Function to shift any atoms that are found to be periodic based on the minimum image convention
def shift_bonded_atoms(mm, angle, lx, ly, lz):
    shifted = {} # {atomID : (x, y, z)}
    id1, id2, id3 = angle
    
    # Do not shift id2 as this is the anchor (center atom)
    atom2 = mm.atoms[id2]
    x2 = atom2.x; y2 = atom2.y; z2 = atom2.z;
    shifted[id2] = (x2, y2, z2)
    
    # Shift id1 and id3 (outer atoms)
    atom1 = mm.atoms[id1]; atom3 = mm.atoms[id3];
    x1 = atom1.x; y1 = atom1.y; z1 = atom1.z;
    x3 = atom3.x; y3 = atom3.y; z3 = atom3.z;
    shifted[id1] = shift_atom_via_anchor(x1, y1, z1, x2, y2, z2, lx, ly, lz)
    shifted[id3] = shift_atom_via_anchor(x3, y3, z3, x2, y2, z2, lx, ly, lz)
    return shifted

  


###############################################
# Function to create hybridization data table #
###############################################
class Data:
    pass #  .size .mass .pmass .psize
    
def hybridization_data(mm):
        
    """
    EXAMPLE DATA STRUCT:
    
    DATA of Sp2 C member rings:
        mm.hybridization['C']['Sp2'].size
        mm.hybridization['C']['Sp2'].mass
        mm.hybridization['C']['Sp2'].pmass
        mm.hybridization['C']['Sp2'].psize
    
    """
    
    # Intialize new attribute of mm as hybridization
    mm.hybridization = {} # { element symbol: {'Sp1': data object, 'Sp2': data object, 'Sp3': data object, 'unknown': data object } }  mass_Sp2_C = data['C']['Sp2'].mass
    
    ###############################
    # Initialization of self.data #
    # structure to hold all info  #
    ###############################
    # All categoreis of sub dictionary
    elems2skip = ['H'] # element types to skip hybridization and stick in 'all' key
    hybrid = ['Sp1', 'Sp2', 'Sp3', 'unknown', 'all', 'Terminal']
    
    # Create data dict in dict with instances
    for i in mm.elements:
        tmp_dict = {};
        for j in hybrid:
            d = Data()
            d.size = 0
            d.mass = 0
            d.pmass = 0
            d.psize = 0
            tmp_dict[j] = d
        
        # Create dict in dict with instances as values of sub dict
        mm.hybridization[i] = tmp_dict   
        
    #####################################################
    # Build self.data object with re-hybridized results #
    #####################################################
    for i in mm.atoms:
        atom = mm.atoms[i]
        element = atom.element
        hybridization = atom.hybridization
        mass = mm.masses[mm.atoms[i].type].coeffs[0]
        
        # add element to specific hybridization
        if element not in elems2skip:
            mm.hybridization[element][hybridization].size += 1
            mm.hybridization[element][hybridization].mass += mass 
        
        # add elememt data to all
        mm.hybridization[element]['all'].size += 1
        mm.hybridization[element]['all'].mass += mass  
            
    # Find percents and update
    for i in mm.hybridization:
        for j in mm.hybridization[i]:
            # Find pmass and psize
            pmass = 0; psize = 0;
            if mm.total_system_mass > 0: pmass = round(100*mm.hybridization[i][j].mass/mm.total_system_mass, 2)
            if mm.total_system_size > 0: psize = round(100*mm.hybridization[i][j].size/mm.total_system_size, 2)
            
            # Update mm.hybridization['element']['hybridization']
            mm.hybridization[i][j].pmass = pmass
            mm.hybridization[i][j].psize = psize
    return mm