# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 18th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

This file contains useful functions to find 
desired information from mm.atoms[ID] instances
to help in finding atom-types for different
force feilds
"""


# Function to write to assumed file when needed
def write2assumed(file, tag, assumed, atom, ff_name, atomid):
    file.write(f'{tag.replace("-", "  ")} FAILED to be found with search criteria of ring size, number of connects, and connected elements.\n')
    file.write(f'Most likely cause is that {ff_name} doesnt have an atom type for this application and assumed atom type was put in its place.\n')
    file.write('You may go into the ouputted *.nta file and give this atom type a different type if you would like.\n')
    file.write('neigh info lst = [[element1, ring1, nb1], [element2, ring2, nb2], ....] sorted by nb -> ring -> element\n')
    try: file.write(f'1st-neigh info : {str(atom.neighbor_info[1])}\n')
    except: pass
    try: file.write(f'2nd-neigh info : {str(atom.neighbor_info[2])}\n')
    except: pass
    try: file.write(f'3rd-neigh info : {str(atom.neighbor_info[3])}\n')
    except: pass
    file.write(f'atom number    : {atomid} could not be characterized\n')
    file.write(f'element        : {atom.element}\n')
    file.write(f'ring size      : {atom.ring}\n')
    file.write(f'# of connects  : {atom.nb}\n')
    file.write(f'atom type      : {assumed} was assumed\n\n\n')
    return


# Function to write to failed file when needed
def write2failed(file, atom, ff_name, atomid):
    file.write('FAILED TO BE CHARACTERIZED with search criteria of ring size, number of connects, and connected elements.\n')
    file.write(f'Most likely cause is that {ff_name} doesnt have an atom type for this application and YOU MUST SET ATOM TYPE MANUALLY.\n')
    file.write('You may go into the ouputted *.nta file and give this atom type a different type based on the following info.\n')
    file.write('neigh info lst = [[element1, ring1, nb1], [element2, ring2, nb2], ....] sorted by nb -> ring -> element\n')
    try: file.write(f'1st-neigh info : {str(atom.neighbor_info[1])}\n')
    except: pass
    try: file.write(f'2nd-neigh info : {str(atom.neighbor_info[2])}\n')
    except: pass
    try: file.write(f'3rd-neigh info : {str(atom.neighbor_info[3])}\n')
    except: pass
    file.write(f'atom number    : {atomid} could not be characterized\n')
    file.write(f'element        : {atom.element}\n')
    file.write(f'ring size      : {atom.ring}\n')
    file.write(f'# of connects  : {atom.nb}\n\n\n')
    return


# Function to return neigh information
def neigh_extract(atom, depth, info):
    lst = []; info_map = {'element':0, 'ring':1, 'nb':2};
    for i in atom.neighbor_info[depth]:
        if i: lst.append(i[info_map[info]])
    return lst


# Function to find neighbors ringIDs (Initiallly built for DREIDING)
def get_neigh_ringIDs(atomID, m, depth, criteria):
    neigh_ringIDs = [] # lst of neigh ringIDs
    
    # Find neighs to find ringIDs for
    neighs = m.atoms[atomID].neighbor_ids[depth]
    for i in neighs:
        ringID = m.atoms[i].ringID
        ringsize = m.atoms[i].ring
        flag = True # Intialize as True and update based on criteria
        
        # ringsize criteria
        if 'minringsize' in criteria:
            if ringsize < criteria['minringsize']:
                flag = False
        
        # log based on flag
        if flag:
            neigh_ringIDs.append(ringID)
    
    return neigh_ringIDs


# Function to count neighbor info (set element, ring, or nb to False if not desired to count)
def count_neigh(neigh_info, element, ring, nb):
    count = 0
    for i in neigh_info:
        
        # If element, ring, and nb use all criteria to count
        if element and ring and nb:
            if i[0] == element and i[1] == ring and i[2] == nb:
                count += 1
                
        # elif only search using element and ring
        elif element and ring:
            if i[0] == element and i[1] == ring:
                count += 1
                
        # elif only search using element and nb
        elif element and nb:
            if i[0] == element and i[2] == nb:
                count += 1
                
        # elif only search using ring and nb
        elif element and nb:
            if i[1] == ring and i[2] == nb:
                count += 1
                
        # else raise Exception
        else:
            raise Exception('count_neigh_info function is being requested for unsupported neigh counting')
    return count


# Function to find percent neighbor info (set element, ring, or nb to False if not desired to count)
def percent_neigh(neigh_info, element, ring, nb):
    count = count_neigh(neigh_info, element, ring, nb)
    return 100*count/len(neigh_info)


# Function to count heavy atoms
def count_heavies(bonded_elements, heavies):
    count = 0;
    for i in bonded_elements:
        if i in heavies:
            count += 1
    return count


# Function to check for PEO (polyethylene oxide) topologies
def check_peo_topo(atom, element_type):
    return_boolean = False # Intialize as False and update as True if topology makes sense
    
    # Find 1st-neighbor info (denoted by nameN; where N=neighbor depth)
    elements1 = neigh_extract(atom, depth=1, info='element') # Example: ['C', 'C', 'N']
    rings1 = neigh_extract(atom, depth=1, info='ring') # Example: [6, 6, 0]
    nbs1 = neigh_extract(atom, depth=1, info='nb') # Example: [3, 3, 3]
    
    # Find 2nd-neighbor info (denoted by nameN; where N=neighbor depth)
    elements2 = neigh_extract(atom, depth=2, info='element') # Example: ['C', 'C', 'N']
    rings2 = neigh_extract(atom, depth=2, info='ring') # Example: [6, 6, 0]
    nbs2 = neigh_extract(atom, depth=2, info='nb') # Example: [3, 3, 3]
    
    # Find 2nd-neighbor info (denoted by nameN; where N=neighbor depth)
    elements3 = neigh_extract(atom, depth=3, info='element') # Example: ['C', 'C', 'N']
    rings3 = neigh_extract(atom, depth=3, info='ring') # Example: [6, 6, 0]
    nbs3 = neigh_extract(atom, depth=3, info='nb') # Example: [3, 3, 3]
    
    # Find neighN lsts with shorter name
    neigh1 = atom.neighbor_info[1]
    neigh2 = atom.neighbor_info[2]
    neigh3 = atom.neighbor_info[3]
    
    #############################################
    # Perform PEO typing for the Carbon element #
    #############################################
    if element_type == 'C':
        # PEO Carbon Backbone check for 1st neighboring elements, then check 2nd neighboring elements
        if elements1[0] == 'C' and nbs1[0] == 4 and elements1.count('H') == 2 and elements1[3] == 'O' and nbs1[3] == 2 and all(rings1) == 0:
            if elements2[0] == 'C' and nbs2[0] == 4 and elements2.count('H') == 2 and elements2[3] == 'O' and nbs2[3] == 2 and all(rings2) == 0:
                return_boolean = True
                
        # PEO Carbon Terminal group check for 1st neighboring elements, then check 2nd neighboring elements
        if elements1[0] == 'C' and nbs1[0] == 4 and elements1.count('H') == 3 and elements2.count('H') == 2 and elements2.count('O') == 2 and len(elements2) == 3:
            return_boolean = True
            
    ###############################################
    # Perform PEO typing for the Hydrogen element #
    ###############################################
    elif element_type == 'H':
        # PEO Hydrogen Backbone -CH check that 1st-neigh looks like PEO, then 2nd-neigh lsts looks like PEO
        if len(neigh1) >= 1 and len(neigh2) >= 3:
            if neigh1[0] == ['C', 0, 4] and neigh2[0] == ['C', 0, 4] and neigh2[1] == ['H', 0, 1] and neigh2[2] == ['O', 0, 2]:
                return_boolean = True
                
        # PEO Hydrogen Terminal -CH group check that 1st-neigh looks like PEO, then 2nd-neigh lsts looks like PEO, then 3rd neigh O-element
        if len(neigh1) >= 1 and len(neigh2) >= 3:
            if neigh1[0] == ['C', 0, 4] and neigh2[0] == ['C', 0, 4] and neigh2[1] == ['H', 0, 1] and neigh2[2] == ['H', 0, 1] and elements3.count('O') == 2:
                return_boolean = True
                
        # PEO Hydrogen Terminal -OH group check that 1st-neigh looks like PEO, then 2nd-neigh lsts looks like PEO and 3rd
        if len(neigh1) >= 1 and len(neigh2) >= 1:
            if neigh1[0] == ['O', 0, 1] and neigh2[0] == ['C', 0, 4] and neigh3[0] == ['C', 0, 4] and neigh3[1] == ['H', 0, 1] and neigh3[2] == ['H', 0, 1]:
                return_boolean = True
        
    #############################################
    # Perform PEO typing for the Oxygen element #
    #############################################
    elif element_type == 'O':
        # PEO Oxygen Backbone C-O-C check that 1st-neigh looks like PEO, then 2nd-neigh lsts looks like PEO
        if len(neigh1) >= 2 and len(neigh2) >= 6:
            if neigh1.count(['C', 0, 4]) == 2 and neigh2.count(['C', 0, 4]) == 2 and neigh2.count(['H', 0, 1]) == 4:
                return_boolean = True
                
        # PEO Oxygen Terminal H-O-C group check that 1st-neigh looks like PEO, then 2nd-neigh lsts looks like PEO
        if len(neigh1) >= 2 and len(neigh2) >= 3:
            if neigh1.count(['C', 0, 4]) == 1 and neigh1.count(['H', 0, 1]) == 1 and neigh2.count(['C', 0, 4]) == 1 and neigh2.count(['H', 0, 1]) == 2:
                return_boolean = True
                
    else:
        raise Exception('check_peo_function does not contain topology checks for elements other then C, H, or O')

    return return_boolean


# Function to update atomtype in supported_types dictionary with flag if flag is needed for typing
def update_supported_types(supported_types, element, atomtype, flag):
    # Find index in list
    typeindex = supported_types[element].index(atomtype)
    
    # Convert flag to string and then set status from 1st character
    status = str(flag)[0]
    
    # Create and insert new type into supported_types lst
    newtype = '{} ({})'.format(atomtype, status)
    supported_types[element][typeindex] = newtype
    return supported_types