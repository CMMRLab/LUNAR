# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
September 30, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

This file contains useful functions to find 
desired information from mm.atoms[ID] instances
to help in finding atom-types for different
force feilds
"""

def is_sp3_dialkyl_hydrazone_nitrogen(mm, nitrogen_id):
    """
    Recognize the pyramidal/sp3-like dialkylamino nitrogen in a
    dialkylhydrazone environment:

        R2N-N=C(R)R

    Expected first-neighbor environment of the candidate nitrogen:

        two saturated/sp3 carbon atoms
        one imine nitrogen atom

    This is a PCFF atom-type analogy for assigning 'na'; chemically,
    the atom belongs to a hydrazone/hydrazine functional group rather
    than an ordinary tertiary amine.

    Call only after charged nitrogen and other specialized nitrogen
    environments have been excluded.

    Explicit hydrogen atoms are required for the carbon coordination
    check.
    """
    atom = mm.atoms[nitrogen_id]
    neighbor_ids = atom.neighbor_ids[1]

    # Candidate must be a neutral-looking, three-coordinate nitrogen.
    if atom.element != 'N' or int(atom.nb) != 3:
        return False

    # Exclude nitrogen that is itself part of an aromatic ring.
    if check_aromaticity(nitrogen_id, mm.atoms, check_rings=True):
        return False

    carbon_neighbor_ids = [
        neighbor_id
        for neighbor_id in neighbor_ids
        if mm.atoms[neighbor_id].element == 'C'
    ]

    nitrogen_neighbor_ids = [
        neighbor_id
        for neighbor_id in neighbor_ids
        if mm.atoms[neighbor_id].element == 'N'
    ]

    # Dialkylhydrazone pattern: C-N(-C)-N
    if len(carbon_neighbor_ids) != 2:
        return False

    if len(nitrogen_neighbor_ids) != 1:
        return False

    # Both alkyl carbons directly bonded to the candidate N
    # must be saturated/sp3.
    for carbon_id in carbon_neighbor_ids:
        carbon = mm.atoms[carbon_id]
        if int(carbon.nb) != 4:
            return False

    # Examine the N directly bonded to the candidate nitrogen.
    imine_nitrogen_id = nitrogen_neighbor_ids[0]
    imine_nitrogen = mm.atoms[imine_nitrogen_id]

    # In R2N-N=C, the imine N normally has two neighbors:
    # the candidate N and the imine carbon.
    if int(imine_nitrogen.nb) != 2:
        return False

    if check_aromaticity(imine_nitrogen_id, mm.atoms, check_rings=True):
        return False

    imine_nitrogen_other_neighbors = [
        neighbor_id
        for neighbor_id in imine_nitrogen.neighbor_ids[1]
        if neighbor_id != nitrogen_id
    ]

    if len(imine_nitrogen_other_neighbors) != 1:
        return False

    # The other neighbor must be a three-coordinate, nonaromatic
    # carbon consistent with an imine carbon.
    imine_carbon_id = imine_nitrogen_other_neighbors[0]
    imine_carbon = mm.atoms[imine_carbon_id]

    if imine_carbon.element != 'C':
        return False

    if int(imine_carbon.nb) != 3:
        return False

    if check_aromaticity(imine_carbon_id, mm.atoms, check_rings=True):
        return False

    return True


def is_aromatic_amine_nitrogen(mm, nitrogen_id):
    """
    Recognize an aromatic-amine candidate.

    Call only after amide, guanidine, charged-N, and other specialized
    nitrogen environments have been excluded.
    """
    amine_class = amine_candidate_class(mm, nitrogen_id)
    if amine_class not in ('primary', 'secondary', 'tertiary'):
        return False

    # Exclude nitrogen that is itself part of an aromatic ring,
    # such as pyrrole or carbazole nitrogen.
    if check_aromaticity(nitrogen_id, mm.atoms, check_rings=True):
        return False

    # Aromatic amines must be directly bonded to an aromatic carbon.
    for neighbor_id in mm.atoms[nitrogen_id].neighbor_ids[1]:
        neighbor = mm.atoms[neighbor_id]
        is_aromatic = check_aromaticity(neighbor_id, mm.atoms, check_rings=True)
        if neighbor.element == 'C' and is_aromatic:
            return True

    return False


def amine_candidate_class(mm, nitrogen_id):
    """
    Classify a three-coordinate nitrogen having only C/H neighbors.

    This identifies a possible primary, secondary, or tertiary amine.
    Amides, guanidines, aromatic-ring nitrogens, charged nitrogens,
    and other specialized environments must be excluded separately.
    """
    atom = mm.atoms[nitrogen_id]
    nb = int(atom.nb)
    elements1 = neigh_extract(atom, depth=1, info='element') # Example: ['C', 'C', 'H']
    if atom.element != 'N' or nb != 3:
        return ''

    # Conventional amines have only C/H directly bonded to nitrogen.
    if any(element not in {'C', 'H'} for element in elements1):
        return ''

    number_of_carbons = elements1.count('C')
    amine_classes = {1: 'primary', 2: 'secondary', 3: 'tertiary'}
    amine_class = amine_classes.get(number_of_carbons, '')
    return amine_class


def is_sp3_amine_nitrogen(mm, nitrogen_id):
    """
    Recognize an ordinary sp3-amine candidate.

    Call only after charged nitrogen, small-ring nitrogen, and other
    specialized nitrogen environments have been excluded.

    Explicit hydrogen atoms are required.
    """
    amine_class = amine_candidate_class(mm, nitrogen_id)
    if amine_class not in ('primary', 'secondary', 'tertiary'):
        return False

    # Exclude nitrogen that is itself part of an aromatic ring.
    if check_aromaticity(nitrogen_id, mm.atoms, check_rings=True):
        return False

    # Every carbon directly bonded to N must be saturated/sp3.
    for neighbor_id in mm.atoms[nitrogen_id].neighbor_ids[1]:
        neighbor = mm.atoms[neighbor_id]

        if neighbor.element == 'C' and int(neighbor.nb) != 4:
            return False

    return True


# Function to check the atom is aromatic
def check_aromaticity(atomid, atoms, check_rings=False):
    aromaticity = True
    atom        = atoms[atomid]
    rings       = atom.rings
    cycles      = atom.cycles
    if not rings: # early break-out check
        return False
    
    # Optional check for ring sizes to continue checking for aromaticity.
    # For let ring size checks occur out of this function
    if check_rings and not all(isinstance(x, int) and x in (5, 6) for x in rings):
        return False
    
    # Ensure every ring this atom is in that all other atoms only have 2 or 3-nbs
    aromatic_cycles = 0
    nonaromatic_cycles = 0
    for cycle in cycles:
        nbs_lst = [atoms[i].nb for i in cycle]
        if max(nbs_lst) > 3:
            nonaromatic_cycles += 1
        else:
            aromatic_cycles += 1
    
    # If an atom is in both an aromatic and non-aromatic 
    # cycle, let its atom type be assigned to aromatic
    if aromatic_cycles >= 1 and nonaromatic_cycles >= 1:
        aromaticity = True
    elif aromatic_cycles == 0 and nonaromatic_cycles >= 1:
        aromaticity = False

    return aromaticity

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
    # add element if not in supported types
    if element not in supported_types:
        supported_types[element] = [atomtype]
    
    # Find index in list
    typeindex = supported_types[element].index(atomtype)
    
    # Convert flag to string and then set status from 1st character
    status = str(flag)[0]
    
    # Create and insert new type into supported_types lst
    newtype = '{} ({})'.format(atomtype, status)
    supported_types[element][typeindex] = newtype
    return supported_types