#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.6
November 13th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


This set of functions and class for finding the bonds, angles, dihedrals,
impropers, atom types, bond types, angle types, dihedral types, improper
types. Bond types will be ordered in alphabetical order of atom types in
bond, angle/dihedral types will be ordered by the alphabetical order of 
the two outer most atoms, and improper types will be order in alphabetical 
ordering of the 3 outer most atoms and the central atom in the 2nd location
of the improper type.
"""
#############################################
# Class to access all sub routine functions #
#############################################
class BADI:
    def __init__(self, m, nta, name, ff_class, log):
        
        # Find info of ff_class for class 1 or 2 or 'd' or 's1' or 's2'
        if ff_class in [0, 1, 2, 'd', 's1', 's2', '0', '1', '2']:            
            # Find BADI            
            self.graph = generate_graph(m)
            self.bonds = get_bonds(m)
            self.angles = GTangles(self.graph, self.bonds)
            self.dihedrals, self.impropers, self.angleangles = GTdihedrals_impropers_angleangle(self.graph, self.angles)
            
            # Option to reduce DREIDING angle, dihedral, and improper types since outer atoms dont have any meaning in DREIDING (True or False)
            reduce_DREIDING_ADI_types = True
                       
            # Find unique types in list and dictionary format
            self.atom_types_lst, self.atom_types_dict = atom_types(nta, name)
            self.bond_types_lst, self.bond_types_dict = bond_types(self.bonds, nta)
            self.angle_types_lst, self.angle_types_dict = angle_types(self.angles, nta, ff_class, reduce_DREIDING_ADI_types, log)
            self.dihedral_types_lst, self.dihedral_types_dict = dihedral_types(self.dihedrals, nta, ff_class, reduce_DREIDING_ADI_types, log)
            self.improper_types_lst, self.improper_types_dict, self.angleangle_types_lst, self.angleangle_types_dict, self.flagged_angleangles = improper_types(self.impropers, self.angleangles, nta, m, ff_class, reduce_DREIDING_ADI_types, log)
            
            # if ff_class is 2 of 'd' or 's2' add impropers and angleangles lists together
            if ff_class in [2, 'd', 's2', '2']:
                self.impropers.extend(self.angleangles)
                self.improper_types_lst.extend(self.angleangle_types_lst)  
                
        # If ff_class is 'ilmp' for interatomic force fields like reaxFF, REBO, AIREBO, SNAP,...reaxFF find atom types only
        # the 'lmp' ending means the topofile comes from LAMMPS and implies keeping the LAMMPS atomTypeID ordering, instead
        # of resetting all atomTypeIDs when using the 'i' ff_class.
        if ff_class == 'ilmp':
            self.atom_types_lst, self.atom_types_dict = atom_types_from_lmp_datafile(nta, name, m)
            
            # Will need to adjust name dictionary
            for i in m.atoms:
                atom = m.atoms[i]
                name[i] = '{}{}'.format(name[i], atom.type)

        # If ff_class is 'i' for interatomic force fields like reaxFF, REBO, AIREBO, SNAP,...reaxFF find atom types only
        if ff_class == 'i': 
            self.atom_types_lst, self.atom_types_dict = atom_types(nta, name)

            
        
##############################
# Function for finding bonds #
##############################
def get_bonds(m):
    # Set to add bonds to
    bonds = set([])
    for i in m.bonds:
        bonds.add(tuple(sorted(m.bonds[i].atomids)))
    return sorted(bonds)


################################
# Function to generate a graph #
################################
def generate_graph(m):
    # Generate/build graph
    graph = {i:[] for i in m.atoms}
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        graph[id1].append(id2)
        graph[id2].append(id1)
    return graph

#######################################################
# NEW Graph Theory (GT) Based approach to find angles #
#######################################################
def GTangles(graph, bonds):
    # Set to add angles to
    angles = set([])

    # Loop through bonds and graph to generate angles
    for (id1, id2) in bonds:        
        # Loop through id1 adjacency list
        for id3 in graph[id1]:
            id123 = set([id1, id2, id3])
            if len(id123) == 3:
                # Sort angle in ascending order
                if id2 < id3:
                    angle = (id2, id1, id3) # id1 must be next to id3
                else:
                    angle = (id3, id1, id2) # id1 must be next to id3
                    
                # log sorted angle
                if angle not in angles and tuple(reversed(angle)) not in angles:
                    angles.add(angle)
                    
        # Loop through id2 adjacency list
        for id3 in graph[id2]:
            id123 = set([id1, id2, id3])
            if len(id123) == 3:
                # Sort angle in ascending order
                if id1 < id3:
                    angle = (id1, id2, id3) # id2 must by next to id3 
                else:
                    angle = (id3, id2, id1) # id2 must by next to id3
                    
                # log sorted angle
                if angle not in angles and tuple(reversed(angle)) not in angles:
                    angles.add(angle)
    return sorted(angles, key=lambda x: x[1])


#####################################################################################
# NEW Graph Theory (GT) Based approach to find dihedrals, impropers and angleangles #
#####################################################################################
def GTdihedrals_impropers_angleangle(graph, angles):
    # Set to add dihedrals to
    dihedrals = set([])
    
    # Set to add impropers to
    impropers = set([])
    
    # Set to add angleangles to
    angleangles = set([])
    
    # Loop through angles and graph to generate dihedrals, impropers, and angleangle
    for (id1, id2, id3) in angles:
        #--------------------------------------#
        # look for dihedrals using id1 and id3 #
        #--------------------------------------#
        # Loop through id1 adjacency list
        for id4 in graph[id1]:
            id1234 = set([id1, id2, id3, id4])
            if len(id1234) == 4:
                # Sort dihedral in ascending order
                if id3 < id4:
                    dihedral = (id3, id2, id1, id4) # id1 must be next to id4
                else:
                    dihedral = (id4, id1, id2, id3) # id1 must be next to id4
                    
                # log sorted dihedral
                if dihedral not in dihedrals and tuple(reversed(dihedral)) not in dihedrals:
                    dihedrals.add(dihedral)
                    
        # Loop through id3 adjacency list
        for id4 in graph[id3]:
            id1234 = set([id1, id2, id3, id4])
            if len(id1234) == 4:
                # Sort dihedral in ascending order
                if id1 < id4:
                    dihedral = (id1, id2, id3, id4) # id3 must be next to id4
                else:
                    dihedral = (id4, id3, id2, id1) # id3 must be next to id4
                
                # log sorted dihedral
                if dihedral not in dihedrals and tuple(reversed(dihedral)) not in dihedrals:
                    dihedrals.add(dihedral)

                    
                    
        #----------------------------------------------#
        # look for impropers and angleangles using id2 #
        #----------------------------------------------#
        # Loop through id2 adjacency list
        for id4 in graph[id2]:
            id1234 = set([id4, id1, id2, id3])
            if len(id1234) == 4:
                outsides = sorted([id1, id3, id4]) # order outside atoms
                oop = (outsides[0], id2, outsides[1], outsides[2]) # id2 must be center
                
                # If central atom has 3-bonded neighs it is an improper
                if len(graph[id2]) == 3:
                    impropers.add(oop)
                    
                # elif central atom has more then 3-bonded neighs it is an angleangle
                elif len(graph[id2]) > 3:
                    angleangles.add(oop)
  
    return sorted(dihedrals, key=lambda x: x[1]), sorted(impropers, key=lambda x: x[1]), sorted(angleangles, key=lambda x: x[1])


###################################
# Function for finding atom types #
###################################
def atom_types(nta, name):
    # Set to add atom types to
    atom_types_set = set([]); atom_types_dict = {} # { atom type letters : atom type number(s) }
    
    # Find unique atom types
    for i in name:
        atom_types_set.add(name[i])   
        
    # Sort atom types
    atom_types_lst = list(atom_types_set)
    atom_types_lst = sorted(atom_types_lst)

    # Number atom types
    for n, i in enumerate(atom_types_lst):
        atom_types_dict[i] = n + 1
        
    return atom_types_lst, atom_types_dict


####################################################
# Function to find most frequent occurance in list #
####################################################
def most_frequent(List):
    return max(set(List), key = List.count)


####################################################
# Function for finding atom types from ReaxFF type #
# of datafiles, while maintaining atomTypeIDs.     #
####################################################
def atom_types_from_lmp_datafile(nta, name, m):
    # Set to add atom types to
    atom_types_lst = []; atom_types_dict = {} # { atom type letters : atom type number(s) }
    
    # Need to find all currently defined atom types
    atom_types = {i:set() for i in m.masses } # { atomTypeID:set(to put atom type names in)}
    
    # Find atom type symbols and how they map back to the LAMMPS atomTypeID
    for i in name:
        atom = m.atoms[i]
        typeid = atom.type
        symbol = '{}{}'.format(name[i], typeid)
        atom_types[typeid].add(symbol)
        
    # Build inputs
    atom_type_ids = sorted(list(atom_types.keys()))
    for i in atom_type_ids:
        if atom_types[i]: symbol = most_frequent(list(atom_types[i]))
        else: symbol = 'NOTusedBYanyATOMS'
        atom_types_lst.append(symbol)
        
        # This will be over written ... This variable however is a dummy
        # vararible that will not be used in "fill_in_parameter.py" 
        # find_interatomic_atom_parameters() function.
        atom_types_dict[symbol] = i
    
    return atom_types_lst, atom_types_dict
    


###################################
# Function for finding bond types #
###################################
def bond_types(bonds, nta):
    # Set to add bond types to
    bond_types_set = set([]); bond_types_dict = {} # { tuple(bond type letters) : bond type number(s) }
    
    # Find unique bond types
    for id1, id2 in bonds:
        bond_type = sorted([nta[id1], nta[id2]])
        bond_types_set.add(tuple(bond_type))
    
    # Sort all angle types sub tuples by last atom type, 2nd to last
    bond_types_lst = list(bond_types_set)
    bond_types_lst = sorted(bond_types_lst, key=lambda x: x[1])
    bond_types_lst = sorted(bond_types_lst, key=lambda x: x[0])
    
    # Number bond types
    for n, i in enumerate(bond_types_lst):
        bond_types_dict[i] = n + 1
        
    return bond_types_lst, bond_types_dict


####################################
# Function for finding angle types #
####################################
def angle_types(angles, nta, ff_class, reduce_DREIDING_ADI_types, log):
    # Set to add angle types to  (order direction will be based on aplhabetical ordering or the two outer most atoms)
    angle_types_set = set([]); angle_types_dict = {} # { tuple(angle type letters) : angle type number(s) }
    
    # Find unique angle types
    for id1, id2, id3 in angles:
        nta1 = nta[id1]
        nta2 = nta[id2]
        nta3 = nta[id3]
 
        #####################
        # If nta_1 != nta_3 # 
        # use for sorting   #
        #####################
        if nta1 != nta3:
        
            # Find sorted end atoms
            pairing = sorted([nta1, nta3])
            
            # Find index of nta_1 and nta_3 in sorted pairing
            index_nta1 = pairing.index(nta1)
            index_nta3 = pairing.index(nta3)
            
            # If index of nta_1 is zero, nta_1 will be the 1st atom type to use to find angle types
            if index_nta1 == 0:
                angle_type = (nta1, nta2, nta3)
                
            # If index of nta_3 is zero, nta_3 will be the 1st atom type to use to find angle types
            elif index_nta3 == 0:
                angle_type = (nta3, nta2, nta1)
            
            # Let user now if there is a problem   
            else: log.error('ERROR Finding angle types')
                
        #############################
        # Safety just in case outer #
        # atoms are the same        #
        #############################
        else:
            angle_type = (nta1, nta2, nta3)
            
        # ff_class == 'd' and reduce_DREIDING_ADI_types swap outer atoms to X
        if ff_class == 'd' and reduce_DREIDING_ADI_types:
            angle_type = ('X', nta2, 'X')
        
        # Add angle type to angle set
        angle_types_set.add(angle_type)
        
        
    # Sort all angle types sub tuples by last atom type, 2nd to last, and 3rd to last
    angle_types_lst = list(angle_types_set)
    angle_types_lst = sorted(angle_types_lst, key=lambda x: x[2])
    angle_types_lst = sorted(angle_types_lst, key=lambda x: x[1])
    angle_types_lst = sorted(angle_types_lst, key=lambda x: x[0])
    
    # Number angle types
    for n, i in enumerate(angle_types_lst):
        angle_types_dict[i] = n + 1
        
    return angle_types_lst, angle_types_dict


#######################################
# Function for finding dihedral types #
#######################################
def dihedral_types(dihedrals, nta, ff_class, reduce_DREIDING_ADI_types, log):
    # Set to add dihedral types to
    dihedral_types_set = set([]); dihedral_types_dict = {} # { tuple(dihedral type letters) : dihedral type number(s) }
    
    # Find unique dihedral types  (order direction will be based on aplhabetical ordering or the two outer most atoms)
    for id1, id2, id3, id4 in dihedrals:
        nta1 = nta[id1]
        nta2 = nta[id2]
        nta3 = nta[id3]
        nta4 = nta[id4]
        
        # Set dihedral type and attempt sorting twice
        dihedral_type = (nta1, nta2, nta3, nta4)
        
        ##################################
        # Use outer atoms to find order  #
        # dihedral typer ordering if     #
        # outer atom types != each other #
        ##################################
        if nta1 != nta4:
            # Find sorted end atoms
            pairing = sorted([nta1, nta4])
            
            # Find index of nta_1 and nta_4 in sorted pairing
            index_nta1 = pairing.index(nta1)
            index_nta4 = pairing.index(nta4)
            
            
            # If index of nta1 is zero, nta1 will be the 1st atom type to use to find angle types
            if index_nta1 == 0:
                dihedral_type = (nta1, nta2, nta3, nta4)
                
            # If index of nta4 is zero, nta1 will be the 1st atom type to use to find angle types
            elif index_nta4 == 0:
                dihedral_type = (nta4, nta3, nta2, nta1)
            
            # Let user now if there is a problem   
            else: log.error('ERROR Finding dihedral types')
        
        ##################################
        # Use inner atoms to find order  #
        # dihedral typer ordering if     #
        # outer atom types == each other #
        ##################################
        elif nta2 != nta3:
            # Find sorted end atoms
            pairing = sorted([nta2, nta3])
            
            # Find index of nta1 and nta4 in sorted pairing
            index_nta2 = pairing.index(nta2)
            index_nta3 = pairing.index(nta3)
            
            
            # If index of nta_1 is zero, nta_1 will be the 1st atom type to use to find angle types
            if index_nta2 == 0:
                dihedral_type = (nta1, nta2, nta3, nta4)
                
            # If index of nta_4 is zero, nta_1 will be the 1st atom type to use to find angle types
            elif index_nta3 == 0:
                dihedral_type = (nta4, nta3, nta2, nta1)
            
            # Let user now if there is a problem   
            else: log.error('ERROR Finding dihedral types')
        
        #########################################
        # Safety just in case outer atoms are   #
        # the same and inner atoms are the same #
        #########################################
        else:
            dihedral_type = (nta1, nta2, nta3, nta4)
            
        # ff_class == 'd' and reduce_DREIDING_ADI_types swap outer atoms to X
        if ff_class == 'd' and reduce_DREIDING_ADI_types:
            dihedral_type = ('X', nta2, nta3, 'X')
        
        # Add dihedral type to angle set
        dihedral_types_set.add(dihedral_type)
    
    # Sort all dihedral types sub tuples by last atom
    # type, 2nd to last, 3rd to last, and 4th to last
    dihedral_types_lst = list(dihedral_types_set)
    dihedral_types_lst = sorted(dihedral_types_lst, key=lambda x: x[3])
    dihedral_types_lst = sorted(dihedral_types_lst, key=lambda x: x[2])
    dihedral_types_lst = sorted(dihedral_types_lst, key=lambda x: x[1])
    dihedral_types_set = sorted(dihedral_types_lst, key=lambda x: x[0])
    
    # Number dihedral types
    for n, i in enumerate(dihedral_types_lst):
        dihedral_types_dict[i] = n + 1
        
    return dihedral_types_lst, dihedral_types_dict


#######################################
# Function for finding improper types #
#######################################
def improper_types(impropers, angleangles, nta, m, ff_class, reduce_DREIDING_ADI_types, log): 
    # Set to add impropers types to
    improper_types_set = set([]);
    improper_types_dict = {} # { tuple(improper type letters) : improper type number(s) }
    
    # Find unique improper types (order direction will be based on aplhabetical ordering or the two outer most atoms)
    for id1, id2, id3, id4 in impropers:
        nta1 = nta[id1]
        nta2 = nta[id2]
        nta3 = nta[id3]
        nta4 = nta[id4]        

        ############################################################
        # Find sorted outer atoms to find unique ordering of outer #
        # atoms of the improper to help find unique improper types #
        ############################################################
        outer_atoms = sorted([nta1, nta3, nta4])
            
        # Find improper type by sorted outer atoms types
        improper_type = (outer_atoms[0], nta2, outer_atoms[1], outer_atoms[2])
        
        # ff_class == 'd' and reduce_DREIDING_ADI_types swap outer atoms to X
        if ff_class == 'd' and reduce_DREIDING_ADI_types:
            improper_type = (nta2, 'X', 'X', 'X')
        
        # Add improper type to improper set
        improper_types_set.add(improper_type)
        
    # Sort all improper types sub tuples by last atom
    # type, 2nd to last, 3rd to last, and 4th to last
    improper_types_lst = list(improper_types_set)
    improper_types_lst = sorted(improper_types_lst, key=lambda x: x[3])
    improper_types_lst = sorted(improper_types_lst, key=lambda x: x[2])
    improper_types_lst = sorted(improper_types_lst, key=lambda x: x[1])
    improper_types_lst = sorted(improper_types_lst, key=lambda x: x[0])
    
    
    # Set to add angleangles types to
    angleangle_types_set = set([]);
    angleangle_types_dict = {} # { tuple(improper type letters) : improper type number(s) }
    
    # Find unique improper types (order direction will be based on aplhabetical ordering or the two outer most atoms)
    for id1, id2, id3, id4 in angleangles:
        nta1 = nta[id1]
        nta2 = nta[id2]
        nta3 = nta[id3]
        nta4 = nta[id4]        
        
        ############################################################
        # Find sorted outer atoms to find unique ordering of outer #
        # atoms of the improper to help find unique improper types #
        ############################################################
        outer_atoms = sorted([nta1, nta3, nta4])
            
        # Find angleangle type by sorted outer atoms types
        angleangle_type = (outer_atoms[0], nta2, outer_atoms[1], outer_atoms[2])
        
        # ff_class == 'd' and reduce_DREIDING_ADI_types swap outer atoms to X
        if ff_class == 'd' and reduce_DREIDING_ADI_types:
            angleangle_type = (nta2, 'X', 'X', 'X')
        
        # Add improper type to improper set
        angleangle_types_set.add(angleangle_type)
        
    # Sort all improper types sub tuples by last atom
    # type, 2nd to last, 3rd to last, and 4th to last
    angleangle_types_lst = list(angleangle_types_set)
    angleangle_types_lst = sorted(angleangle_types_lst, key=lambda x: x[3])
    angleangle_types_lst = sorted(angleangle_types_lst, key=lambda x: x[2])
    angleangle_types_lst = sorted(angleangle_types_lst, key=lambda x: x[1])
    angleangle_types_lst = sorted(angleangle_types_lst, key=lambda x: x[0])
    
    # Loop through improper_types_lst and warn if improper type and angleangle type are the same
    if ff_class in [2, 'd', 's2', '2']:
        for i in improper_types_lst:
            if i in angleangle_types_lst:
                log.out('OOP Improper Coeff {} {} {} {} has same angleangle coeff - Duplicates of improper Improper coeffs will exists (one for OOP and one for angleangle)'.format(i[0], i[1], i[2], i[3]))
    
    
    # combine impropers and angleangles together based on ff_class and assign numbers
    if ff_class in [2, 'd', 's2', '2']:
        total_impropers = improper_types_lst + angleangle_types_lst 
    elif ff_class in [1, '1']:
        total_impropers = improper_types_lst
    else:
        total_impropers = improper_types_lst
    
    # Number improper/angleangle types and flag angleangle types
    flagged_angleangles = []
    for n, i in enumerate(total_impropers):
        # if n < len(improper_types_lst) it is an improper type set. Create improper coeffs ordering.
        if n < len(improper_types_lst):
            improper_types_dict[i] = n + 1 
        
        # if n >= len(improper_types_lst) flag the number id as an angleangle set. Create angleangle coeffs ordering.
        if n >= len(improper_types_lst):
            angleangle_types_dict[i] = (n + 1)
            flagged_angleangles.append(n + 1)

    return improper_types_lst, improper_types_dict, angleangle_types_lst, angleangle_types_dict, flagged_angleangles