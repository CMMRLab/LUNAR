# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 2.4
November 13th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

Will preform bond typing by assigning dissasciation energy to bonds
that are longer then min_bond_length angstroms in length. If a 
dissasciation energy can not be assigned, the coeff will be left as
a harmonic coeff.
"""
##############################
# Import Necessary Libraries #
##############################
import os

# Function to find most frequent occurance in list
def most_frequent(List):
    return max(set(List), key = List.count)

# Function to check if variable is a float
def check_float(variable):
    try: float(variable); return_boolean = True
    except: return_boolean = False
    return return_boolean

# Open and read the file into a dictionary
class Parameters: pass # .d .bo .use 
def read(filename):
    parameters = {} # { (itype, jtype) : Parameters object }; where i/jtype = element,ring,nb
    with open(filename, 'r') as f:   
        parms_flag = False
        for line in f:
            # Strip comment's, white space and split line
            line = line.split('!')[0]; line = line.rstrip();
            line = line.strip(); line = line.split();
            if len(line) > 0 and line[0][0] == '#':
                parms_flag = False
            if '#parameters' in line:
                parms_flag = True
            if parms_flag and len(line) >= 5 and check_float(line[2]) and check_float(line[3]):
                i = line[0]; j = line[1]; use = line[4]; alpha = line[5];
                if check_float(alpha): alpha = float(line[5])
                if use == 'T':
                    P = Parameters()
                    P.d = float(line[2])
                    P.bo = float(line[3])
                    P.alpha = alpha
                    parameters[tuple([i, j])] = P  
    return parameters

# Class to find bond types and set info that will be passed onto alpha parameters code
class topology:
    def __init__(self, m, min_bond_length, coeffs2skip, ff_class, morsefile, class2xe_update, log):
        self.data = {} # {bond coeff id: []}
        self.types = {} # { bond coeff id: string type}
        self.messages = [] # [lst of messages]
        self.effected_r0s = [] # [lst of r0's that have been converted to morse bond]
        self.updated = {i:False for i in m.bond_coeffs}
        
        # elements known to have one bonded atom that are usualy terminal
        elements2skip = ['Br', 'Ca', 'Cl', 'F', 'H'] 
        
        # Read parameters from file
        if os.path.isfile(morsefile):
            parameters = read(morsefile)
            log.out('Typing bonds ....')
            log.out('    Bond typing rules and parameters read from: {}'.format(morsefile))
        else: log.error(f'ERROR bond type rules and Morse bond parameters file: {morsefile} does not exist')
        
        # Find atoms info for each bond coeff
        for i in m.bond_coeffs:
            coeff = m.bond_coeffs[i].coeffs
            if ff_class in [1, '1']: 
                k, r0 = coeff
                self.data[i] = [k, r0] # {bond coeff id: []}
            if ff_class in [2, '2']:
                r0, k2, k3, k4 = coeff
                self.data[i] = [r0, k2, k3, k4] # {bond coeff id: []}
            
            # Initialize flag and types
            found_flag = False; self.types[i] = 'Not Typed (likely not used by any bonds)';
            
            # Try finding bonds with atomids to find bond info
            atomids = [] # intialize as empty and append atomids
            for j in m.bonds:
                bond = m.bonds[j]
                if bond.type == i:
                    atomids.append(bond.atomids)
                    
                    
            # Warn User if no atomids currently use the bond-type
            if not atomids and r0 > min_bond_length and i not in coeffs2skip:
                log.warn('WARNING bond-type {} currently does not have any atomids that use this bond-type so the coeff will get left as class{}'.format(i, ff_class))
                

            # If atomids could be found find bond type and set dissociation energy
            if atomids:
                # Loop through atomids and find each atomids structure info. We will then take the most frequently found set of [[elem1,  ring1,  nb1],
                # [elem2,  ring2,  nb2]], since if conversion comes from ReaxFF each bond type may have slightly different bonding configs thru atomids
                total_info = []; IDs = []
                for id1, id2 in atomids:
                    atom1 = m.atoms[id1];  atom2 = m.atoms[id2];
                    IDs.append([id1, id2])
                    
                    # Find element, ring, and nb for atom1 and atom2
                    info1 = (atom1.element, atom1.ring, atom1.nb) # (elem,  ring,  nb)
                    info2 = (atom2.element, atom2.ring, atom2.nb) # (elem,  ring,  nb)
                    
                    # Find and sort info such that it is order by nb, ring, and then element to order as best as possible.
                    info = sorted(tuple([info1, info2]), key=lambda x: x[2])  
                    info = sorted(info, key=lambda x: x[1])  
                    info = sorted(info, key=lambda x: x[0])
                    total_info.append(tuple(info))
                    
                # Find most frequetly associated bonding pair structure info with bond coeff
                # i to use for setting bond type based on atomids found in the bond coeff
                info1, info2 = most_frequent(total_info)
                
                # Find index of most frequent atom types info to find corresponding atomids 1st neighbors. Example 1st neigh list:
                # [['C', 6, 3], ['C', 6, 3], ['O', 0, 2]] -> [ [neigh1], [neigh2], [neigh3]] -> [neigh1] = [element, ring, nb]
                index_most_frequent = total_info.index(tuple([info1, info2]))
                id1, id2 = atomids[index_most_frequent]
                id1_neigh1 = m.atoms[id1].neighbor_info[1]
                id2_neigh1 = m.atoms[id2].neighbor_info[1]

                
                # Then back out element1, ring1, nb1, element2, ring2, and nb2. This will allow to assume element/ring/nb ordering and reduce
                # the rest of the logic knowing that it is sorted from low# to high for numbers like nb/ring and alphabetized for elements.
                element1 = info1[0]; ring1 = info1[1]; nb1 = info1[2];
                element2 = info2[0]; ring2 = info2[1]; nb2 = info2[2];
                self.types[i] = '{},{},{} - {},{},{}  {}  {}'.format(element1, ring1, nb1, element2, ring2, nb2, 'BO: N/A', 'parms: class'+str(ff_class))
                
                # Generate types and equivs to check against the parameters dict for setting dissociation energy. Possible combinations:
                #    - Permutation1-2:   (itype, jtype)     or (jtype, itype)
                #    - Permutation3-4:   (itype, jequiv1)   or (jequiv1, itype)
                #    - Permutation5-6:   (itype, jequiv2)   or (jequiv2, itype)
                #    - Permutation7-8:   (iequiv1, jtype)   or (jtype, iequiv1)
                #    - Permutation9-10:  (iequiv2, jtype)   or (jtype, iequiv2)
                #    - Permutation11-12: (iequiv1, jequiv1) or (jequiv1, iequiv1)
                #    - Permutation13-14: (iequiv1, jequiv2) or (jequiv1, iequiv2)
                #    - Permutation15-16: (iequiv2, jequiv1) or (jequiv1, iequiv2)
                #    - Permutation17-18: (iequiv2, jequiv2) or (jequiv2, iequiv2)
                #    - Permutation19-20: (iequiv1, jequiv2) or (jequiv2, iequiv1)   
                #    - Permutation21-22: (iequiv1, jequiv3) or (jequiv3, iequiv1)
                #    - Permutation23-24: (iequiv2, jequiv3) or (jequiv3, iequiv2)
                #    - Permutation25-26: (iequiv3, jequiv1) or (jequiv1, iequiv3)
                #    - Permutation27-28: (iequiv3, jequiv2) or (jequiv2, iequiv3)
                #    - Permutation29-30: (iequiv3, jequiv3) or (jequiv3, iequiv3)
                #    - Permutation31-32: (iequiv3, jtype)   or (jtype, iequiv3)
                #    - Permutation33-34: (itype, jequiv3)   or (jequiv3, itype)
                itype = '{},{},{}'.format(element1, ring1, nb1)
                jtype = '{},{},{}'.format(element2, ring2, nb2)
                iequiv1 = '{},{},{}'.format(element1, '*', nb1)
                jequiv1 = '{},{},{}'.format(element2, '*', nb2)
                iequiv2 = '{},{},{}'.format(element1, ring1, '*')
                jequiv2 = '{},{},{}'.format(element2, ring2, '*')
                iequiv3 = '{},{},{}'.format(element1, '*', '*')
                jequiv3 = '{},{},{}'.format(element2, '*', '*')
                
                # If atomids could be found and r0 is large enough and bond coeff is not listed as one
                # to skip over try finding bond type info based on atomids and set dissociation energy
                if r0 > min_bond_length and i not in coeffs2skip or class2xe_update:
                    
                    # Graphene/CNT with pi-electrons check of IFF/IFF-R
                    if itype == 'C,6,5' and jtype == 'C,6,5' and id1_neigh1.count(['C', 6, 5]) == 3 or id2_neigh1.count(['C', 6, 5]) == 3 and ('GRAPHENE', 'GRAPHENE') in parameters:
                        parms = parameters['GRAPHENE', 'GRAPHENE']
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    # Graphene/CNT check
                    elif itype == 'C,6,3' and jtype == 'C,6,3' and  id1_neigh1.count(['C', 6, 3]) == 3 or id2_neigh1.count(['C', 6, 3]) == 3 and ('GRAPHENE', 'GRAPHENE') in parameters:
                        parms = parameters['GRAPHENE', 'GRAPHENE']
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation1-2 check
                    elif (itype, jtype) in parameters:
                        parms = parameters[(itype, jtype)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jtype, itype) in parameters:
                        parms = parameters[(jtype, itype)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation3-4 check
                    elif (itype, jequiv1) in parameters:
                        parms = parameters[(itype, jequiv1)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv1, itype) in parameters:
                        parms = parameters[(jequiv1, itype)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation5-6 check
                    elif (itype, jequiv2) in parameters:
                        parms = parameters[(itype, jequiv2)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv2, itype) in parameters:
                        parms = parameters[(jequiv2, itype)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation7-8 check
                    elif (iequiv1, jtype) in parameters:
                        parms = parameters[(iequiv1, jtype)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jtype, iequiv1) in parameters:
                        parms = parameters[(jtype, iequiv1)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation9-10 check
                    elif (iequiv2, jtype) in parameters:
                        parms = parameters[(iequiv2, jtype)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jtype, iequiv2) in parameters:
                        parms = parameters[(jtype, iequiv2)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation11-12 check
                    elif (iequiv1, jequiv1) in parameters:
                        parms = parameters[(iequiv1, jequiv1)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv1, iequiv1) in parameters:
                        parms = parameters[(jequiv1, iequiv1)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation13-14 check
                    elif (iequiv1, jequiv2) in parameters:
                        parms = parameters[(iequiv1, jequiv2)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv1, iequiv2) in parameters:
                        parms = parameters[(jequiv1, iequiv2)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation15-16 check
                    elif (iequiv2, jequiv1) in parameters:
                        parms = parameters[(iequiv2, jequiv1)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv1, iequiv2) in parameters:
                        parms = parameters[(jequiv1, iequiv2)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation17-18 check
                    elif (iequiv2, jequiv2) in parameters:
                        parms = parameters[(iequiv2, jequiv2)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv2, iequiv2) in parameters:
                        parms = parameters[(jequiv2, iequiv2)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation19-20 check
                    elif (iequiv1, jequiv2) in parameters:
                        parms = parameters[(iequiv1, jequiv2)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv2, iequiv1) in parameters:
                        parms = parameters[(jequiv2, iequiv1)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation21-22 check
                    elif (iequiv1, jequiv3) in parameters:
                        parms = parameters[(iequiv1, jequiv3)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv3, iequiv1) in parameters:
                        parms = parameters[(jequiv3, iequiv1)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation23-24 check
                    elif (iequiv2, jequiv3) in parameters:
                        parms = parameters[(iequiv2, jequiv3)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv3, iequiv2) in parameters:
                        parms = parameters[(jequiv3, iequiv2)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation25-26 check
                    elif (iequiv3, jequiv1) in parameters:
                        parms = parameters[(iequiv3, jequiv1)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv1, iequiv3) in parameters:
                        parms = parameters[(jequiv1, iequiv3)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation27-28 check
                    elif (iequiv3, jequiv2) in parameters:
                        parms = parameters[(iequiv3, jequiv2)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv2, iequiv3) in parameters:
                        parms = parameters[(jequiv2, iequiv3)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation29-30 check
                    elif (iequiv3, jequiv3) in parameters:
                        parms = parameters[(iequiv3, jequiv3)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv3, iequiv3) in parameters:
                        parms = parameters[(jequiv3, iequiv3)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation31-32 check
                    elif (iequiv3, jtype) in parameters:
                        parms = parameters[(iequiv3, jtype)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jtype, iequiv3) in parameters:
                        parms = parameters[(jtype, iequiv3)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN  
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                        
                    # Permutation33-34 check
                    elif (itype, jequiv3) in parameters:
                        parms = parameters[(itype, jequiv3)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                    elif (jequiv3, itype) in parameters:
                        parms = parameters[(jequiv3, itype)]
                        self.data[i] = [r0, parms.d]; found_flag = True;
                        self.types[i] = self.types[i].replace('N/A', str(parms.bo)) # Update BO with same spacing as N/A
                        self.types[i] = self.types[i].replace('class'+str(ff_class), 'file  ') # Update parms with same spacing as classN
                        if check_float(parms.alpha): self.data[i] = [r0, parms.d, parms.alpha] # Add alpha if in file
                
                    
                    #############################################################################
                    # if not found_flag let user know bonding pair not coded yet or not updated #
                    #############################################################################
                    elif not found_flag and element1 not in elements2skip and element2 not in elements2skip:
                        message_part1 = 'Bond Coeff type {} not currently supported for bond typing and dissociation energy definition.\n'.format(i)
                        message_part2 = '{}-{} bond with {}-atom ring: {} nb: {} and {}-atom ring: {} nb: {}\n'.format(element1, element2, element1, ring1, nb1, element2, ring2, nb2)
                        message = message_part1 + message_part2; self.messages.append(message)
                        
                    
                    ##########################################
                    # If found add r0's to self.effected_r0s #
                    ##########################################
                    if found_flag: self.effected_r0s.append(r0); self.updated[i] = True;