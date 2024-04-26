# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
October 28th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to find most frequent occurance in list
def most_frequent(List):
    return max(set(List), key = List.count)

# Function for stringing together parameter types
def string_parameter_type(parameter_type):
    string = ''; str_buffer = 6 # Set string buffer size
    for n, i in enumerate(parameter_type):
        if n == 0: string += '{:<{str_buffer}} '.format(i, str_buffer=str_buffer)
        elif n < len(parameter_type)-1: string += ' {:^{str_buffer}} '.format(i, str_buffer=str_buffer+2)
        else: string += ' {:>{str_buffer}} {:^2}'.format(i, '', str_buffer=str_buffer)
    return string

# Function to add comments to m class to write to datafile later on
def comment(m, cta, rm_unused_coeffs, log):
    
    ################################
    # Build graph to user later on #
    ################################
    # Initialize graph { atomID : [list of 1st neighs]}
    graph = {i:[] for i in m.atoms}
    
    # Add adjacent atoms to graph
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        graph[id1].append(id2)
        graph[id2].append(id1)
        
        
    ##################################################################
    # Add to the following attributes to each atomID of m.atoms[ID]  #
    #     .symbol                                                    # 
    #     .element                                                   #
    # symbol is the atom type and element is the element type.       # 
    # This info will come from cta = { atomID : atomtype/element }   #
    ##################################################################
    for i in m.atoms:
        comment = cta[i].split('/')
        try:
            m.atoms[i].symbol = comment[0] # 1st-index should be atom type
            m.atoms[i].element = comment[1] # 2nd-index should be element
        except: log.error('ERROR cta file not formated properly. Format ID atomtype/element')
    
    
    ###############################################################
    # Find new comments for atom-types/Masses/Pair Coeffs section #
    ###############################################################
    if m.natomtypes > 0:
        # Find all atom type symbols associated with the LAMMPS integer type
        commentID = {} # { TypeID : [lst of comments from atom-types] }
        for i in m.atoms:
            atom = m.atoms[i]
            if atom.type in commentID:
                commentID[atom.type].append(atom.symbol)
            else:
                commentID[atom.type] = [atom.symbol]
                
        # Loop through m.masses and add/update the .type attribute with
        # the most common symbol found in each lst value from commentID
        for i in m.masses:
            try: newtype = most_frequent(commentID[i])
            except: newtype = 'N/A'
            m.masses[i].type = newtype
            m.pair_coeffs[i].type = newtype
            
    ########################################################
    # Find new comments for bond-types/Bond Coeffs section #
    ########################################################
    if m.nbondtypes > 0:
        # Find all atom type symbols associated with the LAMMPS integer type
        commentID = {} # { TypeID : [ (atomTypeID1, atomTypeID2), (), nbonds using TypeID ] }
        for i in m.bonds:
            bond = m.bonds[i]
            
            # Find bonding atomIDs and bond type (Leaving type ordering as ordering of atomIDs)
            id1, id2 = bond.atomids
            TypeID1 = m.atoms[id1].symbol
            TypeID2 = m.atoms[id2].symbol
            bondsymbols = (TypeID1, TypeID2)
            
            # Add bondsymbols tuple to commentID
            if bond.type in commentID:
                commentID[bond.type].append(bondsymbols)
            else:
                commentID[bond.type] = [bondsymbols]
                
        # Loop through m.bond_coeffs and add/update the .type attribute with
        # the most common symbol found in each lst value from commentID
        for i in m.bond_coeffs:
            try: newtype = most_frequent(commentID[i])
            except: newtype = ('N/A', 'N/A')
            m.bond_coeffs[i].type = string_parameter_type(newtype)
            m.bond_coeffs[i].tuple_type = newtype
            
            
    ##########################################################
    # Find new comments for angle-types/Angle Coeffs section #
    ##########################################################
    if m.nangletypes > 0:
        # Find all atom type symbols associated with the LAMMPS integer type
        commentID = {} # { TypeID : [ (atomTypeID1, atomTypeID2, atomTypeID3), (), nangles using TypeID ] }
        for i in m.angles:
            angle = m.angles[i]
            
            # Find bonding atomIDs and bond type (Leaving type ordering as ordering of atomIDs)
            id1, id2, id3 = angle.atomids
            TypeID1 = m.atoms[id1].symbol
            TypeID2 = m.atoms[id2].symbol
            TypeID3 = m.atoms[id3].symbol
            anglesymbols = (TypeID1, TypeID2, TypeID3)
            
            # Add anglesymbols tuple to commentID
            if angle.type in commentID:
                commentID[angle.type].append(anglesymbols)
            else:
                commentID[angle.type] = [anglesymbols]
                
        # Loop through m.angle_coeffs and add/update the .type attribute with
        # the most common symbol found in each lst value from commentID
        for i in m.angle_coeffs:
            try: newtype = most_frequent(commentID[i])
            except: newtype = ('N/A', 'N/A', 'N/A')
            m.angle_coeffs[i].type = string_parameter_type(newtype)
            m.angle_coeffs[i].tuple_type = newtype
                

    ################################################################
    # Find new comments for dihedral-types/Dihedral Coeffs section #
    ################################################################
    if m.ndihedraltypes > 0:
        # Find all atom type symbols associated with the LAMMPS integer type
        commentID = {} # { TypeID : [ (atomTypeID1, atomTypeID2, atomTypeID3, atomTypeID4), (), ndihedrals using TypeID ] }
        for i in m.dihedrals:
            dihedral = m.dihedrals[i]
            
            # Find bonding atomIDs and bond type (Leaving type ordering as ordering of atomIDs)
            id1, id2, id3, id4 = dihedral.atomids
            TypeID1 = m.atoms[id1].symbol
            TypeID2 = m.atoms[id2].symbol
            TypeID3 = m.atoms[id3].symbol
            TypeID4 = m.atoms[id4].symbol
            dihedralsymbols = (TypeID1, TypeID2, TypeID3, TypeID4)
            
            # Add dihedralsymbols tuple to commentID
            if dihedral.type in commentID:
                commentID[dihedral.type].append(dihedralsymbols)
            else:
                commentID[dihedral.type] = [dihedralsymbols]
                
        # Loop through m.dihedral_coeffs and add/update the .type attribute with
        # the most common symbol found in each lst value from commentID
        for i in m.dihedral_coeffs:
            try: newtype = most_frequent(commentID[i])
            except: newtype = ('N/A', 'N/A', 'N/A', 'N/A')
            m.dihedral_coeffs[i].type = string_parameter_type(newtype)
            m.dihedral_coeffs[i].tuple_type = newtype
            

    ################################################################
    # Find new comments for improper-types/Improper Coeffs section #
    ################################################################
    if m.nimpropertypes > 0:
        # Find all atom type symbols associated with the LAMMPS integer type
        nbcentralID = {} # { TypeID : [ nb on central atom, nimpropers using TypeID ] }
        commentID = {} # { TypeID : [ (atomTypeID1, atomTypeID2, atomTypeID3, atomTypeID4), (), nimpropers using TypeID ] }
        for i in m.impropers:
            improper = m.impropers[i]
            
            # Find bonding atomIDs and bond type (Leaving type ordering as ordering of atomIDs)
            id1, id2, id3, id4 = improper.atomids
            TypeID1 = m.atoms[id1].symbol
            TypeID2 = m.atoms[id2].symbol
            TypeID3 = m.atoms[id3].symbol
            TypeID4 = m.atoms[id4].symbol
            impropersymbols = (TypeID1, TypeID2, TypeID3, TypeID4)
            
            # Find nb on central atom (id2) of the defined improper
            nb_on_central = len(graph[id2])
            
            # Add impropersymbols tuple to commentID and nb_on_central to nbcentralID
            if improper.type in commentID:
                commentID[improper.type].append(impropersymbols)
                nbcentralID[improper.type].append(nb_on_central)
            else:
                commentID[improper.type] = [impropersymbols]
                nbcentralID[improper.type] = [nb_on_central]
                
        # Loop through m.improper_coeffs and add/update the .type and .nb 
        # attribute with the most common symbol found in each lst value
        # from commentID and nbcentralID
        for i in m.improper_coeffs:
            try: 
                newtype = most_frequent(commentID[i])
                nb = most_frequent(nbcentralID[i])
            except:
                newtype = ('N/A', 'N/A', 'N/A', 'N/A')
                nb = 0
                
            # Build nb_comment
            if nb == 3: nb_comment = 'nb==3'
            else: nb_comment = 'nb!=3'
            
            # Add to impropers
            full_comment = list(newtype)
            full_comment.append(nb_comment)
            full_comment = tuple(full_comment)
            m.improper_coeffs[i].type = string_parameter_type(full_comment)
            m.improper_coeffs[i].nb = nb_comment
            m.improper_coeffs[i].tuple_type = newtype
            

    ########################################################################
    # Update crossterm .type / .nb (when applicable) from classical coeffs #
    ########################################################################
    # BondBond and BondAngle Coeffs section
    if m.nangletypes > 0:
        try:
            update_crossterms_from_classical(m.angle_coeffs, m.bondbond_coeffs, oop_flag=False)
            update_crossterms_from_classical(m.angle_coeffs, m.bondangle_coeffs, oop_flag=False)
        except: pass
        
    # AngleAngleTorsion, EndBondTorsion, MiddleBondTorsion, BondBond13, and AngleTorsion Coeffs section
    if m.ndihedraltypes > 0:
        try:
            update_crossterms_from_classical(m.dihedral_coeffs, m.angleangletorsion_coeffs, oop_flag=False)
            update_crossterms_from_classical(m.dihedral_coeffs, m.endbondtorsion_coeffs, oop_flag=False)
            update_crossterms_from_classical(m.dihedral_coeffs, m.middlebondtorsion_coeffs, oop_flag=False)
            update_crossterms_from_classical(m.dihedral_coeffs, m.bondbond13_coeffs, oop_flag=False)
            update_crossterms_from_classical(m.dihedral_coeffs, m.angletorsion_coeffs, oop_flag=False)
        except: pass
           
    # AngleAngle Coeffs section
    if m.nimpropertypes > 0:
        try : update_crossterms_from_classical(m.improper_coeffs, m.angleangle_coeffs, oop_flag=True)
        except: pass
        
    
    #########################################################
    # if rm_unused_coeffs find new coeffIDs and coeffID map #
    #########################################################
    if rm_unused_coeffs:
        log.out('\n\nRemoving any unused CoeffTypeIDs:')
        
        #--------------#
        # Atom TypeIDs #
        #--------------#
        # Find new atomID mapping
        new_atomIDs = []
        for i in m.masses:
            if m.masses[i].type != 'N/A':
                new_atomIDs.append(i)
        new_atomIDs = sorted(new_atomIDs); norg_atomIDs = len(m.masses);
        atom_typeID_map = {ID:n for n, ID in enumerate(new_atomIDs, 1)} # {oldID : newID}
        
        # Apply new atomTypeID mapping to masses, pair_coeffs, and atoms
        masses = {}; pair_coeffs = {};
        for i in new_atomIDs:
            ID = atom_typeID_map[i]
            masses[ID] = m.masses[i]
            pair_coeffs[ID] = m.pair_coeffs[i]
        m.masses = masses; m.pair_coeffs = pair_coeffs; m.natomtypes = len(new_atomIDs);
        for i in m.atoms:
            atom = m.atoms[i]
            typeID = atom_typeID_map[atom.type]
            atom.type = typeID
        log.out('{:>25}, {:>25} = {:>5}, {:>5}'.format('N-orginal atomTypeIDs', 'N-new atomTypeIDs', norg_atomIDs, len(new_atomIDs)))
        
        #--------------#
        # Bond TypeIDs #
        #--------------#
        # Find new bondID mapping
        new_bondIDs = []
        for i in m.bond_coeffs:
            if m.bond_coeffs[i].tuple_type != ('N/A', 'N/A'):
                new_bondIDs.append(i)
        new_bondIDs = sorted(new_bondIDs); norg_bondIDs = len(m.bond_coeffs);
        bond_typeID_map = {ID:n for n, ID in enumerate(new_bondIDs, 1)} # {oldID : newID}
        
        # Apply new bondTypeID mapping to bond_coeffs and bonds
        bond_coeffs = {};
        for i in new_bondIDs:
            ID = bond_typeID_map[i]
            bond_coeffs[ID] = m.bond_coeffs[i]
        m.bond_coeffs = bond_coeffs; m.nbondtypes = len(new_bondIDs);
        for i in m.bonds:
            bond = m.bonds[i]
            typeID = bond_typeID_map[bond.type]
            bond.type = typeID
        log.out('{:>25}, {:>25} = {:>5}, {:>5}'.format('N-orignal bondTypeIDs', 'N-new bondTypeIDs', norg_bondIDs, len(new_bondIDs)))
        
        #---------------#
        # Angle TypeIDs #
        #---------------#
        # Find new angleID mapping
        new_angleIDs = []
        for i in m.angle_coeffs:
            if m.angle_coeffs[i].tuple_type != ('N/A', 'N/A', 'N/A'):
                new_angleIDs.append(i)
        new_angleIDs = sorted(new_angleIDs); norg_angleIDs = len(m.angle_coeffs);
        angle_typeID_map = {ID:n for n, ID in enumerate(new_angleIDs, 1)} # {oldID : newID}
        
        # Apply new angleTypeID mapping to angle_coeffs, bondbond_coeffs, bondangle_coeffs and angles
        angle_coeffs = {}; bondbond_coeffs = {}; bondangle_coeffs = {};
        for i in new_angleIDs:
            ID = angle_typeID_map[i]
            angle_coeffs[ID] = m.angle_coeffs[i]
            bondbond_coeffs[ID] = m.bondbond_coeffs[i]
            bondangle_coeffs[ID] = m.bondangle_coeffs[i]
        m.angle_coeffs = angle_coeffs; m.bondbond_coeffs = bondbond_coeffs;
        m.bondangle_coeffs = bondangle_coeffs; m.nangletypes = len(new_angleIDs);
        for i in m.angles:
            angle = m.angles[i]
            typeID = angle_typeID_map[angle.type]
            angle.type = typeID
        log.out('{:>25}, {:>25} = {:>5}, {:>5}'.format('N-orignal angleTypeIDs', 'N-new angleTypeIDs', norg_angleIDs, len(new_angleIDs)))
        
        #------------------#
        # Dihedral TypeIDs #
        #------------------#
        # Find new dihedralID mapping
        new_dihedralIDs = []
        for i in m.dihedral_coeffs:
            if m.dihedral_coeffs[i].tuple_type != ('N/A', 'N/A', 'N/A', 'N/A'):
                new_dihedralIDs.append(i)
        dihedralIDs = sorted(new_dihedralIDs); norg_dihedralIDs = len(m.dihedral_coeffs);
        dihedral_typeID_map = {ID:n for n, ID in enumerate(dihedralIDs, 1)} # {oldID : newID}
        
        # Apply new dihedralTypeID mapping to dihedral_coeffs, angleangletorsion_coeffs,
        # endbondtorsion_coeffs, middlebondtorsion_coeffs, bondbond13_coeffs, 
        # angletorsion_coeffs, and dihedrals
        dihedral_coeffs = {}; angleangletorsion_coeffs = {}; endbondtorsion_coeffs = {};
        middlebondtorsion_coeffs = {}; bondbond13_coeffs = {}; angletorsion_coeffs = {};
        for i in new_dihedralIDs:
            ID = dihedral_typeID_map[i]
            dihedral_coeffs[ID] = m.dihedral_coeffs[i]
            try:
                angleangletorsion_coeffs[ID] = m.angleangletorsion_coeffs[i]
                endbondtorsion_coeffs[ID] = m.endbondtorsion_coeffs[i]
                middlebondtorsion_coeffs[ID] = m.middlebondtorsion_coeffs[i]
                bondbond13_coeffs[ID] = m.bondbond13_coeffs[i]
                angletorsion_coeffs[ID] = m.angletorsion_coeffs[i]
            except: pass
        m.dihedral_coeffs = dihedral_coeffs; m.angleangletorsion_coeffs = angleangletorsion_coeffs;
        m.endbondtorsion_coeffs = endbondtorsion_coeffs; m.middlebondtorsion_coeffs = middlebondtorsion_coeffs;
        m.bondbond13_coeffs = bondbond13_coeffs; m.angletorsion_coeffs = angletorsion_coeffs;
        m.ndihedraltypes = len(new_dihedralIDs);
        for i in m.dihedrals:
            dihedral = m.dihedrals[i]
            typeID = dihedral_typeID_map[dihedral.type]
            dihedral.type = typeID
        log.out('{:>25}, {:>25} = {:>5}, {:>5}'.format('N-orignal dihedralTypeIDs', 'N-new dihedralTypeIDs', norg_dihedralIDs, len(new_dihedralIDs)))
        
        
        #------------------#
        # Improper TypeIDs #
        #------------------#
        # Find new improperID mapping
        new_improperIDs = []
        for i in m.improper_coeffs:
            if m.improper_coeffs[i].tuple_type != ('N/A', 'N/A', 'N/A', 'N/A'):
                new_improperIDs.append(i)
        improperIDs = sorted(new_improperIDs); norg_improperIDs = len(m.improper_coeffs);
        improper_typeID_map = {ID:n for n, ID in enumerate(improperIDs, 1)} # {oldID : newID}
        
        # Apply new improperTypeID mapping to improper_coeffs and angleangle_coeffs
        improper_coeffs = {}; angleangle_coeffs = {};
        for i in new_improperIDs:
            ID = improper_typeID_map[i]
            improper_coeffs[ID] = m.improper_coeffs[i]
            try: angleangle_coeffs[ID] = m.angleangle_coeffs[i]
            except: pass
        m.improper_coeffs = improper_coeffs; m.angleangle_coeffs = angleangle_coeffs;
        m.nimpropertypes = len(new_improperIDs);
        for i in m.impropers:
            improper = m.impropers[i]
            typeID = improper_typeID_map[improper.type]
            improper.type = typeID
        log.out('{:>25}, {:>25} = {:>5}, {:>5}\n\n'.format('N-orignal improperTypeIDs', 'N-new improperTypeIDs', norg_improperIDs, len(new_improperIDs)))
    return m



############################################################
# Function to update crossterms info from classical coeffs #
############################################################
def update_crossterms_from_classical(classical_coeffs, crossterm_coeffs, oop_flag):
    for i in classical_coeffs:
        # Update crossterm coeff info
        crossterm_coeffs[i].type = classical_coeffs[i].type    
        
        # Update oop_flag update nb as well
        if oop_flag:
            crossterm_coeffs[i].nb = classical_coeffs[i].nb  
    return