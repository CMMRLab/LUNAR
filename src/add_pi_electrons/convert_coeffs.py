# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
April 10th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function for stringing together parameter types
def string_parameter_type(parameter_type):
    string = ''; str_buffer = 6 # Set string buffer size
    for n, i in enumerate(parameter_type):
        if n == 0: string += '{:<{str_buffer}} '.format(i, str_buffer=str_buffer)
        elif n < len(parameter_type)-1: string += ' {:^{str_buffer}} '.format(i, str_buffer=str_buffer+2)
        else: string += ' {:>{str_buffer}} {:^2}'.format(i, '', str_buffer=str_buffer)
    return string

# Function to update dictionary of coeffs
def update_coeff_dict(dict2update, typeID, coeffs, typestring):
    coeff = Coeff_class()
    coeff.coeffs = coeffs
    coeff.type = typestring
    dict2update.update({typeID:coeff})
    return dict2update


# Function to convert coeffs to cg1 coeffs
class Coeff_class:
    pass  # .type .coeffs = []
def cg1(m, types2convert, graph):
    
    ################################
    # Find pair coeffs that are cp #
    # cp  = [0.064, 4.01]          #
    # cg1 = [0.0620, 3.9320]       #
    ################################
    # Finding coeff TypeID to convert
    pair_types = set()
    for i in m.pair_coeffs:
        if i in types2convert:
            pair_types.add(i)
            
    # Update pair coeffs
    for i in pair_types:
        m.pair_coeffs = update_coeff_dict(m.pair_coeffs, i, [0.06200, 3.93200],  'cp->cg1') # cg1 pair
        m.masses = update_coeff_dict(m.masses, i, m.masses[i].coeffs,  'cp->cg1') # cg1 masses
        
    ##################################################
    # Find bond coeffs that are cp-cp                #
    # cp  = [1.4170, 470.8361, -627.6179, 1327.6345] #
    # cg1 = [1.4200, 480.0000, 0.0000, 0.0000]       #
    ##################################################
    bond_types = set()
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids

        # Finding type from atom ids
        type1 = m.atoms[id1].type 
        type2 = m.atoms[id2].type
        
        # Finding coeff number that has the cp-cp bond coeff
        if type1 in types2convert and type2 in types2convert:
            bond_types.add(m.bonds[i].type)
            
    # Update bond coeffs
    for i in bond_types:
        coeff = [1.4200, 480.0000, 0.0000, 0.0000]
        typestring = string_parameter_type(('cg1', 'cg1'))
        m.bond_coeffs = update_coeff_dict(m.bond_coeffs, i, coeff, typestring) # cg1-cg1 bond
        
    ################################################
    # Find angle coeffs that are cp-cp-cp          #
    # cp  = [118.9000, 61.0226, -34.9931,  0.0000] #
    # cg1 = [120.0000, 90.0000, 0.0000, 0.0000]    #
    ################################################
    angle_types = set()
    for i in m.angles:
        id1, id2, id3 = m.angles[i].atomids
        
        # Finding type from atom ids
        type1 = m.atoms[id1].type 
        type2 = m.atoms[id2].type
        type3 = m.atoms[id3].type
        
        # Finding coeff number that has the cp-cp-cp angle coeff
        if type1 in types2convert and type2 in types2convert and type3 in types2convert:
            angle_types.add(m.angles[i].type)
            
    # Update angle, bondbond, and bondangle coeffs
    for i in angle_types:
        coeff = [120.0000, 90.0000, 0.0000, 0.0000] 
        typestring = string_parameter_type(('cg1', 'cg1', 'cg1'))
        m.angle_coeffs = update_coeff_dict(m.angle_coeffs, i, coeff, typestring) # cg1-cg1-cg1 angle
        m.bondbond_coeffs = update_coeff_dict(m.bondbond_coeffs, i, [0.0, 0.0, 0.0], typestring) # cg1-cg1-cg1 bondbond
        m.bondangle_coeffs = update_coeff_dict(m.bondangle_coeffs, i, [0.0, 0.0, 0.0, 0.0], typestring) # cg1-cg1-cg1 bondangle
        
    #################################################
    # Find dihedral coeffs that are cp-cp-cp-cp     #
    # cp  = [8.3667, 0.0, 1.1932, 0.0, 0.0000, 0.0] #
    # cg1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]          #
    #################################################
    dihedral_types = set()
    for i in m.dihedrals:
        id1, id2, id3, id4 = m.dihedrals[i].atomids
        
        # Finding type from atom ids
        type1 = m.atoms[id1].type 
        type2 = m.atoms[id2].type
        type3 = m.atoms[id3].type
        type4 = m.atoms[id4].type
        
        # Finding coeff number that has the cp-cp-cp-cp dihedral coeff
        if type1 in types2convert and type2 in types2convert and type3 in types2convert and type4 in types2convert:
            dihedral_types.add(m.dihedrals[i].type)
            
    # Update dihedral, AngleAngleTorsion, EndBondTorsion, MiddleBondTorsion, BondBond13, and AngleTorsion coeffs
    for i in dihedral_types:
        coeff = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  
        typestring = string_parameter_type(('cg1', 'cg1', 'cg1', 'cg1'))
        m.dihedral_coeffs = update_coeff_dict(m.dihedral_coeffs, i, coeff, typestring) # cg1-cg1-cg1-cg1 dihedral
        m.angleangletorsion_coeffs = update_coeff_dict(m.angleangletorsion_coeffs, i, [0.0, 0.0, 0.0], typestring) # cg1-cg1-cg1-cg1 angleangletorsion
        m.endbondtorsion_coeffs = update_coeff_dict(m.endbondtorsion_coeffs, i, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], typestring) # cg1-cg1-cg1-cg1 endbondtorsion
        m.middlebondtorsion_coeffs = update_coeff_dict(m.middlebondtorsion_coeffs, i, [0.0, 0.0, 0.0, 0.0], typestring) # cg1-cg1-cg1-cg1 middlebondtorsion
        m.bondbond13_coeffs = update_coeff_dict(m.bondbond13_coeffs, i, [0.0, 0.0, 0.0], typestring) # cg1-cg1-cg1-cg1 bondbond13
        m.angletorsion_coeffs = update_coeff_dict(m.angletorsion_coeffs, i, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], typestring) # cg1-cg1-cg1-cg1 angletorsion
        
        
    #################################################
    # Find improper coeffs that are cp-cp-cp-cp     #
    # cp  = [7.1794, 0.0000]                        #
    # cg1 = [0.0, 0.0]                              #
    #################################################
    improper_types = set(); angleangle_types = set();
    for i in m.impropers:
        id1, id2, id3, id4 = m.impropers[i].atomids
        
        # Finding type from atom ids
        type1 = m.atoms[id1].type 
        type2 = m.atoms[id2].type
        type3 = m.atoms[id3].type
        type4 = m.atoms[id4].type
        
        # Finding coeff number that has the cp-cp-cp-cp dihedral coeff
        if type1 in types2convert and type2 in types2convert and type3 in types2convert and type4 in types2convert:
            if len(graph[id2]) == 3:
                improper_types.add(m.impropers[i].type)
            else:
                angleangle_types.add(m.impropers[i].type)
    
    # Group improper and angleange types together
    oop = list(improper_types) + list(angleangle_types)
    
    # Update improper and angleangle coeffs
    for i in oop:
        if i in improper_types:
            typestring = string_parameter_type(('cg1', 'cg1', 'cg1', 'cg1', 'nb==3'))
        else:
            typestring = string_parameter_type(('cg1', 'cg1', 'cg1', 'cg1', 'nb!=3'))
        m.improper_coeffs = update_coeff_dict(m.improper_coeffs, i, [0.0, 0.0], typestring) # cg1-cg1-cg1-cg1 improper
        m.angleangle_coeffs = update_coeff_dict(m.angleangle_coeffs, i, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], typestring) # cg1-cg1-cg1-cg1 angleangle            
    return m