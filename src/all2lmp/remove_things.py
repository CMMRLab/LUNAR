# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
October 11th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to modify parameters based on info found in the remove dictionary
def post_processor(parameters, remove, ff_class, BADI, log):
    #--------------------------------#
    # Print to user the .nta options #
    #--------------------------------#
    log.out('\n\nRemoving topology and parameters (based on .nta file settings)')
    for i in remove:
        if remove[i]: log.debug('{} = {}'.format(i, remove[i]))
            
    #---------------#
    # Angle section #
    #---------------#
    # Find which angle TypeIDs to remove based on info from remove dict
    angleCoeffIDs2remove = set(); rm_angle = {'types':0, 'angle':0}
    if remove['zero']['angle']:
        for i in parameters.angle_coeffs:
            coeffs = parameters.angle_coeffs[i].coeffs
            if all(p == 0 for p in coeffs):
                angleCoeffIDs2remove.add(i)
    rm_nta = remove['angle-nta']
    for i in parameters.angles:
        angle = parameters.angles[i]
        forward = angle.symbol; reverse = tuple(reversed(forward));
        wildcardf = ('*', forward[1], '*'); wildcardr = ('*', reverse[1], '*');
        if forward in rm_nta or reverse in rm_nta or wildcardf in rm_nta or wildcardr in rm_nta:
            angleCoeffIDs2remove.add(angle.type)
    rm_angle['types'] = len(angleCoeffIDs2remove)
    
    # Reset angles, angle coeffs and crossterm coeffs
    angles = {}; angle_coeffs = {}; bondbond_coeffs = {}; bondangle_coeffs = {};
    coeffIDmap = generate_coeffID_map(parameters.angle_coeffs, angleCoeffIDs2remove); ID = 0;
    rm_id = remove['angle-ID']
    for i in parameters.angles:
        forward = tuple(parameters.angles[i].atomids); reverse = tuple(reversed(forward));
        wildcardf = ('*', forward[1], '*'); wildcardr = ('*', reverse[1], '*');
        if forward in rm_id or reverse in rm_id: rm_angle['angle'] += 1; continue
        if wildcardf in rm_id or wildcardr in rm_id: rm_angle['angle'] += 1; continue
        if parameters.angles[i].type in coeffIDmap:
            ID += 1
            angles[ID] = parameters.angles[i]
            angles[ID].type = coeffIDmap[parameters.angles[i].type]
        else: rm_angle['angle'] += 1
    for i in parameters.angle_coeffs:
        if i in coeffIDmap:
            angle_coeffs[coeffIDmap[i]] = parameters.angle_coeffs[i]
            if ff_class == 2:
                bondbond_coeffs[coeffIDmap[i]] = parameters.bondbond_coeffs[i]
                bondangle_coeffs[coeffIDmap[i]] = parameters.bondangle_coeffs[i]
                
    # Update parameters object dicts w/new reduced dicts
    parameters.angles = angles
    parameters.angle_coeffs = angle_coeffs
    parameters.bondbond_coeffs = bondbond_coeffs
    parameters.bondangle_coeffs = bondangle_coeffs
    parameters.nangles = len(angles)
    parameters.nangletypes = len(angle_coeffs)


    #------------------#
    # Dihedral section #
    #------------------#    
    # Find which dihedral TypeIDs to remove based on info from remove dict
    dihedralCoeffIDs2remove = set(); rm_dihedral = {'types':0, 'dihedral':0}
    if remove['zero']['dihedral']:
        for i in parameters.dihedral_coeffs:
            coeffs = parameters.dihedral_coeffs[i].coeffs
            if all(p == 0 for p in coeffs):
                dihedralCoeffIDs2remove.add(i)
    rm_nta = remove['dihedral-nta']
    for i in parameters.dihedrals:
        dihedral = parameters.dihedrals[i]
        forward = dihedral.symbol; reverse = tuple(reversed(forward));
        wildcardf = ('*', forward[1], forward[2], '*'); wildcardr = ('*', reverse[1], reverse[2], '*');
        if forward in rm_nta or reverse in rm_nta or wildcardf in rm_nta or wildcardr in rm_nta:
            dihedralCoeffIDs2remove.add(dihedral.type)  
    rm_dihedral['types'] = len(dihedralCoeffIDs2remove)
    
    # Reset dihedrals, dihedral coeffs and crossterm coeffs
    dihedrals = {}; dihedral_coeffs = {}; bondbond13_coeffs = {};
    angleangletorsion_coeffs = {}; endbondtorsion_coeffs = {};
    middlebondtorsion_coeffs = {}; angletorsion_coeffs = {}; 
    coeffIDmap = generate_coeffID_map(parameters.dihedral_coeffs, dihedralCoeffIDs2remove); ID = 0;
    rm_id = remove['dihedral-ID']
    for i in parameters.dihedrals:
        forward = tuple(parameters.dihedrals[i].atomids); reverse = tuple(reversed(forward));
        wildcardf = ('*', forward[1], forward[2], '*'); wildcardr = ('*', reverse[1], reverse[2], '*');
        if forward in rm_id or reverse in rm_id: rm_dihedral['dihedral'] += 1; continue
        if wildcardf in rm_id or wildcardr in rm_id: rm_dihedral['dihedral'] += 1; continue
        if parameters.dihedrals[i].type in coeffIDmap:
            ID += 1
            dihedrals[ID] = parameters.dihedrals[i]
            dihedrals[ID].type = coeffIDmap[parameters.dihedrals[i].type]
        else: rm_dihedral['dihedral'] += 1
    for i in parameters.dihedral_coeffs:
        if i in coeffIDmap:
            dihedral_coeffs[coeffIDmap[i]] = parameters.dihedral_coeffs[i]
            if ff_class == 2:
                bondbond13_coeffs[coeffIDmap[i]] = parameters.bondbond13_coeffs[i]
                angleangletorsion_coeffs[coeffIDmap[i]] = parameters.angleangletorsion_coeffs[i]
                endbondtorsion_coeffs[coeffIDmap[i]] = parameters.endbondtorsion_coeffs[i]
                middlebondtorsion_coeffs[coeffIDmap[i]] = parameters.middlebondtorsion_coeffs[i]
                angletorsion_coeffs[coeffIDmap[i]] = parameters.angletorsion_coeffs[i]
    
    # Update parameters object dicts w/new reduced dicts
    parameters.dihedrals = dihedrals
    parameters.dihedral_coeffs = dihedral_coeffs
    parameters.bondbond13_coeffs = bondbond13_coeffs
    parameters.angleangletorsion_coeffs = angleangletorsion_coeffs
    parameters.endbondtorsion_coeffs = endbondtorsion_coeffs
    parameters.middlebondtorsion_coeffs = middlebondtorsion_coeffs
    parameters.angletorsion_coeffs = angletorsion_coeffs
    parameters.ndihedrals = len(dihedrals)
    parameters.ndihedraltypes = len(dihedral_coeffs)
    
    
    #------------------#
    # Improper section #
    #------------------#
    # Find which improper TypeIDs to remove based on info from remove dict
    improperCoeffIDs2remove = set(); rm_improper = {'types':0, 'improper':0}
    if remove['zero']['improper']:
        for i in parameters.improper_coeffs:
            coeffs = parameters.improper_coeffs[i].coeffs
            if all(p == 0 for p in coeffs):
                improperCoeffIDs2remove.add(i)
    rm_nta = remove['improper-nta']
    flagged_angleangles = BADI.flagged_angleangles
    for i in parameters.impropers:
        improper = parameters.impropers[i]
        if improper.type in flagged_angleangles:
            nb = 'nb!=3'
        else: nb = 'nb==3'
        sym = improper.symbol;
        const12_order34 = (sym[0], sym[1], sym[2], sym[3], nb)
        const12_order43 = (sym[0], sym[1], sym[3], sym[2], nb)
        const32_order14 = (sym[2], sym[1], sym[0], sym[3], nb)
        const32_order41 = (sym[2], sym[1], sym[3], sym[0], nb)
        const42_order13 = (sym[3], sym[1], sym[0], sym[2], nb)
        const42_order31 = (sym[3], sym[1], sym[2], sym[1], nb)
        wildard_nb_eql3 = ('*', sym[1], '*', '*', 'nb==3')
        wildard_nb_grt3 = ('*', sym[1], '*', '*', 'nb!=3')
        if const12_order34 in rm_nta or const12_order43 in rm_nta or const32_order14 in rm_nta or const32_order41 in rm_nta or const42_order13 in rm_nta or const42_order31 in rm_nta or wildard_nb_eql3 in rm_nta or wildard_nb_grt3 in rm_nta:
            improperCoeffIDs2remove.add(improper.type)
    rm_improper['types'] = len(improperCoeffIDs2remove)
    
    # Reset impropers, improper coeffs and crossterm coeffs
    impropers = {}; improper_coeffs = {}; angleangle_coeffs = {};
    coeffIDmap = generate_coeffID_map(parameters.improper_coeffs, improperCoeffIDs2remove); ID = 0;
    rm_id = remove['improper-ID']
    for i in parameters.impropers:
        sym = parameters.impropers[i].atomids;
        const12_order34 = (sym[0], sym[1], sym[2], sym[3], nb)
        const12_order43 = (sym[0], sym[1], sym[3], sym[2], nb)
        const32_order14 = (sym[2], sym[1], sym[0], sym[3], nb)
        const32_order41 = (sym[2], sym[1], sym[3], sym[0], nb)
        const42_order13 = (sym[3], sym[1], sym[0], sym[2], nb)
        const42_order31 = (sym[3], sym[1], sym[2], sym[1], nb)
        wildard_nb_eql3 = ('*', sym[1], '*', '*', 'nb==3')
        wildard_nb_grt3 = ('*', sym[1], '*', '*', 'nb!=3')
        if const12_order34 in rm_id or const12_order43 in rm_id: rm_improper['improper'] += 1; continue
        if const32_order14 in rm_id or const32_order41 in rm_id: rm_improper['improper'] += 1; continue
        if const42_order13 in rm_id or const42_order31 in rm_id: rm_improper['improper'] += 1; continue
        if wildard_nb_eql3 in rm_id or wildard_nb_grt3 in rm_id: rm_improper['improper'] += 1; continue
        if parameters.impropers[i].type in coeffIDmap:
            ID += 1
            impropers[ID] = parameters.impropers[i]
            impropers[ID].type = coeffIDmap[parameters.impropers[i].type]
        else: rm_improper['improper'] += 1
    for i in parameters.improper_coeffs:
        if i in coeffIDmap:
            improper_coeffs[coeffIDmap[i]] = parameters.improper_coeffs[i]
            if ff_class == 2:
                angleangle_coeffs[coeffIDmap[i]] = parameters.angleangle_coeffs[i]

    # Update parameters object dicts w/new reduced dicts
    parameters.impropers = impropers
    parameters.improper_coeffs = improper_coeffs
    parameters.angleangle_coeffs = angleangle_coeffs
    parameters.nimpropers = len(impropers)
    parameters.nimpropertypes = len(improper_coeffs)
    
    # Print updates of how many things were remove
    log.out('   Number of removed angles, types = {:^10} {:^5}'.format(rm_angle['angle'], rm_angle['types']))
    log.out('   Number of removed dihedrals, types = {:^10} {:^5}'.format(rm_dihedral['dihedral'], rm_dihedral['types']))
    log.out('   Number of removed impropers, types = {:^10} {:^5}'.format(rm_improper['improper'], rm_improper['types']))
    log.out('')
    return parameters


# Function to generate coeffID map
def generate_coeffID_map(coeff_dict, CoeffIDs2remove):
    newID = 0;  coeffIDmap = {} # { oldID : newID }
    for i in sorted(list(coeff_dict.keys())):
        if i not in CoeffIDs2remove:
            newID += 1
            coeffIDmap[i] = newID
    return coeffIDmap          