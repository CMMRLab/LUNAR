# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
December 1st, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

ff_functions is a secondary script to fill_in_parameters that
contains various functions that fill_in_parameters will call
and use to match atom-types in coeff and find corresponding
parameters. ff_functions also contains other miscellaneous
functions that all2lmp will call and there was no logical
place to stick the function but in this script.
"""
##############################
# Import Necessary Libraries #
##############################
import math
import os


#####################
# Rotational matrix #
#####################
cos = math.cos
sin = math.sin

def Rx(phi):
    phi = math.radians(phi)
    return [[ 1, 0,        0        ],
            [ 0, cos(phi), -sin(phi)],
            [ 0, sin(phi), cos(phi)]]
  
def Ry(theta):
    theta = math.radians(theta)
    return [[ cos(theta), 0, sin(theta)],
            [ 0,          1, 0         ],
            [-sin(theta), 0, cos(theta)]]
  
def Rz(psi):
    psi = math.radians(psi)
    return [[ cos(psi),  -sin(psi), 0 ],
            [ sin(psi),  cos(psi) , 0 ],
            [ 0,         0,         1 ]]

# Function to multiply two matrices together
def matrix_by_matrix(m1, m2):
    # Intialize matrix
    ncolumns = max( [len(m1[0]), len(m2[0])] ) 
    nrows = max( [len(m1), len(m2)] )
    result = [ncolumns*[0] for row in range(nrows)]
    
    # Iterate through rows of m1
    for i in range(len(m1)):
       # Iterate through columns of m2
       for j in range(len(m2[0])):
           # Iterate through rows of m2
           for k in range(len(m2)):
               result[i][j] += m1[i][k] * m2[k][j]
    return result

# Function to multiply a vector and a matrix together
def vector_by_matrix(m1, v1):
    return [sum([x*y for x, y in zip(v1, v2)]) for v2 in m1]


#############################################################
# Function to rotate a molecule at centered about (0, 0, 0) #
#############################################################
def rotate_system(m, phi, theta, psi):
    # Update simulation cell
    xline = m.xbox_line.split(); yline = m.ybox_line.split(); zline = m.zbox_line.split();
    xlo = float(xline[0]); xhi = float(xline[1]);
    ylo = float(yline[0]); yhi = float(yline[1]);
    zlo = float(zline[0]); zhi = float(zline[1]);
        
    # Find rotational matrix
    RzRy = matrix_by_matrix(Rz(psi), Ry(theta))
    RzRyRx = matrix_by_matrix(RzRy, Rx(phi))
    old_pos = {'x':set(), 'y':set(), 'z':set()}
    new_pos = {'x':set(), 'y':set(), 'z':set()}
    for i in m.atoms:
        atom = m.atoms[i]
        old_pos['x'].add(atom.x)
        old_pos['y'].add(atom.y)
        old_pos['z'].add(atom.z)
        
        
        # Rotate around atoms
        x, y, z = vector_by_matrix(RzRyRx, [atom.x, atom.y, atom.z])
        
        # Save data into m.atoms
        atom.x = x
        atom.y = y
        atom.z = z
        new_pos['x'].add(atom.x)
        new_pos['y'].add(atom.y)
        new_pos['z'].add(atom.z)
        
    # Find dilo and diho buffers from orginal system
    if xlo < 0: dxlo = xlo - min(old_pos['x'])
    else: dxlo = min(old_pos['x']) - xlo
    if xhi > 0: dxhi = xhi - max(old_pos['x'])
    else: dxhi = xhi - max(old_pos['x'])
    
    if ylo < 0: dylo = ylo - min(old_pos['y'])
    else: dylo = min(old_pos['y']) - ylo
    if yhi > 0: dyhi = yhi - max(old_pos['y'])
    else: dyhi = yhi - max(old_pos['y'])
    
    if zlo < 0: dzlo = zlo - min(old_pos['z'])
    else: dzlo = min(old_pos['z']) - zlo
    if zhi > 0: dzhi = zhi - max(old_pos['z'])
    else: dzhi = zhi - max(old_pos['z'])
    
    # reset box dims  
    #print(dxlo, dxhi)
    #print(dylo, dyhi)
    #print(dzlo, dzhi)
    xlo = min(new_pos['x']) + dxlo; xhi = max(new_pos['x']) + dxhi;
    ylo = min(new_pos['y']) + dylo; yhi = max(new_pos['y']) + dyhi;
    zlo = min(new_pos['z']) + dzlo; zhi = max(new_pos['z']) + dzhi;
    m.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(xlo, xhi, 'xlo', 'xhi')
    m.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(ylo, yhi, 'ylo', 'yhi')
    m.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(zlo, zhi, 'zlo', 'zhi')
    return m

############################################
# Function to apply shift to atoms and box #
############################################
def shift_system(m, Shifts, log):
    try:
        # Update simulation cell
        xline = m.xbox_line.split(); yline = m.ybox_line.split(); zline = m.zbox_line.split();
        xlo = float(xline[0]) + Shifts['x']; xhi = float(xline[1]) + Shifts['x'];
        ylo = float(yline[0]) + Shifts['y']; yhi = float(yline[1]) + Shifts['y'];
        zlo = float(zline[0]) + Shifts['z']; zhi = float(zline[1]) + Shifts['z'];
        m.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(xlo, xhi, 'xlo', 'xhi')
        m.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(ylo, yhi, 'ylo', 'yhi')
        m.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(zlo, zhi, 'zlo', 'zhi')
        
        # Update atom positions
        for i in m.atoms:
            atom = m.atoms[i]
            atom.x += Shifts['x']
            atom.y += Shifts['y']
            atom.z += Shifts['z']
    except: log.error(f'ERROR failed to apply shifts {Shifts} to {os.path.basename(m.filename)}')
    return m


############################################
# Function for finding molecules in system #
############################################
def clusters(m, log, pflag=True):
    # Generate graph
    graph = {i:[] for i in m.atoms}
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        graph[id1].append(id2)
        graph[id2].append(id1)
    
    # Initialize clusters and molids
    molids = {} # { atom id : molid }
    clusters = set([]); checked = {ID:False for ID in m.atoms};
    for ID in graph:
        molids[ID] = 1; # initialize as 1 and update later
        if checked[ID]: continue
        visited=set([ID]); queue=[ID];
        while queue:
            s = queue.pop(0) 
            for neighbor in graph[s]:
                if checked[neighbor]: continue
                visited.add(neighbor)
                queue.append(neighbor)
                checked[neighbor]=True
        clusters.add( tuple(sorted(visited)) )
    
    # Sort clusters in a very unique fashion
    clusters = sorted(clusters, key=lambda x: x[0]) # Sort all clusters based on 1st atomID in cluster
    clusters = sorted(clusters, key=len, reverse=True) # Sort all clusters by number of atoms
    
    # Find molids after sorting clusters. molid will be count + 1 where the largest
    # cluster will be molid = 1 and smallet cluster will be molid = nclusters in system
    for molID, cluster in enumerate(clusters, 1):
        for ID in cluster:
            molids[ID] = molID

    # analyze clusters
    natoms = [len(c) for c in clusters];                                                                                                                                                              
    size_tot = sum(natoms);                                                                                    
    
    # Print if pflag
    if pflag:
        maxID = 50
        log.out('-----------------Cluster Analysis-----------------')
        log.out('{:^10} {:^15} {:^25}'.format('molID', 'Size of Fragment', '%Size of Fragment'))
        log.out('--------------------------------------------------')                                    
        for n, cnatoms in enumerate(natoms):                                                                                                                                    
            natoms = '{: >6}'.format(cnatoms)                                                                                                         
            size_per = '{:.2f}'.format(round(100*cnatoms/size_tot, 2))   
        
            # Printing info and saving for log file
            if n+1 <= maxID:
                log.out('{:^10} {:^15} {:^25}'.format(n+1, natoms, size_per)) 
        if n+1 > maxID:
            for i in range(3):
                log.out('{:^10} {:^15} {:^25}'.format('.', '.', '.'))                                                                          
    return clusters, molids


##################################################
# Function to sum bond increments and set charge #
##################################################
def charge_from_bond_increments(id1, bonded, nta, frc, use_auto_equivalence, skip_printouts, log):
    new_charge = 0; bond_incs = {} # { tuple(id1, id2): lst[(nta1_a, nta2_a), (b), ...] }
    for id2 in bonded[id1]:
        # Find new types
        type1 = nta[id1]
        type2 = nta[id2]
        
        # Flag things
        types_flag = False
        equivalent_flag = False
        auto_equivalent_flag = False
        
        # add to bond_incs dict
        bond_incs[(id1, id2)] = []
        
        
        # Try getting bond-increments without equivalences in forward ordering
        if (type1, type2) in frc.bond_increments:
            new_charge += frc.bond_increments[(type1, type2)].ij
            bond_incs[(id1, id2)].append((type1, type2))
            bond_incs[(id1, id2)].append(frc.bond_increments[(type1, type2)].ij)
            types_flag = True
        
        # Try getting bond-increments without equivalences in reverse ordering
        elif (type2, type1) in frc.bond_increments:
            new_charge += frc.bond_increments[(type2, type1)].ji
            bond_incs[(id1, id2)].append((type2, type1))
            bond_incs[(id1, id2)].append(frc.bond_increments[(type2, type1)].ji)
            types_flag = True
        
        
        # Try finding equivalences for bond to sum bond increments
        elif not types_flag:
            # Try finding auto_equivalent types
            if type1 in frc.equivalences and type2 in frc.equivalences:
                equiv1 = frc.equivalences[type1].bond
                equiv2 = frc.equivalences[type2].bond
            
                # Try forward ordering
                if (equiv1, equiv2) in frc.bond_increments:
                    new_charge += frc.bond_increments[(equiv1, equiv2)].ij
                    bond_incs[(id1, id2)].append((equiv1, equiv2))
                    bond_incs[(id1, id2)].append(frc.bond_increments[(equiv1, equiv2)].ij)
                    equivalent_flag = True
                
                # Try reverse ordering
                elif (equiv2, equiv1) in frc.bond_increments:
                    new_charge += frc.bond_increments[(equiv2, equiv1)].ji
                    bond_incs[(id1, id2)].append((equiv2, equiv1))
                    bond_incs[(id1, id2)].append(frc.bond_increments[(equiv2, equiv1)].ji)
                    equivalent_flag = True
                    
            # Warn if equivalences not found
            else:
                failed = ', '.join([i for i in (type1, type2) if i not in frc.equivalences])
                log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed)) 
                
        # Try finding auto-equivalences for bond to sum bond increments
        elif not types_flag and not equivalent_flag and use_auto_equivalence:
            # Try finding auto_equivalent types
            if type1 in frc.auto_equivalences and type2 in frc.auto_equivalences:
                equiv1 = frc.auto_equivalences[type1].bond_inct
                equiv2 = frc.auto_equivalences[type2].bond_inct
            
                # Try forward ordering
                if (equiv1, equiv2) in frc.bond_increments:
                    new_charge += frc.bond_increments[(equiv1, equiv2)].ij
                    bond_incs[(id1, id2)].append((equiv1, equiv2))
                    bond_incs[(id1, id2)].append(frc.bond_increments[(equiv1, equiv2)].ij)
                    auto_equivalent_flag = True
                
                # Try reverse ordering
                elif (equiv2, equiv1) in frc.bond_increments:
                    new_charge += frc.bond_increments[(equiv2, equiv1)].ji
                    bond_incs[(id1, id2)].append((equiv2, equiv1))
                    bond_incs[(id1, id2)].append(frc.bond_increments[(equiv2, equiv1)].ji)
                    auto_equivalent_flag = True
                    
            # Warn if equivalences not found
            else:
                failed = ', '.join([i for i in (type1, type2) if i not in frc.auto_equivalences])
                log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed)) 
            
            
        # Else sum a zero and print failure
        if not types_flag and not equivalent_flag and not auto_equivalent_flag:
            new_charge += 0
            if not skip_printouts:
                log.warn('WARNING No Bond Increment for atom {} and atom {}, bond pair {} {} '.format(id1, id2, type1, type2))
    return new_charge, bond_incs


###################################################################################
# Function to sum bond increments and set charge if edge id was used in .nta file #
###################################################################################
def nta_edge_id_bond_incs_extend(id1, type1, edges, frc, use_auto_equivalence, log):
    extended_charge = 0;
    for type2 in edges:
        
        # Flag things
        types_flag = False
        equivalent_flag = False
        auto_equivalent_flag = False
        
        # Try getting bond-increments without equivalences in forward ordering
        if (type1, type2) in frc.bond_increments:
            extended_charge += frc.bond_increments[(type1, type2)].ij
            types_flag = True
        
        # Try getting bond-increments without equivalences in reverse ordering
        elif (type2, type1) in frc.bond_increments:
            extended_charge += frc.bond_increments[(type2, type1)].ji
            types_flag = True
        
        
        # Try finding equivalences for bond to sum bond increments
        elif not types_flag:
            # Try finding auto_equivalent types
            if type1 in frc.equivalences and type2 in frc.equivalences:
                equiv1 = frc.equivalences[type1].bond
                equiv2 = frc.equivalences[type2].bond
            
                # Try forward ordering
                if (equiv1, equiv2) in frc.bond_increments:
                    extended_charge += frc.bond_increments[(equiv1, equiv2)].ij
                    equivalent_flag = True
                
                # Try reverse ordering
                elif (equiv2, equiv1) in frc.bond_increments:
                    extended_charge += frc.bond_increments[(equiv2, equiv1)].ji
                    equivalent_flag = True
                    
            # Warn if equivalences not found
            else:
                failed = ', '.join([i for i in (type1, type2) if i not in frc.equivalences])
                log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed)) 
                
        # Try finding auto-equivalences for bond to sum bond increments
        elif not types_flag and not equivalent_flag and use_auto_equivalence:
            # Try finding auto_equivalent types
            if type1 in frc.auto_equivalences and type2 in frc.auto_equivalences:
                equiv1 = frc.auto_equivalences[type1].bond_inct
                equiv2 = frc.auto_equivalences[type2].bond_inct
            
                # Try forward ordering
                if (equiv1, equiv2) in frc.bond_increments:
                    extended_charge += frc.bond_increments[(equiv1, equiv2)].ij
                    auto_equivalent_flag = True
                
                # Try reverse ordering
                elif (equiv2, equiv1) in frc.bond_increments:
                    extended_charge += frc.bond_increments[(equiv2, equiv1)].ji
                    auto_equivalent_flag = True
                    
            # Warn if equivalences not found
            else:
                failed = ', '.join([i for i in (type1, type2) if i not in frc.auto_equivalences])
                log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed)) 
            
            
        # Else sum a zero and print failure
        if not types_flag and not equivalent_flag and not auto_equivalent_flag:
            extended_charge += 0
            log.warn('WARNING No Bond Increment for atom {} and atom {}, bond pair {} {} '.format(id1, 'edge extend', type1, type2)) 
    return extended_charge


##################################################
# Function to compute system mass/volume/density #
##################################################
def compute_mass_volume_density(parameters, BADI, ff_class, remove_booleans, reset_charges, ignore_missing_parameters, frc_file, log):
    system_mass_sum = 0; amu2grams = 1/6.02214076e+23;
    for i in parameters.atoms:
        atom = parameters.atoms[i]
        mass = parameters.masses[atom.type]
        system_mass_sum += mass.coeffs
        
    # convert system mass in amu to grams
    mass = system_mass_sum*amu2grams
    
    # Find box dimensions to compute volume
    angstromcubed2cmcubed = 1e-24
    lx = parameters.box.lx
    ly = parameters.box.ly
    lz = parameters.box.lz
    volume = lx*ly*lz*angstromcubed2cmcubed
    
    # Compute density
    density = mass/volume
    
    ###########################################
    # print density and box results to screen #
    ###########################################
    log.out('\n---------------------------------------')
    log.out('System Unit cell, volume, mass, density')
    log.out('---------------------------------------')
    log.out('{:<10} {:^10.5f} {:<10}'.format('Lx:', lx, 'angstrom'))
    log.out('{:<10} {:^10.5f} {:<10}'.format('Ly:', ly, 'angstrom'))
    log.out('{:<10} {:^10.5f} {:<10}'.format('Lz:', lz, 'angstrom'))
    log.out('{:<10} {:^10.4E} {:<10}'.format('volume:', volume, 'cm^3'))
    log.out('{:<10} {:^10.4E} {:<10}'.format('mass:', mass, 'grams'))
    log.out('{:<10} {:^10.5f} {:<10}' .format('density:', density, 'g/cm^3'))
    
    #######################################################
    # Print out number of found topologies and parameters #
    #######################################################
    log.out('\n\n\n--------------------------------------')
    log.out('Atoms/Bonds/Angles/Dihedrals/Impropers')
    log.out('--------------------------------------')
    log.out('{:<18} {:^5}'.format('natoms: ', parameters.natoms))
    log.out('{:<18} {:^5}'.format('nbonds: ', parameters.nbonds))
    log.out('{:<18} {:^5}'.format('nangles: ', parameters.nangles))
    log.out('{:<18} {:^5}'.format('ndihedrals: ', parameters.ndihedrals))
    log.out('{:<18} {:^5}'.format('nimpropers: ', parameters.nimpropers))
    if ff_class == 2 and not any(remove_booleans):
        log.out('{:<2}{:<16} {:^5}'.format('-', 'oop: ', len(BADI.impropers) - len(BADI.angleangles) ))
        log.out('{:<2}{:<16} {:^5}'.format('-', 'angleangle: ', len(BADI.angleangles) ))
    
    log.out('\n--------------------------------------------')
    log.out('Atoms/Bonds/Angles/Dihedrals/Impropers Types')
    log.out('--------------------------------------------')
    log.out('{:<18} {:^5}'.format('natomtypes:', parameters.natomtypes))
    log.out('{:<18} {:^5}'.format('nbondtypes:', parameters.nbondtypes))
    log.out('{:<18} {:^5}'.format('nangletypes:', parameters.nangletypes))
    log.out('{:<18} {:^5}'.format('ndihedraltypes:', parameters.ndihedraltypes))
    log.out('{:<18} {:^5}'.format('nimpropertypes:', parameters.nimpropertypes))
    if ff_class == 2 and not any(remove_booleans):
        log.out('{:<2}{:<16} {:^5}'.format('-', 'oop: ', len(BADI.improper_types_dict) ))
        log.out('{:<2}{:<16} {:^5}'.format('-', 'angleangle: ', len(BADI.angleangle_types_dict) ))
    
    ##########################################
    # Print out parameterization percentages #
    ##########################################
    def calculate_percent_found(dict2search):
        percent = 100; total = len(dict2search);
        failed = sum([1 for i in dict2search if 'UNABLE' in dict2search[i].comments])
        if total > 0 and failed > 0:
            found = total - failed
            percent = 100*found/total
        return percent
    log.out('\n-------------------------------------')
    log.out('Parameterization percentage breakdown')
    log.out('-------------------------------------')
    percents = {} # { 'Coeff name' : float value of percent parameterized }
    if ff_class in [0, 1, 2, 'd']:
        if parameters.masses: percents['Masses'] = calculate_percent_found(parameters.masses)
        if parameters.pair_coeffs: percents['Pair Coeffs'] =calculate_percent_found(parameters.pair_coeffs)
        if parameters.bond_coeffs: percents['Bond Coeffs'] = calculate_percent_found(parameters.bond_coeffs)
        if parameters.angle_coeffs: percents['Angle Coeffs'] = calculate_percent_found(parameters.angle_coeffs)
        if parameters.dihedral_coeffs: percents['Dihedral Coeffs'] = calculate_percent_found(parameters.dihedral_coeffs)
        if parameters.improper_coeffs: percents['Improper Coeffs'] = calculate_percent_found(parameters.improper_coeffs)
        percents['Average MPBADI Coeffs'] = sum(list(percents.values()))/len(percents)
        if ff_class == 2:
            if parameters.bondbond_coeffs: percents['BondBond Coeffs'] = calculate_percent_found(parameters.bondbond_coeffs)
            if parameters.bondangle_coeffs: percents['BondAngle Coeffs'] = calculate_percent_found(parameters.bondangle_coeffs)
            if parameters.angleangletorsion_coeffs: percents['AngleAngleTorsion Coeffs'] = calculate_percent_found(parameters.angleangletorsion_coeffs)
            if parameters.endbondtorsion_coeffs: percents['EndBondTorsion Coeffs'] = calculate_percent_found(parameters.endbondtorsion_coeffs)
            if parameters.middlebondtorsion_coeffs: percents['MiddleBondTorsion Coeffs'] = calculate_percent_found(parameters.middlebondtorsion_coeffs)
            if parameters.bondbond13_coeffs: percents['BondBond13 Coeffs'] = calculate_percent_found(parameters.bondbond13_coeffs)
            if parameters.angletorsion_coeffs: percents['AngleTorsion Coeffs'] = calculate_percent_found(parameters.angletorsion_coeffs)
            if parameters.angleangle_coeffs: percents['AngleAngle Coeffs'] = calculate_percent_found(parameters.angleangle_coeffs)
        percents['Average (without Bond-incs)'] = sum(list(percents.values()))/len(percents)
        if reset_charges:
            total = 0; failed = 0; percent = 100;
            for i in parameters.atoms:
                bond_incs = parameters.atoms[i].bond_incs
                total += len(bond_incs)
                failed += list(bond_incs.values()).count([])
            if total > 0 and failed > 0:
                found = total - failed
                percent = 100*found/total
            percents['Bond-incs'] = percent
        for i in percents:
            log.out('{:<28} {:^6.2f} %'.format(i+str(':'), percents[i]))
    parameters.percents = percents
     
    ########################################################################################
    # Zero all coeffs if percents['Average (without Bond-incs)'] < 100 and not ignore flag #
    ########################################################################################
    def zero_coeffs(dict2zero, mass_flag=False):
        # Special mass zeroing
        if mass_flag:
            for i in dict2zero:
                dict2zero[i].coeffs = float(0.0)
                dict2zero[i].equivalent = 'N/A'
                dict2zero[i].comments = 'ZEROED'
        else: # multiple coeffs zeroing
            ncoeff_types = [len(dict2zero[i].coeffs) for i in dict2zero]
            ncoeff_most_frequent = max(set(ncoeff_types), key = ncoeff_types.count)
            for i in dict2zero:
                dict2zero[i].coeffs = ncoeff_most_frequent*[float(0.0)]
                dict2zero[i].equivalent = 'N/A'
                dict2zero[i].match = 'N/A'
                dict2zero[i].comments = 'ZEROED'
        return None
    if ff_class in [0, 1, 2, 'd'] and percents['Average (without Bond-incs)'] < 100:
        ignore_printout = True
        if not ignore_missing_parameters:
            ignore_printout = False
            log.out('\nignore_missig_parameters is False')
            log.out('  - zeroing all Coeffs making them non operative')
            log.out('  - set ignore_missig_parameters as True to get partially parmaterized datafile')
            if parameters.masses: zero_coeffs(parameters.masses, mass_flag=True)
            if parameters.pair_coeffs: zero_coeffs(parameters.pair_coeffs, mass_flag=False)
            if parameters.bond_coeffs: zero_coeffs(parameters.bond_coeffs, mass_flag=False)
            if parameters.angle_coeffs: zero_coeffs(parameters.angle_coeffs, mass_flag=False)
            if parameters.dihedral_coeffs: zero_coeffs(parameters.dihedral_coeffs, mass_flag=False)
            if parameters.improper_coeffs: zero_coeffs(parameters.improper_coeffs, mass_flag=False)
            if ff_class == 2:
                if parameters.bondbond_coeffs: zero_coeffs(parameters.bondbond_coeffs, mass_flag=False)
                if parameters.bondangle_coeffs: zero_coeffs(parameters.bondangle_coeffs, mass_flag=False)
                if parameters.angleangletorsion_coeffs: zero_coeffs(parameters.angleangletorsion_coeffs, mass_flag=False)
                if parameters.endbondtorsion_coeffs: zero_coeffs(parameters.endbondtorsion_coeffs, mass_flag=False)
                if parameters.middlebondtorsion_coeffs: zero_coeffs(parameters.middlebondtorsion_coeffs, mass_flag=False)
                if parameters.bondbond13_coeffs: zero_coeffs(parameters.bondbond13_coeffs, mass_flag=False)
                if parameters.angletorsion_coeffs: zero_coeffs(parameters.angletorsion_coeffs, mass_flag=False)
                if parameters.angleangle_coeffs: zero_coeffs(parameters.angleangle_coeffs, mass_flag=False)
        if ignore_printout:
            log.out('\nUser is acknowledging that all2lmp.py can not supply a fully parameterized datafile as')
            log.out(f'the {frc_file} doesnt contain all parameters')
            log.out('to parameterize every interaction. The user is responsible for deciding on how to handle')
            log.out('missing parameters in the written datafile. Please note that it is possible that the inputs')
            log.out('to all2lmp.py are incorrect causing this message and that a certain atom type maybe wrong, ')
            log.out('where it might be easier to correct the inputs to all2lmp.py rather then correct the outputs')
            log.out('of all2lmp.py. Alternativley missing parameters are set to ZEROs, by all2lmp.py and it is')
            log.out('possible that interaction doesnt require parameters. For instance not all systems require')
            log.out('dihedral interactions.')
    return


###################################
# Function to match 2-body coeffs #
###################################
def match_2_body(type1, type2, log, dict2search, equivalences, form):
    """
    match_2_body has the following inputs  and is for matching bond coeffs/crossterms:
        type1: atom-type1 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        type2: atom-type2 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        
        dict2search: coeff dictionary to search from read_frc
        
        equivalences: coeff equivalence dictionary to find equivalences of if needed
            - can be a dictionary or False. The dictionary will tell the code to use
              the equivalences and False will skip the use of the equivalences
    """
    
    # Set match/order as (type1, type2) to return if not found
    boolean = False; match = (type1, type2); order = (type1, type2); equiv = (type1, type2);
    
    # Try matching in forward/reverse order without equivalences
    if (type1, type2) in dict2search:
        match = (type1, type2)
        order = (type1, type2)
        equiv = (type1, type2)
        boolean = True
    elif (type2, type1) in dict2search:
        match = (type2, type1)
        order = (type2, type1)
        equiv = (type2, type1)
        boolean = True
        
    # Next try using equivalences if equivalences not False
    if not boolean and equivalences:
        equiv_flag = False # set flag and if found update as True
        
        # Find equivalences for quartic or quadratic form
        if form == 'equiv':
            if type1 in equivalences and type2 in equivalences:
                equiv1 = equivalences[type1].bond
                equiv2 = equivalences[type2].bond
                equiv_flag = True
            else:
                failed = ', '.join([i for i in (type1, type2) if i not in equivalences])
                log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed))      
        elif form == 'auto-equiv':
            if type1 in equivalences and type2 in equivalences:
                equiv1 = equivalences[type1].bond
                equiv2 = equivalences[type2].bond
                equiv_flag = True
            else:
                failed = ', '.join([i for i in (type1, type2) if i not in equivalences])
                log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed))
        else:
            log.error('ERROR unsupported equivalences are for either equiv or auto-equiv')
            
        # If equivalences exists start using equivalences to try matching
        if equiv_flag:            
            # Try matching in forward/reverse order with equivalences
            if (equiv1, equiv2) in dict2search:
                match = (equiv1, equiv2)
                order = (type1, type2)
                equiv = (equiv1, equiv2)
                boolean = True     
            elif (equiv2, equiv1) in dict2search:
                match = (equiv2, equiv1)
                order = (type2, type1)
                equiv = (equiv2, equiv1)
                boolean = True
    return boolean, match, order, equiv


###################################
# Function to get cross-term r0's #
###################################
def get_crossterm_r0(type1, type2, log, use_auto_equivalence, frc):
    """
    Function to get r0's for cross terms. Will search
    both qurtic and quadratic bond coeffs to find r0.
    
    Will return r0 = 0, equiv = (N/A, N/A), and match =
    (N/A, N/A) if r0 bond type can not be found.
    """
    
    # Set flags and intial return values
    quartic_flag = False
    quadratic_flag = False
    equiv = ('N/A', 'N/A')
    match = ('N/A', 'N/A')
    r0 = 0.0
    
    
    # Try matching equiv bonds in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
    if match_2_body(type1, type2, log, frc.quartic_bonds, equivalences=frc.equivalences, form='equiv')[0]:
        boolean, match, order, equiv = match_2_body(type1, type2, log, frc.quartic_bonds, equivalences=frc.equivalences, form='equiv')
        bond_coeff = frc.quartic_bonds[match]
        quartic_flag = True
        
    # Try matching auto-equivc bonds in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
    elif match_2_body(type1, type2, log, frc.quadratic_bonds, equivalences=frc.auto_equivalences, form='auto-equiv')[0] and use_auto_equivalence:
        boolean, match, order, equiv = match_2_body(type1, type2, log, frc.quadratic_bonds, equivalences=frc.auto_equivalences, form='auto-equiv')
        bond_coeff = frc.quadratic_bonds[match]
        quadratic_flag = True
        
    # Find r0 and equivalent/match used
    if quartic_flag and not quadratic_flag:
        r0 = bond_coeff.r0
        equiv = equiv
        match = match
    elif quadratic_flag and not quartic_flag:
        r0 = bond_coeff.r0
        equiv = equiv
        match = match
    return r0, match, equiv


###################################
# Function to match 3-body coeffs #
###################################
def match_angle_wildcards(coeff, dict2search, test, log):
    # Set match as coeff and update later
    match = coeff;  return_boolean = False;
    
    # If test = [1, 1, 0] look for exact match of first two atoms, but last atom can be a wild card
    if test == [1, 1, 0]:
        for i, j, k in dict2search:
            if i == coeff[0] and j == coeff[1] and k[0] == '*':
                match = (i, j, k); return_boolean = True; break;
    
    # If test = [0, 1, 1] look for exact match of last two atoms, but first atom can be a wild card
    elif test == [0, 1, 1]:
        for i, j, k in dict2search:
            if i[0] == '*' and j == coeff[1] and k == coeff[2]:
                match = (i, j, k); return_boolean = True; break;
                
    # If test = [0, 1, 0] look for exact match of center atoms, but outer atom can be a wild card
    elif test == [0, 1, 0]:
        for i, j, k in dict2search:
            if i[0] == '*' and j == coeff[1] and k[0] == '*':
                match = (i, j, k); return_boolean = True; break;
                
    # If user tries accessing unsupported wild card check raise exception
    else:
        log.error(f'ERROR {str(test)} wild card test not supported')
    return return_boolean, match


def match_3_body(type1, type2, type3, log, dict2search, equivalences, wildcard_search, form):
    """
    match_3_body has the following inputs and is for matching angle coeffs/crossterms:
        type1: atom-type1 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        type2: atom-type2 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        type3: atom-type3 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        
        dict2search: coeff dictionary to search from read_frc
        
        equivalences: coeff equivalence dictionary to find equivalences of if needed
            - can be a dictionary or False. The dictionary will tell the code to use
              the equivalences and False will skip the use of the equivalences
              
        wildcard_search: True or False to tell the function to search for wild card atom-types in coeff
    """
    
    # Set match/order as (type1, type2, type3) to return if not found
    boolean = False; match = (type1, type2, type3); order = (type1, type2, type3); equiv = (type1, type2, type3); 
    
    # Try matching in forward/reverse order without equivalences
    if (type1, type2, type3) in dict2search:
        match = (type1, type2, type3)
        order = (type1, type2, type3)
        equiv = (type1, type2, type3)
        boolean = True
    elif (type3, type2, type1) in dict2search:
        match = (type3, type2, type1)
        order = (type3, type2, type1)
        equiv = (type3, type2, type1)
        boolean = True
        
    # Next try using equivalences if equivalences not False
    if not boolean and equivalences:
        equiv_flag = False # set flag and if found update as True
        
        # Find equivalences for quartic or quadratic form
        if form == 'equiv':
            if type1 in equivalences and type2 in equivalences and type3 in equivalences:
                equiv1 = equivalences[type1].angle
                equiv2 = equivalences[type2].angle
                equiv3 = equivalences[type3].angle
                equiv_flag = True
            else:
                failed = ', '.join([i for i in (type1, type2, type3) if i not in equivalences])
                log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed))      
        elif form == 'auto-equiv':
            if type1 in equivalences and type2 in equivalences and type3 in equivalences:
                equiv1 = equivalences[type1].angle_end
                equiv2 = equivalences[type2].angle_apex
                equiv3 = equivalences[type3].angle_end
                equiv_flag = True
            else:
                failed = ', '.join([i for i in (type1, type2, type3) if i not in equivalences])
                log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed))
        else:
            log.error('ERROR unsupported equivalences are for either equiv or auto-equiv')
            
        # If equivalences exists start using equivalences to try matching
        if equiv_flag:            
            # Try matching in forward/reverse order with equivalences
            if (equiv1, equiv2, equiv3) in dict2search:
                match = (equiv1, equiv2, equiv3)
                order = (type1, type2, type3)
                equiv = (equiv1, equiv2, equiv3)
                boolean = True     
            elif (equiv3, equiv2, equiv1) in dict2search:
                match = (equiv3, equiv2, equiv1)
                order = (type3, type2, type1)
                equiv = (equiv3, equiv2, equiv1)
                boolean = True
                
            # Try matching using wildcards in forward and reverse
            elif wildcard_search:
                # If test = [1, 1, 0] look for exact match of first two atoms, but last atom can be a wild card
                if match_angle_wildcards((equiv1, equiv2, equiv3), dict2search, [1, 1, 0], log)[0]:
                    match = match_angle_wildcards((equiv1, equiv2, equiv3), dict2search, [1, 1, 0], log)[1]
                    order = (type1, type2, type3)
                    equiv = (equiv1, equiv2, equiv3)
                    boolean = True 
                elif match_angle_wildcards((equiv3, equiv2, equiv1), dict2search, [1, 1, 0], log)[0]:
                    match = match_angle_wildcards((equiv3, equiv2, equiv1), dict2search, [1, 1, 0], log)[1]
                    order = (type3, type2, type1)
                    equiv = (equiv3, equiv2, equiv1)
                    boolean = True 
                    
                # If test = [0, 1, 1] look for exact match of last two atoms, but first atom can be a wild card
                elif match_angle_wildcards((equiv1, equiv2, equiv3), dict2search, [0, 1, 1], log)[0]:
                    match = match_angle_wildcards((equiv1, equiv2, equiv3), dict2search, [0, 1, 1], log)[1]
                    order = (type1, type2, type3)
                    equiv = (equiv1, equiv2, equiv3)
                    boolean = True 
                elif match_angle_wildcards((equiv3, equiv2, equiv1), dict2search, [0, 1, 1], log)[0]:
                    match = match_angle_wildcards((equiv3, equiv2, equiv1), dict2search, [0, 1, 1], log)[1]
                    order = (type3, type2, type1)
                    equiv = (equiv3, equiv2, equiv1)
                    boolean = True 
    
                # If test = [0, 1, 0] look for exact match of center atoms, but outer atom can be a wild card
                elif match_angle_wildcards((equiv1, equiv2, equiv3), dict2search, [0, 1, 0], log)[0]:
                    match = match_angle_wildcards((equiv1, equiv2, equiv3), dict2search, [0, 1, 0], log)[1]
                    order = (type1, type2, type3)
                    equiv = (equiv1, equiv2, equiv3)
                    boolean = True 
    return boolean, match, order, equiv


#######################################
# Function to get cross-term theta0's #
#######################################
def get_crossterm_theta0(type1, type2, type3, log, use_auto_equivalence, frc):
    """
    Function to get theta0's for cross terms. Will search
    both qurtic and quadratic angle coeffs to find theta0.
    
    Will return theta0 = 0, equiv = (N/A, N/A, N/A), and
    match = (N/A, N/A, N/A) if theta0 angle type can not be found.
    """
    
    # Set flags and intial return values
    quartic_flag = False
    quadratic_flag = False
    equiv = ('N/A', 'N/A', 'N/A')
    match = ('N/A', 'N/A', 'N/A')
    theta0 = 0.0
    
    
    # Try matching equiv angles in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences
    if match_3_body(type1, type2, type3, log, frc.quartic_angles, equivalences=frc.equivalences, wildcard_search=False, form='equiv')[0]:
        boolean, match, order, equiv = match_3_body(type1, type2, type3, log, frc.quartic_angles, equivalences=frc.equivalences, wildcard_search=False, form='equiv')
        angle_coeff = frc.quartic_angles[match]
        quartic_flag = True
        
    # Try matching equiv angles in forward/reverse. Attempts: 1) without equivalences, 2) with equivalences, 3) with equivalences and wild cards
    elif match_3_body(type1, type2, type3, log, frc.quadratic_angles, equivalences=frc.auto_equivalences, wildcard_search=True, form='auto-equiv')[0] and use_auto_equivalence:
        boolean, match, order, equiv = match_3_body(type1, type2, type3, log, frc.quadratic_angles, equivalences=frc.auto_equivalences, wildcard_search=True, form='auto-equiv')
        angle_coeff = frc.quadratic_angles[match]
        quadratic_flag = True
        
    # Find r0 and equivalent/match used
    if quartic_flag and not quadratic_flag:
        theta0 = angle_coeff.theta0
        equiv = equiv
        match = match
    elif quadratic_flag and not quartic_flag:
        theta0 = angle_coeff.theta0
        equiv = equiv
        match = match
    return theta0, match, equiv


############################################
# Function to match 4-body Dihedral coeffs #
############################################
def match_dihedral_wildcards(coeff, dict2search, test, log):
    # Set match as coeff and update later
    match = coeff;  return_boolean = False;
    
    # If test = [1, 1, 1, 0] look for exact match of first three atoms, but last atom can be a wild card
    if test == [1, 1, 1, 0]:
        for i, j, k, l in dict2search:
            if i == coeff[0] and j == coeff[1] and k == coeff[2] and l[0] == '*':
                match = (i, j, k, l); return_boolean = True; break;
                
    # If test = [0, 1, 1, 1] look for exact match of last three atoms, but first atom can be a wild card
    elif test == [0, 1, 1, 1]:
        for i, j, k, l in dict2search:
            if i[0] == '*' and j == coeff[1] and k == coeff[2] and l == coeff[3]:
                match = (i, j, k, l); return_boolean = True; break;
                
    # If test = [0, 1, 1, 0] look for exact match of inner two atoms, but outer atoms can be a wild card
    elif test == [0, 1, 1, 0]:
        for i, j, k, l in dict2search:
            if i[0] == '*' and j == coeff[1] and k == coeff[2] and l[0] == '*':
                match = (i, j, k, l); return_boolean = True; break;
                
    # If test = [0, 1, 0, 0] look for exact match of inner two atoms, but outer atoms can be a wild card
    elif test == [0, 1, 0, 0]:
        for i, j, k, l in dict2search:
            if i[0] == '*' and j == coeff[1] and k[0] == '*' and l[0] == '*':
                match = (i, j, k, l); return_boolean = True; break;
                
    # If test = [0, 0, 1, 0] look for exact match of inner two atoms, but outer atoms can be a wild card
    elif test == [0, 0,  1, 0]:
        for i, j, k, l in dict2search:
            if i[0] == '*' and j[0] == '*' and k[0] == coeff[2] and l[0] == '*':
                match = (i, j, k, l); return_boolean = True; break;
                
    # If user tries accessing unsupported wild card check raise exception
    else:
        log.error(f'ERROR {str(test)} wild card test not supported')
    return return_boolean, match


def match_4_body(type1, type2, type3, type4, log, dict2search, equivalences, wildcard_search, form):
    """
    match_4_body has the following inputs and is for matching dihedral coeffs/crossterms:
        type1: atom-type1 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        type2: atom-type2 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        type3: atom-type3 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        type4: atom-type4 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        
        dict2search: coeff dictionary to search from read_frc
        
        equivalences: coeff equivalence dictionary to find equivalences of if needed
            - can be a dictionary or False. The dictionary will tell the code to use
              the equivalences and False will skip the use of the equivalences
              
        wildcard_search: True or False to tell the function to search for wild card atom-types in coeff
    """
    
    # Set match/order as (type1, type2, type3, type4) to return if not found
    boolean = False; match = (type1, type2, type3, type4); 
    order = (type1, type2, type3, type4); equiv = (type1, type2, type3, type4); 
    
    # Try matching in forward/reverse order without equivalences
    if (type1, type2, type3, type4) in dict2search:
        match = (type1, type2, type3, type4)
        order = (type1, type2, type3, type4)
        equiv = (type1, type2, type3, type4)
        boolean = True
    elif (type4, type3, type2, type1) in dict2search:
        match = (type4, type3, type2, type1)
        order = (type4, type3, type2, type1)
        equiv = (type4, type3, type2, type1)
        boolean = True
        
    # Next try using equivalences if equivalences not False
    if not boolean and equivalences:
        equiv_flag = False # set flag and if found update as True
        
        # Find equivalences for quartic or quadratic form
        if form == 'equiv':
            if type1 in equivalences and type2 in equivalences and type3 in equivalences and type4 in equivalences:
                equiv1 = equivalences[type1].torsion
                equiv2 = equivalences[type2].torsion
                equiv3 = equivalences[type3].torsion
                equiv4 = equivalences[type4].torsion
                equiv_flag = True
            else:
                failed = ', '.join([i for i in (type1, type2, type3, type4) if i not in equivalences])
                log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed))      
        elif form == 'auto-equiv':
            if type1 in equivalences and type2 in equivalences and type3 in equivalences and type4 in equivalences:
                equiv1 = equivalences[type1].torsion_end
                equiv2 = equivalences[type2].torsion_center
                equiv3 = equivalences[type3].torsion_center
                equiv4 = equivalences[type4].torsion_end
                equiv_flag = True
            else:
                failed = ', '.join([i for i in (type1, type2, type3, type4) if i not in equivalences])
                log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed))
        else:
            log.error('ERROR unsupported equivalences are for either quartic or quadratic')
            
        # If equivalences exists start using equivalences to try matching
        if equiv_flag:            
            # Try matching in forward/reverse order with equivalences
            if (equiv1, equiv2, equiv3, equiv4) in dict2search:
                match = (equiv1, equiv2, equiv3, equiv4)
                order = (type1, type2, type3, type4)
                equiv = (equiv1, equiv2, equiv3, equiv4)
                boolean = True     
            elif (equiv4, equiv3, equiv2, equiv1) in dict2search:
                match = (equiv4, equiv3, equiv2, equiv1)
                order = (type4, type3, type2, type1)
                equiv = (equiv4, equiv3, equiv2, equiv1)
                boolean = True
                
            # Try matching using wildcards in forward and reverse
            elif wildcard_search:
                # If test = [1, 1, 1, 0] look for exact match of first three atoms, but last atom can be a wild card
                if match_dihedral_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [1, 1, 1, 0], log)[0]:
                    match = match_dihedral_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [1, 1, 1, 0], log)[1]
                    order = (type1, type2, type3, type4)
                    equiv = (equiv1, equiv2, equiv3, equiv4)
                    boolean = True 
                elif match_dihedral_wildcards((equiv4, equiv3, equiv2, equiv1), dict2search, [1, 1, 1, 0], log)[0]:
                    match = match_dihedral_wildcards((equiv4, equiv3, equiv2, equiv1), dict2search, [1, 1, 1, 0], log)[1]
                    order = (type4, type3, type2, type1)
                    equiv = (equiv4, equiv3, equiv2, equiv1)
                    boolean = True 
                    
                # If test = [0, 1, 1, 1] look for exact match of last three atoms, but first atom can be a wild card
                elif match_dihedral_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [0, 1, 1, 1], log)[0]:
                    match = match_dihedral_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [0, 1, 1, 1], log)[1]
                    order = (type1, type2, type3, type4)
                    equiv = (equiv1, equiv2, equiv3, equiv4)
                    boolean = True 
                elif match_dihedral_wildcards((equiv4, equiv3, equiv2, equiv1), dict2search, [0, 1, 1, 1], log)[0]:
                    match = match_dihedral_wildcards((equiv4, equiv3, equiv2, equiv1), dict2search, [0, 1, 1, 1], log)[1]
                    order = (type4, type3, type2, type1)
                    equiv = (equiv4, equiv3, equiv2, equiv1)
                    boolean = True 
                
                # If test = [0, 1, 1, 0] look for exact match of inner two atoms, but outer atoms can be a wild card
                elif match_dihedral_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [0, 1, 1, 0], log)[0]:
                    match = match_dihedral_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [0, 1, 1, 0], log)[1]
                    order = (type1, type2, type3, type4)
                    equiv = (equiv1, equiv2, equiv3, equiv4)
                    boolean = True 
                    
                # If test = [0, 1, 1, 0] look for exact match of inner two atoms in reverse, but outer atoms can be a wild card
                elif match_dihedral_wildcards((equiv4, equiv3, equiv2, equiv1), dict2search, [0, 1, 1, 0], log)[0]:
                    match = match_dihedral_wildcards((equiv4, equiv3, equiv2, equiv1), dict2search, [0, 1, 1, 0], log)[1]
                    order = (type4, type3, type2, type1)
                    equiv = (equiv4, equiv3, equiv2, equiv1)
                    boolean = True 
                    
                #-----------------------------------------------------------------#
                # LIKELY DREIDING ONLY, but still try finding these cases as well #
                #-----------------------------------------------------------------#
                    
                # If test = [0, 1, 0, 0] look for exact match of inner two atomID2 in forward, but outer atoms can be a wild card
                elif match_dihedral_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [0, 1, 0, 0], log)[0]:
                    match = match_dihedral_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [0, 1, 0, 0], log)[1]
                    order = (type1, type2, type3, type4)
                    equiv = (equiv1, equiv2, equiv3, equiv4)
                    boolean = True 
                    
                # If test = [0, 1, 0, 0] look for exact match of inner two atomID2 in forward, but outer atoms can be a wild card
                elif match_dihedral_wildcards((equiv4, equiv3, equiv2, equiv1), dict2search, [0, 1, 0, 0], log)[0]:
                    match = match_dihedral_wildcards((equiv4, equiv3, equiv2, equiv1), dict2search, [0, 1, 0, 0], log)[1]
                    order = (type4, type3, type2, type1)
                    equiv = (equiv4, equiv3, equiv2, equiv1)
                    boolean = True 
    return boolean, match, order, equiv


###################################################################################
# Function for fitting autoequivalent dihedral coeff to equivalent dihedral coeff #
###################################################################################
def dihedral_autoequiv2equiv(dihedral_coeff):
    # Get phi0, v, and n and intialize coeff_list
    phi0 = dihedral_coeff.phi0; v = dihedral_coeff.kphi; n = dihedral_coeff.n;
    coeff_list = 6*[0.0]
    
    # Build coeff list
    if int(n) != 0:

        # Adjust phi0 accordingly and set as phi
        if phi0 == 180: phi = 0.0; # reset phi if phi0 is 180 back to 0
        #elif phi0 == 0.0: phi = 180; # reset phi if phi0 is 0 back to 180 (Keep commented b/c it is how PCFF_parsing was written)
        else: phi = phi0; # else keep phi as is
    
        # Rebuild coeffs list
        if int(n) == 1:
            coeff_list = [v, phi, 0.0, 0.0, 0.0, 0.0]
        elif int(n) == 2:
            coeff_list = [0.0, 0.0, v, phi, 0.0, 0.0]
        elif int(n) == 3:
            coeff_list = [0.0, 0.0, 0.0, 0.0, v, phi]
    return coeff_list


#######################################################
# Function for fitting autoequivalent dihedral coeff  #
#######################################################
def dihedral_autoequiv(dihedral_coeff):
    # Get phi0, k and n and intialize coeff_list
    phi0 = dihedral_coeff.phi0; k = dihedral_coeff.kphi; n = dihedral_coeff.n;
    coeff_list = [0.0, 0, 0]
    
    # Rebuild coeff list
    if phi0 == 0:
        d = 1;
        coeff_list = [k, d, n]
    elif phi0 == 180:
        d = -1
        coeff_list = [k, d, n]   
    else: # For Drieding support 2/10/2023
        d = 0
        coeff_list = [k, d, n]
    return coeff_list


#######################################
# Function to match 4-body oop coeffs #
#######################################
def match_improper_wildcards(coeff, dict2search, test, log):
    # Set match as coeff and update later
    match = coeff;  return_boolean = False;
    
    # If test = [1, 1, 1, 0] look for exact match of first three atoms, but last atom can be a wild card
    if test == [1, 1, 1, 0]:
        for i, j, k, l in dict2search:
            # Try as 1-2 (indexes 0, 1) be constant and varying 3-4 (indexes 2, 3) ordering
            if i == coeff[0] and j == coeff[1] and k == coeff[2] and l[0] == '*':
                match = (i, j, k, l); return_boolean = True; break;
            elif i == coeff[0] and j == coeff[1] and k == coeff[3] and l[0] == '*':
                match = (i, j, k, l); return_boolean = True; break;
                
            # Try as 3-2 (indexes 2, 1) be constant and varying 1-4 (indexes 0, 3) sordering
            elif i == coeff[2] and j == coeff[1] and k == coeff[0] and l[0] == '*':
                match = (i, j, k, l); return_boolean = True; break;
            elif i == coeff[2] and j == coeff[1] and k == coeff[3] and l[0] == '*':
                match = (i, j, k, l); return_boolean = True; break;
                
            # Try as 4-2 (indexes 3, 1) be constant and varying 1-3 (indexes 0, 2) ordering
            elif i == coeff[3] and j == coeff[1] and k == coeff[0] and l[0] == '*':
                match = (i, j, k, l); return_boolean = True; break;
            elif i == coeff[3] and j == coeff[1] and k == coeff[2] and l[0] == '*':
                match = (i, j, k, l); return_boolean = True; break;
                
    # If test = [0, 1, 1, 1] look for exact match of last three atoms, but first atom can be a wild card
    elif test == [0, 1, 1, 1]:
        for i, j, k, l in dict2search:
            # Try as 1-2 (indexes 0, 1) be constant and varying 3-4 (indexes 2, 3) ordering
            if i[0] == '*' and j == coeff[1] and k == coeff[0] and l == coeff[2]:
                match = (i, j, k, l); return_boolean = True; break;
            elif i[0] == '*' and j == coeff[1] and k == coeff[0] and l == coeff[3]:
                match = (i, j, k, l); return_boolean = True; break;
                
            # Try as 3-2 (indexes 2, 1) be constant and varying 1-4 (indexes 0, 3) sordering
            elif i[0] == '*' and j == coeff[1] and k == coeff[2] and l == coeff[0]:
                match = (i, j, k, l); return_boolean = True; break;
            elif i[0] == '*' and j == coeff[1] and k == coeff[2] and l == coeff[3]:
                match = (i, j, k, l); return_boolean = True; break;
                
            # Try as 4-2 (indexes 3, 1) be constant and varying 1-3 (indexes 0, 2) ordering
            elif i[0] == '*' and j == coeff[1] and k == coeff[3] and l == coeff[0]:
                match = (i, j, k, l); return_boolean = True; break;
            elif i[0] == '*' and j == coeff[1] and k == coeff[3] and l == coeff[2]:
                match = (i, j, k, l); return_boolean = True; break;
                
    # If test = [0, 1, 1, 0] look for exact match of inner two atoms, but outer atoms can be a wild card
    elif test == [0, 1, 1, 0]:
        for i, j, k, l in dict2search:
            if i[0] == '*' and j == coeff[1] and k == coeff[2] and l[0] == '*':
                match = (i, j, k, l); return_boolean = True; break;
                
    # If test = [0, 1, 0, 0] look for exact match of central, but outer atoms can be a wild card
    elif test == [0, 1, 0, 0]:
        for i, j, k, l in dict2search:
            if i[0] == '*' and j == coeff[1] and k[0] == '*' and l[0] == '*':
                match = (i, j, k, l); return_boolean = True; break;

    # If user tries accessing unsupported wild card check raise exception
    else:
        log.error(f'ERROR {str(test)} wild card test not supported')
    return return_boolean, match

def match_4_body_oop(type1, type2, type3, type4, log, dict2search, equivalences, wildcard_search, form):
    """
    match_4_body has the following inputs and is for matching improper coeffs/crossterms:
        type1: atom-type1 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        type2: atom-type2 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        type3: atom-type3 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        type4: atom-type4 (MUST BE A STRING VALUE AND NOT AN EQUIVALENT FORM)
        
        dict2search: coeff dictionary to search from read_frc
        
        equivalences: coeff equivalence dictionary to find equivalences of if need
            - can be a dictionary or False. The dictionary will tell the code to use
              the equivalences and False will skip the use of the equivalences
              
        wildcard_search: True or False to tell the function to search for wild card atom-types in coeff
    """    
    # Set match/order as (type1, type2, type3, type4) to return if not found
    boolean = False; match = (type1, type2, type3, type4); 
    order = (type1, type2, type3, type4); equiv = (type1, type2, type3, type4); 
    
    # Try matching in 6-permuation orders without equivalences
    # Try as 1-2 be constant and varying 3-4 ordering
    if (type1, type2, type3, type4) in dict2search:
        match = (type1, type2, type3, type4)
        order = (type1, type2, type3, type4)
        equiv = (type1, type2, type3, type4)
        boolean = True
    elif (type1, type2, type4, type3) in dict2search:
        match = (type1, type2, type4, type3)
        order = (type1, type2, type4, type3)
        equiv = (type1, type2, type4, type3)
        boolean = True
        
    # Try as 3-2 be constant and varying 1-4 ordering
    elif (type3, type2, type1, type4) in dict2search:
        match = (type3, type2, type1, type4)
        order = (type3, type2, type1, type4)
        equiv = (type3, type2, type1, type4)
        boolean = True
    elif (type3, type2, type4, type1) in dict2search:
        match = (type3, type2, type4, type1)
        order = (type3, type2, type4, type1)
        equiv = (type3, type2, type4, type1)
        boolean = True
        
    # Try as 4-2 be constant and varying 1-3 ordering
    elif (type4, type2, type1, type3) in dict2search:
        match = (type4, type2, type1, type3)
        order = (type4, type2, type1, type3)
        equiv = (type4, type2, type1, type3)
        boolean = True
    elif (type4, type2, type3, type1) in dict2search:
        match = (type4, type2, type3, type1)
        order = (type4, type2, type3, type1)
        equiv = (type4, type2, type3, type1)
        boolean = True   
        
    # Next try using equivalences if equivalences not False
    if not boolean and equivalences:
        equiv_flag = False # set flag and if found update as True
        
        # Find equivalences for quartic or quadratic form
        if form == 'equiv':
            if type1 in equivalences and type2 in equivalences and type3 in equivalences and type4 in equivalences:
                equiv1 = equivalences[type1].oop
                equiv2 = equivalences[type2].oop
                equiv3 = equivalences[type3].oop
                equiv4 = equivalences[type4].oop
                equiv_flag = True
            else:
                failed = ', '.join([i for i in (type1, type2, type3, type4) if i not in equivalences])
                log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed))      
        elif form == 'auto-equiv':
            if type1 in equivalences and type2 in equivalences and type3 in equivalences and type4 in equivalences:
                equiv1 = equivalences[type1].oop_end
                equiv2 = equivalences[type2].oop_center
                equiv3 = equivalences[type3].oop_end
                equiv4 = equivalences[type4].oop_end
                equiv_flag = True
            else:
                failed = ', '.join([i for i in (type1, type2, type3, type4) if i not in equivalences])
                log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed))
        else:
            log.error('Unsupported equivalences are for either quartic or quadratic')
            
        # If equivalences exists start using equivalences to try matching
        if equiv_flag:    
            # Try as 1-2 be constant and varying 3-4 ordering
            if (equiv1, equiv2, equiv3, equiv4) in dict2search:
                match = (equiv1, equiv2, equiv3, equiv4)
                order = (type1, type2, type3, type4)
                equiv = (equiv1, equiv2, equiv3, equiv4)
                boolean = True
            elif (equiv1, equiv2, equiv4, equiv3) in dict2search:
                match = (equiv1, equiv2, equiv4, equiv3)
                order = (type1, type2, type4, type3)
                equiv = (equiv1, equiv2, equiv4, equiv3)
                boolean = True
                
            # Try as 3-2 be constant and varying 1-4 ordering
            elif (equiv3, equiv2, equiv1, equiv4) in dict2search:
                match = (equiv3, equiv2, equiv1, equiv4)
                order = (type3, type2, type1, type4)
                equiv = (equiv3, equiv2, equiv1, equiv4)
                boolean = True
            elif (equiv3, equiv2, equiv4, equiv1) in dict2search:
                match = (equiv3, equiv2, equiv4, equiv1)
                order = (type3, type2, type4, type1)
                equiv = (equiv3, equiv2, equiv4, equiv1)
                boolean = True
                
            # Try as 4-2 be constant and varying 1-3 ordering
            elif (equiv4, equiv2, equiv1, equiv3) in dict2search:
                match = (equiv4, equiv2, equiv1, equiv3)
                order = (type4, type2, type1, type3)
                equiv = (equiv4, equiv2, equiv1, equiv3)
                boolean = True
            elif (equiv4, equiv2, equiv3, equiv1) in dict2search:
                match = (equiv4, equiv2, equiv3, equiv1)
                order = (type4, type2, type3, type1)
                equiv = (equiv4, equiv2, equiv3, equiv1)
                boolean = True  
                
            # Try matching using wildcards in forward and reverse
            elif wildcard_search:
                # look for exact match of 3 atoms including center/central atom
                if match_improper_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [1, 1, 1, 0], log)[0]:
                   match = match_improper_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [1, 1, 1, 0], log)[1]
                   order = (type1, type2, type3, type4)
                   equiv = (equiv1, equiv2, equiv3, equiv4)
                   boolean = True 
                    
                # look for exact match of 3 atoms including center/central atom
                elif match_improper_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [0, 1, 1, 1], log)[0]:
                     match = match_improper_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [0, 1, 1, 1], log)[1]
                     order = (type1, type2, type3, type4)
                     equiv = (equiv1, equiv2, equiv3, equiv4)
                     boolean = True 
                    
                # look for exact match of 2 atoms including center/central atom
                elif match_improper_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [0, 1, 1, 0], log)[0]:
                     match = match_improper_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [0, 1, 1, 0], log)[1]
                     order = (type1, type2, type3, type4)
                     equiv = (equiv1, equiv2, equiv3, equiv4)
                     boolean = True 
                
                # Look for exact match of only center/central atom
                elif match_improper_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [0, 1, 0, 0], log)[0]:
                     match = match_improper_wildcards((equiv1, equiv2, equiv3, equiv4), dict2search, [0, 1, 0, 0], log)[1]
                     order = (type1, type2, type3, type4)
                     equiv = (equiv1, equiv2, equiv3, equiv4)
                     boolean = True  
    return boolean, match, order, equiv


##################################################
# Function to find angleangle coeff 1, 2, 3 data #
##################################################
def get_angleangle_data(type1, type2, type3, type4, log, frc):
    
    # Function to find match of angleangle data
    def match_angleangle(type1, type2, type3, type4):
        matches = [];
        
        # Try 1 2 3 4 and 4 2 3 1 ordering 1st
        if (type1, type2, type3, type4) in frc.angleangle:
            k_theta_thetap = frc.angleangle[(type1, type2, type3, type4)].k_theta_thetap
            matches.append(k_theta_thetap)
        elif (type4, type2, type3, type1) in frc.angleangle:
            k_theta_thetap = frc.angleangle[(type4, type2, type3, type1)].k_theta_thetap
            matches.append(k_theta_thetap)
            
        # Try 4 2 1 3 and 3 2 1 4 ordering 2nd
        if (type4, type2, type1, type3) in frc.angleangle:
            k_theta_thetap = frc.angleangle[(type4, type2, type1, type3)].k_theta_thetap
            matches.append(k_theta_thetap)
        elif (type3, type2, type1, type4) in frc.angleangle:
            k_theta_thetap = frc.angleangle[(type3, type2, type1, type4)].k_theta_thetap
            matches.append(k_theta_thetap)
            
        # Try 1 2 4 3 and 3 2 4 1 ordering 3rd
        if (type1, type2, type4, type3) in frc.angleangle:
            k_theta_thetap = frc.angleangle[(type1, type2, type4, type3)].k_theta_thetap
            matches.append(k_theta_thetap)
        elif (type3, type2, type4, type1) in frc.angleangle:
            k_theta_thetap = frc.angleangle[(type3, type2, type4, type1)].k_theta_thetap
            matches.append(k_theta_thetap)
        return matches
    
    # set flags (types_flag is for without equivs and equivs flag is for with equivalence)
    types_flag = False
    equivs_flag = False
    found_flag = False
        
    # Try matching without equivalences
    matches = match_angleangle(type1, type2, type3, type4)
    
    # If len(matches) >= 3: set types_flag
    if len(matches) >= 3:
        types_flag = True
    
    # If not types_flag try matching with equivalences
    if not types_flag:
        
        # Set intial flags and data structures up
        equiv_existance_flag = False # set flag and if found update as True
        
        # Find equivalences for quartic form to try getting angleangle data for
        if type1 in frc.equivalences and type2 in frc.equivalences and type3 in frc.equivalences and type4 in frc.equivalences:
            equiv1 = frc.equivalences[type1].oop
            equiv2 = frc.equivalences[type2].oop
            equiv3 = frc.equivalences[type3].oop
            equiv4 = frc.equivalences[type4].oop
            equiv_existance_flag = True
        else:
            failed = ', '.join([i for i in (type1, type2, type3, type4) if i not in frc.equivalences])
            log.warn('{} {}'.format('WARNING requesting equivalences of unsupported type:', failed)) 
            
            
        # If equivalent flag try matching
        if equiv_existance_flag:
            matches = match_angleangle(equiv1, equiv2, equiv3, equiv4)
            
            # If len(matches) >= 3: set types_flag
            if len(matches) >= 3:
                equivs_flag = True
                
    # Find matchs else set as zeros if they do not exists
    if types_flag:
        matches = [matches[i] for i in range(3)]
        found_flag = True
        equivalent = (type1, type2, type3, type4)
    elif equivs_flag:
        matches = [matches[i] for i in range(3)]
        found_flag = True
        equivalent = (equiv1, equiv2, equiv3, equiv4)
    else:
        matches = [0.0, 0.0, 0.0]; found_flag = False
        equivalent = (type1, type2, type3, type4)
    return matches, found_flag, equivalent