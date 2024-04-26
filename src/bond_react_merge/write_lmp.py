# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
March 27th, 2024
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

# Function to string type labels (No white spacing)
def string_type_labels(parameter_type):
    string = '';
    for n, i in enumerate(parameter_type):
        string += i
        if n+1 < len(parameter_type):
            string += '-'
    return string

# Function for stringing together float values for parameters
def string_parameters(coeff):
    string = ''
    for i in coeff:
        if isinstance(i, float):
            string += '{:>16.8f}'.format(i)
        elif isinstance(i, int):
            string += '{:>16}'.format(i)
        else:
            string += '{:>16}'.format(i)
    return string

# Function to split coeff types into tuple
def split_coeff(types):
    # Find coeff types and split
    types = types.strip()
    types = types.split()
    types = tuple(types)
    return types

    

# Function for writing lammps datafile
def file(newname, atom_style, m, coeff_flag, coeffs_only, new, ff_class, file, version, include_type_labels):
        
    
    # Writing new file with new atoms
    with open(newname,'w') as f: 
        # Write header
        header = '{} > bond_react_merge: {} datafile (filetag: {})'.format(m.header, version, file)
        f.write(f'{header[-220:len(header)]}\n\n') # Make max header length of 220 characters 
        
        # Write structure quantities
        if not coeffs_only:
            f.write(f'{m.natoms} atoms\n')
            f.write(f'{m.nbonds} bonds\n')
            if m.nangles > 0:
                f.write(f'{m.nangles} angles\n')
            if m.ndihedrals > 0:
                f.write(f'{m.ndihedrals} dihedrals\n')   
            if m.nimpropers > 0:
                f.write(f'{m.nimpropers} impropers\n')
            f.write('\n')
        
        # Write structure quantity types
        if coeff_flag:
            f.write(f'{new.natomtypes} atom types\n')
            f.write(f'{new.nbondtypes} bond types\n')
            f.write(f'{new.nangletypes} angle types\n')
            f.write(f'{new.ndihedraltypes} dihedral types\n')
            f.write(f'{new.nimpropertypes} improper types\n')
            f.write('\n')
        
        # write box size
        if not coeffs_only: f.write('{0}\n{1}\n{2}\n'.format(m.xbox_line, m.ybox_line, m.zbox_line))
        
        
        # Write coeffs if coeff flag
        if coeff_flag:
            
            # Write Atom Type Labels if user wants
            if include_type_labels:
                f.write('\nAtom Type Labels\n\n')
                for i in new.masses:
                    mass = new.masses[i]
                    f.write('{:^3} {:^2}\n'.format(i, mass.type))
                    
            # Write Bond Type Labels if user wants
            if include_type_labels and new.bond_coeffs:
                f.write('\nBond Type Labels\n\n')
                for i in new.bond_coeffs:
                    bond = new.bond_coeffs[i]
                    f.write('{:^3} {:^2}\n'.format(i, string_type_labels(bond.type)))
                    
            # Write Angle Type Labels if user wants
            if include_type_labels and new.angle_coeffs:
                f.write('\nAngle Type Labels\n\n')
                for i in new.angle_coeffs: 
                    angle = new.angle_coeffs[i]
                    f.write('{:^3} {:^2}\n'.format(i, string_type_labels(angle.type)))    
    
            # Write Dihedral Type Labels if user wants
            if include_type_labels and new.dihedral_coeffs:
                f.write('\nDihedral Type Labels\n\n')
                for i in new.dihedral_coeffs: 
                    dihedral = new.dihedral_coeffs[i]
                    f.write('{:^3} {:^2}\n'.format(i, string_type_labels(dihedral.type))) 
    
            # Write Improper Type Labels if user wants
            if include_type_labels and new.improper_coeffs:
                f.write('\nImproper Type Labels\n\n')
                for i in new.improper_coeffs: 
                    improper = new.improper_coeffs[i]
                    f.write('{:^3} {:^2}\n'.format(i, string_type_labels(improper.type))) 


            # Write massses
            #f.write(f'\nMasses # {new.mass_coeffs_style_hint}\n\n')
            f.write('\nMasses\n\n')
            for i in new.masses: 
                coeff = new.masses[i]
                parms = coeff.coeffs
                comment = '{:^2} {:5}'.format('#', coeff.type)
                f.write('{:^3} {:^10.5f} {:^2}\n'.format(i, parms, comment))
                
            # Write pair coeffs
            if new.pair_coeffs:
                if m.pair_coeffs_style_hint != 'N/A':
                    f.write(f'\nPair Coeffs  # {m.pair_coeffs_style_hint}\n\n')
                else:
                    f.write('\nPair Coeffs\n\n')
                for i in new.pair_coeffs: 
                    pair = new.pair_coeffs[i]
                    comment = '{:^2} {:5}'.format('#', pair.type)
                    f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(pair.coeffs), comment))
    
                
            # Write bond coeffs
            if new.bond_coeffs:
                if m.bond_coeffs_style_hint != 'N/A':
                    f.write(f'\nBond Coeffs  # {m.bond_coeffs_style_hint}\n\n')
                else:
                    f.write('\nBond Coeffs\n\n')
                for i in new.bond_coeffs: 
                    bond = new.bond_coeffs[i]                
                    coeff_type = string_parameter_type(bond.type)
                    comment = '{:^2} {:10}'.format('#', coeff_type)
                    f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(bond.coeffs), comment))
                
                
            # Write angle coeffs
            if new.angle_coeffs:
                if m.angle_coeffs_style_hint != 'N/A':
                    f.write(f'\nAngle Coeffs  # {m.angle_coeffs_style_hint}\n\n')
                else:
                    f.write('\nAngle Coeffs\n\n')
                for i in new.angle_coeffs: 
                    angle = new.angle_coeffs[i]
                    coeff_type = string_parameter_type(angle.type)
                    comment = '{:^2} {:15}'.format('#', coeff_type)
                    f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(angle.coeffs), comment))
    
                
            # Write dihedral coeffs
            if new.dihedral_coeffs:
                if m.dihedral_coeffs_style_hint != 'N/A':
                    f.write(f'\nDihedral Coeffs  # {m.dihedral_coeffs_style_hint}\n\n')
                else:
                    f.write('\nDihedral Coeffs\n\n')
                for i in new.dihedral_coeffs: 
                    dihedral = new.dihedral_coeffs[i]
                    coeff_type = string_parameter_type(dihedral.type)
                    comment = '{:^2} {:15}'.format('#', coeff_type)
                    f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(dihedral.coeffs), comment))
               
                
            # Write Improper coeffs
            if new.improper_coeffs:
                if m.improper_coeffs_style_hint != 'N/A':
                    f.write(f'\nImproper Coeffs  # {m.improper_coeffs_style_hint}\n\n')
                else:
                    f.write('\nImproper Coeffs\n\n')
                for i in new.improper_coeffs: 
                    improper = new.improper_coeffs[i]                    
                    coeff_type = string_parameter_type(improper.type)
                    comment = '{:^2} {:15}'.format('#', coeff_type)
                    f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(improper.coeffs), comment))
    
            # write crossterms if class2 ff
            if ff_class == 2:                

                # Write bondbond coeffs
                if new.angle_coeffs and new.bondbond_coeffs:
                    if m.bondbond_coeffs_style_hint != 'N/A':
                        f.write(f'\nBondBond Coeffs  # {m.bondbond_coeffs_style_hint}\n\n')
                    else:
                        f.write('\nBondBond Coeffs\n\n')
                    for i in new.bondbond_coeffs: 
                        bondbond = new.bondbond_coeffs[i]
                        coeff_type = string_parameter_type(bondbond.type)
                        comment = '{:^2} {:15}'.format('#', coeff_type)
                        f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(bondbond.coeffs), comment))
        
                    
                # Write bondangle coeffs
                if new.angle_coeffs and new.bondangle_coeffs:
                    if m.bondbond_coeffs_style_hint != 'N/A':
                        f.write(f'\nBondAngle Coeffs  # {m.bondangle_coeffs_style_hint}\n\n')
                    else:
                        f.write('\nBondAngle Coeffs\n\n')
                    for i in new.bondangle_coeffs: 
                        bondangle = new.bondangle_coeffs[i]
                        coeff_type = string_parameter_type(bondangle.type)
                        comment = '{:^2} {:15}'.format('#', coeff_type)
                        f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(bondangle.coeffs), comment))
                    
                    
                # Write angleangletorsion coeffs
                if new.dihedral_coeffs and new.angleangletorsion_coeffs:
                    if m.angleangletorsion_coeffs_style_hint != 'N/A':
                        f.write(f'\nAngleAngleTorsion Coeffs  # {m.angleangletorsion_coeffs_style_hint}\n\n')
                    else:
                        f.write('\nAngleAngleTorsion Coeffs\n\n')
                    for i in new.angleangletorsion_coeffs: 
                        angleangletorsion = new.angleangletorsion_coeffs[i]
                        coeff_type = string_parameter_type(angleangletorsion.type)
                        comment = '{:^2} {:15}'.format('#', coeff_type)
                        f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(angleangletorsion.coeffs), comment))
        
                    
                # Write endbondtorsion coeffs
                if new.dihedral_coeffs and new.endbondtorsion_coeffs:
                    if m.endbondtorsion_coeffs_style_hint != 'N/A':
                        f.write(f'\nEndBondTorsion Coeffs  # {m.endbondtorsion_coeffs_style_hint}\n\n')
                    else:
                        f.write('\nEndBondTorsion Coeffs\n\n')
                    for i in new.endbondtorsion_coeffs: 
                        endbondtorsion = new.endbondtorsion_coeffs[i]
                        coeff_type = string_parameter_type(endbondtorsion.type)
                        comment = '{:^2} {:15}'.format('#', coeff_type)
                        f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(endbondtorsion.coeffs), comment))
        
        
                # Write middlebondtorsion coeffs
                if new.dihedral_coeffs and new.middlebondtorsion_coeffs:
                    if m.middlebondtorsion_coeffs_style_hint != 'N/A':
                        f.write(f'\nMiddleBondTorsion Coeffs  # {m.middlebondtorsion_coeffs_style_hint}\n\n')
                    else:
                        f.write('\nMiddleBondTorsion Coeffs\n\n')
                    for i in new.middlebondtorsion_coeffs: 
                        middlebondtorsion  = new.middlebondtorsion_coeffs[i]
                        coeff_type = string_parameter_type(middlebondtorsion.type)
                        comment = '{:^2} {:15}'.format('#', coeff_type)
                        f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(middlebondtorsion.coeffs), comment))
        
                    
                # Write bondbond13 coeffs
                if new.dihedral_coeffs and new.bondbond13_coeffs:
                    if m.bondbond13_coeffs_style_hint != 'N/A':
                        f.write(f'\nBondBond13 Coeffs  # {m.bondbond13_coeffs_style_hint}\n\n')
                    else:
                        f.write('\nBondBond13 Coeffs\n\n')
                    for i in new.bondbond13_coeffs: 
                        bondbond13 = new.bondbond13_coeffs[i]
                        coeff_type = string_parameter_type(bondbond13.type)
                        comment = '{:^2} {:15}'.format('#', coeff_type)
                        f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(bondbond13.coeffs), comment))
        
        
                # Write angletorsion coeffs
                if new.dihedral_coeffs and new.angletorsion_coeffs:
                    if m.angletorsion_coeffs_style_hint != 'N/A':
                        f.write(f'\nAngleTorsion Coeffs  # {m.angletorsion_coeffs_style_hint}\n\n')
                    else:
                        f.write('\nAngleTorsion Coeffs\n\n')
                    for i in new.angletorsion_coeffs: 
                        angletorsion = new.angletorsion_coeffs[i]
                        coeff_type = string_parameter_type(angletorsion.type)
                        comment = '{:^2} {:15}'.format('#', coeff_type)
                        f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(angletorsion.coeffs), comment))
    
                    
                # Write angleangle coeffs
                if new.improper_coeffs and new.angleangle_coeffs:
                    if m.angleangle_coeffs_style_hint != 'N/A':
                        f.write(f'\nAngleAngle Coeffs  # {m.angleangle_coeffs_style_hint}\n\n')
                    else:
                        f.write('\nAngleAngle Coeffs\n\n')
                    for i in new.angleangle_coeffs: 
                        angleangle = new.angleangle_coeffs[i]
                        coeff_type = string_parameter_type(angleangle.type)
                        comment = '{:^2} {:15}'.format('#', coeff_type)
                        f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(angleangle.coeffs), comment))

        #################################################
        # write all other quantities if not coeffs_only #
        #################################################
        if not coeffs_only:
        
            # Write atoms    
            f.write(f'\nAtoms # {atom_style}\n\n')            
            for i in m.atoms:
                # Set new type and update later on
                new_type = 0
                    
                atom = m.atoms[i]
    
                # Find str_type by compairing to pair_coeffs
                coeff = m.pair_coeffs[atom.type]
                str_type = coeff.type
                
                # Find new type via mapping from atom_types_map
                new_type = int(new.atom_types_map[str_type])
                if include_type_labels: new_type = '{t:<{s}}'.format(t=str_type, s=5)
                
                comment = '{:^2} {:5}'.format('#', atom.comment)
                
                # write charge atom style
                if atom_style == 'charge':
                    f.write('{:^6} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2} {:^2}\n'.format(i, new_type, atom.charge, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment)) 
                
                # write molecular atom style
                elif atom_style == 'molecular':
                    f.write('{:^6} {:^4} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2} {:^2}\n'.format(i, atom.molid, new_type, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment))                          
                
                # write full atom style
                elif atom_style == 'full':
                    f.write('{:^6} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2} {:^2}\n'.format(i, atom.molid, new_type, atom.charge, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment)) 
                
                # write angle or bond atom style
                elif atom_style == 'angle' or atom_style == 'bond':
                    f.write('{:^6} {:^4} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.molid, new_type, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment)) 
                    
                # write atomic atom style
                elif atom_style == 'atomic':
                    f.write('{:^6} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, new_type, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment)) 
                    
                # write dipole atom style
                elif atom_style == 'dipole':
                    try: mux = atom.mux; muy = atom.muy; muz = atom.muz;
                    except:
                        mux = 0; muy = 0; muz = 0;
                        print(f'WARNING writing LAMMPS datafile Atom section in dipole style and initializing mux = {mux}, muy = {muy}, muz = {muz}')
                    f.write('{:^6} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, new_type, atom.charge, atom.x, atom.y, atom.z, mux, muy, muz, atom.ix, atom.iy, atom.iz, comment)) 
                        
                # write charge atom style
                elif atom_style == 'dpd':
                    try: theta = atom.theta
                    except:
                        theta = 0
                        print(f'WARNING writing LAMMPS datafile Atom section in dpd style and initializing theta = {theta}')
                    f.write('{:^6} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, new_type, theta, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment))  
                    
                # write line style
                elif atom_style == 'line':
                    try: lineflag = atom.lineflag; density = atom.density;
                    except:
                        lineflag = 0; density = 0;
                        print(f'WARNING writing LAMMPS datafile Atom section in line style and initializing lineflag = {lineflag}, density = {density}')
                    f.write('{:^6} {:^4} {:^2} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.molid, new_type, lineflag, density, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment))    
        
                # else raise and exception
                else:
                    raise Exception('Atom Style Error - write_lmp.py does not support atom style {atom_style}')
                
            # Write bonds
            if m.nbonds > 0:
                f.write('\nBonds\n\n')
                for i in m.bonds:
                    # Set new type and update later on
                    new_type = 0
                        
                    bond = m.bonds[i]
                    id1, id2 = bond.atomids
                    
                    # Find str_type by compairing to bond_coeffs
                    coeff = m.bond_coeffs[bond.type]
                    str_type = coeff.type
                    str_type_tuple = split_coeff(str_type)
                    
                    # Find new type via mapping from atom_types_map
                    new_type = int(new.bond_types_map[str_type_tuple])
                    if include_type_labels: new_type = '{t:<{s}}'.format(t=string_type_labels(str_type_tuple), s=10)
                    
                    # build comment and write results
                    comment = '{:^2} {:5}'.format('#', str_type)
                    f.write('{:^2} {:^2} {:^2} {:^2} {:^6}\n'.format(i, new_type, id1, id2, comment))
                
            # Write angles
            if m.nangles > 0:
                f.write('\nAngles\n\n')
                for i in m.angles:
                    # Set new type and update later on
                    new_type = 0
                        
                    angle = m.angles[i]
                    id1, id2, id3 = angle.atomids
                    
                    # Find str_type by compairing to angle_coeffs
                    coeff = m.angle_coeffs[angle.type]
                    str_type = coeff.type
                    str_type_tuple = split_coeff(str_type)
                    
                    # Find new type via mapping from atom_types_map
                    new_type = int(new.angle_types_map[str_type_tuple])
                    if include_type_labels: new_type = '{t:<{s}}'.format(t=string_type_labels(str_type_tuple), s=15)
                    
                    # build comment and write results
                    comment = '{:^2} {:5}'.format('#', str_type)
                    f.write('{:^2} {:^2} {:^2} {:^2} {:^2} {:^6}\n'.format(i, new_type, id1, id2, id3, comment))
                
            # Write dihedrals
            if m.ndihedrals > 0:
                f.write('\nDihedrals\n\n')
                for i in m.dihedrals:
                    # Set new type and update later on
                    new_type = 0
                    
                    dihedral = m.dihedrals[i]
                    id1, id2, id3, id4 = dihedral.atomids
                    
                    # Find str_type by compairing to dihedral_coeffs
                    coeff = m.dihedral_coeffs[dihedral.type]
                    str_type = coeff.type
                    str_type_tuple = split_coeff(str_type)
                    
                    # Find new type via mapping from atom_types_map
                    new_type = int(new.dihedral_types_map[str_type_tuple])
                    if include_type_labels: new_type = '{t:<{s}}'.format(t=string_type_labels(str_type_tuple), s=20)
                    
                    # build comment and write results
                    comment = '{:^2} {:5}'.format('#', str_type)
                    f.write('{:^2} {:^2} {:^2} {:^2} {:^2} {:^2} {:^6}\n'.format(i, new_type, id1, id2, id3, id4, comment))
                        
            # Write impropers
            if m.nimpropers > 0:
                f.write('\nImpropers\n\n')
                for i in m.impropers:
                    # Set new type and update later on
                    new_type = 0
                    
                    improper = m.impropers[i]
                    id1, id2, id3, id4 = improper.atomids
                    
                    # Find str_type by compairing to improper_coeffs
                    coeff = m.improper_coeffs[improper.type]
                    str_type = coeff.type
                    str_type_tuple = split_coeff(str_type)
                    
                    # Find new type via mapping from atom_types_map
                    new_type = int(new.improper_types_map[str_type_tuple])
                    if include_type_labels: new_type = '{t:<{s}}'.format(t=string_type_labels(str_type_tuple), s=25)
                    
                    # build comment and write results
                    comment = '{:^2} {:5}'.format('#', str_type)
                    f.write('{:^2} {:^2} {:^2} {:^2} {:^2} {:^2} {:^6}\n'.format(i, new_type, id1, id2, id3, id4, comment))   
    return
