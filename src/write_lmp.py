# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
April 1st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
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

# Function to string type labels (No white spacing)
def string_type_labels(parameter_type):
    string = '';
    for n, i in enumerate(parameter_type):
        string += i
        if n+1 < len(parameter_type):
            string += '-'
    return string

# Function to split coeff types into tuple
def split_coeff(types):
    # Find coeff types and split
    types = types.strip()
    types = types.split()
    types = tuple(types)
    return types    

# Function for writing lammps datafile
def file(m, filename, header, atom_style, include_type_labels, log, force_field_only=False):
        
    # Writing m file with m atoms
    with open(filename,'w') as f: 
        # Write header
        f.write(f'{header[-220:len(header)]}\n\n') # Make max header length of 220 characters 
        
        # Write structure quantities
        if not force_field_only:
            f.write(f'{m.natoms} atoms\n')
            if m.nbonds > 0: f.write(f'{m.nbonds} bonds\n')
            if m.nangles > 0: f.write(f'{m.nangles} angles\n')
            if m.ndihedrals > 0: f.write(f'{m.ndihedrals} dihedrals\n')   
            if m.nimpropers > 0: f.write(f'{m.nimpropers} impropers\n')
            f.write('\n')
        
        # Write structure quantity types
        f.write(f'{m.natomtypes} atom types\n')
        if m.nbondtypes > 0: f.write(f'{m.nbondtypes} bond types\n')
        if m.nangletypes > 0: f.write(f'{m.nangletypes} angle types\n')
        if m.ndihedraltypes > 0: f.write(f'{m.ndihedraltypes} dihedral types\n')
        if m.nimpropertypes > 0: f.write(f'{m.nimpropertypes} improper types\n')
        
        # write box size
        if not force_field_only:
            f.write('\n{0}\n{1}\n{2}\n'.format(m.xbox_line, m.ybox_line, m.zbox_line))
            if m.xy != 0 or m.xz != 0 or m.yz != 0:
                f.write('{:>12.9f} {:^9.9f} {:^9.9f} {} {} {}\n'.format(m.xy, m.xz, m.yz, 'xy', 'xz', 'yz'))


        # Write Atom Type Labels if user wants
        if include_type_labels:
            f.write('\nAtom Type Labels\n\n')
            used = set()
            for i in m.masses:
                mass = m.masses[i]
                if i in m.atom_type_labels_reverse: typelabel = m.atom_type_labels_reverse[i]
                elif mass.type != 'N/A': typelabel = mass.type
                else: log.error('\nERROR datafile(s) do NOT have all2lmp style comments for Masses or Type Labels are not defined')
                f.write('{:^3} {:^2}\n'.format(i, typelabel))
                if typelabel in used:
                    log.warn('WARNING Atom Type Label: {} has been defined multiple times.'.format(typelabel))
                used.add(typelabel)
                
        # Write Bond Type Labels if user wants
        if include_type_labels and m.bond_coeffs:
            f.write('\nBond Type Labels\n\n')
            used = set()
            for i in m.bond_coeffs:
                bond = m.bond_coeffs[i]
                if i in m.bond_type_labels_reverse: typelabel = m.bond_type_labels_reverse[i]
                elif bond.type != 'N/A' and len(split_coeff(bond.type)) == 2: typelabel = string_type_labels( split_coeff(bond.type) )
                else: log.error('\nERROR datafile(s) do NOT have all2lmp style comments for Bond Coeffs or Type Labels are not defined')
                f.write('{:^3} {:^2}\n'.format(i, typelabel))
                if typelabel in used:
                    log.warn('WARNING Bond Type Label: {} has been defined multiple times.'.format(typelabel))
                used.add(typelabel)
                
        # Write Angle Type Labels if user wants
        if include_type_labels and m.angle_coeffs:
            f.write('\nAngle Type Labels\n\n')
            used = set()
            for i in m.angle_coeffs:
                angle = m.angle_coeffs[i]
                if i in m.angle_type_labels_reverse: typelabel = m.angle_type_labels_reverse[i]
                elif angle.type != 'N/A' and len(split_coeff(angle.type)) == 3: typelabel = string_type_labels( split_coeff(angle.type) )
                else: log.error('\nERROR datafile(s) do NOT have all2lmp style comments for Angle Coeffs or Type Labels are not defined')
                f.write('{:^3} {:^2}\n'.format(i, typelabel))
                if typelabel in used:
                    log.warn('WARNING Angle Type Label: {} has been defined multiple times.'.format(typelabel))
                used.add(typelabel)
                
        # Write Dihedral Type Labels if user wants
        if include_type_labels and m.dihedral_coeffs:
            f.write('\nDihedral Type Labels\n\n')
            used = set()
            for i in m.dihedral_coeffs:
                dihedral = m.dihedral_coeffs[i]
                if i in m.dihedral_type_labels_reverse: typelabel = m.dihedral_type_labels_reverse[i]
                elif dihedral.type != 'N/A' and len(split_coeff(dihedral.type)) == 4: typelabel = string_type_labels( split_coeff(dihedral.type) )
                else: log.error('\nERROR datafile(s) do NOT have all2lmp style comments for Dihedral Coeffs or Type Labels are not defined')
                f.write('{:^3} {:^2}\n'.format(i, typelabel))
                if typelabel in used:
                    log.warn('WARNING Dihedral Type Label: {} has been defined multiple times.'.format(typelabel))
                used.add(typelabel)
                
        # Write Improper Type Labels if user wants
        if include_type_labels and m.improper_coeffs:
            f.write('\nImproper Type Labels\n\n')
            used = set()
            for i in m.improper_coeffs:
                improper = m.improper_coeffs[i]
                if i in m.improper_type_labels_reverse: typelabel = m.improper_type_labels_reverse[i]
                elif improper.type != 'N/A' and len(split_coeff(improper.type)) == 5: typelabel = string_type_labels( split_coeff(improper.type) )
                else: log.error('\nERROR datafile(s) do NOT have all2lmp style comments for Improper Coeffs or Type Labels are not defined')
                f.write('{:^3} {:^2}\n'.format(i, typelabel))
                if typelabel in used:
                    log.warn('WARNING Improper Type Label: {} has been defined multiple times.'.format(typelabel))
                used.add(typelabel)

        # Write massses
        #f.write(f'Masses  # {m.mass_coeffs_style_hint}\n\n')
        f.write('\nMasses\n\n')
        for i in m.masses: 
            coeff = m.masses[i]
            parms = coeff.coeffs
            if coeff.type != 'N/A':
                f.write('{:^3} {:^10.5f} # {}\n'.format(i, parms[0], coeff.type))
            else:
                f.write('{:^3} {:^10.5f}\n'.format(i, parms[0]))
            
        # Write pair coeffs
        if m.pair_coeffs:
            if m.pair_coeffs_style_hint != 'N/A':
                f.write(f'\nPair Coeffs  # {m.pair_coeffs_style_hint}\n\n')
            else:
                f.write('\nPair Coeffs\n\n')
            for i in m.pair_coeffs: 
                pair = m.pair_coeffs[i]
                if pair.type != 'N/A':
                    f.write('{:^3} {} # {}\n'.format(i, string_parameters(pair.coeffs), pair.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(pair.coeffs)))

            
        # Write bond coeffs
        if m.bond_coeffs:
            if m.bond_coeffs_style_hint != 'N/A':
                f.write(f'\nBond Coeffs  # {m.bond_coeffs_style_hint}\n\n')
            else:
                f.write('\nBond Coeffs\n\n')
            for i in m.bond_coeffs: 
                bond = m.bond_coeffs[i]
                if bond.type != 'N/A':
                    f.write('{:^3} {:<80} # {}\n'.format(i, string_parameters(bond.coeffs), bond.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(bond.coeffs)))
            
            
        # Write angle coeffs
        if m.angle_coeffs:
            if m.angle_coeffs_style_hint != 'N/A':
                f.write(f'\nAngle Coeffs  # {m.angle_coeffs_style_hint}\n\n')
            else:
                f.write('\nAngle Coeffs\n\n')
            for i in m.angle_coeffs: 
                angle = m.angle_coeffs[i]
                if angle.type != 'N/A':
                    f.write('{:^3} {} # {}\n'.format(i, string_parameters(angle.coeffs), angle.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(angle.coeffs)))

            
        # Write dihedral coeffs
        if m.dihedral_coeffs:
            if m.dihedral_coeffs_style_hint != 'N/A':
                f.write(f'\nDihedral Coeffs  # {m.dihedral_coeffs_style_hint}\n\n')
            else:
                f.write('\nDihedral Coeffs\n\n')
            for i in m.dihedral_coeffs: 
                dihedral = m.dihedral_coeffs[i]
                if dihedral.type != 'N/A':
                    f.write('{:^3} {} # {}\n'.format(i, string_parameters(dihedral.coeffs), dihedral.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(dihedral.coeffs)))
           
            
        # Write Improper coeffs
        if m.improper_coeffs:
            if m.improper_coeffs_style_hint != 'N/A':
                f.write(f'\nImproper Coeffs  # {m.improper_coeffs_style_hint}\n\n')
            else:
                f.write('\nImproper Coeffs\n\n')
            for i in m.improper_coeffs: 
                improper = m.improper_coeffs[i] 
                if improper.type != 'N/A':                   
                    f.write('{:^3} {} # {}\n'.format(i, string_parameters(improper.coeffs), improper.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(improper.coeffs)))


        # Write bondbond coeffs
        if m.angle_coeffs and m.bondbond_coeffs:
            if m.bondbond_coeffs_style_hint != 'N/A':
                f.write(f'\nBondBond Coeffs  # {m.bondbond_coeffs_style_hint}\n\n')
            else:
                f.write('\nBondBond Coeffs\n\n')
            for i in m.bondbond_coeffs: 
                bondbond = m.bondbond_coeffs[i]
                if bondbond.type != 'N/A':
                    f.write('{:^3} {} # {}\n'.format(i, string_parameters(bondbond.coeffs), bondbond.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(bondbond.coeffs)))

            
        # Write bondangle coeffs
        if m.angle_coeffs and m.bondangle_coeffs:
            if m.bondangle_coeffs_style_hint != 'N/A':
                f.write(f'\nBondAngle Coeffs  # {m.bondangle_coeffs_style_hint}\n\n')
            else:
                f.write('\nBondAngle Coeffs\n\n')
            for i in m.bondangle_coeffs: 
                bondangle = m.bondangle_coeffs[i]
                if bondangle.type != 'N/A':
                    f.write('{:^3} {} # {}\n'.format(i, string_parameters(bondangle.coeffs), bondangle.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(bondangle.coeffs)))
            
        # Write angleangletorsion coeffs
        if m.dihedral_coeffs and m.angleangletorsion_coeffs:
            if m.angleangletorsion_coeffs_style_hint != 'N/A':
                f.write(f'\nAngleAngleTorsion Coeffs  # {m.angleangletorsion_coeffs_style_hint}\n\n')
            else:
                f.write('\nAngleAngleTorsion Coeffs\n\n')
            for i in m.angleangletorsion_coeffs: 
                angleangletorsion = m.angleangletorsion_coeffs[i]
                if angleangletorsion.type != 'N/A':
                    f.write('{:^3} {} # {}\n'.format(i, string_parameters(angleangletorsion.coeffs), angleangletorsion.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(angleangletorsion.coeffs)))

            
        # Write endbondtorsion coeffs
        if m.dihedral_coeffs and m.endbondtorsion_coeffs:
            if m.endbondtorsion_coeffs_style_hint != 'N/A':
                f.write(f'\nEndBondTorsion Coeffs  # {m.endbondtorsion_coeffs_style_hint}\n\n')
            else:
                f.write('\nEndBondTorsion Coeffs\n\n')
        for i in m.endbondtorsion_coeffs: 
            endbondtorsion = m.endbondtorsion_coeffs[i]
            if endbondtorsion.type != 'N/A':
                f.write('{:^3} {} # {}\n'.format(i, string_parameters(endbondtorsion.coeffs), endbondtorsion.type))
            else:
                f.write('{:^3} {}\n'.format(i, string_parameters(endbondtorsion.coeffs)))


        # Write middlebondtorsion coeffs
        if m.dihedral_coeffs and m.middlebondtorsion_coeffs:
            if m.middlebondtorsion_coeffs_style_hint != 'N/A':
                f.write(f'\nMiddleBondTorsion Coeffs  # {m.middlebondtorsion_coeffs_style_hint}\n\n')
            else:
                f.write('\nMiddleBondTorsion Coeffs\n\n')
            for i in m.middlebondtorsion_coeffs: 
                middlebondtorsion  = m.middlebondtorsion_coeffs[i]
                if middlebondtorsion.type != 'N/A':
                    f.write('{:^3} {} # {}\n'.format(i, string_parameters(middlebondtorsion.coeffs), middlebondtorsion.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(middlebondtorsion.coeffs)))

            
        # Write bondbond13 coeffs
        if m.dihedral_coeffs and m.bondbond13_coeffs:
            if m.bondbond13_coeffs_style_hint != 'N/A':
                f.write(f'\nBondBond13 Coeffs  # {m.bondbond13_coeffs_style_hint}\n\n')
            else:
                f.write('\nBondBond13 Coeffs\n\n')
            for i in m.bondbond13_coeffs: 
                bondbond13 = m.bondbond13_coeffs[i]
                if bondbond13.type != 'N/A':
                    f.write('{:^3} {} # {}\n'.format(i, string_parameters(bondbond13.coeffs), bondbond13.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(bondbond13.coeffs)))


        # Write angletorsion coeffs
        if m.dihedral_coeffs and m.angletorsion_coeffs:
            if m.angletorsion_coeffs_style_hint != 'N/A':
                f.write(f'\nAngleTorsion Coeffs  # {m.angletorsion_coeffs_style_hint}\n\n')
            else:
                f.write('\nAngleTorsion Coeffs\n\n')
            for i in m.angletorsion_coeffs: 
                angletorsion = m.angletorsion_coeffs[i]
                if angletorsion.type != 'N/A':
                    f.write('{:^3} {} # {}\n'.format(i, string_parameters(angletorsion.coeffs), angletorsion.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(angletorsion.coeffs)))

            
        # Write angleangle coeffs
        if m.improper_coeffs and m.angleangle_coeffs:
            if m.angleangle_coeffs_style_hint != 'N/A':
                f.write(f'\nAngleAngle Coeffs  # {m.angleangle_coeffs_style_hint}\n\n')
            else:
                f.write('\nAngleAngle Coeffs\n\n')
            for i in m.angleangle_coeffs: 
                angleangle = m.angleangle_coeffs[i]
                if angleangle.type != 'N/A':
                    f.write('{:^3} {} # {}\n'.format(i, string_parameters(angleangle.coeffs), angleangle.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(angleangle.coeffs)))

        ###############################################
        # Cut out and atoms, bonds, angles, dihedrals #
        # and impropers force_field_only              #
        ###############################################
        if not force_field_only:
            # Write atoms    
            f.write(f'\nAtoms # {atom_style}\n\n')            
            for i in m.atoms:
                atom = m.atoms[i]
                
                if atom.comment != 'N/A':
                    comment = '# {}'.format(atom.comment)
                else: comment = ''
                
                # if include_type_labels and atom.comment != 'N/A':
                #     atomtype = '{t:<{s}}'.format(t=m.masses[atom.type].type, s=5)
                # else: atomtype = atom.type
                
                # write charge atom style
                if atom_style == 'charge':
                    f.write('{:^6} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2} {}\n'.format(i, atom.type, atom.charge, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment)) 
                
                # write molecular atom style
                elif atom_style == 'molecular':
                    f.write('{:^6} {:^4} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2} {}\n'.format(i, atom.molid, atom.type, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment))                          
                
                # write full atom style
                elif atom_style == 'full':
                    f.write('{:^6} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2} {}\n'.format(i, atom.molid, atom.type, atom.charge, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment)) 
                
                # write angle or bond atom style
                elif atom_style == 'angle' or atom_style == 'bond':
                    f.write('{:^6} {:^4} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.molid, atom.type, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment)) 
                    
                # write atomic atom style
                elif atom_style == 'atomic':
                    f.write('{:^6} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.type, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment)) 
                    
                # write dipole atom style
                elif atom_style == 'dipole':
                    try: mux = atom.mux; muy = atom.muy; muz = atom.muz;
                    except:
                        mux = 0; muy = 0; muz = 0;
                        log.warn(f'WARNING writing LAMMPS datafile Atom section in dipole style and initializing mux = {mux}, muy = {muy}, muz = {muz}')
                    f.write('{:^6} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.type, atom.charge, atom.x, atom.y, mux, muy, muz, atom.z, atom.ix, atom.iy, atom.iz, comment)) 
                        
                # write charge atom style
                elif atom_style == 'dpd':
                    try: theta = atom.theta
                    except:
                        theta = 0
                        log.warn(f'WARNING writing LAMMPS datafile Atom section in dpd style and initializing theta = {theta}')
                    f.write('{:^6} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.type, theta, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment))   
                    
                # write line style
                elif atom_style == 'line':
                    try: lineflag = atom.lineflag; density = atom.density;
                    except:
                        lineflag = 0; density = 0;
                        log.warn(f'WARNING writing LAMMPS datafile Atom section in line style and initializing lineflag = {lineflag}, density = {density}')
                    f.write('{:^6} {:^4} {:^2} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.molid, atom.type, lineflag, density, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment))
                
                # else raise and exception
                else:
                    log.error('Atom Style Error - write_lmp.py does not support atom style {atom_style}')
                    
            # Write velocities
            try:
                if m.velocities:
                    velocities = {m.velocities[i] for i in m.velocities}
                    atoms_in_velocities = True
                    for i in m.atoms:
                        try: vel = m.velocities[i]
                        except: atoms_in_velocities = False; break;
                    if atoms_in_velocities and len(velocities) > 1 and (0, 0, 0) not in velocities and len(m.velocities) == len(m.atoms):
                        f.write('\nVelocities\n\n')
                        for i in m.velocities:
                            vx, vy, vz = m.velocities[i]
                            f.write('{:^6} {:^15.9f} {:^15.9f} {:^15.9f}\n'.format(i, vx, vy, vz))
            except: pass
    
    
                
            # Write bonds
            if m.nbonds > 0:
                f.write("\nBonds\n\n")
                for i in m.bonds:
                    bond = m.bonds[i]
                    id1, id2 = bond.atomids
                    # if include_type_labels and m.bond_coeffs[bond.type].type != 'N/A':
                    #     bondtype = '{t:<{s}}'.format(t=string_type_labels(split_coeff(m.bond_coeffs[bond.type].type)), s=10)
                    # else: bondtype = bond.type
                    f.write('{:^2} {:^2} {:^2} {:^2}\n'.format(i, bond.type, id1, id2))
                
            # Write angles
            if m.nangles > 0:
                f.write("\nAngles\n\n")
                for i in m.angles:
                    angle = m.angles[i]
                    id1, id2, id3 = angle.atomids
                    # if include_type_labels and m.angle_coeffs[angle.type].type != 'N/A':
                    #     angletype = '{t:<{s}}'.format(t=string_type_labels(split_coeff(m.angle_coeffs[angle.type].type)), s=15)
                    # else: angletype = angle.type
                    f.write('{:^2} {:^2} {:^2} {:^2} {:^2}\n'.format(i, angle.type, id1, id2, id3))
                
            # Write dihedrals
            if m.ndihedrals > 0:
                f.write("\nDihedrals\n\n")
                for i in m.dihedrals:
                    dihedral = m.dihedrals[i]
                    id1, id2, id3, id4 = dihedral.atomids
                    # if include_type_labels and m.dihedral_coeffs[dihedral.type].type != 'N/A':
                    #     dihedraltype = '{t:<{s}}'.format(t=string_type_labels(split_coeff(m.dihedral_coeffs[dihedral.type].type)), s=20)
                    # else: dihedraltype = dihedral.type
                    f.write('{:^2} {:^2} {:^2} {:^2} {:^2} {:^2}\n'.format(i, dihedral.type, id1, id2, id3, id4))
                        
            # Write impropers
            if m.nimpropers > 0:
                f.write("\nImpropers\n\n")
                for i in m.impropers:
                    improper = m.impropers[i]
                    id1, id2, id3, id4 = improper.atomids
                    # if include_type_labels and m.improper_coeffs[improper.type].type != 'N/A':
                    #     impropertype = '{t:<{s}}'.format(t=string_type_labels(split_coeff(m.improper_coeffs[improper.type].type)), s=20)
                    # else: impropertype = improper.type
                    f.write('{:^2} {:^2} {:^2} {:^2} {:^2} {:^2}\n'.format(i, improper.type, id1, id2, id3, id4))