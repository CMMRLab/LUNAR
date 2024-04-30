#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.7
December 1st, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""


# Function for stringing together parameter types
def string_parameter_type(parameter_type):
    if parameter_type == 'N/A':
        string = 'N/A'
    else:
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
   

# Function for writing lammps datafile's
def datafile(doc_title, atom_style, parameters, ff_class, version, include_type_labels, log):
        
    
    # Writing new file with new atoms
    with open(doc_title+'.data','w') as f: 
        # Write header
        header = '{} > all2lmp: {}  Class: {}'.format(parameters.header, version, ff_class)
        f.write(f'{header[-220:len(header)]}\n\n') # Make max header length of 220 characters
        
        # Write structure quantities
        f.write(f'{parameters.natoms} atoms\n')
        if parameters.nbonds > 0:
            f.write(f'{parameters.nbonds} bonds\n')
        if parameters.nangles > 0:
            f.write(f'{parameters.nangles} angles\n')
        if parameters.ndihedrals > 0:
            f.write(f'{parameters.ndihedrals} dihedrals\n')          
        if parameters.nimpropers > 0:
            f.write(f'{parameters.nimpropers} impropers\n')
        f.write('\n')
        
        # Write structure quantity types
        f.write(f'{parameters.natomtypes} atom types\n')
        if parameters.nbondtypes > 0:
            f.write(f'{parameters.nbondtypes} bond types\n')
        if parameters.nangletypes > 0:
            f.write(f'{parameters.nangletypes} angle types\n')
        if parameters.ndihedraltypes > 0:
            f.write(f'{parameters.ndihedraltypes} dihedral types\n')
        if parameters.nimpropertypes > 0:
            f.write(f'{parameters.nimpropertypes} improper types\n')
        f.write('\n')
        
        # write box size
        f.write('{:>12.9f} {:^9.9f} {} {}\n'.format(parameters.box.xlo, parameters.box.xhi, 'xlo', 'xhi'))
        f.write('{:>12.9f} {:^9.9f} {} {}\n'.format(parameters.box.ylo, parameters.box.yhi, 'ylo', 'yhi'))
        f.write('{:>12.9f} {:^9.9f} {} {}\n'.format(parameters.box.zlo, parameters.box.zhi, 'zlo', 'zhi'))
        if parameters.box.xy != 0 or parameters.box.xz != 0 or parameters.box.yz != 0:
            f.write('{:>12.9f} {:^9.9f} {:^9.9f} {} {} {}\n'.format(parameters.box.xy, parameters.box.xz, parameters.box.yz, 'xy', 'xz', 'yz'))
        
        
        
        # Write Atom Type Labels if user wants
        if include_type_labels:
            f.write('\nAtom Type Labels\n\n')
            for i in parameters.masses:
                mass = parameters.masses[i]
                f.write('{:^3} {:^2}\n'.format(i, mass.type))
 
        # Write Bond Type Labels if user wants
        if include_type_labels and parameters.nbonds > 0:
            f.write('\nBond Type Labels\n\n')
            for i in parameters.bond_coeffs:
                bond = parameters.bond_coeffs[i]
                f.write('{:^3} {:^2}\n'.format(i, string_type_labels(bond.type)))
                
        # Write Angle Type Labels if user wants
        if include_type_labels and parameters.nangles > 0:
            f.write('\nAngle Type Labels\n\n')
            for i in parameters.angle_coeffs: 
                angle = parameters.angle_coeffs[i]
                f.write('{:^3} {:^2}\n'.format(i, string_type_labels(angle.type)))    

        # Write Dihedral Type Labels if user wants
        if include_type_labels and parameters.ndihedraltypes > 0:
            f.write('\nDihedral Type Labels\n\n')
            for i in parameters.dihedral_coeffs: 
                dihedral = parameters.dihedral_coeffs[i]
                f.write('{:^3} {:^2}\n'.format(i, string_type_labels(dihedral.type))) 

        # Write Improper Type Labels if user wants
        if include_type_labels and parameters.nimpropertypes > 0:
            f.write('\nImproper Type Labels\n\n')
            for i in parameters.improper_coeffs: 
                improper = parameters.improper_coeffs[i]
                f.write('{:^3} {:^2}\n'.format(i, string_type_labels(list(improper.type) + [improper.nb]))) 
            
        
        # Write massses
        #f.write(f'\nMasses # {parameters.mass_comment}\n\n')
        f.write('\nMasses\n\n')
        for i in parameters.masses: 
            mass = parameters.masses[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
            #if include_type_labels: ID = '{t:<{s}}'.format(t=mass.type, s=5)            
            comment = '{:^2} {:5}'.format('#', mass.type)
            f.write('{} {:^12.8f} {:^2}\n'.format(ID, mass.coeffs, comment))

        # Write pair coeffs
        if ff_class in [0, 1, 2, 'd', 's1', 's2']:
            f.write(f'\nPair Coeffs  # {parameters.pair_comment}\n\n')
            for i in parameters.pair_coeffs: 
                pair = parameters.pair_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                #if include_type_labels: ID = '{t:<{s}}'.format(t=pair.type, s=5)  
                comment = '{:^2} {:5}'.format('#', pair.type)
                f.write('{} {} {:^2}\n'.format(ID, string_parameters(pair.coeffs), comment))
            
        # Write bond coeffs
        if parameters.nbondtypes > 0:
            f.write(f'\nBond Coeffs  # {parameters.bond_comment}\n\n')
            for i in parameters.bond_coeffs: 
                bond = parameters.bond_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(bond.type), s=10)
                comment = '{:^2} {:10}'.format('#', string_parameter_type(bond.type))
                f.write('{} {} {:^2}\n'.format(ID, string_parameters(bond.coeffs), comment))

        # Write angle coeffs
        if parameters.nangletypes > 0:
            f.write(f'\nAngle Coeffs  # {parameters.angle_comment}\n\n')
            for i in parameters.angle_coeffs: 
                angle = parameters.angle_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(angle.type), s=15)
                comment = '{:^2} {:15}'.format('#', string_parameter_type(angle.type))
                f.write('{} {} {:^2}\n'.format(ID, string_parameters(angle.coeffs), comment))
            
        # Write dihedral coeffs
        if parameters.ndihedraltypes > 0:
            f.write(f'\nDihedral Coeffs  # {parameters.dihedral_comment}\n\n')
            for i in parameters.dihedral_coeffs: 
                dihedral = parameters.dihedral_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(dihedral.type), s=20)
                comment = '{:^2} {:15}'.format('#', string_parameter_type(dihedral.type))
                f.write('{} {} {:^2}\n'.format(ID, string_parameters(dihedral.coeffs), comment))
                    
        # Write Improper coeffs
        if parameters.nimpropertypes > 0:
            f.write(f'\nImproper Coeffs  # {parameters.improper_comment}\n\n')
            for i in parameters.improper_coeffs: 
                improper = parameters.improper_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(list(improper.type) + [improper.nb]), s=25)
                comment = '{:^2} {:15}'.format('#', string_parameter_type(list(improper.type) + [improper.nb]))
                f.write('{} {} {:^2}\n'.format(ID, string_parameters(improper.coeffs), comment))


        # Write class2 parameters if class II FF or 's1'
        if ff_class in [2, 's2']:
            
            # Write bondbond coeffs
            if parameters.nangletypes > 0:
                f.write('\nBondBond Coeffs  # class2\n\n')
                for i in parameters.bondbond_coeffs: 
                    bondbond = parameters.bondbond_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(bondbond.type), s=15)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(bondbond.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(bondbond.coeffs), comment))
                
            # Write bondangle coeffs
            if parameters.nangletypes > 0:
                f.write('\nBondAngle Coeffs  # class2\n\n')
                for i in parameters.bondangle_coeffs: 
                    bondangle = parameters.bondangle_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(bondangle.type), s=15)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(bondangle.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(bondangle.coeffs), comment))
                
            # Write angleangletorsion coeffs
            if parameters.ndihedraltypes > 0:
                f.write('\nAngleAngleTorsion Coeffs  # class2\n\n')
                for i in parameters.angleangletorsion_coeffs: 
                    angleangletorsion = parameters.angleangletorsion_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(angleangletorsion.type), s=20)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(angleangletorsion.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(angleangletorsion.coeffs), comment))
    
            # Write endbondtorsion coeffs
            if parameters.ndihedraltypes > 0:
                f.write('\nEndBondTorsion Coeffs  # class2\n\n')
                for i in parameters.endbondtorsion_coeffs: 
                    endbondtorsion = parameters.endbondtorsion_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(endbondtorsion.type), s=20)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(endbondtorsion.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(endbondtorsion.coeffs), comment))
    
            # Write middlebondtorsion coeffs
            if parameters.ndihedraltypes > 0:
                f.write('\nMiddleBondTorsion Coeffs  # class2\n\n')
                for i in parameters.middlebondtorsion_coeffs: 
                    middlebondtorsion = parameters.middlebondtorsion_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(middlebondtorsion.type), s=20)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(middlebondtorsion.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(middlebondtorsion.coeffs), comment))
    
            # Write bondbond13 coeffs
            if parameters.ndihedraltypes > 0:
                f.write('\nBondBond13 Coeffs  # class2\n\n')
                for i in parameters.bondbond13_coeffs: 
                    bondbond13 = parameters.bondbond13_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(bondbond13.type), s=20)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(bondbond13.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(bondbond13.coeffs), comment))
    
            # Write angletorsion coeffs
            if parameters.ndihedraltypes > 0:
                f.write('\nAngleTorsion Coeffs  # class2\n\n')
                for i in parameters.angletorsion_coeffs: 
                    angletorsion = parameters.angletorsion_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(angletorsion.type), s=20)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(angletorsion.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(angletorsion.coeffs), comment))
                
            # Write angleangle coeffs
            if parameters.nimpropertypes > 0:
                f.write('\nAngleAngle Coeffs  # class2\n\n')
                for i in parameters.angleangle_coeffs: 
                    angleangle = parameters.angleangle_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(list(angleangle.type) + [angleangle.nb]), s=25)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(list(angleangle.type) + [angleangle.nb]))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(angleangle.coeffs), comment))

        # Write atoms    
        f.write(f'\nAtoms # {atom_style}\n\n')            
        for i in parameters.atoms:
            atom = parameters.atoms[i]; atomtype = '{t:<{s}}'.format(t=atom.type, s=3)
            #if include_type_labels: atomtype = '{t:<{s}}'.format(t=atom.name, s=5)
            comment = '{:^2} {}/{}'.format('#', atom.symbol, atom.element)
            
            # write charge atom style
            if atom_style == 'charge':
                f.write('{:^6} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atomtype, atom.charge, atom.x, atom.y,\
                                                                                                                atom.z, atom.ix, atom.iy, atom.iz, comment)) 
            # write molecular atom style
            elif atom_style == 'molecular':
                f.write('{:^6} {:^4} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.molid, atomtype, atom.x, atom.y,\
                                                                                                                atom.z, atom.ix, atom.iy, atom.iz, comment))                          
            # write full atom style
            elif atom_style == 'full':
                f.write('{:^6} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.molid, atomtype, atom.charge, atom.x, atom.y,\
                                                                                                                atom.z, atom.ix, atom.iy, atom.iz, comment)) 
            # write angle or bond atom style
            elif atom_style == 'angle' or atom_style == 'bond':
                f.write('{:^6} {:^4} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.molid, atomtype, atom.x, atom.y,\
                                                                                                         atom.z, atom.ix, atom.iy, atom.iz, comment))                 
            # write atomic atom style
            elif atom_style == 'atomic':
                f.write('{:^6} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atomtype, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment)) 
            
            # write dipole atom style
            elif atom_style == 'dipole':
                mux = 0; muy = 0; muz = 0;
                log.warn(f'WARNING writing LAMMPS datafile Atom section in dipole style and initializing mux = {mux}, muy = {muy}, muz = {muz}')
                f.write('{:^6} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.type, atom.charge, atom.x,
                                                                                                                                             atom.y, atom.z, mux, muy, muz, atom.ix,
                                                                                                                                             atom.iy, atom.iz, comment)) 
            # write charge atom style
            elif atom_style == 'dpd':
                theta = 0
                log.warn(f'WARNING writing LAMMPS datafile Atom section in dpd style and initializing theta = {theta}')
                f.write('{:^6} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atomtype, theta, atom.x, atom.y,\
                                                                                                                atom.z, atom.ix, atom.iy, atom.iz, comment)) 
                    
            # write angle or bond atom style
            elif atom_style == 'line':
                lineflag = 0; density = 0;
                log.warn(f'WARNING writing LAMMPS datafile Atom section in line style and initializing lineflag = {lineflag}, density = {density}')
                f.write('{:^6} {:^4} {:^2} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.molid, atomtype, lineflag, density, 
                                                                                                                           atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment))    
            # else raise and exception
            else: log.warn('Atom Style Error - write_lmp.py does not support atom style {atom_style}')
                
        # Write bonds
        if parameters.nbonds > 0:
            f.write("\nBonds\n\n")
            for i in parameters.bonds:
                bond = parameters.bonds[i]; bondtype = '{t:<{s}}'.format(t=str(bond.type), s=3)
                #if include_type_labels: bondtype = '{t:<{s}}'.format(t=string_type_labels(bond.symbol), s=10)
                id1, id2 = bond.atomids
                f.write('{:^5} {} {:^5} {:^5}\n'.format(i, bondtype, id1, id2))
            
        # Write angles
        if parameters.nangles > 0:
            f.write("\nAngles\n\n")
            for i in parameters.angles:
                angle = parameters.angles[i]; angletype = '{t:<{s}}'.format(t=str(angle.type), s=3)
                #if include_type_labels: angletype = '{t:<{s}}'.format(t=string_type_labels(angle.symbol), s=15)
                id1, id2, id3 = angle.atomids
                f.write('{:^5} {} {:^5} {:^5} {:^5}\n'.format(i, angletype, id1, id2, id3))
            
        # Write dihedrals
        if parameters.ndihedrals > 0:
            f.write("\nDihedrals\n\n")
            for i in parameters.dihedrals:
                dihedral = parameters.dihedrals[i]; dihedraltype = '{t:<{s}}'.format(t=str(dihedral.type), s=3)
                #if include_type_labels: dihedraltype = '{t:<{s}}'.format(t=string_type_labels(dihedral.symbol), s=20)
                id1, id2, id3, id4 = dihedral.atomids
                f.write('{:^5} {} {:^5} {:^5} {:^5} {:^5}\n'.format(i, dihedraltype, id1, id2, id3, id4))
                    
        # Write impropers
        if parameters.nimpropers > 0:
            f.write("\nImpropers\n\n")
            for i in parameters.impropers:
                improper = parameters.impropers[i]; impropertype = '{t:<{s}}'.format(t=str(improper.type), s=3)
                #if include_type_labels: impropertype = '{t:<{s}}'.format(t=string_type_labels(list(improper.symbol) + [improper.nb]), s=25)
                id1, id2, id3, id4 = improper.atomids
                f.write('{:^5} {} {:^5} {:^5} {:^5} {:^5}\n'.format(i, impropertype, id1, id2, id3, id4))
            

            
            
            
def comments(doc_title, atom_style, parameters, ff_class, version, log):
    
    def bondics2string(bond_incs):
        # bond_incs is dict formated as:
        # { (atom-ids) : [ (type1, type2), bond-inc)]   (          ): [                    ], ....}
        # {(2755, 2734): [('c=',   'cp'),  0.0],        (2755, 2776): [('c=', 'h'), -0.1268], ....}
        
        # intialize string
        string = ''
        for i in bond_incs:
            try:
                types, bond_inc = bond_incs[i]
                string += '{:<5}'.format('[')
                string += '{:^10} {:^6} - {:^6}'.format('atom-ids: ', i[0], i[1])
                string += '{:^10} {:^5} - {:^5}'.format('equivs: ', types[0], types[1])
                string += '{:^15} {:8.5f}'.format('bond-inc: ', bond_inc)
                string += '{:>5}'.format('],   ')
            except: 
                string += '{:<5}'.format('[')
                string += '{:^10} {:^6} - {:^6}'.format('atom-ids: ', i[0], i[1])
                string += '{:^10} {:^5} - {:^5}'.format('equivs: ', 'N/A', 'N/A')
                string += '{:^15} {:^8.5f}'.format('bond-inc: ', 0)
                string += '{:>5}'.format('],   ')
        return string
    
    # Writing comments file
    with open(doc_title +'.txt','w') as f: 
        # Write header
        f.write(f'{parameters.header} > all2lmp: {version}  Class: {ff_class}\n\n') 
        
        # Write percetage breakdown
        f.write('# -------------------------------------\n')
        f.write('# Parameterization percentage breakdown\n')
        f.write('# -------------------------------------\n')
        for i in parameters.percents:
            f.write('# {:<28} {:^6.2f} %\n'.format(i+str(':'), parameters.percents[i]))
    
        # Write massses
        f.write(f'\nMasses # {parameters.mass_comment}\n\n')
        for i in parameters.masses: 
            mass = parameters.masses[i]
            comment = '{:^2} {:5} {:^10} {:5} {:25}'.format('#', mass.type, 'equivalent:', mass.equivalent, mass.comments)
            f.write('{:^3} {:^12.8f} {:^2}\n'.format(i, mass.coeffs, comment))
            
        # Write pair coeffs
        if ff_class in [0, 1, 2, 'd', 's1', 's2']:
            f.write(f'\nPair Coeffs  # {parameters.pair_comment}\n\n')
            for i in parameters.pair_coeffs: 
                pair = parameters.pair_coeffs[i] 
                comment = '{:^2} {:5} {:^15} {:5} {:25}'.format('#', pair.type, 'equivalent:', pair.equivalent, pair.comments)
                f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(pair.coeffs), comment))
            
        # Write bond coeffs
        if parameters.nbonds > 0:
            f.write(f'\nBond Coeffs  # {parameters.bond_comment}\n\n')
            for i in parameters.bond_coeffs: 
                bond = parameters.bond_coeffs[i]
                coeff_type = string_parameter_type(bond.type)
                equiv_type = string_parameter_type(bond.equivalent)
                try: match_type = string_parameter_type(bond.match)
                except:  match_type = 'N/A'
                comment = '{:^2} {:15} {:^10} {:^25} {:^10} {:^25} {:15}'.format('#', coeff_type, 'equivalent:', equiv_type, 'match:', match_type, bond.comments)
                f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(bond.coeffs), comment))
            
        # Write angle coeffs
        if parameters.nangles > 0:
            f.write(f'\nAngle Coeffs  # {parameters.angle_comment}\n\n')
            for i in parameters.angle_coeffs: 
                angle = parameters.angle_coeffs[i]
                coeff_type = string_parameter_type(angle.type)
                equiv_type = string_parameter_type(angle.equivalent)
                try: match_type = string_parameter_type(angle.match)
                except:  match_type = 'N/A'
                comment = '{:^2} {:15} {:^10} {:^30} {:^10} {:^30} {:15}'.format('#', coeff_type, 'equivalent:', equiv_type, 'match:', match_type, angle.comments)
                f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(angle.coeffs), comment))
                
        # Write dihedral coeffs
        if parameters.ndihedraltypes > 0:
            f.write(f'\nDihedral Coeffs  # {parameters.dihedral_comment}\n\n')
            for i in parameters.dihedral_coeffs: 
                dihedral = parameters.dihedral_coeffs[i]
                coeff_type = string_parameter_type(dihedral.type)
                equiv_type = string_parameter_type(dihedral.equivalent)
                try: match_type = string_parameter_type(dihedral.match)
                except:  match_type = 'N/A'
                comment = '{:^2} {:20} {:^15} {:^35} {:^15} {:^35} {:20}'.format('#', coeff_type, 'equivalent:', equiv_type, 'match:', match_type, dihedral.comments)
                f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(dihedral.coeffs), comment))
            
        # Write Improper coeffs
        if parameters.nimpropertypes > 0:
            f.write(f'\nImproper Coeffs  # {parameters.improper_comment}\n\n')
            for i in parameters.improper_coeffs: 
                improper = parameters.improper_coeffs[i]
                coeff_type = string_parameter_type(list(improper.type) + [improper.nb])
                equiv_type = string_parameter_type(improper.equivalent)
                try: match_type = string_parameter_type(improper.match)
                except:  match_type = 'N/A'
                comment = '{:^2} {:25} {:^15} {:^35} {:^15} {:^35} {:20}'.format('#', coeff_type, 'equivalent:', equiv_type, 'match:', match_type, improper.comments)
                f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(improper.coeffs), comment))

            
        # Write class2 parameters if class II FF or 's1'
        if ff_class in [2, 's2']:
            
            # Write bondbond coeffs
            if parameters.nangles > 0:
                f.write('\nBondBond Coeffs  # class2\n\n')
                for i in parameters.bondbond_coeffs: 
                    bondbond = parameters.bondbond_coeffs[i]
                    coeff_type = string_parameter_type(bondbond.type)
                    equiv_type = string_parameter_type(bondbond.equivalent)
                    comment = '{:^2} {:15} {:^10} {:^30} {:15}'.format('#', coeff_type, 'equivalent:', equiv_type, bondbond.comments)
                    f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(bondbond.coeffs), comment))
                    
            # Write bondangle coeffs
            if parameters.nangles > 0:
                f.write('\nBondAngle Coeffs  # class2\n\n')
                for i in parameters.bondangle_coeffs: 
                    bondangle = parameters.bondangle_coeffs[i]
                    coeff_type = string_parameter_type(bondangle.type)
                    equiv_type = string_parameter_type(bondangle.equivalent)
                    comment = '{:^2} {:25} {:^15} {:^30} {:15}'.format('#', coeff_type, 'equivalent:', equiv_type, bondangle.comments)
                    f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(bondangle.coeffs), comment))
                
            # Write angleangletorsion coeffs
            f.write('\nAngleAngleTorsion Coeffs  # class2\n\n')
            for i in parameters.angleangletorsion_coeffs: 
                angleangletorsion = parameters.angleangletorsion_coeffs[i]
                coeff_type = string_parameter_type(angleangletorsion.type)
                equiv_type = string_parameter_type(angleangletorsion.equivalent)
                comment = '{:^2} {:25} {:^15} {:25} {:15}'.format('#', coeff_type, 'equivalent:', equiv_type, angleangletorsion.comments)
                f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(angleangletorsion.coeffs), comment))
                
            # Write endbondtorsion coeffs
            f.write('\nEndBondTorsion Coeffs  # class2\n\n')
            for i in parameters.endbondtorsion_coeffs: 
                endbondtorsion = parameters.endbondtorsion_coeffs[i]
                coeff_type = string_parameter_type(endbondtorsion.type)
                equiv_type = string_parameter_type(endbondtorsion.equivalent)
                comment = '{:^2} {:25} {:^15} {:25} {:15}'.format('#', coeff_type, 'equivalent:', equiv_type, endbondtorsion.comments)
                f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(endbondtorsion.coeffs), comment))
    
            # Write middlebondtorsion coeffs
            f.write('\nMiddleBondTorsion Coeffs  # class2\n\n')
            for i in parameters.middlebondtorsion_coeffs: 
                middlebondtorsion = parameters.middlebondtorsion_coeffs[i]
                coeff_type = string_parameter_type(middlebondtorsion.type)
                equiv_type = string_parameter_type(middlebondtorsion.equivalent)
                comment = '{:^2} {:25} {:^15} {:25} {:15}'.format('#', coeff_type, 'equivalent:', equiv_type, middlebondtorsion.comments)
                f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(middlebondtorsion.coeffs), comment))
                
            # Write bondbond13 coeffs
            f.write('\nBondBond13 Coeffs  # class2\n\n')
            for i in parameters.bondbond13_coeffs: 
                bondbond13 = parameters.bondbond13_coeffs[i]
                coeff_type = string_parameter_type(bondbond13.type)
                equiv_type = string_parameter_type(bondbond13.equivalent)
                comment = '{:^2} {:25} {:^15} {:25} {:15}'.format('#', coeff_type, 'equivalent:', equiv_type, bondbond13.comments)
                f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(bondbond13.coeffs), comment))
    
            # Write angletorsion coeffs
            f.write('\nAngleTorsion Coeffs  # class2\n\n')
            for i in parameters.angletorsion_coeffs: 
                angletorsion = parameters.angletorsion_coeffs[i]
                coeff_type = string_parameter_type(angletorsion.type)
                equiv_type = string_parameter_type(angletorsion.equivalent)
                comment = '{:^2} {:25} {:^15} {:25} {:15}'.format('#', coeff_type, 'equivalent:', equiv_type, angletorsion.comments)
                f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(angletorsion.coeffs), comment))  
                
            # Write angleangle coeffs
            f.write('\nAngleAngle Coeffs  # class2\n\n')
            for i in parameters.angleangle_coeffs: 
                angleangle = parameters.angleangle_coeffs[i]
                coeff_type = string_parameter_type(list(angleangle.type) + [angleangle.nb])
                equiv_type = string_parameter_type(angleangle.equivalent)
                comment = '{:^2} {:25} {:^15} {:25} {:15}'.format('#', coeff_type, 'equivalent:', equiv_type, angleangle.comments)
                f.write('{:^3} {} {:^2}\n'.format(i, string_parameters(angleangle.coeffs), comment))

        # Write atoms    
        f.write(f'\nAtoms # {atom_style}\n\n')            
        for i in parameters.atoms:
            atom = parameters.atoms[i]
            
            # Find formatted bond-incs string
            bond_incs = bondics2string(atom.bond_incs)
            
            # Generate comment
            comment = '{:^2} {:^6} {:^72} {}'.format('#', atom.symbol, atom.comments, bond_incs)
        
            
            # write charge atom style
            if atom_style == 'charge':
                f.write('{:^3} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {}\n'.format(i, atom.type, atom.charge, atom.x, atom.y,\
                                                                                                            atom.z, atom.ix, atom.iy, atom.iz, comment)) 
            # write molecular atom style
            elif atom_style == 'molecular':
                f.write('{:^3} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {}\n'.format(i, atom.type, atom.charge, atom.x, atom.y,\
                                                                                                            atom.z, atom.ix, atom.iy, atom.iz, comment))                          
            # write full atom style
            elif atom_style == 'full':
                f.write('{:^3} {:^2} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {}\n'.format(i, atom.molid, atom.type, atom.charge, atom.x, atom.y,\
                                                                                                                  atom.z, atom.ix, atom.iy, atom.iz, comment)) 
                    
            # write angle or bond atom style
            elif atom_style == 'angle' or atom_style == 'bond':
                f.write('{:^6} {:^4} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.molid, atom.type, atom.x, atom.y,\
                                                                                                         atom.z, atom.ix, atom.iy, atom.iz, comment))                 
            # write atomic atom style
            elif atom_style == 'atomic':
                f.write('{:^6} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.type, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment)) 
            
            # write dipole atom style
            elif atom_style == 'dipole':
                mux = 0; muy = 0; muz = 0;
                log.warn(f'WARNING writing LAMMPS datafile Atom section in dipole style and initializing mux = {mux}, muy = {muy}, muz = {muz}')
                f.write('{:^6} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.type, atom.charge, atom.x,
                                                                                                                                             atom.y, atom.z, mux, muy, muz, atom.ix,
                                                                                                                                             atom.iy, atom.iz, comment)) 
            # write charge atom style
            elif atom_style == 'dpd':
                theta = 0
                log.warn(f'WARNING writing LAMMPS datafile Atom section in dpd style and initializing theta = {theta}')
                f.write('{:^6} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.type, theta, atom.x, atom.y,\
                                                                                                                atom.z, atom.ix, atom.iy, atom.iz, comment)) 
                    
            # write angle or bond atom style
            elif atom_style == 'line':
                lineflag = 0; density = 0;
                log.warn(f'WARNING writing LAMMPS datafile Atom section in line style and initializing lineflag = {lineflag}, density = {density}')
                f.write('{:^6} {:^4} {:^2} {:^2} {:^15.9f} {:^15.9f} {:^15.9f} {:^15.9f} {:>4} {:>4} {:>4} {:^5}\n'.format(i, atom.molid, atom.type, lineflag, density, 
                                                                                                                           atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz, comment))    
                
            # else raise and exception
            else: log.error('Atom Style Error - write_lmp.py does not support atom style {atom_style}')

        # Write bonds
        if parameters.nbonds > 0:
            f.write("\nBonds\n\n")
            for i in parameters.bonds:
                bond = parameters.bonds[i]
                id1, id2 = bond.atomids
                comment = '{:^6} {:10}'.format('#', string_parameter_type(bond.symbol))
                f.write('{:^5} {:^5} {:^5} {:^5} {:^6}\n'.format(i, bond.type, id1, id2, comment))
            
        # Write angles
        if parameters.nangles > 0:
            f.write("\nAngles\n\n")
            for i in parameters.angles:
                angle = parameters.angles[i]
                id1, id2, id3 = angle.atomids
                comment = '{:^6} {:10}'.format('#', string_parameter_type(angle.symbol))
                f.write('{:^5} {:^5} {:^5} {:^5} {:^5} {:^6}\n'.format(i, angle.type, id1, id2, id3, comment))
            
        # Write dihedrals
        if parameters.ndihedrals > 0:
            f.write("\nDihedrals\n\n")
            for i in parameters.dihedrals:
                dihedral = parameters.dihedrals[i]
                id1, id2, id3, id4 = dihedral.atomids
                comment = '{:^6} {:10}'.format('#', string_parameter_type(dihedral.symbol))
                f.write('{:^5} {:^5} {:^5} {:^5} {:^5} {:^5} {:^6}\n'.format(i, dihedral.type, id1, id2, id3, id4, comment))
                    
        # Write impropers
        if parameters.nimpropers > 0:
            f.write("\nImpropers\n\n")
            for i in parameters.impropers:
                improper = parameters.impropers[i]
                id1, id2, id3, id4 = improper.atomids
                comment = '{:^6} {:10}'.format('#', string_parameter_type(list(improper.symbol) + [improper.nb]))
                f.write('{:^5} {:^5} {:^5} {:^5} {:^5} {:^5} {:^6}\n'.format(i, improper.type, id1, id2, id3, id4, comment))            
        
            
