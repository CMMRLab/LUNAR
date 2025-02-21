#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
November 13th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function for stringing together parameter types
def string_parameter_type(parameter_type):
    string = ''; str_buffer = 6 # Set string buffer size
    for n, i in enumerate(parameter_type):
        if n == 0: string += '{:<{str_buffer}}'.format(i, str_buffer=str_buffer)
        elif n < len(parameter_type)-1: string += '{:^{str_buffer}}'.format(i, str_buffer=str_buffer+2)
        else: string += '{:>{str_buffer}} {:^2}'.format(i, '', str_buffer=str_buffer)
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
        string += '{:^16.8f}'.format(i)
    return string


# Function to write LAMMPS molecule file
def write_moltemp(doc_title, parameters, ff_class, include_type_labels, version):
    
    # Optional molecule sections not need for bond/react (hidden from user's since
    # its assumed they wont need these, but can be activated easily when needed)
    write_molids = False
    write_shake_flags = False
    
    # Writing new file with new atoms
    with open(doc_title+'.lmpmol','w') as f: 
        # Write header
        f.write(f'{parameters.header} > all2lmp: {version}  Class: {ff_class}\n\n') 
        
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
        
        # Write types    
        f.write('Types\n\n')            
        for i in parameters.atoms:
            atom = parameters.atoms[i]; atomtype = '{t:<{s}}'.format(t=atom.type, s=3)
            if include_type_labels: atomtype = '{t:<{s}}'.format(t=atom.symbol, s=5)
            f.write('{:^3} {} {:^2} {:^2}\n'.format(i, atomtype, '#', atom.symbol)) 
        
        # Write charges    
        f.write('\nCharges\n\n')            
        for i in parameters.atoms:
            atom = parameters.atoms[i]
            f.write('{:^3} {:^15.6f} {:^2} {:^2}\n'.format(i, atom.charge, '#', atom.symbol)) 
            
        # Write coordinates   
        f.write('\nCoords\n\n')            
        for i in parameters.atoms:
            atom = parameters.atoms[i]
            f.write('{:^3} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2}\n'.format(i, atom.x, atom.y, atom.z, '#', atom.symbol)) 
            
        # Write molecule ids if desired
        if write_molids:
            f.write('\nMolecules\n\n')            
            for i in parameters.atoms:
                atom = parameters.atoms[i]
                f.write('{:^3} {:^2} {:^2} {:^2}\n'.format(i, atom.molid, '#', atom.symbol)) 
                
        # Write Shake Flags if desired (intialize as zero)
        if write_shake_flags:
            f.write('\nShake Flags\n\n')            
            for i in parameters.atoms:
                atom = parameters.atoms[i]
                f.write('{:^3} {:^2} {:^2} {:^2}\n'.format(i, 0, '#', atom.symbol)) 
            
        # Write bonds
        if parameters.nbonds > 0:
            f.write("\nBonds\n\n")
            for i in parameters.bonds:
                bond = parameters.bonds[i]; bondtype = '{t:<{s}}'.format(t=str(bond.type), s=3)
                if include_type_labels: bondtype = '{t:<{s}}'.format(t=string_type_labels(bond.symbol), s=10)
                id1, id2 = bond.atomids
                comment = '{:^6} {:4}'.format('#', string_parameter_type(parameters.bond_coeffs[bond.type].type))
                f.write('{:^2} {} {:^5} {:^5} {:^6}\n'.format(i, bondtype, id1, id2, comment))
            
        # Write angles
        if parameters.nangles > 0:
            f.write("\nAngles\n\n")
            for i in parameters.angles:
                angle = parameters.angles[i]; angletype = '{t:<{s}}'.format(t=str(angle.type), s=3)
                if include_type_labels: angletype = '{t:<{s}}'.format(t=string_type_labels(angle.symbol), s=15)
                id1, id2, id3 = angle.atomids
                comment = '{:^6} {:4}'.format('#', string_parameter_type(parameters.angle_coeffs[angle.type].type))
                f.write('{:^2} {} {:^5} {:^5} {:^5} {:^6}\n'.format(i, angletype, id1, id2, id3, comment))
            
        # Write dihedrals
        if parameters.ndihedrals > 0:
            f.write("\nDihedrals\n\n")
            for i in parameters.dihedrals:
                dihedral = parameters.dihedrals[i]; dihedraltype = '{t:<{s}}'.format(t=str(dihedral.type), s=3)
                if include_type_labels: dihedraltype = '{t:<{s}}'.format(t=string_type_labels(dihedral.symbol), s=20)
                id1, id2, id3, id4 = dihedral.atomids
                comment = '{:^6} {:4}'.format('#', string_parameter_type(parameters.dihedral_coeffs[dihedral.type].type))
                f.write('{:^2} {} {:^5} {:^5} {:^5} {:^5} {:^6}\n'.format(i, dihedraltype, id1, id2, id3, id4, comment))
                    
        # Write impropers
        if parameters.nimpropers > 0:
            f.write("\nImpropers\n\n")
            for i in parameters.impropers:
                improper = parameters.impropers[i]; impropertype = '{t:<{s}}'.format(t=str(improper.type), s=3)
                if include_type_labels: impropertype = '{t:<{s}}'.format(t=string_type_labels(list(improper.symbol) + [improper.nb]), s=25)
                id1, id2, id3, id4 = improper.atomids
                comment = '{:^6} {:4}'.format('#', string_parameter_type(parameters.improper_coeffs[improper.type].type))
                f.write('{:^2} {} {:^5} {:^5} {:^5} {:^5} {:^6}\n'.format(i, impropertype, id1, id2, id3, id4, comment))         
    return

# Function to write a LAMMPS data file, but only writes energy coeffs               
def write_ecoeffs(doc_title, parameters, ff_class, version, include_type_labels):
    
    # Writing new file with new atoms
    with open(doc_title+'.ecoeffs','w') as f:         
        # Write header
        header = '{} > all2lmp: {}  Class: {}'.format(parameters.header, version, ff_class)
        f.write(f'{header[-220:len(header)]}\n\n') # Make max header length of 220 characters
        
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

        
        # Write Atom Type Labels if user wants
        if include_type_labels:
            f.write('\nAtom Type Labels\n\n')
            for i in parameters.masses:
                mass = parameters.masses[i]
                f.write('{:^3} {:^2}\n'.format(i, mass.type))
        
        # Write massses
        f.write(f'\nMasses # {parameters.mass_comment}\n\n')
        for i in parameters.masses: 
            mass = parameters.masses[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
            #if include_type_labels: ID = '{t:<{s}}'.format(t=mass.type, s=5)            
            comment = '{:^2} {:5}'.format('#', mass.type)
            f.write('{} {:^12.8f} {:^2}\n'.format(ID, mass.coeffs, comment))
        
        # Write pair coeffs
        if ff_class in [0, 1, 2, 'd', 's1', 's2', '0', '1', '2']:
            f.write(f'\nPair Coeffs  # {parameters.pair_comment}\n\n')
            for i in parameters.pair_coeffs: 
                pair = parameters.pair_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                #if include_type_labels: ID = '{t:<{s}}'.format(t=pair.type, s=5)  
                comment = '{:^2} {:5}'.format('#', pair.type)
                f.write('{} {} {:^2}\n'.format(ID, string_parameters(pair.coeffs), comment))
                
                
                
        # Write Bond Type Labels if user wants
        if include_type_labels and parameters.nbonds > 0:
            f.write('\nBond Type Labels\n\n')
            for i in parameters.bond_coeffs:
                bond = parameters.bond_coeffs[i]
                f.write('{:^3} {:^2}\n'.format(i, string_type_labels(bond.type)))
            
        # Write bond coeffs
        if parameters.nbonds > 0:
            f.write(f'\nBond Coeffs  # {parameters.bond_comment}\n\n')
            for i in parameters.bond_coeffs: 
                bond = parameters.bond_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(bond.type), s=10)
                comment = '{:^2} {:10}'.format('#', string_parameter_type(bond.type))
                f.write('{} {} {:^2}\n'.format(ID, string_parameters(bond.coeffs), comment))
            
                        
            
        # Write Angle Type Labels if user wants
        if include_type_labels and parameters.nangles > 0:
            f.write('\nAngle Type Labels\n\n')
            for i in parameters.angle_coeffs: 
                angle = parameters.angle_coeffs[i]
                f.write('{:^3} {:^2}\n'.format(i, string_type_labels(angle.type)))    
            
        # Write angle coeffs
        if parameters.nangles > 0:
            f.write(f'\nAngle Coeffs  # {parameters.angle_comment}\n\n')
            for i in parameters.angle_coeffs: 
                angle = parameters.angle_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(angle.type), s=15)
                comment = '{:^2} {:15}'.format('#', string_parameter_type(angle.type))
                f.write('{} {} {:^2}\n'.format(ID, string_parameters(angle.coeffs), comment))
            
            
            
        # Write Dihedral Type Labels if user wants
        if include_type_labels and parameters.ndihedraltypes > 0:
            f.write('\nDihedral Type Labels\n\n')
            for i in parameters.dihedral_coeffs: 
                dihedral = parameters.dihedral_coeffs[i]
                f.write('{:^3} {:^2}\n'.format(i, string_type_labels(dihedral.type))) 
            
        # Write dihedral coeffs
        if parameters.ndihedraltypes > 0:
            f.write(f'\nDihedral Coeffs  # {parameters.dihedral_comment}\n\n')
            for i in parameters.dihedral_coeffs: 
                dihedral = parameters.dihedral_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(dihedral.type), s=20)
                comment = '{:^2} {:15}'.format('#', string_parameter_type(dihedral.type))
                f.write('{} {} {:^2}\n'.format(ID, string_parameters(dihedral.coeffs), comment))
           
            
           
        # Write Improper Type Labels if user wants
        if include_type_labels and parameters.nimpropertypes > 0:
            f.write('\nImproper Type Labels\n\n')
            for i in parameters.improper_coeffs: 
                improper = parameters.improper_coeffs[i]
                f.write('{:^3} {:^2}\n'.format(i, string_type_labels(list(improper.type) + [improper.nb]))) 
                    
        # Write Improper coeffs
        if parameters.nimpropertypes > 0:
            f.write(f'\nImproper Coeffs  # {parameters.improper_comment}\n\n')
            for i in parameters.improper_coeffs: 
                improper = parameters.improper_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(list(improper.type) + [improper.nb]), s=25)
                comment = '{:^2} {:15}'.format('#', string_parameter_type(list(improper.type) + [improper.nb]))
                f.write('{} {} {:^2}\n'.format(ID, string_parameters(improper.coeffs), comment))
        
        
        # Write class2 parameters if class II FF or 's1'
        if ff_class in [2, '2', 's2']:
            
            # Write bondbond coeffs
            if parameters.nangles > 0:
                f.write('\nBondBond Coeffs  # class2\n\n')
                for i in parameters.bondbond_coeffs: 
                    bondbond = parameters.bondbond_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(bondbond.type), s=15)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(bondbond.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(bondbond.coeffs), comment))
                
            # Write bondangle coeffs
            if parameters.nangles > 0:
                f.write('\nBondAngle Coeffs  # class2\n\n')
                for i in parameters.bondangle_coeffs: 
                    bondangle = parameters.bondangle_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(bondangle.type), s=15)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(bondangle.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(bondangle.coeffs), comment))
                
            # Write angleangletorsion coeffs
            f.write('\nAngleAngleTorsion Coeffs  # class2\n\n')
            if parameters.ndihedraltypes > 0:
                for i in parameters.angleangletorsion_coeffs: 
                    angleangletorsion = parameters.angleangletorsion_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(angleangletorsion.type), s=20)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(angleangletorsion.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(angleangletorsion.coeffs), comment))
        
            # Write endbondtorsion coeffs
            f.write('\nEndBondTorsion Coeffs  # class2\n\n')
            if parameters.ndihedraltypes > 0:
                for i in parameters.endbondtorsion_coeffs: 
                    endbondtorsion = parameters.endbondtorsion_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(endbondtorsion.type), s=20)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(endbondtorsion.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(endbondtorsion.coeffs), comment))
        
            # Write middlebondtorsion coeffs
            f.write('\nMiddleBondTorsion Coeffs  # class2\n\n')
            if parameters.ndihedraltypes > 0:
                for i in parameters.middlebondtorsion_coeffs: 
                    middlebondtorsion = parameters.middlebondtorsion_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(middlebondtorsion.type), s=20)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(middlebondtorsion.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(middlebondtorsion.coeffs), comment))
        
            # Write bondbond13 coeffs
            f.write('\nBondBond13 Coeffs  # class2\n\n')
            if parameters.ndihedraltypes > 0:
                for i in parameters.bondbond13_coeffs: 
                    bondbond13 = parameters.bondbond13_coeffs[i]; ID = '{t:<{s}}'.format(t=str(i), s=3)
                    #if include_type_labels: ID = '{t:<{s}}'.format(t=string_type_labels(bondbond13.type), s=20)
                    comment = '{:^2} {:15}'.format('#', string_parameter_type(bondbond13.type))
                    f.write('{} {} {:^2}\n'.format(ID, string_parameters(bondbond13.coeffs), comment))
        
            # Write angletorsion coeffs
            f.write('\nAngleTorsion Coeffs  # class2\n\n')
            if parameters.ndihedraltypes > 0:
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
    return