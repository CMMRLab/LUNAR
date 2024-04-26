# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
December, 2nd 2022
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

# Function to split coeff types into tuple
def split_coeff(types):
    # Find coeff types and split
    types = types.strip()
    types = types.split()
    types = tuple(types)
    return types

    

# Function for writing lammps molecule file
def file(newname, m, new, file, version, include_type_labels):
    # Writing new file with new atoms
    with open(newname+'.lmpmol','w') as f: 
        # Write header
        header = '{} > bond_react_merge: {} molecule file (filetag: {})'.format(m.header, version, file)
        f.write(f'{header[-220:len(header)]}\n\n') # Make max header length of 220 characters 
        
        # Write structure quantities
        f.write(f'{m.natoms} atoms\n')
        f.write(f'{m.nbonds} bonds\n')
        if m.nangles > 0:
            f.write(f'{m.nangles} angles\n')
        if m.ndihedrals > 0:
            f.write(f'{m.ndihedrals} dihedrals\n')          
        if m.nimpropers > 0:
            f.write(f'{m.nimpropers} impropers\n')  
        if hasattr(m,'fragments'):
            f.write(f'{len(m.fragments)} fragments\n')
        
        # Write fragments if they exist
        if hasattr(m,'fragments'):
            f.write('\nFragments\n\n')
            for i in m.fragments:
                frag = m.fragments[i]
                fragments = [str(i) for i in frag.atoms] # make all a string
                fragments = ' '.join(fragments) # Join string atomIDs together
                comment = '{:^2} option: {:5}'.format('#', frag.comment)
                f.write('{:^3}   {} {:^2}\n'.format(i, fragments, comment)) 
                
        # # Write molecules if they exist
        if hasattr(m,'molecules'):
            f.write('\nMolecules\n\n')
            for i in m.molecules.atoms:
                atom = m.molecules.atoms[i]
                comment = '{:^2} {:5} {:5}'.format('#', atom.type, atom.comment)
                f.write('{:^3}   {} {:^2}\n'.format(i, atom.molid, comment))
        
        # # Write Types   
        f.write('\nTypes\n\n')            
        for i in m.atoms:
            atom = m.atoms[i]
            coeff = m.pair_coeffs[atom.type]
            str_type = coeff.type
            
            # Find new type via mapping from atom_types_map
            new_type = atom.type
            if include_type_labels: new_type = '{t:<{s}}'.format(t=str_type, s=5)
            
            comment = '{:^2} {:5}'.format('#', str_type)
            f.write('{:^3} {:^2} {:^2}\n'.format(i, new_type, comment)) 
            
        # # Write Charges  
        f.write('\nCharges\n\n')            
        for i in m.atoms:   
            atom = m.atoms[i]
            coeff = m.pair_coeffs[atom.type]
            str_type = coeff.type
            new_type = atom.type
            
            # Find if atom has hasattr mapped_charge and mapped_comment if see user wanted to
            # map charges near and around edge atoms from their respective datafiles since edge
            # atoms may not have correct charge on them depending on pre-processing or charge method
            if hasattr(atom,'mapped_charge'):
                charge = atom.mapped_charge
                comment = '{:^2} {:5} {}'.format('#', str_type, atom.mapped_comment)
                
            else:
                charge = atom.charge
                comment = '{:^2} {:5}'.format('#', str_type)
            f.write('{:^3} {:>10.6f} {:^2}\n'.format(i, charge, comment)) 
            
        # # Write Coordinates 
        f.write('\nCoords\n\n')            
        for i in m.atoms:
            atom = m.atoms[i]
            coeff = m.pair_coeffs[atom.type]
            str_type = coeff.type
            new_type = atom.type
            
            comment = '{:^2} {:5}'.format('#', str_type)
            f.write('{:^3} {:^15.6f} {:^15.6f} {:^15.6f} {:^2}\n'.format(i, atom.x, atom.y, atom.z, comment)) 

            
        # Write bonds
        if m.nbonds > 0:
            f.write('\nBonds\n\n')
            for i in m.bonds:
                bond = m.bonds[i]
                id1, id2 = bond.atomids
                coeff = m.bond_coeffs[bond.type]
                str_type = coeff.type
                str_type_tuple = split_coeff(str_type)
                
                # Find new type via mapping from atom_types_map
                new_type = bond.type
                if include_type_labels: new_type = '{t:<{s}}'.format(t=string_type_labels(str_type_tuple), s=10)
                
                # build comment and write results
                comment = '{:^2} {:5}'.format('#', str_type)
                f.write('{:^2} {:^2} {:^2} {:^2} {:^6}\n'.format(i, new_type, id1, id2, comment))
            
        # Write angles
        if m.nangles > 0:
            f.write('\nAngles\n\n')
            for i in m.angles:
                angle = m.angles[i]
                id1, id2, id3 = angle.atomids
                coeff = m.angle_coeffs[angle.type]
                str_type = coeff.type
                str_type_tuple = split_coeff(str_type)
                
                # Find new type via mapping from atom_types_map
                new_type = angle.type
                if include_type_labels: new_type = '{t:<{s}}'.format(t=string_type_labels(str_type_tuple), s=15)
                
                # build comment and write results
                comment = '{:^2} {:5}'.format('#', str_type)
                f.write('{:^2} {:^2} {:^2} {:^2} {:^2} {:^6}\n'.format(i, new_type, id1, id2, id3, comment))
            
        # Write dihedrals
        if m.ndihedrals > 0:
            f.write('\nDihedrals\n\n')
            for i in m.dihedrals:
                dihedral = m.dihedrals[i]
                id1, id2, id3, id4 = dihedral.atomids
                coeff = m.dihedral_coeffs[dihedral.type]
                str_type = coeff.type
                str_type_tuple = split_coeff(str_type)
                
                # Find new type via mapping from atom_types_map
                new_type = dihedral.type
                if include_type_labels: new_type = '{t:<{s}}'.format(t=string_type_labels(str_type_tuple), s=20)
                
                # build comment and write results
                comment = '{:^2} {:5}'.format('#', str_type)
                f.write('{:^2} {:^2} {:^2} {:^2} {:^2} {:^2} {:^6}\n'.format(i, new_type, id1, id2, id3, id4, comment))
                    
        # Write impropers
        if m.nimpropers > 0:
            f.write('\nImpropers\n\n')
            for i in m.impropers:
                improper = m.impropers[i]
                id1, id2, id3, id4 = improper.atomids
                coeff = m.improper_coeffs[improper.type]
                str_type = coeff.type
                str_type_tuple = split_coeff(str_type)
                
                # Find new type via mapping from atom_types_map
                new_type = improper.type
                if include_type_labels: new_type = '{t:<{s}}'.format(t=string_type_labels(str_type_tuple), s=25)
                
                # build comment and write results
                comment = '{:^2} {:5}'.format('#', str_type)
                f.write('{:^2} {:^2} {:^2} {:^2} {:^2} {:^2} {:^6}\n'.format(i, new_type, id1, id2, id3, id4, comment))