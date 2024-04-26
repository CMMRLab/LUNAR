# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.5
January 5th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.read_lmp as read_lmp
import os



# Function to read -> merge -> add missing data from each file that atom_typing will be able to read
def merge(topofile, mass_map, ff_class, log):
        
    ################################################
    # read topofile and bondfile (when applicable) #
    ################################################
    if os.path.isfile(topofile):
        m = read_lmp.Molecule_File(topofile, method='forward', sections=['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers'])
        log.out(f'Read in {m.filename} LAMMPS datafile')
    else: log.error(f'ERROR lammps datafile: {topofile} does not exist')
            
    # If read in a LAMMPS .data file it will not have an element 
    # attribute in m.atoms so find element symbol from mass_map
    elements = set(); total_system_mass = 0; total_system_size = 0;
    for i in m.atoms:
        # Find element from mass_map and add to atom instance
        atom = m.atoms[i]; mass = m.masses[atom.type].coeffs[0];
        total_system_mass += mass; total_system_size += 1;
        try:
            element = [i for i in mass_map if mass in mass_map[i]][0];
            atom.element = element; atom.comment = element; elements.add(element)
        except: log.error(f'ERROR Not all masses in {topofile} are in the mass_map dictionary. AtomID: {i} atomTypeID: {m.atoms[i].type}')
            
    # Add element symbol to m.masses[ID].type
    for i in m.masses:
        mass = m.masses[i]
        mass.element = [i for i in mass_map if mass.coeffs[0] in mass_map[i]][0]
        
    # add elements, total_system_mass, and total_system_size to m
    m.elements = sorted(elements);
    m.total_system_mass = total_system_mass
    m.total_system_size = total_system_size
    
    # Check that bond coeffs exists else exit code
    if len(m.bond_coeffs) == 0:
        log.error(f'ERROR No bond coeffs currently exist in the read {topofile}')
                
    # Check for harmonic or class2 coeffs
    r0_range = [0.0, 3.0] # Set min/max r0s to check that r0 is at index 0
    for i in m.bond_coeffs:
        coeffs = m.bond_coeffs[i].coeffs
        if ff_class == 1 and len(coeffs) != 2: log.error(f'ERROR read in {topofile} Bond Coeff {i} is not harmonic (class1)')
        if ff_class == 2 and len(coeffs) != 4: log.error(f'ERROR read in {topofile} Bond Coeff {i} is not class2')
        
        # Get r0 based on harmonic vs class2
        try:
            if ff_class == 1: k, r0 = coeffs
            if ff_class == 2: r0, k2, k3, k4 = coeffs
                
            # check that r0 makes sense
            if r0 < min(r0_range): log.error(f'ERROR read in {topofile} Bond Coeff {i} r0 does not make sense: {r0}')
            if r0 > max(r0_range): log.error(f'ERROR read in {topofile} Bond Coeff {i} r0 does not make sense: {r0}')
        except: log.error(f'ERROR Read in topofile FF class is not supported or not understood. Please make sure ff_class: {ff_class} is consistent with {topofile} setup')
    return m
