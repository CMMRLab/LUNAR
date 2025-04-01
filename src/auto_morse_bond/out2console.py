# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
June 31st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to write log file
def out(mm, log, version, min_bond_length, coeffs2skip, zero_effected_xterms, alpha_specs, alpha_scale, include_type_labels):
    
    ##################
    # Write settings #
    ##################
    log.out('\n\n\nauto_morse_bond settings:')
    log.out(f'- min_bond_length:      {min_bond_length}')
    log.out(f'- coeffs2skip:          {str(coeffs2skip)}')
    log.out(f'- alpha_specs:          {str(alpha_specs)}')
    log.out(f'- alpha_scale:          {alpha_scale}')
    log.out(f'- zero_effected_xterms: {zero_effected_xterms}')
    log.out(f'- include_type_labels:  {include_type_labels}')
    
    ##########################################
    # Write the elements found in the system #
    ##########################################
    log.out('\n\nElements found in system:')
    for i in mm.elements:
        log.out('{} {}'.format('-', i))
     
    ####################################
    # Write molecules/cluster findings #
    ####################################
    # Write molecules table          
    log.out('\n\n\n-----------------------------------------------------Cluster Analysis----------------------------------------------')
    log.out('{:^10} {:^16} {:^25} {:^20} {:^20} {:^20}'.format('##', 'Size of Cluster', 'Mass', '%Mass Cluster', '%Size Cluster', 'Cluster Formula'))
    log.out('-------------------------------------------------------------------------------------------------------------------')  
    for i in mm.molecules.data:
        data = mm.molecules.data[i]
        size = '{: >6}'.format(data.size)
        mass = '{:.2f}'.format(data.mass)
        pmass = '{:.2f}'.format(data.pmass)
        psize = '{:.2f}'.format(data.psize)
        formula = '{:^10}'.format(data.formula)
        log.out('{:^10} {:^16} {:^25} {:^20} {:^20} {:^20}'.format(i, size, mass, pmass, psize, formula))
        
        
    ######################
    # Write bond lengths #
    ######################
    log.out('\n\n\nBond length and bond type statistics and info:')
    log.out('Bond order (BO) will let you know what type of bond the code determine to coeff was from atom types in bond. BO options:')
    log.out('- single:   1')
    log.out('- double:   2')
    log.out('- triple:   3')
    log.out('- aromatic: 1.5\n')
    log.out('Parameters (parms) will let you know how the morse parms were found. parms options:')
    log.out('- classN: coeff left as class1 or class2')
    log.out('- file:   dissociation energy was set by morsefile and and alpha was minized by the code')
    log.out('----------------------------------------------------Bond type bond length statistics----------------------------------------------------')
    log.out('{:<5} {:^45} {:^10} {:^10} {:^10} {:^15} {:^15} {:^15}'.format('Bond', 'Bond Type: atom 1-2 (element,ring,nb)', 'Bond', 'Bond', 'Bond Length', 'Bond Length', 'Bond Length', 'Bond Length'))
    log.out('{:<5} {:^45} {:^10} {:^10} {:^10} {:^15} {:^15} {:^15}'.format('id','BO: 1, 2, 3, 0; parms: file, classN;', 'Count', 'r0', 'Average', 'Minimum', 'Maximum', 'Standard Deviation'))
    log.out('----------------------------------------------------------------------------------------------------------------------------------------')
    for i in mm.bonds_lengths.types:
        bond = mm.bonds_lengths.types[i]
        r0 = mm.bond_coeffs[i].coeffs[0]
        btype = mm.bond_info.types[i]
        log.out('{:<5} {:^45} {:^10} {:^10} {:^10} {:^15} {:^15} {:^15}'.format(i, btype, bond.count, r0, bond.avg, bond.min, bond.max, bond.std))
        
        
    ###########################
    # Write ringed atoms data #
    ###########################
    # Write all inputs of find_pyramidalization
    log.out('\n\n\n----Inputs used for find_rings----')
    for i in mm.find_rings:
        log.out('{}:  {}'.format(i, mm.find_rings[i]))
    
    log.out('\n')
    log.out('{} {}'.format('Checked for ring sizes : ', mm.find_rings['rings2check']))
    log.out('{} {}'.format('Total rings found : ', mm.rings.total))
    
    log.out(f'{mm.rings.partitioned_count} atoms along a ring seam had to be partitioned amoung the ring') 
    log.out('types (To find accurate %Mass of ring type to entire system).')
    log.out('Giving preference of partionioning in this order:')
    log.out('- 6 member ring')
    log.out('- 5 member ring')
    log.out('- 7 member ring')
    log.out('- 4 member ring')
    log.out('- 3 member ring')
    log.out('- 8 member ring')
    log.out('- minimum ring size')
    log.out('*NOTE: If count of rings exists, but no atoms exist for that ring, this means the')
    log.out('atoms for that ring were partionted to other rings that the atoms belong to.*\n')                                                                                                     
    for i in mm.rings.data:
        data = mm.rings.data[i]
        count = '{:^5}'.format(data.count)
        pcount = '{:.2f}'.format(data.pcount)
         
        # If count is greater then zero write
        if data.count > 0:
            log.out('---------------------------------------------------------------------------------------------')
            log.out('|{:^25} | {:^25} | {:^35}|'.format('Ring', 'Count', '%Ring count'))
            log.out('|{:^25} | {:^25} | {:^35}|'.format('Type', 'of Rings', 'per all rings'))
            log.out('---------------------------------------------------------------------------------------------')
            log.out('|{:^25} | {:^25} | {:^35}|'.format(i, count, pcount))
            log.out('---------------------------------------------------------------------------------------------')
            log.out('|{:^16} | {:^15} | {:^16} | {:^16} | {:^16}|'.format('Element', 'natoms', 'Mass', '%Mass', '%natoms'))
            log.out('---------------------------------------------------------------------------------------------')
            for j in mm.rings.data[i].partitioned:
                data1 = mm.rings.data[i].partitioned[j]
                size = '{:^5}'.format(data1.size)
                mass = '{:.2f}'.format(data1.mass)
                pmass = '{:.2f}'.format(data1.pmass)
                mass = '{:.2f}'.format(data1.mass)
                psize = '{:.2f}'.format(data1.psize)
                log.out('|{:^16} | {:^15} | {:^16} | {:^16} | {:^16}|'.format(j, size, mass, pmass, psize))
            log.out('---------------------------------------------------------------------------------------------\n\n')

    ########################
    # Write skipped coeffs #
    ########################
    for i in mm.bond_info.messages:
        log.out(f'{i}')
    log.out('\n\n')
        
    ###################################################################
    # Write recommendations to starting a simulation after conversion #
    ###################################################################
    log.out('\n\nRecommendation for Morse bond conversion process:') 
    log.out(' Step1:')
    log.out('  Fully equilibrate the system in the harmonic form of the force field in NVT or NPT (preferably')
    log.out('  NPT) as the harmonic form can generate large restoring forces to ensure the system is at the')
    log.out('  lowest energy state as possible. Then convert the harmonic bonds to morse bonds and apply and')
    log.out('  constraint to class2 crossterms.')
    
    log.out('\n Step2:')
    log.out('  # Morse bond simulation initialization')
    log.out('  timestep       0.5 # may need to be changed to 0.1')
    log.out('  fix            1 all nve/limit 0.1')
    log.out('  run            20000 # (10000/0.5 = 20000) = 10ps with 0.5 dt')
    log.out('  write_data     morse_bond_initialization.data')
    log.out('  unfix          1\n')
    log.out('  # Start simulation that is desired')
    log.out('  :\n')
    
    log.out('  # This is because the conversion from a harmonic form to a Morse bond form of a force field')
    log.out('  # may result in discontinuities in the energy and gradients (forces) during the first few')
    log.out('  # timesteps due the geometries resting in a different potential energy enviroment.')
    
    ########################################################
    # Write out notes about commands like "fix bond/react" #
    ########################################################
    log.out('\n\n*NOTE this code can only update coeff types if atoms are used by the bond types, thus')
    log.out('it can be dangerous to use the updated force field parameters from one system to the next')
    log.out('if using LAMMPS commands like "fix bond/react" as those type of commands are constantly')
    log.out('switching out which part of the force field is being used by the current system.')
    
    log.out('\nThus for use with the "fix bond/react" command you SHOULD update each force field for each')
    log.out('system and not "mix-and-match" morse bond updates.*')
    
    return