# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
October 8th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

def update(neutralize, parameters, nta_file, charges, log):
    #############################################################
    # Update charges from charges if any exists in charges dict #
    #############################################################
    if charges:
        log.out(f'\n\nThe read in nta file {nta_file}')
        log.out('contains charge section(s) of either id or type or nta or combinations.')
        log.out('These charges are now being applied after possible bond-incs search.')
        total_charge = 0
        for ID in charges:
            atom = parameters.atoms[ID]
            atom.charge = charges[ID]
            atom.comments = '{} [ {} ]'.format(atom.comments, 'charge set by nta file')
            total_charge += atom.charge
        log.out(f'  {len(charges)} atomIDs had charges updated via flags set in nta file')
        log.out(f'  Total system charge: {round(total_charge, 6)}')
        
    ###########################################################
    # neutralize system charge based on user specified method #
    ###########################################################
    # if "neutralize system charge all" add a fixed value to all atoms in the system to neutralize the system charge
    if neutralize['all']:
        log.out(f'\n\nApplying "neutralize system charge all" option which was set in {nta_file} file')
        natoms = len(parameters.atoms)
        if natoms > 0:
            sys_charge = sum([parameters.atoms[ID].charge for ID in parameters.atoms]) # find current system charge
            fixed_charge = -sys_charge/natoms # Find fixed charge value
            total_charge = 0
            for ID in parameters.atoms:
                atom = parameters.atoms[ID]
                atom.charge += fixed_charge
                atom.comments = '{} [ {} {} {} {} ]'.format(atom.comments, 'fixed charged added: ', fixed_charge, 'via', '"neutralize system charge all"')
                total_charge += atom.charge
            log.out(f'  Total system charge: {round(total_charge, 6)} after adding {round(fixed_charge, 6)} to "all" atomIDs')
        else:
            log.out('number of atoms that meet the "all" criteria is ZERO. NO fixed charge was added to any atomID')
            
    # if "neutralize system charge bond-inc" add a fixed value to all atoms in the system to neutralize the system charge
    if neutralize['bond-inc']:
        log.out(f'\n\nApplying "neutralize system charge bond-inc" option which was set in {nta_file} file')
        natoms = len(parameters.atoms) - len(charges)
        if natoms > 0:
            sys_charge = sum([parameters.atoms[ID].charge for ID in parameters.atoms if ID not in charges]) # find current system charge
            fixed_charge = -sys_charge/natoms # Find fixed charge value
            total_charge = 0
            for ID in parameters.atoms:
                if ID not in charges:
                    atom = parameters.atoms[ID]
                    atom.charge += fixed_charge
                    atom.comments = '{} [ {} {} {} {} ]'.format(atom.comments, 'fixed charged added: ', fixed_charge, 'via', '"neutralize system charge bond-inc"')
                total_charge += atom.charge
            log.out(f'  Total system charge: {round(total_charge, 6)} after adding {round(fixed_charge, 6)} to "bond-inc" atomIDs')
        else:
            log.out('number of atoms that meet the "bond-inc" criteria is ZERO. NO fixed charge was added to any atomID')
        
    # if "neutralize system charge user-defined" add a fixed value to all atoms in the system to neutralize the system charge
    if neutralize['user-defined']:
        log.out(f'\n\nApplying "neutralize system charge user-defined" option which was set in {nta_file} file')
        natoms = len(charges)
        if natoms > 0:
            sys_charge = sum([parameters.atoms[ID].charge for ID in parameters.atoms if ID in charges]) # find current system charge
            fixed_charge = -sys_charge/natoms # Find fixed charge value
            total_charge = 0
            for ID in parameters.atoms:
                if ID in charges:
                    atom = parameters.atoms[ID]
                    atom.charge += fixed_charge
                    atom.comments = '{} [ {} {} {} {} ]'.format(atom.comments, 'fixed charged added: ', fixed_charge, 'via', '"neutralize system charge user-defined"')
                total_charge += atom.charge
            log.out(f'  Total system charge: {round(total_charge, 6)} after adding {round(fixed_charge, 6)} to "user-defined" atomIDs')
        else:
            log.out('number of atoms that meet the "user-defined" criteria is ZERO. NO fixed charge was added to any atomID')
            
    # if "neutralize system charge zero" add a fixed value to all atoms in the system to neutralize the system charge
    if neutralize['zero']:
        log.out(f'\n\nApplying "neutralize system charge all" option which was set in {nta_file} file')
        natoms = len([parameters.atoms[ID].charge for ID in parameters.atoms if abs(parameters.atoms[ID].charge) == 0.0])
        if natoms > 0:
            sys_charge = sum([parameters.atoms[ID].charge for ID in parameters.atoms]) # find current system charge
            fixed_charge = -sys_charge/natoms # Find fixed charge value
            total_charge = 0
            for ID in parameters.atoms:
                atom = parameters.atoms[ID]
                if abs(atom.charge) == 0.0:
                    atom.charge += fixed_charge
                    atom.comments = '{} [ {} {} {} {} ]'.format(atom.comments, 'fixed charged added: ', fixed_charge, 'via', '"neutralize system charge zero"')
                total_charge += atom.charge
            log.out(f'  Total system charge: {round(total_charge, 6)} after adding {round(fixed_charge, 6)} to "all" atomIDs')
        else:
            log.out('number of atoms that meet the "zero" criteria is ZERO. NO fixed charge was added to any atomID')
    return parameters