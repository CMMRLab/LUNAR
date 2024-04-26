# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.4
April 2nd, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
# Function to check if xterm r0 was updated to a morse bond
def check_morse_r0s(xterm_r0s, effected_r0s, zero_all_xterm_r0s):
    checked = []; return_boolean = False
    for i in xterm_r0s:
        if i in effected_r0s:
            checked.append(i)
    
    # If checked update return_boolean to True
    if checked: return_boolean = True
    
    # Zero all xterm r0s if carbon based system
    if zero_all_xterm_r0s: return_boolean = True
    return return_boolean

# Function for writing lammps datafile
def zero(m, log):
    # True or False to zero all xterms with r0s (True) or to check if r0 exists in a morse bond (False)
    zero_all_xterm_r0s = True
    
    # IFF-RX0 hc/hpan matching of pair coeffs to update to new pair coeff
    hc_pair_coeff_match = [[0.02, 2.3], [0.02, 2.7], [0.02, 2.995]] # replace any of these pair coeffs
    hc_new_pair_coeff = [0.02, 2.0] # with this pair coeff for IFF-RX0 implementation
    #hc_pair_coeff_match = [];
    
    # xterms to zero that have r0s. Put acronym in xterms2zero list to zero specifc xterm. Acronyms and meanings:
    #    'bb'    = bondbond (computed with angles)
    #    'ba'    = bondangle (computed with angles)
    #    'ebt'   = endbondtorion (computed with dihedrals)
    #    'mbt'   = middlebondtorion (computed with dihedrals)
    #    'bb13'  = bondbond13 (computed with dihedrals)
    #    'all'   = shortcut to zero all xterms with r0s
    #    'every' = to zero all xterms with r0s, even if the r0 was found not to be used by a morse bond
    xterms2zero = ['bb', 'ba', 'ebt', 'mbt', 'bb13']
    xterms2zero = ['all']
    if 'every' in xterms2zero: zero_all_xterm_r0s = True
    log.out('\n\n\nZero crossterms that were effected by the morse bond update have been ZEROED. ZEROED crossterm statistics:')
        
    
    # Update pair_coeff for hc/hpan if class2 force field (IFF-RX0)
    for i in m.pair_coeffs: 
        pair = m.pair_coeffs[i]
        coeffs = pair.coeffs
        comment = '{:^3} {:^3}'.format('#', m.masses[i].element)
        if coeffs in hc_pair_coeff_match and m.masses[i].element == 'H':
            log.out('Pair Coeffs {} was found to be either hc/hpan and was switched to {} due to'.format(i, str(coeffs)))
            log.out('zeroing of xterms leading to usage of IFF-RX0 requiring smaller hydrogen atom to maintain proper density\n\n')
            comment += ' vdw sigma changed due to usage of IFF-RX0. Orginal coeff: {}'.format(str(coeffs))
            pair.type = comment
            pair.coeffs = hc_new_pair_coeff
            
    # update bondbond coeffs
    if m.bondbond_coeffs:
        tally = 0
        for i in m.bondbond_coeffs: 
            bondbond = m.bondbond_coeffs[i]
            comment = ''
            
            # Zero effect xterms if user wants
            if 'bb' in xterms2zero or 'all' in xterms2zero:
                r12 = bondbond.coeffs[1]; r23 = bondbond.coeffs[2]; r0s = [r12, r23];
                if check_morse_r0s(r0s, m.bond_info.effected_r0s, zero_all_xterm_r0s):
                    tally += 1
                    bondbond.type = 'effected by morse bond update - ZEROED'
                    bondbond.coeffs = [0, r12, r23]
        nxterms = len(m.bondbond_coeffs);  percent = 100*tally/nxterms;
        log.out('{} of {} ({:.2f}%) BondBond coeffs had r0s that were effected by the morse bond update'.format(tally, nxterms, percent))

    # Write bondangle coeffs
    if m.bondangle_coeffs:
        tally = 0
        for i in m.bondangle_coeffs: 
            bondangle = m.bondangle_coeffs[i]
            comment = ''
            
            # Zero effect xterms if user wants
            if 'ba' in xterms2zero or 'all' in xterms2zero:
                r12 = bondangle.coeffs[2]; r23 = bondangle.coeffs[3]; r0s = [r12, r23];
                if check_morse_r0s(r0s, m.bond_info.effected_r0s, zero_all_xterm_r0s):
                    tally += 1  
                    bondangle.type = 'effected by morse bond update - ZEROED'
                    bondangle.coeffs = [0, 0, r12, r23]
        nxterms = len(m.bondangle_coeffs);  percent = 100*tally/nxterms;
        log.out('{} of {} ({:.2f}%) BondAngle coeffs had r0s that were effected by the morse bond update'.format(tally, nxterms, percent))
        
    # Write endbondtorsion coeffs
    if m.endbondtorsion_coeffs:
        tally = 0
        for i in m.endbondtorsion_coeffs: 
            endbondtorsion = m.endbondtorsion_coeffs[i]
            comment = '';
            
            # Zero effect xterms if user wants
            if 'ebt' in xterms2zero or 'all' in xterms2zero:
                r12 = endbondtorsion.coeffs[6]; r34 = endbondtorsion.coeffs[7]; r0s = [r12, r34];
                if check_morse_r0s(r0s, m.bond_info.effected_r0s, zero_all_xterm_r0s):
                    tally += 1     
                    endbondtorsion.type = 'effected by morse bond update - ZEROED'
                    endbondtorsion.coeffs = [0, 0, 0, 0, 0, 0, r12, r34]
        nxterms = len(m.endbondtorsion_coeffs);  percent = 100*tally/nxterms;
        log.out('{} of {} ({:.2f}%) EndBondTorsion coeffs had r0s that were effected by the morse bond update'.format(tally, nxterms, percent))

    # Write middlebondtorsion coeffs
    if m.middlebondtorsion_coeffs:
        tally = 0
        for i in m.middlebondtorsion_coeffs: 
            middlebondtorsion  = m.middlebondtorsion_coeffs[i]
            comment = ''
            
            # Zero effect xterms if user wants
            if 'mbt' in xterms2zero or 'all' in xterms2zero:
                r23 = middlebondtorsion.coeffs[3]; r0s = [r23];
                if check_morse_r0s(r0s, m.bond_info.effected_r0s, zero_all_xterm_r0s):
                    tally += 1
                    middlebondtorsion.type = 'effected by morse bond update - ZEROED'
                    middlebondtorsion.coeffs = [0, 0, 0, r23]
        nxterms = len(m.middlebondtorsion_coeffs);  percent = 100*tally/nxterms;
        log.out('{} of {} ({:.2f}%) MiddleBondTorsion coeffs had r0s that were effected by the morse bond update'.format(tally, nxterms, percent))

    # Write bondbond13 coeffs
    if m.bondbond13_coeffs:
        tally = 0
        for i in m.bondbond13_coeffs: 
            bondbond13 = m.bondbond13_coeffs[i]
            comment = ''
            
            # Zero effect xterms if user wants
            if 'bb13' in xterms2zero or 'all' in xterms2zero:
                r12 = bondbond13.coeffs[1]; r34 = bondbond13.coeffs[2]; r0s = [r12, r34];
                if check_morse_r0s(r0s, m.bond_info.effected_r0s, zero_all_xterm_r0s):
                    tally += 1
                    bondbond13.type = 'effected by morse bond update - ZEROED'
                    bondbond13.coeffs = [0, r12, r34]
        nxterms = len(m.bondbond13_coeffs);  percent = 100*tally/nxterms;
        log.out('{} of {} ({:.2f}%) BondBond13 coeffs had r0s that were effected by the morse bond update'.format(tally, nxterms, percent))
    return m
