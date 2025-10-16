# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
July 31st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.auto_morse_bond.class2xe_plotting as plotting
import traceback
import math


# Set which xterms are plotted
plot_xterm = {'bb':   True,
              'ba':   True,
              'ebt':  True,
              'mbt':  True,
              'bb13': True,
              }

# plot_xterm['bb'] =   False
# plot_xterm['ba'] =   False
# plot_xterm['ebt'] =  False
# plot_xterm['mbt'] =  False
# plot_xterm['bb13'] = False


# Main linker function for updating crossterms that use r0's
def update(m, morsefile, potential_styles, include_rcut, basename, log):
    # Option to force unused coeffs (such as when using "fix bond/react" to map
    # coeffs to different topologies, leaving unused coeffs) to zero's. This is
    # useful as it will make all bond coeffs be 'morse', all angle and dihedral
    # coeffs be 'class2xe' as LAMMPS "write_data" command currently doesnt write
    # hybrid forms of potentials.
    force_unused_to_zeros = True # True or False
    force_unused_comments = True # True or False to update comments of coeffs
    
    # Option to use rcut in convert crossterms (True or False). NOTE the LAMMPS
    # 'class2xe' angle and dihedral potentials do not currently have the rcut
    # shifts implemented. So for now leave this as False ...
    use_rcut_in_xterms = False
    
    # Set old/new syle names
    style_old = 'class2' # orginal LAMMPS style
    style_new = 'class2xe' # new LAMMPS style
    log.out('\n\nUpdating class2 force field to {} force field ....'.format(style_new))
    
    # Check that all bond coeffs are morse
    log.out('  - checking morse bond updates for consistency')
    forced_to_zeros = False
    for i in m.bond_coeffs:
        coeffs = m.bond_coeffs[i].coeffs
        btype = m.bond_info.types[i]
        if coeffs[0] != 'morse' and 'Not Typed' not in btype:
            parms = ' '.join([str(i) for i in coeffs])
            log.out('')
            log.out(f'ERROR bond coeff {i} "{parms}" has not been')
            log.out('  converted to a Morse bond. All bond coeffs must be Morse bond')
            log.out('  for class2 Morse FF. Please update Morse file')
            log.out(f'  "{morsefile}" with bond type')
            log.error(f'  "{btype}"')
        elif coeffs[0] != 'morse' and 'Not Typed' in btype and force_unused_to_zeros:
            style, r0, k2, k3, k4 = coeffs
            d = 0; alpha = 0; rcut = 2*r0;
            if include_rcut:
                m.bond_coeffs[i].coeffs = ['morse', d, alpha, r0, rcut]
            else:
                m.bond_coeffs[i].coeffs = ['morse', d, alpha, r0]
            if force_unused_comments: m.bond_coeffs[i].type += '   {}-forced-to-zeros/NO-atoms-use-coeff'.format(style_new)
            log.warn( '    WARNING forcing bond coeff {} to "morse {} {} {}" as it is not used by any bonding atoms and using {} FF'.format(i, d, alpha, r0, style_new) )
            forced_to_zeros = True
    if forced_to_zeros:
        log.out("")
        log.out("    Forcing zero's was used to avoid hybrid styles, as LAMMPS 'write_data' command currently does not")
        log.out("    write out hybrid style coeffs and '{}' requires most of the force field to be adjusted. This".format(style_new))
        log.out("    is a convenience option. If you would like zero's to not be forced set 'force_unused_to_zeros' to")
        log.out("    False in LUNAR/src/auto_morse_bond/class2xe_conversion.py in the 'update' function. Additionally you")
        log.out("    may shut off the comment appending in the write datafile by setting 'force_unused_comments' to False.")
            
    # Set map for bond atomIDs to bond type
    bond2type = {} # {tuple(sorted[id1, id2]): bond-type}
    for i in m.bonds:
        bond = m.bonds[i]
        atomids = tuple(sorted(list(bond.atomids)))
        bond2type[atomids] = bond.type
            
    # Update crossterm coeffs that are effected by a Morse bond inclusion
    m = angle_crossterm(m, bond2type, style_old, style_new, potential_styles, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log)
    m = dihedral_crossterm(m, bond2type, style_old, style_new, potential_styles, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log)
    
    return m

# Function to find most frequent occurance in list
def most_frequent(List):
    return max(set(List), key = List.count)

# Function to find bond type from atomIDs in bond
def find_bond_type(bond2type, id1, id2):
    atomids = tuple(sorted([id1, id2]))
    try: bondtype = bond2type[atomids]
    except: bondtype = None
    return bondtype

# Function to get sign of a value
def sign(value):
    sign = 1
    try:
        if value >= 0: sign = 1
        else: sign = -1
    except: pass
    return sign

# Function to update angle crossterm
def angle_crossterm(m, bond2type, style_old, style_new, potential_styles, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log):
    # Find bond 12 and bond 23 type from angles
    angle_bond_types_master = {i:([], []) for i in m.angle_coeffs} # { angle-type: ([bondtype1, ...], [bondtype2, ...]) }
    for i in m.angles:
        angle = m.angles[i]
        id1, id2, id3 = angle.atomids
        bondtype_12 = find_bond_type(bond2type, id1, id2)
        bondtype_23 = find_bond_type(bond2type, id2, id3)
        if bondtype_12 is not None:
            angle_bond_types_master[angle.type][0].append(bondtype_12)
        if bondtype_23 is not None:
            angle_bond_types_master[angle.type][1].append(bondtype_23)
    
    # Find most frequently used bond types
    angle_bond_types = {i:() for i in m.angle_coeffs} # { angle-type: (bondtype1, bondtype2) }
    for i in angle_bond_types_master:
        types1, types2 = angle_bond_types_master[i]
        if types1 and types2:
            type1 = most_frequent(types1)
            type2 = most_frequent(types2)
            angle_bond_types[i] = (type1, type2)
            
    # Update bondbond, and bondangle coeffs
    m = update_bondbond_coeffs(m, angle_bond_types, style_old, style_new, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log)
    m = update_bondangle_coeffs(m, angle_bond_types, style_old, style_new, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log)
    
    # Check that angle, bondbond, and bondangle use the same potential for each coeffID
    angle_style_check = {} # {angle-CoeffID : style }
    for i in m.angle_coeffs:
        angle = m.angle_coeffs[i].coeffs
        bondbond = m.bondbond_coeffs[i].coeffs
        bondangle = m.bondangle_coeffs[i].coeffs
        
        # Default to old style if there is a mismatch
        if bondbond[0] != bondangle[0]: style = style_old
        else: style = bondbond[0]
        angle_style_check[i] = style
        
    # If there is a mismatch between bondbond and bondangle adjust accordingly
    angle_styles = set(); 
    for i in angle_style_check:
        style = angle_style_check[i]
        angle_styles.add(style)
        
        # Get current coeffs
        a_coeffs = m.angle_coeffs[i].coeffs
        bb_coeffs = m.bondbond_coeffs[i].coeffs
        ba_coeffs = m.bondangle_coeffs[i].coeffs
        
        # Ensure angle, bondbond, and bondangle have consistent style types amoungst the CoeffIDs
        a_coeffs.insert(0, style)
        if style not in bb_coeffs:
            bb_coeffs = m.bondbond_coeffs[i].harmonic
            bb_coeffs.insert(0, style)
        if style not in ba_coeffs:
            ba_coeffs = m.bondangle_coeffs[i].harmonic
            ba_coeffs.insert(0, style)
        
        # Update new coeff setup
        m.angle_coeffs[i].coeffs = a_coeffs
        m.bondbond_coeffs[i].coeffs = bb_coeffs
        m.bondangle_coeffs[i].coeffs = ba_coeffs
    
    # If only one style is used remove style name from coeffs list
    if len(angle_styles) == 1:
        for i in angle_style_check:
            style = angle_style_check[i]
            a_coeffs = m.angle_coeffs[i].coeffs
            bb_coeffs = m.bondbond_coeffs[i].coeffs
            ba_coeffs = m.bondangle_coeffs[i].coeffs
            a_coeffs.pop(0); bb_coeffs.pop(0); ba_coeffs.pop(0);
            m.angle_coeffs[i].coeffs = a_coeffs
            m.bondbond_coeffs[i].coeffs = bb_coeffs
            m.bondangle_coeffs[i].coeffs = ba_coeffs
        angle_style = list(angle_styles)[0]
        m.angle_coeffs_style_hint = angle_style
        m.bondbond_coeffs_style_hint = angle_style
        m.bondangle_coeffs_style_hint = angle_style
        potential_styles.append(('angle', angle_style))
    else: # else use hybrid style setup
        angle_style = ' '.join(angle_styles)
        m.angle_coeffs_style_hint = angle_style
        m.bondbond_coeffs_style_hint = angle_style
        m.bondangle_coeffs_style_hint = angle_style
        potential_styles.append(('angle', 'hybrid '+angle_style))
    return m

# Function to update bondbond coeffs
def update_bondbond_coeffs(m, angle_bond_types, style_old, style_new, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log):
    """
    output coeffs format:
      class2xe D alpha r1 r2 <r1cut r2cut>
      class2 M r1 r2
      
      where:
          class2xe - is the new angle style (if coeff could be updated)
              D - dissociation energy (Kcal/mol)
              alpha - unitless alpha parameter from Morse bond fit
              r1 - bond1 length
              r2 - bond2 length
              r1cut - bond1 length for rcut option (optional)
              r2cut - bond2 length for rcut option (optional)
          
          class2 - is the orginal class2 angle style (if no bonds are using that coeff)
              M - bondbond stiffness parameter  ({Kcal/mol}/{angstrom^2})
              r1 - bond1 length
              r2 - bond2 length
    """
    log.out('  - updating bondbond coeffs')
    
    # Find bond 1 and 2 for each bondbond coeff
    plotting_coeffs = {} # {TypeID:{'class2':[], 'class2xe':[], 'harmonic1':[], 'harmonic2':[], 'morse1':[], 'morse2':[], 'bond_types':[]}
    for i in m.bondbond_coeffs:
        coeff = m.bondbond_coeffs[i]
        coeff.harmonic = coeff.coeffs
        M, r1, r2 = coeff.coeffs
        use_standard = True
        use_zeros = False
        plot_coeffs = {'class2':coeff.coeffs}
        if len(angle_bond_types[i]) == 2:
            type1, type2 = angle_bond_types[i]
            plot_coeffs['bond_types'] = [type1, type2]
            plot_coeffs['harmonic1'] = m.bond_coeffs[type1].harmonic
            plot_coeffs['harmonic2'] = m.bond_coeffs[type2].harmonic
            plot_coeffs['morse1'] = m.bond_coeffs[type1].coeffs
            plot_coeffs['morse2'] = m.bond_coeffs[type2].coeffs
            if len(m.bond_coeffs[type1].coeffs) == 4:
                style1, D1, alpha1, r01  = m.bond_coeffs[type1].coeffs
                r1cut = 0
            else: style1, D1, alpha1, r01, r1cut  = m.bond_coeffs[type1].coeffs
            if len(m.bond_coeffs[type2].coeffs) == 4:
                style2, D2, alpha2, r02 = m.bond_coeffs[type2].coeffs
                r2cut = 0
            else: style2, D2, alpha2, r02, r2cut = m.bond_coeffs[type2].coeffs
            if style1 == 'morse' and style2 == 'morse':
                D = (D1 + D2)/2
                alpha = math.sqrt( abs(M)/ (2*D) )
                #if M < 0: alpha = 0
                D = sign(M)*D # TODO: Checking curvature direction
                if include_rcut and use_rcut_in_xterms:
                    coeff.coeffs = [style_new, D, alpha, r01, r02, r1cut, r2cut]
                else:
                    coeff.coeffs = [style_new, D, alpha, r01, r02]
                use_standard = False
            elif force_unused_to_zeros: use_zeros = True
        elif force_unused_to_zeros or use_zeros:
            D = 0; alpha = 0;
            if include_rcut and use_rcut_in_xterms:
                coeff.coeffs = [style_new, D, alpha, r1, r2, 2*r1, 2*r2]
            else:
                coeff.coeffs = [style_new, D, alpha, r1, r2]
            if force_unused_comments: m.bondbond_coeffs[i].type += '   {}-forced-to-zeros/NO-atoms-use-coeff'.format(style_new)
            log.warn( '    WARNING forcing bondbond coeff {} to "{} {} {} {} {}" as it is not used by any bonding atoms'.format(i, style_new, D, alpha, r1, r2) )
            use_standard = False
        if use_standard:
            coeff.coeffs = [style_old, M, r1, r2]
        
        # Save coeffs for the plot
        plot_coeffs['class2xe'] = coeff.coeffs
        plotting_coeffs[i] = plot_coeffs
        
    # Plot the results
    #basename = ''
    if basename and plot_xterm['bb']:
        try:
            basename = basename + '_bb'
            coeff_name = 'Bond/bond'
            plotting.bondbond(m, plotting_coeffs, basename, coeff_name)
        except:
            stack_trace_string = traceback.format_exc()
            log.out(f'ERROR {stack_trace_string}')
    return m

# Function to update bondangle coeffs
def update_bondangle_coeffs(m, angle_bond_types, style_old, style_new, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log):
    """
    output coeffs format:
      class2xe D1 D2 alpha1 alpha2 r1 r2  <r1cut r2cut thetacut>
      class2 N1 N2 r1 r2
      
      where:
          class2-morse - is the new angle style (if coeff could be updated)
              D1 - dissociation energy (Kcal/mol) for bond1
              D2 - dissociation energy (Kcal/mol) for bond2
              alpha1 - unitless alpha parameter from Morse bond fit for bond1
              alpha2 - unitless alpha parameter from Morse bond fit for bond2
              r1 - bond1 length
              r2 - bond2 length
              r1cut - bond1 length for rcut option (optional)
              r2cut - bond2 length for rcut option (optional)
              thetacut - angle for shift when bond "breaks" option (optional - default will be theta0 from angle)
          
          class2 - is the orginal class2 angle style (if no bonds are using that coeff)
              N1 - bondangle stiffness parameter ({Kcal/mol}/{angstrom}) for bond1
              N2 - bondangle stiffness parameter ({Kcal/mol}/{angstrom}) for bond1
              r1 - bond1 length
              r2 - bond2 length
    """
    log.out('  - updating bondangle coeffs')
    
    # Find bond 1 and 2 for each bondangle coeff
    plotting_coeffs = {} # {TypeID:{'class2':[], 'class2xe':[], 'harmonic1':[], 'harmonic2':[], 'morse1':[], 'morse2':[], 'bond_types':[]}
    for i in m.bondangle_coeffs:
        coeff = m.bondangle_coeffs[i]
        coeff.harmonic = coeff.coeffs
        N1, N2, r1, r2 = coeff.coeffs
        theta0, k2, k3, k4 = m.angle_coeffs[i].coeffs
        use_standard = True
        use_zeros = False
        plot_coeffs = {'class2':coeff.coeffs}
        plot_coeffs['angle'] = m.angle_coeffs[i].coeffs
        if len(angle_bond_types[i]) == 2:
            type1, type2 = angle_bond_types[i]
            plot_coeffs['bond_types'] = [type1, type2]
            plot_coeffs['harmonic1'] = m.bond_coeffs[type1].harmonic
            plot_coeffs['harmonic2'] = m.bond_coeffs[type2].harmonic
            plot_coeffs['morse1'] = m.bond_coeffs[type1].coeffs
            plot_coeffs['morse2'] = m.bond_coeffs[type2].coeffs
            if len(m.bond_coeffs[type1].coeffs) == 4:
                style1, D1, alpha1, r01  = m.bond_coeffs[type1].coeffs
                r1cut = 0
            else: style1, D1, alpha1, r01, r1cut  = m.bond_coeffs[type1].coeffs
            if len(m.bond_coeffs[type2].coeffs) == 4:
                style2, D2, alpha2, r02 = m.bond_coeffs[type2].coeffs
                r2cut = 0
            else: style2, D2, alpha2, r02, r2cut = m.bond_coeffs[type2].coeffs
            if style1 == 'morse' and style2 == 'morse':
                alpha1 = math.sqrt( abs(N1)/(2*D1) )
                alpha2 = math.sqrt( abs(N2)/(2*D2) )
                #if N1 < 0: alpha1 = 0
                #if N2 < 0: alpha2 = 0
                D1 = sign(N1)*D1 # TODO: Checking curvature direction
                D2 = sign(N2)*D2 # TODO: Checking curvature direction
                if include_rcut and use_rcut_in_xterms:
                    coeff.coeffs = [style_new, D1, D2, alpha1, alpha2, r01, r02, r1cut, r2cut, theta0]
                else:
                    coeff.coeffs = [style_new, D1, D2, alpha1, alpha2, r01, r02]
                use_standard = False
            elif force_unused_to_zeros: use_zeros = True
        elif force_unused_to_zeros or use_zeros:
            D1 = 0; alpha1 = 0; D2 = 0; alpha2 = 0;
            if include_rcut and use_rcut_in_xterms:
                coeff.coeffs = [style_new, D1, D2, alpha1, alpha2, r1, r2, 2*r1, 2*r2, theta0]
            else:
                coeff.coeffs = [style_new, D1, D2, alpha1, alpha2, r1, r2]
            if force_unused_comments: m.bondangle_coeffs[i].type += '   {}-forced-to-zeros/NO-atoms-use-coeff'.format(style_new)
            log.warn( '    WARNING forcing bondangle coeff {} to "{} {} {} {} {} {} {}" as it is not used by any bonding atoms'.format(i, style_new, D1, D2, alpha1, alpha2, r1, r2) )
            use_standard = False
                
        if use_standard:
            coeff.coeffs = [style_old, N1, N2, r1, r2]
            
        # Save coeffs for the plot
        plot_coeffs['class2xe'] = coeff.coeffs
        plotting_coeffs[i] = plot_coeffs

    # Plot the results
    #basename = ''
    if basename and plot_xterm['ba']:
        try:
            basename = basename + '_ba'
            plotting.bondangle(m, plotting_coeffs, basename)
        except:
            stack_trace_string = traceback.format_exc()
            log.out(f'ERROR {stack_trace_string}')
    return m

# Function to update dihedral crossterm
def dihedral_crossterm(m, bond2type, style_old, style_new, potential_styles, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log):
    # Find bond 12, bond 23, and bond 34 type from dihedrals
    dihedral_bond_types_master = {i:([], [], []) for i in m.dihedral_coeffs} # { dihedral-type: ([bondtype1, ...], [bondtype2, ...], [bondtype3, ...]) }
    for i in m.dihedrals:
        dihedral = m.dihedrals[i]
        id1, id2, id3, id4 = dihedral.atomids
        bondtype_12 = find_bond_type(bond2type, id1, id2)
        bondtype_23 = find_bond_type(bond2type, id2, id3)
        bondtype_34 = find_bond_type(bond2type, id3, id4)
        if bondtype_12 is not None:
            dihedral_bond_types_master[dihedral.type][0].append(bondtype_12)
        if bondtype_23 is not None:
            dihedral_bond_types_master[dihedral.type][1].append(bondtype_23)
        if bondtype_34 is not None:
            dihedral_bond_types_master[dihedral.type][2].append(bondtype_34)
    
    # Find most frequently used bond types
    dihedral_bond_types = {i:() for i in m.dihedral_coeffs} # { dihedral-type: (bondtype1, bondtype2, bondtype3) }
    for i in dihedral_bond_types_master:
        types1, types2, types3 = dihedral_bond_types_master[i]
        if types1 and types2 and types3:
            type1 = most_frequent(types1)
            type2 = most_frequent(types2)
            type3 = most_frequent(types3)
            dihedral_bond_types[i] = (type1, type2, type3)
            
    # Update bondbond13, ... and .. coeffs
    m = update_middlebondtorsion_coeffs(m, dihedral_bond_types, style_old, style_new, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log)
    m = update_endbondtorsion_coeffs(m, dihedral_bond_types, style_old, style_new, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log)
    m = update_bondbond13_coeffs(m, dihedral_bond_types, style_old, style_new, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log)
    
    # Check that dihedral, middlebondtorsion, endbondtorsion, and bondbond13 use the same potential for each coeffID
    dihedral_style_check = {} # {dihedral-CoeffID : style }
    for i in m.dihedral_coeffs:
        dihedral = m.dihedral_coeffs[i].coeffs
        middlebondtorsion = m.middlebondtorsion_coeffs[i].coeffs
        endbondtorsion = m.endbondtorsion_coeffs[i].coeffs
        bondbond13 = m.bondbond13_coeffs[i].coeffs
        
        # Default to old style if there is a mismatch
        if middlebondtorsion[0] != endbondtorsion[0] and endbondtorsion[0] != bondbond13[0]: style = style_old
        else: style = middlebondtorsion[0]
        dihedral_style_check[i] = style
        
    # If there is a mismatch between middlebondtorsion, endbondtorsion, and bondbond13 adjust accordingly
    dihedral_styles = set(); styles_reset = {}
    for i in dihedral_style_check:
        style = dihedral_style_check[i]
        dihedral_styles.add(style)
        
        # Get current coeffs
        d_coeffs = m.dihedral_coeffs[i].coeffs
        mbt_coeffs = m.middlebondtorsion_coeffs[i].coeffs
        ebt_coeffs = m.endbondtorsion_coeffs[i].coeffs
        bb13_coeffs = m.bondbond13_coeffs[i].coeffs
        
        # Ensure angle, bondbond, and bondangle have consistent style types amoungst the CoeffIDs
        d_coeffs.insert(0, style)
        styles_reset[i] = style
        if style not in mbt_coeffs:
            mbt_coeffs = m.middlebondtorsion_coeffs[i].harmonic
            mbt_coeffs.insert(0, style) 
        if style not in ebt_coeffs:
            ebt_coeffs = m.endbondtorsion_coeffs[i].harmonic
            ebt_coeffs.insert(0, style)
        if style not in bb13_coeffs:
            bb13_coeffs = m.bondbond13_coeffs[i].harmonic
            bb13_coeffs.insert(0, style)
        
        # Update new coeff setup
        m.dihedral_coeffs[i].coeffs = d_coeffs
        m.middlebondtorsion_coeffs[i].coeffs = mbt_coeffs
        m.endbondtorsion_coeffs[i].coeffs = ebt_coeffs
        m.bondbond13_coeffs[i].coeffs = bb13_coeffs
        
    # If only one style is used remove style name from coeffs list
    if len(dihedral_styles) == 1:
        for i in dihedral_style_check:
            style = dihedral_style_check[i]
            d_coeffs = m.dihedral_coeffs[i].coeffs
            mbt_coeffs = m.middlebondtorsion_coeffs[i].coeffs
            ebt_coeffs = m.endbondtorsion_coeffs[i].coeffs
            bb13_coeffs = m.bondbond13_coeffs[i].coeffs
            d_coeffs.pop(0); mbt_coeffs.pop(0); ebt_coeffs.pop(0); bb13_coeffs.pop(0);
            m.dihedral_coeffs[i].coeffs = d_coeffs
            m.middlebondtorsion_coeffs[i].coeffs = mbt_coeffs
            m.endbondtorsion_coeffs[i].coeffs = ebt_coeffs
            m.bondbond13_coeffs[i].coeffs = bb13_coeffs
        dihedral_style = list(dihedral_styles)[0]
        m.dihedral_coeffs_style_hint = dihedral_style
        m.middlebondtorsion_coeffs_style_hint = dihedral_style
        m.endbondtorsion_coeffs_style_hint = dihedral_style
        m.bondbond13_coeffs_style_hint = dihedral_style
        m.angleangletorsion_coeffs_style_hint = dihedral_style
        m.angletorsion_coeffs_style_hint = dihedral_style
        potential_styles.append(('dihedral', dihedral_style))
    else: # else use hybrid style setup
        dihedral_style = ' '.join(dihedral_styles)
        m.dihedral_coeffs_style_hint = dihedral_style
        m.middlebondtorsion_coeffs_style_hint = dihedral_style
        m.endbondtorsion_coeffs_style_hint = dihedral_style
        m.bondbond13_coeffs_style_hint = dihedral_style  
        m.angleangletorsion_coeffs_style_hint = dihedral_style
        m.angletorsion_coeffs_style_hint = dihedral_style
        potential_styles.append(('dihedral', 'hybrid '+dihedral_style))
        
        # Update crossterms that have not been updated to a hybrid format 
        for i in styles_reset:
            style = styles_reset[i]
            m.angleangletorsion_coeffs[i].coeffs.insert(0, style)
            m.angletorsion_coeffs[i].coeffs.insert(0, style)     
    return m

# Function to update middlebondtorsion coeffs
def update_middlebondtorsion_coeffs(m, dihedral_bond_types, style_old, style_new, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log):
    """
    output coeffs format:
      class2xe alpha A1 A2 A3 r2 <r2cut phi1cut phi2cut phi3cut>
      class2 A1 A2 A3 r2
      
      where:
          class2xe - is the new angle style (if coeff could be updated)
              alpha - unitless alpha parameter from Morse bond fit
              A1 - middlebondtorsion stiffness parameter for cos(1*phi) ({Kcal/mol}/{angstrom^2})
              A2 - middlebondtorsion stiffness parameter for cos(2*phi) ({Kcal/mol}/{angstrom^2})
              A3 - middlebondtorsion stiffness parameter for cos(3*phi) ({Kcal/mol}/{angstrom^2})
              r2 - bond2 length
              r2cut - bond2 length for rcut option (optional)
              phi1cut - phi for shift when bond "breaks" option (optional - default will be phi1 from dihedral)
              phi2cut - phi for shift when bond "breaks" option (optional - default will be phi2 from dihedral)
              phi3cut - phi for shift when bond "breaks" option (optional - default will be phi3 from dihedral)
          
          class2 - is the orginal class2 angle style (if no bonds are using that coeff)
              A1 - middlebondtorsion stiffness parameter for cos(1*phi) ({Kcal/mol}/{angstrom^2})
              A2 - middlebondtorsion stiffness parameter for cos(2*phi) ({Kcal/mol}/{angstrom^2})
              A3 - middlebondtorsion stiffness parameter for cos(3*phi) ({Kcal/mol}/{angstrom^2})
              r2 - bond2 length
    """
    log.out('  - updating middlebondtorsion coeffs')
    
    # Find bond 2 for each middlebondtorsion coeff
    plotting_coeffs = {} # {TypeID:{'class2':[], 'class2xe':[], 'harmonic':[], 'morse':[], 'bond_type1':[]}
    for i in m.middlebondtorsion_coeffs:
        coeff = m.middlebondtorsion_coeffs[i]
        coeff.harmonic = coeff.coeffs
        A1, A2, A3, r2, = coeff.coeffs
        k1, phi1, k2, phi2, k3, phi3 = m.dihedral_coeffs[i].coeffs
        use_standard = True
        use_zeros = False
        plot_coeffs = {'class2':coeff.coeffs}
        plot_coeffs['dihedral'] = m.dihedral_coeffs[i].coeffs
        if len(dihedral_bond_types[i]) == 3:
            type1, type2, type3 = dihedral_bond_types[i]
            plot_coeffs['bond_type'] = [type2]
            plot_coeffs['harmonic'] = m.bond_coeffs[type2].harmonic
            plot_coeffs['morse'] = m.bond_coeffs[type2].coeffs
            if len(m.bond_coeffs[type2].coeffs) == 4:
                style2, D2, alpha2, r02  = m.bond_coeffs[type2].coeffs
                r2cut = 0
            else: style2, D2, alpha2, r02, r2cut  = m.bond_coeffs[type2].coeffs
            if style2 == 'morse':
                if include_rcut and use_rcut_in_xterms:
                    coeff.coeffs = [style_new, alpha2, A1, A2, A3, r02, r2cut, phi1, phi2, phi3]
                else:
                    coeff.coeffs = [style_new, alpha2, A1, A2, A3, r02]
                use_standard = False
            elif force_unused_to_zeros: use_zeros = True
        elif force_unused_to_zeros or use_zeros:
            alpha2 = 0;
            if include_rcut and use_rcut_in_xterms:
                coeff.coeffs = [style_new, alpha2, A1, A2, A3, r2, 2*r2, phi1, phi2, phi3]
            else:
                coeff.coeffs = [style_new, alpha2, A1, A2, A3, r2]
            if force_unused_comments: m.middlebondtorsion_coeffs[i].type += '   {}-forced-to-zeros/NO-atoms-use-coeff'.format(style_new)
            log.warn( '    WARNING forcing middlebondtorsion coeff {} to "{} {} {} {} {} {}" as it is not used by any bonding atoms'.format(i, style_new, alpha2, A1, A2, A3, r2) )
            use_standard = False
        if use_standard:
            coeff.coeffs = [style_old, A1, A2, A3, r2]
            
        # Save coeffs for the plot
        plot_coeffs['class2xe'] = coeff.coeffs
        plotting_coeffs[i] = plot_coeffs
            
    # Plot the results
    #basename = ''
    if basename and plot_xterm['mbt']:
        try:
            basename = basename + '_mbt'
            plotting.middlebondtorsion(m, plotting_coeffs, basename)
        except:
            stack_trace_string = traceback.format_exc()
            log.out(f'ERROR {stack_trace_string}')
    return m

# Function to update endbondtorsion coeffs
def update_endbondtorsion_coeffs(m, dihedral_bond_types, style_old, style_new, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log):
    """
    output coeffs format:
      class2xe alpha1 alpha2 B1 B2 B3 C1 C2 C3 r1 r3  <r1cut r2cut phi1cut phi2cut phi3cut>
      class2 B1 B2 B3 C1 C2 C3 r1 r3
      
      where:
          class2xe - is the new angle style (if coeff could be updated)
              alpha1 - unitless alpha for bond1
              alpha2 - unitless alpha for bond2
              B1 - endbondtorsion stiffness parameter for cos(1*phi) ({Kcal/mol}/{angstrom^2})
              B2 - endbondtorsion stiffness parameter for cos(2*phi) ({Kcal/mol}/{angstrom^2})
              B3 - endbondtorsion stiffness parameter for cos(3*phi) ({Kcal/mol}/{angstrom^2})
              C1 - endbondtorsion stiffness parameter for cos(1*phi) ({Kcal/mol}/{angstrom^2})
              C2 - endbondtorsion stiffness parameter for cos(2*phi) ({Kcal/mol}/{angstrom^2})
              C3 - endbondtorsion stiffness parameter for cos(3*phi) ({Kcal/mol}/{angstrom^2})
              r1 - bond1 length
              r3 - bond3 length
              r1cut - bond1 length for rcut option (optional)
              r3cut - bond3 length for rcut option (optional)
              phi1cut - phi for shift when bond "breaks" option (optional - default will be phi1 from dihedral)
              phi2cut - phi for shift when bond "breaks" option (optional - default will be phi2 from dihedral)
              phi3cut - phi for shift when bond "breaks" option (optional - default will be phi3 from dihedral)
          
          class2 - is the orginal class2 angle style (if no bonds are using that coeff)
              B1 - endbondtorsion stiffness parameter for cos(1*phi) ({Kcal/mol}/{angstrom^2})
              B2 - endbondtorsion stiffness parameter for cos(2*phi) ({Kcal/mol}/{angstrom^2})
              B3 - endbondtorsion stiffness parameter for cos(3*phi) ({Kcal/mol}/{angstrom^2})
              C1 - endbondtorsion stiffness parameter for cos(1*phi) ({Kcal/mol}/{angstrom^2})
              C2 - endbondtorsion stiffness parameter for cos(2*phi) ({Kcal/mol}/{angstrom^2})
              C3 - endbondtorsion stiffness parameter for cos(3*phi) ({Kcal/mol}/{angstrom^2})
              r1 - bond1 length
              r3 - bond3 length
    """
    log.out('  - updating endbondtorsion coeffs')
    
    # Find bond 1 and 3 for each endbondtorsion coeff
    plotting_coeffs = {} # {TypeID:{'class2':[], 'class2xe':[], 'harmonic1':[], 'harmonic2':[], 'morse1':[], 'morse2':[], 'bond_types':[]}
    for i in m.endbondtorsion_coeffs:
        coeff = m.endbondtorsion_coeffs[i]
        coeff.harmonic = coeff.coeffs
        B1, B2, B3, C1, C2, C3, r1, r3 = coeff.coeffs
        k1, phi1, k2, phi2, k3, phi3 = m.dihedral_coeffs[i].coeffs
        use_standard = True
        use_zeros = False
        plot_coeffs = {'class2':coeff.coeffs}
        plot_coeffs['dihedral'] = m.dihedral_coeffs[i].coeffs
        if len(dihedral_bond_types[i]) == 3:
            type1, type2, type3 = dihedral_bond_types[i]
            plot_coeffs['bond_types'] = [type1, type3]
            plot_coeffs['harmonic1'] = m.bond_coeffs[type1].harmonic
            plot_coeffs['harmonic2'] = m.bond_coeffs[type3].harmonic
            plot_coeffs['morse1'] = m.bond_coeffs[type1].coeffs
            plot_coeffs['morse2'] = m.bond_coeffs[type3].coeffs
            if len(m.bond_coeffs[type1].coeffs) == 4:
                style1, D1, alpha1, r01  = m.bond_coeffs[type1].coeffs
                r1cut = 0
            else: style1, D1, alpha1, r01, r1cut  = m.bond_coeffs[type1].coeffs
            if len(m.bond_coeffs[type3].coeffs) == 4:
                style3, D3, alpha3, r03  = m.bond_coeffs[type3].coeffs
                r3cut = 0
            else: style3, D3, alpha3, r03, r3cut  = m.bond_coeffs[type3].coeffs
            if style1 == 'morse' and style3 =='morse':
                if include_rcut and use_rcut_in_xterms:
                    coeff.coeffs = [style_new, alpha1, alpha3, B1, B2, B3, C1, C2, C3, r01, r03, r1cut, r3cut, phi1, phi2, phi3]
                else:
                    coeff.coeffs = [style_new, alpha1, alpha3, B1, B2, B3, C1, C2, C3, r01, r03]
                use_standard = False
            elif force_unused_to_zeros: use_zeros = True
        elif force_unused_to_zeros or use_zeros:
            alpha1 = 0; alpha3 = 0;
            if include_rcut and use_rcut_in_xterms:
                coeff.coeffs = [style_new, alpha1, alpha3, B1, B2, B3, C1, C2, C3, r1, r3, 2*r1, 2*r3, phi1, phi2, phi3]
            else:
                coeff.coeffs = [style_new, alpha1, alpha3, B1, B2, B3, C1, C2, C3, r1, r3]
            if force_unused_comments: m.endbondtorsion_coeffs[i].type += '   {}-forced-to-zeros/NO-atoms-use-coeff'.format(style_new)
            log.warn( '    WARNING forcing edbondtorsion coeff {} to "{} {} {} {} {} {} {} {} {} {} {}" as it is not used by any bonding atoms'.format(i, style_new, alpha1, alpha3, B1, B2, B3, C1, C2, C3, r1, r3) )
            use_standard = False
        if use_standard:
            coeff.coeffs = [style_old, B1, B2, B3, C1, C2, C3, r1, r3]
            
        # Save coeffs for the plot
        plot_coeffs['class2xe'] = coeff.coeffs
        plotting_coeffs[i] = plot_coeffs
            
    # Plot the results
    #basename = ''
    if basename and plot_xterm['ebt']:
        try:
            basename = basename + '_ebt'
            plotting.endbondtorsion(m, plotting_coeffs, basename)
        except:
            stack_trace_string = traceback.format_exc()
            log.out(f'ERROR {stack_trace_string}')
    return m

# Function to update bondbond13 coeffs
def update_bondbond13_coeffs(m, dihedral_bond_types, style_old, style_new, force_unused_to_zeros, force_unused_comments, include_rcut, use_rcut_in_xterms, basename, log):
    """
    output coeffs format:
      class2xe D alpha r1 r3 <r1cut r3cut>
      class2 N r1 r3
      
      where:
          class2xe - is the new angle style (if coeff could be updated)
              D - dissociation energy (Kcal/mol)
              alpha - unitless alpha parameter from Morse bond fit
              r1 - bond1 length
              r3 - bond3 length
              r1cut - bond1 length for rcut option (optional)
              r3cut - bond3 length for rcut option (optional)
          
          class2 - is the orginal class2 angle style (if no bonds are using that coeff)
              N - bondbond13 stiffness parameter  ({Kcal/mol}/{angstrom^2})
              r1 - bond1 length
              r3 - bond3 length
    """
    log.out('  - updating bondbond13 coeffs')
    
    # Find bond 1 and 3 for each bondbond13 coeff
    plotting_coeffs = {} # {TypeID:{'class2':[], 'class2xe':[], 'harmonic1':[], 'harmonic2':[], 'morse1':[], 'morse2':[], 'bond_types':[]}
    for i in m.bondbond13_coeffs:
        coeff = m.bondbond13_coeffs[i]
        coeff.harmonic = coeff.coeffs
        N, r1, r3 = coeff.coeffs
        use_standard = True
        use_zeros = False
        plot_coeffs = {'class2':coeff.coeffs}
        if len(dihedral_bond_types[i]) == 3:
            type1, type2, type3 = dihedral_bond_types[i]
            plot_coeffs['bond_types'] = [type1, type3]
            plot_coeffs['harmonic1'] = m.bond_coeffs[type1].harmonic
            plot_coeffs['harmonic2'] = m.bond_coeffs[type3].harmonic
            plot_coeffs['morse1'] = m.bond_coeffs[type1].coeffs
            plot_coeffs['morse2'] = m.bond_coeffs[type3].coeffs
            if len(m.bond_coeffs[type1].coeffs) == 4:
                style1, D1, alpha1, r01  = m.bond_coeffs[type1].coeffs
                r1cut = 0
            else: style1, D1, alpha1, r01, r1cut  = m.bond_coeffs[type1].coeffs
            if len(m.bond_coeffs[type3].coeffs) == 4:
                style3, D3, alpha3, r03  = m.bond_coeffs[type3].coeffs
                r3cut = 0
            else: style3, D3, alpha3, r03, r3cut  = m.bond_coeffs[type3].coeffs
            if style1 == 'morse' and style3 == 'morse':
                D = (D1 + D3)/2
                alpha = math.sqrt( abs(N)/ (2*D) )
                #if N < 0: alpha = 0
                D = sign(N)*D # TODO: Checking curvature direction
                if include_rcut and use_rcut_in_xterms:
                    coeff.coeffs = [style_new, D, alpha, r01, r03, r1cut, r3cut]
                else:
                    coeff.coeffs = [style_new, D, alpha, r01, r03]
                use_standard = False
            elif force_unused_to_zeros: use_zeros = True
        elif force_unused_to_zeros or use_zeros:
            D = 0; alpha = 0;
            if include_rcut and use_rcut_in_xterms:
                coeff.coeffs = [style_new, D, alpha, r1, r3, 2*r1, 2*r3]
            else:
                coeff.coeffs = [style_new, D, alpha, r1, r3]
            if force_unused_comments: m.bondbond13_coeffs[i].type += '   {}-forced-to-zeros/NO-atoms-use-coeff'.format(style_new)
            log.warn( '    WARNING forcing bondbond13 coeff {} to "{} {} {} {} {}" as it is not used by any bonding atoms'.format(i, style_new, D, alpha, r1, r3) )
            use_standard = False
        if use_standard:
            coeff.coeffs = [style_old, N, r1, r3]
            
        # Save coeffs for the plot
        plot_coeffs['class2xe'] = coeff.coeffs
        plotting_coeffs[i] = plot_coeffs
            
    # Plot the results
    #basename = ''
    if basename and plot_xterm['bb13']:
        try:
            basename = basename + '_bb13'
            coeff_name = 'Bond/bond 13'
            plotting.bondbond(m, plotting_coeffs, basename, coeff_name)
        except:
            stack_trace_string = traceback.format_exc()
            log.out(f'ERROR {stack_trace_string}')
    return m