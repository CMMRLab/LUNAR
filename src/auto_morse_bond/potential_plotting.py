# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.6
November 13th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Import Necessary Libraries  
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0}) # Will suppress warnings about number of figures
from matplotlib.backends.backend_pdf import PdfPages

# Function for stringing together float values for parameters
def string_parameters(coeff):
    string = ''
    for i in coeff:
        if isinstance(i, float):
            string += '{:>16.8f}'.format(i)
        elif isinstance(i, int):
            string += '{:>6}'.format(i)
        else:
            string += '{:>6}'.format(i)
    return string

# https://docs.lammps.org/bond_class2.html
def figure(m, radius, basename, files2write, version, bondbreak_scale, ff_class, bond2style, class2xe_update):
    ############################################################################################
    ### Plotting multiple figures to .pdf                                                    ###
    # https://www.delftstack.com/howto/matplotlib/how-to-save-plots-as-pdf-file-in-matplotlib/ #
    ############################################################################################
    
    # Finding title
    pdf_title = basename + '.pdf' 
    fbb_title = 'in.' + basename + '.script'
    
    def retFig(m, x1, y1, y2, bondbreak_scale, i, morse_harmonic, ylo, yhi, morse_flag, ff_class, rcut, d0):
        name = 'Harmonic vs Morse Bond fit for bond coeff ' + str(i)
        
        # Plot morse and harmonic curves
        fig = plt.figure()
        plt.plot(x1, y1, color='tab:green', lw=2.5)
        if morse_flag:
            plt.plot(x1, y2, color='tab:cyan', lw=2.5)
            if rcut > 0:
                plt.axvline(x=rcut, color='tab:red', ls = ':', lw=2.0)
                plt.axhline(y=d0, color='black', ls = ':', lw=2.0)
        
        plt.xlim((0, max(radius)+0.1))
        plt.ylim((ylo, yhi))        
        
        plt.xlabel('Bond Radius (A)')
        plt.ylabel('Energy (Kcal/mol)')
        
        plt.suptitle(name)
        plt.title('{} (element,ring,nb)'.format(m.bond_info.types[i]))
        if morse_flag and ff_class in [1, '1']: plt.legend(['Harmonic Bond', 'Morse Bond', 'Bond Break scale: '+str(bondbreak_scale), 'Dissociation'])
        elif morse_flag and ff_class in [2, '2']: plt.legend(['Quartic Bond', 'Morse Bond', 'Bond Break scale: '+str(bondbreak_scale), 'Dissociation'])
        else: 
            if ff_class in [1, '1']: plt.legend(['Harmonic Bond'])
            if ff_class in [2, '2']: plt.legend(['Quartic Bond'])
        return fig
    
    
    # For loop plotting and saving figure to new page in .pdf file
    if files2write['write_pdffile']:
        pdf = PdfPages(pdf_title)
    
        # Loop through bond coeffs and plot/save info for writing LMP script
        break_points = {} # { coeffID : zero break point (if class2 = False)}
        bond_break = {} # { coeffID : bondbreak_scale*r0 (if class2 = False)}
        for i in m.alpha_parameter.e_harmonic:
            # Plot figure
            morse_flag = False; rcut = 0; d0 = 0;
            break_points[i] = False; bond_break[i] = False; 
            if ff_class in [1, '1']: 
                k, r0 = m.bond_coeffs[i].harmonic
            if ff_class in [2, '2']:
                r0, k2, k3, k4 = m.bond_coeffs[i].harmonic 
            ylo = -( abs(min( m.alpha_parameter.e_harmonic[i] + m.alpha_parameter.e_morse[i])) + 10 )   
            yhi = 300
            if bond2style[i] == 'morse': 
                morse_flag = True
                bond_break[i] = m.alpha_parameter.rcut[i]
                d0 = m.alpha_parameter.e_morse[i][-1]
                yhi = 2.5*d0
            
            # Create and save figure if desired
            if files2write['write_pdffile']:
                fig_i = retFig(m, radius, m.alpha_parameter.e_harmonic[i], m.alpha_parameter.e_morse[i], bondbreak_scale,
                               i, m.alpha_parameter.morse_harmonic, ylo, yhi, morse_flag, ff_class, rcut, d0)
                pdf.savefig(fig_i)
            
        # Close pdf
        pdf.close()
    
    ###############################
    # write fix bond/break script #
    ###############################
    # Function to create chunks and create print string
    def divide_chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i + n]
    if files2write['write_bondbreak']:
        with open(fbb_title, 'w') as f:
            # Write header
            f.write('# LAMMPS bond break settings for using fix bond/break with the Rmax raduis set at the dissociation point and\n')
            f.write(f'# based on bondbreak_scale (use which one is desired). Generated by: auto_morse_bond {version}\n')
            
            f.write('\n# NOTE that any quartic bond by default has is not allowed to break thus\n')
            f.write('# only bonds that were updated to morse will have a bond/break fix ID. Fix ID:\n')
            f.write('#     fix      bbBondID  # bb=bond/break; BondID=Bond CoeffID of morse bond;\n')
            
            # Write bond/break section
            f.write('\n\n#------------------------------------------------bond/break settings------------------------------------------------\n')
            f.write('variable      nsteps equal  1      # attempt breaking bonds every timestep (integer)\n')
            f.write('variable      pbreak equal  1.0    # probabilty to break bonds (float)\n')
            f.write('variable      sbreak equal  49829  # seed for bond/break (integer)\n')
            
            # # Write bond/break section for Dissociation
            # fixID = 'fix           bb'; groupID = 'all'; bbIDs = [];
            # f.write('\n\n#------------------------------------------------bond/break Dissociation----------------------------------------------\n')
            # for i in break_points:
            #     ID = '{}{}'.format(fixID, i); r = break_points[i];
            #     if r: 
            #         bbIDs.append(i); comment = '{:^3} {:<10}'.format('#', bond_info.types[i])
            #         f.write('{:<20} {:^4} {} {} {:^4} {:^4.3f} {} {} {} {}\n'.format(ID, groupID, 'bond/break', '${nsteps}', i, r, 'prob', '${pbreak}', '${sbreak}', comment))
             
            # Write bond/break section for bondbreak_scale
            fixID = 'fix           bb'; groupID = 'all'; bbIDs = [];
            f.write(f'\n\n#------------------------------------------------bond/break bondbreak_scale={bondbreak_scale}----------------------------------------\n')
            for i in bond_break:
                ID = '{}{}'.format(fixID, i); r = bond_break[i];
                if r: 
                    bbIDs.append(i); comment = '{:^3} {:<10}'.format('#', m.bond_info.types[i])
                    f.write('{:<20} {:^4} {} {} {:^4} {:^4.3f} {} {} {} {}\n'.format(ID, groupID, 'bond/break', '${nsteps}', i, r, 'prob', '${pbreak}', '${sbreak}', comment))
                    
            # Write possible thermo settings Option: 1
            thermo = 'thermo_style      custom' # 18 spaces
            f.write('\n\n#------------------------------------------------Thermo Settings Opt:1-----------------------------------------------\n')
            chunks = list(divide_chunks(bbIDs, 9))
            for n, chunk in enumerate(chunks):
                chunk = ['f_bb'+str(i)+'[1] ' for i in chunk]
                if n == 0: lst = [thermo] + chunk 
                else: lst = ['                 '] + chunk
                if n+1 < len(chunks): lst = lst + ['&']
                f.write('#{}\n'.format(' '.join(lst)))
                
            # Write possible thermo settings Option: 2
            thermo = 'thermo_style      custom' # 18 spaces
            f.write('\n\n#------------------------------------------------Thermo Settings Opt:2-----------------------------------------------\n')
            chunks = list(divide_chunks(bbIDs, 9))
            for n, chunk in enumerate(chunks):
                chunk = ['f_bb'+str(i)+'[2] ' for i in chunk]
                if n == 0: lst = [thermo] + chunk 
                else: lst = ['                 '] + chunk
                if n+1 < len(chunks): lst = lst + ['&']
                f.write('#{}\n'.format(' '.join(lst)))
                
            # Write possible thermo settings Option: 2
            thermo = 'thermo_style      custom' # 18 spaces
            f.write('\n\n#------------------------------------------------Thermo Settings Lists-----------------------------------------------\n')
            lst1 = ['f_bb'+str(i)+'[1] ' for i in bbIDs]
            lst2 = ['f_bb'+str(i)+'[2] ' for i in bbIDs]
            f.write('#{}\n'.format(' '.join(lst1)))
            f.write('#{}\n'.format(' '.join(lst2)))
            
            # This will control whether or not to comment out the bond_coeff, angle_coeff, .... sections
            comment_character = '#' # '' = will leave no comment and '#' will create comment
            
            # This will control whether or not to include angle_dihedral_coeffs
            include_angle_dihedral_coeffs = False # True or False
            
            # Write bond coeffs
            f.write('\n\n#------------------------------------------------Bond Coeffs to include----------------------------------------------\n')
            if m.bond_coeffs_style_hint != 'N/A':
                f.write(f'# Bond Coeffs  # {m.bond_coeffs_style_hint}\n')
            else:
                f.write('# Bond Coeffs\n')
            for i in m.bond_coeffs: 
                bond = m.bond_coeffs[i]
                if bond.type != 'N/A':
                    f.write('{}bond_coeff {:^3} {:<80} # {}\n'.format(comment_character, i, string_parameters(bond.coeffs), bond.type))
                else:
                    f.write('{}bond_coeff {:^3} {}\n'.format(comment_character, i, string_parameters(bond.coeffs)))
            
            # Write angle and bond coeffs if class2xe
            if class2xe_update and include_angle_dihedral_coeffs:
                # Write angle coeffs
                if m.angle_coeffs:
                    f.write('\n\n#------------------------------------------------Angle Coeffs to include----------------------------------------------\n')
                    #  Write angle coeffs
                    if m.angle_coeffs_style_hint != 'N/A':
                        f.write(f'# Angle Coeffs  # {m.angle_coeffs_style_hint}\n')
                    else:
                        f.write('# Angle Coeffs\n')
                    for i in m.angle_coeffs: 
                        angle = m.angle_coeffs[i]
                        if angle.type != 'N/A':
                            f.write('{}angle_coeff {:^3} {} # {}\n'.format(comment_character, i, string_parameters(angle.coeffs), angle.type))
                        else:
                            f.write('{}angle_coeff {:^3} {}\n'.format(comment_character, i, string_parameters(angle.coeffs)))
                    
                    #  Write Bondbond coeffs
                    if m.bondbond_coeffs:
                        if m.bondbond_coeffs_style_hint != 'N/A':
                            f.write(f'\n# BondBond Coeffs  # {m.bondbond_coeffs_style_hint}\n')
                        else:
                            f.write('\n# BondBond Coeffs\n')
                        for i in m.bondbond_coeffs: 
                            bondbond = m.bondbond_coeffs[i]
                            if bondbond.type != 'N/A':
                                f.write('{}angle_coeff {:^3} bb {} # {}\n'.format(comment_character, i, string_parameters(bondbond.coeffs), bondbond.type))
                            else:
                                f.write('{}angle_coeff {:^3} bb {}\n'.format(comment_character, i, string_parameters(bondbond.coeffs)))
                    
                    #  Write Bondangle paramerers
                    if m.bondangle_coeffs_style_hint != 'N/A':
                        f.write(f'\n# BondAngle Coeffs  # {m.bondangle_coeffs_style_hint}\n')
                    else:
                        f.write('\n# BondAngle Coeffs\n')
                    for i in m.bondangle_coeffs: 
                        bondangle = m.bondangle_coeffs[i]
                        if bondangle.type != 'N/A':
                            f.write('{}angle_coeff {:^3} ba {} # {}\n'.format(comment_character, i, string_parameters(bondangle.coeffs), bondangle.type))
                        else:
                            f.write('{}angle_coeff {:^3} ba {}\n'.format(comment_character, i, string_parameters(bondangle.coeffs)))
                                
                # Write dihedral coeffs
                if m.dihedral_coeffs:
                    f.write('\n\n#------------------------------------------------Dihedral Coeffs to include----------------------------------------------\n')
                    # Write dihedral coeffs
                    if m.dihedral_coeffs_style_hint != 'N/A':
                        f.write(f'# Dihedral Coeffs  # {m.dihedral_coeffs_style_hint}\n')
                    else:
                        f.write('# Dihedral Coeffs\n')
                    for i in m.dihedral_coeffs: 
                        dihedral = m.dihedral_coeffs[i]
                        if dihedral.type != 'N/A':
                            f.write('{}dihedral_coeff {:^3} {} # {}\n'.format(comment_character, i, string_parameters(dihedral.coeffs), dihedral.type))
                        else:
                            f.write('{}dihedral_coeff {:^3} {}\n'.format(comment_character, i, string_parameters(dihedral.coeffs)))
                            
                # Write angleangletorsion coeffs
                if m.angleangletorsion_coeffs:
                    if m.angleangletorsion_coeffs_style_hint != 'N/A':
                        f.write(f'\n# AngleAngleTorsion Coeffs  # {m.angleangletorsion_coeffs_style_hint}\n')
                    else:
                        f.write('\n# AngleAngleTorsion Coeffs\n')
                    for i in m.angleangletorsion_coeffs: 
                        angleangletorsion = m.angleangletorsion_coeffs[i]
                        if angleangletorsion.type != 'N/A':
                            f.write('{}dihedral_coeff {:^3} aat {} # {}\n'.format(comment_character, i, string_parameters(angleangletorsion.coeffs), angleangletorsion.type))
                        else:
                            f.write('{}dihedral_coeff {:^3} aat {}\n'.format(comment_character, i, string_parameters(angleangletorsion.coeffs)))
                            
                # Write endbondtorsion coeffs
                if m.endbondtorsion_coeffs:
                    if m.endbondtorsion_coeffs_style_hint != 'N/A':
                        f.write(f'\n# EndBondTorsion Coeffs  # {m.endbondtorsion_coeffs_style_hint}\n')
                    else:
                        f.write('\n# EndBondTorsion Coeffs\n')
                for i in m.endbondtorsion_coeffs: 
                    endbondtorsion = m.endbondtorsion_coeffs[i]
                    if endbondtorsion.type != 'N/A':
                        f.write('{}dihedral_coeff {:^3} ebt {} # {}\n'.format(comment_character, i, string_parameters(endbondtorsion.coeffs), endbondtorsion.type))
                    else:
                        f.write('{}dihedral_coeff {:^3} ebt {}\n'.format(comment_character, i, string_parameters(endbondtorsion.coeffs)))
                        
                # Write middlebondtorsion coeffs
                if m.middlebondtorsion_coeffs:
                    if m.middlebondtorsion_coeffs_style_hint != 'N/A':
                        f.write(f'\n# MiddleBondTorsion Coeffs  # {m.middlebondtorsion_coeffs_style_hint}\n')
                    else:
                        f.write('\n# MiddleBondTorsion Coeffs\n')
                    for i in m.middlebondtorsion_coeffs: 
                        middlebondtorsion  = m.middlebondtorsion_coeffs[i]
                        if middlebondtorsion.type != 'N/A':
                            f.write('{}dihedral_coeff {:^3} mbt {} # {}\n'.format(comment_character, i, string_parameters(middlebondtorsion.coeffs), middlebondtorsion.type))
                        else:
                            f.write('{}dihedral_coeff {:^3} mbt {}\n'.format(comment_character, i, string_parameters(middlebondtorsion.coeffs)))
                            
                # Write bondbond13 coeffs
                if m.bondbond13_coeffs:
                    if m.bondbond13_coeffs_style_hint != 'N/A':
                        f.write(f'\n# BondBond13 Coeffs  # {m.bondbond13_coeffs_style_hint}\n')
                    else:
                        f.write('\n# BondBond13 Coeffs\n')
                    for i in m.bondbond13_coeffs: 
                        bondbond13 = m.bondbond13_coeffs[i]
                        if bondbond13.type != 'N/A':
                            f.write('{}dihedral_coeff {:^3} bb13 {} # {}\n'.format(comment_character, i, string_parameters(bondbond13.coeffs), bondbond13.type))
                        else:
                            f.write('{}dihedral_coeff {:^3} bb13 {}\n'.format(comment_character, i, string_parameters(bondbond13.coeffs)))
                            
                # Write angletorsion coeffs
                if m.angletorsion_coeffs:
                    if m.angletorsion_coeffs_style_hint != 'N/A':
                        f.write(f'\n# AngleTorsion Coeffs  # {m.angletorsion_coeffs_style_hint}\n')
                    else:
                        f.write('\n# AngleTorsion Coeffs\n')
                    for i in m.angletorsion_coeffs: 
                        angletorsion = m.angletorsion_coeffs[i]
                        if angletorsion.type != 'N/A':
                            f.write('{}dihedral_coeff {:^3} at {} # {}\n'.format(comment_character, i, string_parameters(angletorsion.coeffs), angletorsion.type))
                        else:
                            f.write('{}dihedral_coeff {:^3} at {}\n'.format(comment_character, i, string_parameters(angletorsion.coeffs)))
    return 