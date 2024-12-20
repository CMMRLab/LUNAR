# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.4
July 14th, 2023
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
def figure(m, radius, bond_info, basename, files2write, version, bondbreak_scale, ff_class):
    ############################################################################################
    ### Plotting multiple figures to .pdf                                                    ###
    # https://www.delftstack.com/howto/matplotlib/how-to-save-plots-as-pdf-file-in-matplotlib/ #
    ############################################################################################
    
    # Finding title
    pdf_title = basename + '.pdf' 
    fbb_title = 'in.' + basename + '.script'
    
    def retFig(x1, y1, y2, y3, zero_point, rmax, bondbreak_scale, i, morse_harmonic, ylo, morse_flag, ff_class):
        name = 'Harmonic vs Morse Bond fit for bond coeff ' + str(i)
        if len(morse_harmonic[i]) == 4:
            yhi = 2.5*morse_harmonic[i][1] # Set yhi based on D value at index 1
        else: yhi = 300
        
        # Plot morse and harmonic curves
        fig = plt.figure()
        plt.plot(x1, y1, color='tab:green', lw=2.5)
        if morse_flag:
            plt.plot(x1, y2, color='tab:cyan', lw=2.5)
            plt.plot(x1, y3, color='tab:blue', lw=2.5)
            plt.axvline(x=zero_point, color='tab:red', ls = ':', lw=2.0)
            plt.axvline(x=rmax, color='tab:orange', ls = ':', lw=2.0)
            plt.axhline(y=0.0, color ='black', ls = ':', lw=2.5)
        
        plt.xlim((0, max(radius)+0.1))
        plt.ylim((ylo, yhi))        
        
        plt.xlabel('Bond Radius (A)')
        plt.ylabel('Energy (Kcal/mol)')
        
        plt.suptitle(name)
        plt.title('{} (element,ring,nb)'.format(bond_info.types[i]))
        if morse_flag and ff_class == 1: plt.legend(['Harmonic Bond', 'Morse Bond', 'Shifted Morse bond', 'Dissociation point', 'Bond Break scale: '+str(bondbreak_scale)])
        elif morse_flag and ff_class == 2: plt.legend(['Quartic Bond', 'Morse Bond', 'Shifted Morse bond', 'Dissociation point', 'Bond Break scale: '+str(bondbreak_scale)])
        else: 
            if ff_class == 1: plt.legend(['Harmonic Bond'])
            if ff_class == 2: plt.legend(['Quartic Bond'])
        return fig
    
    
    # For loop plotting and saving figure to new page in .pdf file
    if files2write['write_pdffile']:
        pdf = PdfPages(pdf_title)
    
    # Loop through bond coeffs and plot/save info for writing LMP script
    break_points = {} # { coeffID : zero break point (if class2 = False)}
    bond_break = {} # { coeffID : bondbreak_scale*r0 (if class2 = False)}
    for i in m.alpha_parameter.e_harmonic:
        # Plot figure
        ylo = -(abs(min(m.alpha_parameter.e_harmonic[i])) + 10)
        bond = m.alpha_parameter.morse_harmonic[i]; zero_point = 0; rmax = 0;
        shifted_morse = m.alpha_parameter.e_harmonic[i]; morse_flag = False; 
        break_points[i] = False; bond_break[i] = False; 
        if ff_class == 1:
            if len(m.bond_coeffs[i].coeffs) == 2:
                k, r0 = m.bond_coeffs[i].coeffs
            elif len(m.bond_coeffs[i].coeffs) == 3:
                style, k, r0 = m.bond_coeffs[i].coeffs
            else:
                r0 = 1.4
        if ff_class == 2:
            if len(m.bond_coeffs[i].coeffs) == 4:
                r0, k2, k3, k4 = m.bond_coeffs[i].coeffs 
            if len(m.bond_coeffs[i].coeffs) == 5:
                style, r0, k2, k3, k4 = m.bond_coeffs[i].coeffs 
            else:
                r0 = 1.4
        if 'morse' in bond:
            morse_flag = True; zero_point = m.alpha_parameter.dpoint[i]
            shifted_morse = [i-bond[1] for i in m.alpha_parameter.e_morse[i]]; 
            break_points[i] = zero_point; rmax = bondbreak_scale*r0;
            bond_break[i] = rmax;
            ylo = -( abs(min( shifted_morse + m.alpha_parameter.e_harmonic[i] + m.alpha_parameter.e_morse[i])) + 10 )
            
        fig_i = retFig(radius, m.alpha_parameter.e_harmonic[i], m.alpha_parameter.e_morse[i], shifted_morse, zero_point, rmax,
                       bondbreak_scale, i, m.alpha_parameter.morse_harmonic, ylo, morse_flag, ff_class)
        
        # Save figure if desired
        if files2write['write_pdffile']:
            pdf.savefig(fig_i)
    
    if files2write['write_pdffile']:
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
                    bbIDs.append(i); comment = '{:^3} {:<10}'.format('#', bond_info.types[i])
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
            
            # Write bond coeffs
            f.write('\n\n#------------------------------------------------Bond Coeffs to include----------------------------------------------\n')
            for i in m.alpha_parameter.morse_harmonic: 
                bond = m.alpha_parameter.morse_harmonic[i]               
                f.write('#bond_coeff {:^3} {:^6} {}\n'.format(i, bond[0], string_parameters(bond[1:])))
    return 