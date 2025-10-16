# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
October 15, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0}) # Will suppress warnings about number of figures
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import time




# Set the colors for each curve
class2_color    = 'tab:olive'
class2xe_color  = 'tab:purple'

morse1_color    = 'tab:cyan'
morse2_color    = 'tab:blue'

harmonic1_color = 'tab:green'
harmonic2_color = 'darkgreen'




def morse_bond(d, a, r0, r):
    return d*( 1 - np.exp(-a*(r-r0) ) )**2

def class2_bond(r0, k2, k3, k4, r):
    return k2*(r - r0)**2 + k3*(r - r0)**3 + k4*(r - r0)**4

def get_digits(coeffs):
    return [i for i in coeffs if isinstance(i, (int, float))]




def bondbond(m, plotting_coeffs, basename, coeff_name):
    b1, b2 = '1', '2'
    if basename.endswith('13'): b2 = '3'
    
    pdf_title = basename + '.pdf' 
    pdf = PdfPages(pdf_title)
        
    # Plot settings
    legend_fs_scale = 0.75
    fs = 10
    lw = 2
    
    # bond radius:    start,  end,    inc  (Angstrom)
    radius = np.arange(-0.5,  10.1,   0.1)
    start_time = time.time()
    print('      * plotting')
    for i in plotting_coeffs:
        iter_time = (time.time() - start_time)
        print('          {} of {} ({:.4f} seconds)'.format(i, len(plotting_coeffs), iter_time))
        
        potentials = plotting_coeffs[i]
        potentials = {i:get_digits(j) for i, j in potentials.items()}
        if 'bond_types' not in potentials: continue
        type1, type2 = potentials['bond_types']
        bond1 = m.bonds_lengths.types[type1]
        bond2 = m.bonds_lengths.types[type2]
        #if i != 2: continue
        
        
        # Compute the potentials
        n, r1, r2 = potentials['class2']
        class2 = n*(radius - r1)*(radius - r2)
        
        d, a, r1, r2 = potentials['class2xe']
        class2xe = d*(( 1 - np.exp(-a*(radius - r1) ) )**1)*(( 1 - np.exp(-a*(radius - r2) ) )**1)
        
        morse1 = morse_bond(*potentials['morse1'], radius)
        morse2 = morse_bond(*potentials['morse2'], radius)
        harmonic1 = class2_bond(*potentials['harmonic1'], radius)
        harmonic2 = class2_bond(*potentials['harmonic2'], radius)
        d = max([potentials['morse1'][0], potentials['morse2'][0]])
        
        # Plot the results
        fig, ax1 = plt.subplots(figsize=(6, 4.5))

        ax1.plot(radius, harmonic1, ls='--', color=harmonic1_color, lw=lw, label='Harmonic (Bond:{})'.format(b1))
        ax1.plot(radius, harmonic2, ls='--', color=harmonic2_color, lw=lw, label='Harmonic (Bond:{})'.format(b2))
        
        ax1.plot(radius, morse1, color=morse1_color, lw=1.5*lw, label='Morse (Bond:{})'.format(b1))
        ax1.plot(radius, morse2, color=morse2_color, lw=1.5*lw, label='Morse (Bond:{})'.format(b2))
        
        ax1.plot(radius, class2, ls='--', color=class2_color, lw=lw, label='{} (classII)'.format(coeff_name))
        ax1.plot(radius, class2xe, color=class2xe_color, lw=1.5*lw, label='{} (classII-xe)'.format(coeff_name))
        
        ax2 = ax1.twinx()
        ax2.boxplot(bond1.lst, vert=False, patch_artist=True, positions=[round(-d/4)], widths=d/5, label='Bond:{} lengths'.format(b1),
                    labels=[''], boxprops=dict(facecolor=morse1_color, alpha=0.75))
        
        ax2.boxplot(bond2.lst, vert=False, patch_artist=True, positions=[round(-d/2)], widths=d/5, label='Bond:{} lengths'.format(b2),
                    labels=[''], boxprops=dict(facecolor=morse2_color, alpha=0.75))
        
        ax1.set_xlabel('Interatomic distance (Å)', fontsize=fs)
        ax1.set_ylabel('Energy ($kcal$ $mol^{-1}$)', fontsize=fs)
        
        ax1.tick_params(axis='x', which='major', labelsize=fs)
        ax1.tick_params(axis='y', which='major', labelsize=fs)
        ax1.set_xlim((min(radius), max(radius)))
        ax1.set_ylim((-d, 4*d))  
        ax2.set_xlim((min(radius), max(radius)))
        ax2.set_ylim((-d, 4*d))  
        ax2.set_yticks([])
        
        ax1.legend(loc='upper right', bbox_to_anchor=(1.0, 1.00), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs)
        ax2.legend(loc='upper right', bbox_to_anchor=(1.0, 0.675), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs)
        
        name = '{} coeff: {}'.format(coeff_name, i)
        type1 = '{} (element,ring,nb - Bond:{})'.format(m.bond_info.types[type1], b1)
        type2 = '{} (element,ring,nb - Bond:{})'.format(m.bond_info.types[type2], b2)
        plt.title('{}\n{}\n{}'.format(name, type1, type2))
        fig.tight_layout()
        pdf.savefig(fig)
        #plt.close('all')

    pdf.close()
    return




def bondangle(m, plotting_coeffs, basename):
    pdf_title = basename + '.pdf' 
    pdf = PdfPages(pdf_title)
    
    # Plot settings
    legend_fs_scale = 0.75
    fs = 9
    t = 0.75
    
    # bond radius:    start    end, num  (Angstrom)
    radius = np.linspace(0.5,  8,   30)
    
    # Find all theta ranges
    thetas = [plotting_coeffs[i]['angle'][0] for i in plotting_coeffs]
    start = int(round(min(thetas) - 25, -1))
    end   = int(round(max(thetas) + 25, -1))
    
    # theta angle:                 start,    end,   num   (Degrees)
    angle = np.deg2rad(np.linspace(start,    end,   15))
    
    # Make mesh grid
    X, Y = np.meshgrid(radius, angle, sparse=False)
    radius = np.ravel(X)
    theta = np.ravel(Y)
    start_time = time.time()
    print('      * plotting')
    for i in plotting_coeffs:
        iter_time = (time.time() - start_time)
        print('          {} of {} ({:.4f} seconds)'.format(i, len(plotting_coeffs), iter_time))
        
        potentials = plotting_coeffs[i]
        potentials = {i:get_digits(j) for i, j in potentials.items()}
        if 'bond_types' not in potentials: continue
        type1, type2 = potentials['bond_types']
        #if i != 2: continue
        
        # Get theta0
        theta0 = potentials['angle'][0]
        t0 = np.radians(theta0)
        
        # Compute the bonding potentials
        d_max = max([potentials['morse1'][0], potentials['morse2'][0]])
        morse1 = ( morse_bond(*potentials['morse1'], radius)*np.ones(len(theta)) ).reshape(X.shape)
        morse2 = ( morse_bond(*potentials['morse2'], radius)*np.ones(len(theta)) ).reshape(X.shape)
        harmonic1 = ( class2_bond(*potentials['harmonic1'], radius)*np.ones(len(theta)) ).reshape(X.shape)
        harmonic2 = ( class2_bond(*potentials['harmonic2'], radius)*np.ones(len(theta)) ).reshape(X.shape)
        harmonic1[harmonic1 > 5*d_max] = 5*d_max # Will set anything > 5*d_max to 5*d_max
        harmonic2[harmonic2 > 5*d_max] = 5*d_max # Will set anything > 5*d_max to 5*d_max
        
        # Compute the potentials
        n1, n2, r1, r2 = potentials['class2']
        class2 = n1*(radius - r1)*(theta - t0) + n2*(radius - r2)*(theta - t0)
        class2 = class2.reshape(X.shape)
        class2[class2 > 3*d_max] = 3*d_max # Will set anything > 3*d_max to 3*d_max
        
        d1, d2, a1, a2, r1, r2 = potentials['class2xe']
        class2xe = d1*( 1 - np.exp(-a1*(radius - r1)) )*(theta - t0) + d2*( 1 - np.exp(-a2*(radius - r2)) )*(theta - t0)
        class2xe = class2xe.reshape(X.shape)
        

        # Plot the results
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        deg = np.rad2deg(Y)
        
        ax.plot_surface(X, deg, harmonic1, color=harmonic1_color, label='Harmonic (Bond:1)')
        ax.plot_surface(X, deg, harmonic2, color=harmonic2_color, label='Harmonic (Bond:2)')
        ax.plot_surface(X, deg, class2, color=class2_color, alpha=t, label='Bond/angle (classII)')
        
        ax.plot_surface(X, deg, morse1, color=morse1_color , alpha=t, label='Morse (Bond:1)')
        ax.plot_surface(X, deg, morse2, color=morse2_color, alpha=t, label='Morse (Bond:2)')
        ax.plot_surface(X, deg, class2xe, color=class2xe_color, alpha=t, label='Bond/angle (classII-xe)')
        
        ax.legend(loc='upper right', bbox_to_anchor=(1.45, 1.00), fancybox=True, ncol=2, fontsize=legend_fs_scale*fs)
        name = 'Bond/angle coeff: {}'.format(i)
        type1 = '{} (element,ring,nb - Bond:1)'.format(m.bond_info.types[type1])
        type2 = '{} (element,ring,nb - Bond:1)'.format(m.bond_info.types[type2])
        plt.title('{}\n{}\n{}'.format(name, type1, type2), fontsize=fs)
        
        # Adjust axis and set label axis
        ax.set_zlim3d(-1*d_max, 2*d_max)
        ax.set_xlabel('\nInteratomic distance (Å)', fontsize=fs)
        ax.set_ylabel('\n\n\nAngle (°)', fontsize=fs)
        ax.set_zlabel('\n\n\n\nEnergy ($kcal$ $mol^{-1}$)', fontsize=fs)
        ax.tick_params(axis='x', labelsize=fs, pad=4)
        ax.tick_params(axis='y', labelsize=fs, pad=8)
        ax.tick_params(axis='z', labelsize=fs, pad=12)
        ax.view_init(elev=15, azim=-65, roll=0) # tweak viewing angle if you like

        fig.tight_layout()
        pdf.savefig(fig)
        #plt.close('all')

    pdf.close()
    return




def middlebondtorsion(m, plotting_coeffs, basename):
    pdf_title = basename + '.pdf' 
    pdf = PdfPages(pdf_title)
    
    # Plot settings
    legend_fs_scale = 0.75
    fs = 9
    t = 0.75
    
    # bond radius:    start    end, num  (Angstrom)
    radius = np.linspace(0.5,  8,   30)
    
    # Dihedral angle:            start,  end,   num   (Degrees)
    phis = np.deg2rad(np.linspace(-180,  180,   30))
    
    
    # Make mesh grid
    X, Y = np.meshgrid(radius, phis, sparse=False)
    radius = np.ravel(X)
    phi = np.ravel(Y)
    start_time = time.time()
    print('      * plotting')
    for i in plotting_coeffs:
        iter_time = (time.time() - start_time)
        print('          {} of {} ({:.4f} seconds)'.format(i, len(plotting_coeffs), iter_time))
        
        potentials = plotting_coeffs[i]
        potentials = {i:get_digits(j) for i, j in potentials.items()}
        if 'bond_type' not in potentials: continue
        type1 = potentials['bond_type'][0]
        #if i != 2: continue
    
        # Compute the bonding potentials
        d_max = potentials['morse'][0]
        morse = ( morse_bond(*potentials['morse'], radius)*np.ones(len(phi)) ).reshape(X.shape)
        harmonic = ( class2_bond(*potentials['harmonic'], radius)*np.ones(len(phi)) ).reshape(X.shape)
        harmonic[harmonic > 5*d_max] = 5*d_max # Will set anything > 5*d_max to 5*d_max
        
        # Compute the potentials
        a1, a2, a3, r2 = potentials['class2']
        class2 = (radius - r2)*( a1*np.cos(phi) + a2*np.cos(2*phi) + a3*np.cos(3*phi) )
        class2 = class2.reshape(X.shape)
        
        
        a, a1, a2, a3, r2 = potentials['class2xe']
        class2xe = (( 1 - np.exp(-a*(radius - r2) ) )**1)*( a1*np.cos(phi) + a2*np.cos(2*phi) + a3*np.cos(3*phi) )
        class2xe = class2xe.reshape(X.shape)
        
        # Plot the results
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        deg = np.rad2deg(Y)
        
        ax.plot_surface(X, deg, harmonic, color=harmonic1_color, label='Harmonic (Bond:2)')
        ax.plot_surface(X, deg, morse, color=morse1_color, alpha=t, label='Morse (Bond:2)')
        
        ax.plot_surface(X, deg, class2, color=class2_color, alpha=t, label='Middle bond/torsion (classII)')
        ax.plot_surface(X, deg, class2xe, color=class2xe_color, alpha=t, label='Middle bond/torsion (classII-xe)')
        
        ax.legend(loc='upper right', bbox_to_anchor=(1.50, 1.00), fancybox=True, ncol=2, fontsize=legend_fs_scale*fs)
        name = 'Middle bond/torsion coeff: {}'.format(i)
        type1 = '{} (element,ring,nb - Bond:2)'.format(m.bond_info.types[type1])
        plt.title('{}\n{}'.format(name, type1), fontsize=fs)
        
        # Adjust axis and set label axis
        ax.set_zlim3d(-1*d_max, 2*d_max)
        ax.set_xlabel('Interatomic distance (Å)', fontsize=fs)
        ax.set_ylabel('\nTorsion angle (°)', fontsize=fs)
        ax.set_zlabel('\nEnergy ($kcal$ $mol^{-1}$)', fontsize=fs)
        ax.tick_params(axis='x', labelsize=fs, pad=2)
        ax.tick_params(axis='y', labelsize=fs, pad=4)
        ax.tick_params(axis='z', labelsize=fs, pad=4)
        ax.view_init(elev=15, azim=-65, roll=0) # tweak viewing angle if you like

        fig.tight_layout()
        pdf.savefig(fig)
        #plt.close('all')

    pdf.close()
    return




def endbondtorsion(m, plotting_coeffs, basename):
    pdf_title = basename + '.pdf' 
    pdf = PdfPages(pdf_title)
    
    # Plot settings
    legend_fs_scale = 0.75
    fs = 9
    t = 0.75
    
    # bond radius:    start    end, num  (Angstrom)
    radius = np.linspace(0.5,  8,   30)
    
    # Dihedral angle:            start,  end,   num   (Degrees)
    phis = np.deg2rad(np.linspace(-180,  180,   30))
    
    
    # Make mesh grid
    X, Y = np.meshgrid(radius, phis, sparse=False)
    radius = np.ravel(X)
    phi = np.ravel(Y)
    start_time = time.time()
    print('      * plotting')
    for i in plotting_coeffs:
        iter_time = (time.time() - start_time)
        print('          {} of {} ({:.4f} seconds)'.format(i, len(plotting_coeffs), iter_time))
        
        potentials = plotting_coeffs[i]
        potentials = {i:get_digits(j) for i, j in potentials.items()}
        if 'bond_types' not in potentials: continue
        type1, type2 = potentials['bond_types']
        #if i != 2: continue
    
        # Compute the bonding potentials
        d_max = max([potentials['morse1'][0], potentials['morse2'][0]])
        morse1 = ( morse_bond(*potentials['morse1'], radius)*np.ones(len(phi)) ).reshape(X.shape)
        morse2 = ( morse_bond(*potentials['morse2'], radius)*np.ones(len(phi)) ).reshape(X.shape)
        harmonic1 = ( class2_bond(*potentials['harmonic1'], radius)*np.ones(len(phi)) ).reshape(X.shape)
        harmonic2 = ( class2_bond(*potentials['harmonic2'], radius)*np.ones(len(phi)) ).reshape(X.shape)
        harmonic1[harmonic1 > 5*d_max] = 5*d_max # Will set anything > 5*d_max to 5*d_max
        harmonic2[harmonic2 > 5*d_max] = 5*d_max # Will set anything > 5*d_max to 5*d_max
        
        # Compute the potentials
        b1, b2, b3, c1, c2, c3, r1, r3 = potentials['class2']
        B = (radius - r1)*( b1*np.cos(phi) + b2*np.cos(2*phi) + b3*np.cos(3*phi) )
        C = (radius - r3)*( c1*np.cos(phi) + c2*np.cos(2*phi) + c3*np.cos(3*phi) )
        class2 = B + C
        class2 = class2.reshape(X.shape)
        
        
        a1, a3, b1, b2, b3, c1, c2, c3, r1, r3 = potentials['class2xe']
        B = (( 1 - np.exp(-a1*(radius - r1) ) )**1)*( b1*np.cos(phi) + b2*np.cos(2*phi) + b3*np.cos(3*phi) )
        C = (( 1 - np.exp(-a3*(radius - r3) ) )**1)*( c1*np.cos(phi) + c2*np.cos(2*phi) + c3*np.cos(3*phi) )
        class2xe = B + C
        class2xe = class2xe.reshape(X.shape)
        
        # Plot the results
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        deg = np.rad2deg(Y)
        
        ax.plot_surface(X, deg, harmonic1, color=harmonic1_color, label='Harmonic (Bond:1)')
        ax.plot_surface(X, deg, harmonic2, color=harmonic2_color, label='Harmonic (Bond:3)')
        
        ax.plot_surface(X, deg, morse1, color=morse1_color, alpha=t, label='Morse (Bond:1)')
        ax.plot_surface(X, deg, morse2, color=morse2_color, alpha=t, label='Morse (Bond:3)')
        
        ax.plot_surface(X, deg, class2, color=class2_color, alpha=t, label='End bond/torsion (classII)')
        ax.plot_surface(X, deg, class2xe, color=class2xe_color, alpha=t, label='End bond/torsion (classII-xe)')
        
        ax.legend(loc='upper right', bbox_to_anchor=(1.50, 1.00), fancybox=True, ncol=2, fontsize=legend_fs_scale*fs)
        name = 'End bond/torsion coeff: {}'.format(i)
        type1 = '{} (element,ring,nb - Bond:1)'.format(m.bond_info.types[type1])
        type2 = '{} (element,ring,nb - Bond:1)'.format(m.bond_info.types[type2])
        plt.title('{}\n{}\n{}'.format(name, type1, type2), fontsize=fs)
        
        # Adjust axis and set label axis
        ax.set_zlim3d(-1*d_max, 2*d_max)
        ax.set_xlabel('Interatomic distance (Å)', fontsize=fs)
        ax.set_ylabel('\nTorsion angle (°)', fontsize=fs)
        ax.set_zlabel('\nEnergy ($kcal$ $mol^{-1}$)', fontsize=fs)
        ax.tick_params(axis='x', labelsize=fs, pad=2)
        ax.tick_params(axis='y', labelsize=fs, pad=4)
        ax.tick_params(axis='z', labelsize=fs, pad=4)
        ax.view_init(elev=15, azim=-65, roll=0) # tweak viewing angle if you like

        fig.tight_layout()
        pdf.savefig(fig)
        #plt.close('all')

    pdf.close()
    return
    