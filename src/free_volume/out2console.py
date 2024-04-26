# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
December 6th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
# Function to write volume vs postions
def write_csv(basename, posVSvol, direction, free_volume):
    csvname = '{}_spatial_distribution_direction_{}.csv'.format(basename, direction)
    with open(csvname, 'w') as f:
        
        # Find titles and data to write
        titles = '{}, {}, {}'.format('Postion', 'Free Volume (A^3)', '% Free Volume')
        f.write('{}\n'.format(titles))
    
        # Write data 
        posVSvol = dict(sorted(posVSvol.items()))
        for i in posVSvol:
            try: data = '{}, {}, {}'.format(i, posVSvol[i], 100*posVSvol[i]/free_volume)
            except: data = '{}, {}, {}'.format(i, 0, 0)
            f.write('{}\n'.format(data))
    return


# Function to write log file
def out(m, v, grid, execution_time, boundary, max_voxel_size, compute_free_volume_distributions, basename, version, run_mode, files2write, probe_diameter, vdw_method, log):
    # Write csv files
    if files2write['write_spat_dis-x']: write_csv(basename, grid.spatial_distributions['x'], 'x', grid.free_volume)
    if files2write['write_spat_dis-y']: write_csv(basename, grid.spatial_distributions['y'], 'y', grid.free_volume)
    if files2write['write_spat_dis-z']: write_csv(basename, grid.spatial_distributions['z'], 'z', grid.free_volume)
    
    # Compute density
    system_mass_sum = 0; amu2grams = 1/6.02214076e+23;
    for i in m.atoms:
        atom = m.atoms[i]
        mass = m.masses[atom.type]
        system_mass_sum += mass.coeffs[0]
        
    # convert system mass in amu to grams
    mass = system_mass_sum*amu2grams
    
    # Find box dimensions to compute volume
    angstromcubed2cmcubed = 1e-24
    xline = m.xbox_line.split(); yline = m.ybox_line.split(); zline = m.zbox_line.split();
    lx = float(xline[1])-float(xline[0])
    ly = float(yline[1])-float(yline[0])
    lz = float(zline[1])-float(zline[0])
    volume = lx*ly*lz*angstromcubed2cmcubed
    density = mass/volume
        
    # Write header info and file purpose 
    log.out(f'This is a saved .txt file for the print outs that appear when running free_volume.py {version}')
    log.out(f'Running in {run_mode} mode\n')
    
    # Write the elements found in the system
    log.out('\n\nElements found in system:')
    for i in sorted(m.elements):
        log.out('{} {}'.format('-', i))
    
    # density and box results to screen
    log.out('\n\nSystem Unit cell, volume, mass, density')
    log.out('{:<14} {:^10.5f} {:<10}'.format('Lx:', lx, 'angstrom'))
    log.out('{:<14} {:^10.5f} {:<10}'.format('Ly:', ly, 'angstrom'))
    log.out('{:<14} {:^10.5f} {:<10}'.format('Lz:', lz, 'angstrom'))
    log.out('{:<14} {:^10.4E} {:<10}'.format('volume:', volume, 'cm^3'))
    log.out('{:<14} {:^10.4E} {:<10}'.format('mass:', mass, 'grams'))
    log.out('{:<14} {:^10.5f} {:<10}' .format('density:', density, 'g/cm^3'))
    
    # Write domain decomposition info
    if 'dd' in run_mode:
        log.out('\n\nDomain decomposition information:')
        log.out(f'max_domain_size = {v.domain_size}')
        log.out(f'Number of x-domains, x-domain dimension = {v.nxx}, {lx/v.nxx}')
        log.out(f'Number of y-domains, y-domain dimension = {v.nyy}, {ly/v.nyy}')
        log.out(f'Number of z-domains, z-domain dimension = {v.nzz}, {lz/v.nzz}')
        
    
    # Write voxel info
    log.out('\n\nVoxel information:')
    log.out(f'max_voxel_size = {max_voxel_size}')
    log.out(f'Number of x-voxels, x-voxel dimension = {v.nx}, {v.dx}')
    log.out(f'Number of y-voxels, y-voxel dimension = {v.ny}, {v.dy}')
    log.out(f'Number of z-voxels, z-voxel dimension = {v.nz}, {v.dz}')
    log.out(f'Number of voxels created {len(v.voxelIDs)}')
    log.out(f'probe_diameter = {probe_diameter}')
    
    # Write PBC info
    log.out('\n\nPBC information:')
    log.out(f'System periodic boundary flags                        : {boundary}')
    log.out(f'Count of atoms that are near cell edge (possibly PBC) : {grid.npossible_pbc}')
    log.out(f'Total Count of periodically determined occupied voxels: {grid.pbc_count}')
    
    # Write results
    log.out('\n\nFree volume analysis results:')
    #log.out(f'PBC flags             : {boundary}')
    log.out(f'vdw radii method      : {vdw_method}')
    log.out(f'Simulation cell volume: {grid.simulation_volume} A^3')
    log.out(f'Atom volume           : {grid.atom_volume} A^3')
    log.out(f'Free volume           : {grid.free_volume} A^3')
    log.out(f'Percent Free volume   : {grid.pfree_volume} %')
    
    # Write FFV (Fractional Free Volume results)
    if probe_diameter == 0 or probe_diameter in [v.nx, v.ny, v.nz]:
        if mass > 0:
            Vsp = grid.simulation_volume/mass
            Voc = 1.3*grid.atom_volume/mass
        else: Vsp = 0; Voc = 0
        if Vsp > 0: ffv = (Vsp - Voc)/Vsp
        else: ffv = 0
        log.out('\n\nFractional Free Volume for polymers (FFV = Vf/Vsp; where Vf = Vsp - Voc & Voc = 1.3*Vw) results:')
        log.out(f'Van der Waals Volume (Vw)        : {grid.atom_volume} cm^3')
        log.out(f'Occupied Volume (Voc)            : {Voc} cm^3/g')
        log.out(f'Specfic Volume (Vsp)             : {Vsp} cm^3/g')
        log.out(f'Fractional Free volume (100*FFV) : {100*ffv} (%)\n')
        log.out('  If you are using this metric it should be for polymers ONLY and FFV should be positive. If FFV is negative')
        log.out('  you should use the "Percent Free Volume" instead. Addtionally you should cite the following:\n')
        
        log.out('    Yampolskii, Y. (2016). Fractional Free Volume (FFV). In: Drioli, E., Giorno, L. (eds) Encyclopedia of')
        log.out('    Membranes. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-662-44324-8_243\n')
        
        log.out('    Wenqing Zhang, Yang Qing, Weihong Zhong, Gang Sui, Xiaoping Yang (2017). In: Reactive and Functional Polymers,')
        log.out('    Vol 111, pages 60-67, Mechanism of modulus improvement for epoxy resin matrices: A molecular dynamics simulation')
        log.out('    https://doi.org/10.1016/j.reactfunctpolym.2016.12.014\n')
        
        log.out('  It is also recommend to read up on the FFV definiton, its limitations/critiques. A link is provided to get')
        log.out('  started in your search (https://link.springer.com/referenceworkentry/10.1007/978-3-662-44324-8_243#citeas).')
        log.out('  Typical FFV for polymers is in the range of 10-25%, however some polymers such as polytrimethylsilyl propyne')
        log.out('  may have FFV as large as 35%.')

    
    # Write free volume table
    if compute_free_volume_distributions:
        log.out('\n\n')
        log.out('--------------------------------------Free Volume Distribution--------------------------------------')
        log.out('{:<10} {:^12} {:^12} {:^12} {:^12} {:^12} {:^12} {:^12}'.format('clusterID', 'size', '%size', 'volume', '%volume', 'X-span', 'Y-span', 'Z-span'))
        log.out('----------------------------------------------------------------------------------------------------')  
        for ID in grid.free_volume_clusters:
            cluster = grid.free_volume_clusters[ID]
            log.out('{:<10} {:^12} {:^12.2f} {:^12.2f} {:^12.2f} {:^12.2f} {:^12.2f} {:^12.2f}'.format(ID, cluster.size, cluster.psize, cluster.volume, cluster.pvolume, cluster.xspan, cluster.yspan, cluster.zspan))
        #log.out('{:<10} {:^12} {:^12} {:^12} {:^12} {:^12} {:^12} {:^12}'.format(ID+1, '', 'lone', 'free', 'volume', 'voxel', 'clusterID', ''))
        
    # Write excutation time
    log.out(f'\n\nExecution time in seconds:    {execution_time}') 
    return
