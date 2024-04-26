# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
December 30th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

full: atom-ID molecule-ID atom-type q x y z
"""
##############################
# Import Necessary Libraries #
##############################
from tqdm import tqdm


###############
# voxels only #
###############
def voxels_only(m, v, basename, version, log):
    filename = basename+'_voxels_only.data'
    log.out(f'\nWriting {filename} ...')
    # Writing new file with bonds information
    with open(filename, 'w') as f: 
        f.write(f'HEADER, free_volume.py {version}\n\n')
        f.write(f'{len(v.voxelIDs)} atoms\n\n')
        f.write('1 atom types\n\n')
        f.write(f'{m.xbox_line}\n{m.ybox_line}\n{m.zbox_line}\n')
        
        # Write masses
        f.write('\nMasses\n\n')
        f.write('{} {:^10.5f} # {}\n'.format(1, 1.0, 'Voxels'))
            
        # Write atoms  
        f.write('\nAtoms # full\n\n')
        data2write = [];
        for ID in tqdm(v.voxelIDs):
            p = v.voxelIDs[ID]
            #data2write.append('{:^3} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2}\n'.format(ID, 1, 1, 0.0, p[0], p[1], p[2], 0, 0, 0))
            data2write.append('{} {} {} {} {} {} {} {} {} {}\n'.format(ID, 1, 1, 0.0, p[0], p[1], p[2], 0, 0, 0))
        string2write = ''.join(data2write)
        f.write(string2write)
    return

##############
# atoms only #
##############
def atoms_only(m, grid, basename, version, log):
    filename = basename+'_atoms_only.data'
    log.out(f'\nWriting {filename} ...')
    # Writing new file with bonds information
    with open(filename, 'w') as f: 
        f.write(f'HEADER, free_volume.py {version}\n\n')
        f.write(f'{len(grid.voxel_atomIDs)} atoms\n\n')
        f.write(f'{len(grid.atomTypeIDs)} atom types\n\n')
        f.write(f'{m.xbox_line}\n{m.ybox_line}\n{m.zbox_line}\n')
        
        # Write masses
        f.write('\nMasses\n\n')
        for i in grid.atomTypeIDs:
            info = grid.atomTypeIDs[i]
            f.write('{} {:^10.5f}  # {}\n'.format(i, info[1], info[0]))
            
        # Write atoms    
        f.write('\nAtoms # full\n\n')
        data2write = [];
        for ID in tqdm(grid.voxel_atomIDs):
            p = grid.voxel_atomIDs[ID]
            #data2write.append('{:^3} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2}\n'.format(ID, 1, p[3], 0.0, p[0], p[1], p[2], 0, 0, 0)) 
            data2write.append('{} {} {} {} {} {} {} {} {} {}\n'.format(ID, 1, int(p[3]), 0.0, p[0], p[1], p[2], 0, 0, 0)) 
        string2write = ''.join(data2write)
        f.write(string2write)
    return

##############
# free only #
##############
def free_only(m, grid, compute_free_volume_distributions, basename, version, log):
    filename = basename+'_free_only.data'
    log.out(f'\nWriting {filename} ...')
    # Writing new file with bonds information
    with open(filename, 'w') as f: 
        f.write(f'HEADER, free_volume.py {version}\n\n')
        f.write(f'{len(grid.voxel_freeIDs)} atoms\n\n')
        f.write('1 atom types\n\n')
        f.write(f'{m.xbox_line}\n{m.ybox_line}\n{m.zbox_line}\n')
        
        # Write masses
        f.write('\nMasses\n\n')
        f.write('{} {:^10.5f}  # {}\n'.format(1, 1.0, 'Free'))
            
        # Write atoms    
        f.write('\nAtoms # full\n\n')
        data2write = [];
        for ID in tqdm(grid.voxel_freeIDs):
            p = grid.voxel_freeIDs[ID]
            if compute_free_volume_distributions:
                molID = grid.voxelclusterID2molID[ID]
            else:
                molID = 1
            #data2write.append('{:^3} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2}\n'.format(ID, molID, 1, 0.0, p[0], p[1], p[2], 0, 0, 0)) 
            data2write.append('{} {} {} {} {} {} {} {} {} {}\n'.format(ID, molID, 1, 0.0, p[0], p[1], p[2], 0, 0, 0)) 
                
        string2write = ''.join(data2write)
        f.write(string2write)
    return

##################
# free and atoms #
##################
def atoms_free(m, grid, compute_free_volume_distributions, basename, version, log):
    filename = basename+'_atoms_free.data'
    log.out(f'\nWriting {filename} ...')
    # Writing new file with bonds information
    with open(filename, 'w') as f: 
        f.write(f'HEADER, free_volume.py {version}\n\n')
        f.write(f'{len(grid.voxel_atomIDs)+len(grid.voxel_freeIDs)} atoms\n\n')
        f.write('2 atom types\n\n')
        f.write(f'{m.xbox_line}\n{m.ybox_line}\n{m.zbox_line}\n')
        
        # Write masses
        f.write('\nMasses\n\n')
        f.write('{} {:^10.5f}  # {}\n'.format(1, 1.0, 'Atoms'))
        f.write('{} {:^10.5f}  # {}\n'.format(2, 1.0, 'Free'))
            
        # Write atoms    
        f.write('\nAtoms # full\n\n')
        molids = set(); data2write = [];
        for ID in tqdm(grid.voxel_atomIDs):
            p = grid.voxel_atomIDs[ID]
            if compute_free_volume_distributions:
                molID = len(grid.free_volume_clusters)+1
            else:
                molID = 1
            molids.add(molID)
            #data2write.append('{:^3} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2}\n'.format(ID, molID, 1, 0.0, p[0], p[1], p[2], 0, 0, 0)) 
            data2write.append('{} {} {} {} {} {} {} {} {} {}\n'.format(ID, molID, 1, 0.0, p[0], p[1], p[2], 0, 0, 0)) 
            
        molID_offset = max(molids)
        for ID in tqdm(grid.voxel_freeIDs):
            p = grid.voxel_freeIDs[ID]
            if compute_free_volume_distributions:
                molID = grid.voxelclusterID2molID[ID]
            else:
                molID = 1
            molID += molID_offset
            #data2write.append('{:^3} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2}\n'.format(ID, molID, 2, 0.0, p[0], p[1], p[2], 0, 0, 0)) 
            data2write.append('{} {} {} {} {} {} {} {} {} {}\n'.format(ID, molID, 2, 0.0, p[0], p[1], p[2], 0, 0, 0)) 
        string2write = ''.join(data2write)
        f.write(string2write)
    return

########################
# free and atoms/bonds #
########################
def bonds_free(m, grid, compute_free_volume_distributions, basename, version, log):
    filename = basename+'_bonds_free.data'
    log.out(f'\nWriting {filename} ...')
    # Writing new file with bonds information
    with open(filename, 'w') as f: 
        f.write(f'HEADER, free_volume.py {version}\n\n')
        
        # write qtys
        f.write(f'{len(m.atoms)+len(grid.voxel_freeIDs)} atoms\n')
        if m.bonds: f.write(f'{len(m.bonds)} bonds\n')
        f.write('\n')
        
        # write types
        f.write(f'{len(m.masses)+1} atom types\n')
        if m.bond_coeffs: f.write(f'{len(m.bond_coeffs)} bond types\n')
        f.write('\n')
        
        # write box
        f.write(f'{m.xbox_line}\n{m.ybox_line}\n{m.zbox_line}\n')
        
        # Write masses
        f.write('\nMasses\n\n')
        for i in m.masses:
            f.write('{} {:^10.5f}  # {}\n'.format(i, m.masses[i].coeffs[0], m.masses[i].element))
        freeTypeID = i+1
        f.write('{} {:^10.5f}  # {}\n'.format(freeTypeID, 1.0, 'Free'))
        
        # write bond ceoffs
        if m.bond_coeffs:
            f.write('\nBond Coeffs\n\n')
            for i in m.bond_coeffs:
                coeffs = [str(i) for i in m.bond_coeffs[i].coeffs]
                f.write('{} {}\n'.format(i, ' '.join(coeffs)))
            
        # Write atoms    
        f.write('\nAtoms # full\n\n')
        molids = set(); data2write = [];
        for ID in tqdm(m.atoms):
            atom = m.atoms[ID]
            try: molID = atom.molid
            except: molID = 1
            molids.add(molID)
            #data2write.append('{:^3} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2}\n'.format(ID, molID, atom.type, atom.charge, atom.x, atom.y, atom.z, 0, 0, 0)) 
            data2write.append('{} {} {} {} {} {} {} {} {} {}\n'.format(ID, molID, atom.type, atom.charge, atom.x, atom.y, atom.z, 0, 0, 0)) 
        
        # write free voxels
        molID_offset = max(molids)
        for ID in tqdm(grid.voxel_freeIDs):
            p = grid.voxel_freeIDs[ID]
            if compute_free_volume_distributions:
                molID = grid.voxelclusterID2molID[ID]
            else:
                molID = 1
            molID += molID_offset
            #data2write.append('{:^3} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2}\n'.format(ID, molID, freeTypeID, 0.0, p[0], p[1], p[2], 0, 0, 0))
            data2write.append('{} {} {} {} {} {} {} {} {} {}\n'.format(ID, molID, freeTypeID, 0.0, p[0], p[1], p[2], 0, 0, 0)) 
        string2write = ''.join(data2write)
        f.write(string2write)

        # write bonds
        if m.bonds:
            f.write('\nBonds\n\n')
            for i in m.bonds:
                bond = m.bonds[i]
                id1, id2 = bond.atomids
                f.write('{} {} {} {}\n'.format(i, bond.type, id1, id2))
    return