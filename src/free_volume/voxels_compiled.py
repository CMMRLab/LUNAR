# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
December 6th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.free_volume.misc_func as misc_func
import numpy as np
import numba
import math
import time

############################
# Class to generate voxels #
############################
class generate:
    def __init__(self, m, max_voxel_size, boundary, vdw_radius, probe_diameter, vdw_method, run_mode, log):        
        # Get box dimensions
        self.lx, self.ly, self.lz, self.cx, self.cy, self.cz, self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi = misc_func.get_box_dimensions(m)
        
        # Find voxel domain decomposition grid size and check vdw radii is provided
        system_vdw_radii = set()
        if vdw_method == 'dict':
            for element in m.elements:
                if element in vdw_radius:
                    system_vdw_radii.add(vdw_radius[element])
                else: log.error(f'ERROR element {element} is not in vdw_radius dictionary')
        elif vdw_method in ['class1', 'class2']:
            for i in m.pair_coeffs:
                sigma = m.pair_coeffs[i].coeffs[1]
                vdw_radii = sigma/2
                system_vdw_radii.add(vdw_radii)
        else: log.error(f'ERROR vdw_method {vdw_method} not supported')
        system_vdw_radii = sorted(system_vdw_radii)
        if not system_vdw_radii:
            log.error('ERROR unable to assign vdw radii. Perhaps change vdw_method')
        
        # Find images and scaled images 
        self.images, self.pflags = misc_func.generate_iflags(m, boundary, log)
        self.scaled_images = [(ix*self.lx, iy*self.ly, iz*self.lz) for (ix, iy, iz) in self.images]
        
        # Compute nx, ny, nz number of voxels required to build 3D voxel grid
        self.nx = math.ceil(self.lx/max_voxel_size)
        self.ny = math.ceil(self.ly/max_voxel_size)
        self.nz = math.ceil(self.lz/max_voxel_size)
        
        # Compute dx, dy, dz size of each voxel
        self.dx = self.lx/self.nx
        self.dy = self.ly/self.ny
        self.dz = self.lz/self.nz
        self.halfdx = self.dx/2
        self.halfdy = self.dy/2
        self.halfdz = self.dz/2
        self.voxel_diagnol = math.sqrt(self.dx**2 + self.dy**2 + self.dz**2)
        
        # Set probe_diameter if shortcut flag exists
        self.probe_diameter = probe_diameter
        if probe_diameter == 'min-voxel':
            log.out("probe_diameter = 'min-voxel', finding minimum dimension of voxel ...")
            self.probe_diameter = min([self.dx, self.dy, self.dz])
            log.out(f'   updated probe_diameter = {self.probe_diameter}')
        
        # Generate voxels (numba method)
        start_time = time.time()
        if run_mode == 'numba':
            self.numpy_voxels, self.numpy_radial, self.numpy_flags, self.numpyx, self.numpyy, self.numpyz = generate_serial(self.nx, self.ny, self.nz, self.dx, self.dy, self.dz, self.xlo, self.ylo, self.zlo, self.cx, self.cy, self.cz)
        elif run_mode == 'numba-p':
            # Serial is quicker for this method so stick with serial here
            self.numpy_voxels, self.numpy_radial, self.numpy_flags, self.numpyx, self.numpyy, self.numpyz = generate_serial(self.nx, self.ny, self.nz, self.dx, self.dy, self.dz, self.xlo, self.ylo, self.zlo, self.cx, self.cy, self.cz)
        else: log.error(f'ERROR input run_mode = {run_mode} not supported')
        self.voxelIDs = dict(enumerate(self.numpy_voxels, 1))  # { voxelID : (x, y, z) }
        execution_time = (time.time() - start_time)
        log.out(f'Voxel generation execution time: {execution_time} (seconds)')


##########################################################
# numba njit compilation for generating voxels in serial #
##########################################################
@numba.njit(parallel=False)
def generate_serial(nx, ny, nz, dx, dy, dz, xlo, ylo, zlo, cx, cy, cz):
    print('\n\nGenerating voxels (compiled - serial) ...')
    nvoxels = nx*ny*nz; progress_increment = 10; count = 0;
    xadd = dx/2 + xlo; yadd =  dy/2 + ylo; zadd = dz/2 + zlo;
    numpy_voxels = np.zeros(3*nvoxels).reshape((nvoxels, 3))
    numpy_radial = np.zeros(nvoxels); numpy_flags = np.ones(nvoxels)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                xpos = i*dx + xadd
                ypos = j*dy + yadd
                zpos = k*dz + zadd
                diffx = cx - xpos
                diffy = cy - ypos
                diffz = cz - zpos
                r = np.sqrt(diffx*diffx + diffy*diffy + diffz*diffz)
                numpy_voxels[count][0] = xpos
                numpy_voxels[count][1] = ypos
                numpy_voxels[count][2] = zpos
                numpy_radial[count] = r
                count += 1
                if 100*count/nvoxels % progress_increment == 0:
                    print('progress: ', int(100*count/nvoxels),'%')
    numpyx = np.unique(numpy_voxels[:,0])
    numpyy = np.unique(numpy_voxels[:,1])
    numpyz = np.unique(numpy_voxels[:,2])
    return numpy_voxels, numpy_radial, numpy_flags, numpyx, numpyy, numpyz


############################################################
# numba njit compilation for generating voxels in parallel #
############################################################
@numba.njit(parallel=True)
def generate_parallel(nx, ny, nz, dx, dy, dz, xlo, ylo, zlo, cx, cy, cz):
    print('\n\nGenerating voxels (compiled - parallel) ...')
    nvoxels = nx*ny*nz;xadd = dx/2 + xlo; yadd =  dy/2 + ylo; zadd = dz/2 + zlo;
    numpy_voxels = np.zeros(3*nvoxels).reshape((nvoxels, 3))
    numpy_radial = np.zeros(nvoxels); numpy_flags = np.ones(nvoxels)
    for i in numba.prange(nx):
        for j in range(ny):
            for k in range(nz):
                xpos = i*dx + xadd
                ypos = j*dy + yadd
                zpos = k*dz + zadd
                diffx = cx - xpos
                diffy = cy - ypos
                diffz = cz - zpos
                r = np.sqrt(diffx*diffx + diffy*diffy + diffz*diffz)
                numpy_voxels[i*ny*nz+j*nz+k][0] = xpos
                numpy_voxels[i*ny*nz+j*nz+k][1] = ypos
                numpy_voxels[i*ny*nz+j*nz+k][2] = zpos
                numpy_radial[i*ny*nz+j*nz+k] = r
    numpyx = np.unique(numpy_voxels[:,0])
    numpyy = np.unique(numpy_voxels[:,1])
    numpyz = np.unique(numpy_voxels[:,2])
    return numpy_voxels, numpy_radial, numpy_flags, numpyx, numpyy, numpyz
                    

######################################################################
# numba njit compilation for finding atom/free volume in serial mode #
######################################################################
@numba.njit(parallel=False)
def atom_vol_calc_serial(atoms, periodic_atoms, radial_voxel, voxels, flags, max_voxel_size, probe_radius):
    print('Finding atom volumes and free volumes (compiled - serial) ...')
    progress_increment = 10; count = 0; natoms = atoms.shape[:1][0]; nvoxels = voxels.shape[:1][0];
    pbc_count = 0; npossible_pbc = 0; voxel_atomIndexs = np.zeros(voxels.shape[:1])
    for i in range(natoms):
        atom = atoms[i]; count += 1;
        x1 = atom[0]; y1 = atom[1]; z1 = atom[2]; radial_atom = atom[3];
        atompos = np.array([x1, y1, z1])
        vdw_radii = atom[4] + probe_radius; atom_edgeflag = atom[5];
        min_radius = radial_atom - vdw_radii - max_voxel_size
        max_radius = radial_atom + vdw_radii + max_voxel_size
        if atom_edgeflag == 1: npossible_pbc += 1
        if 100*count/natoms % progress_increment == 0:
            print('progress: ', int(100*count/natoms),'%')
        for j in range(nvoxels):
            if flags[j] == 0: continue
            if min_radius < radial_voxel[j] < max_radius:
                if atom_edgeflag == 1:
                    voxel = voxels[j]
                    for x1i, y1i, z1i in periodic_atoms[i]:
                        if abs(voxel[0] - x1i) > vdw_radii: continue
                        elif abs(voxel[1] - y1i) > vdw_radii: continue
                        elif abs(voxel[2] - z1i) > vdw_radii: continue
                        atompos = np.array([x1i, y1i, z1i])
                        distance_ppp = np.sqrt(np.sum(np.square(atompos-voxels[j])))
                        if distance_ppp < vdw_radii:
                            flags[j] = 0; voxel_atomIndexs[j] = atom[6]; pbc_count += 1;
                else:
                    voxel = voxels[j]
                    if abs(voxel[0] - x1) > vdw_radii: continue
                    elif abs(voxel[1] - y1) > vdw_radii: continue
                    elif abs(voxel[2] - z1) > vdw_radii: continue
                    distance_fff = np.sqrt(np.sum(np.square(atompos-voxels[j])))
                    if distance_fff < vdw_radii:
                        flags[j] = 0; voxel_atomIndexs[j] = atom[6];
    return pbc_count, npossible_pbc, voxel_atomIndexs


########################################################################
# numba njit compilation for finding atom/free volume in parallel mode #
########################################################################
@numba.njit(parallel=True)
def atom_vol_calc_parallel(atoms, periodic_atoms, radial_voxel, voxels, flags, max_voxel_size, probe_radius):
    print('Finding atom volumes and free volumes (compiled - parallel) ...')
    progress_increment = 10; checked = np.zeros(len(atoms)); natoms = atoms.shape[:1][0]; nvoxels = voxels.shape[:1][0];
    pbc_count = 0; npossible_pbc = 0; voxel_atomIndexs = np.zeros(voxels.shape[:1])
    for i in numba.prange(natoms):
        atom = atoms[i]; checked[i] = 1;
        x1 = atom[0]; y1 = atom[1]; z1 = atom[2]; radial_atom = atom[3];
        atompos = np.array([x1, y1, z1])
        vdw_radii = atom[4] + probe_radius; atom_edgeflag = atom[5];
        min_radius = radial_atom - vdw_radii - max_voxel_size
        max_radius = radial_atom + vdw_radii + max_voxel_size
        if atom_edgeflag == 1: npossible_pbc += 1
        if 100*np.sum(checked)/natoms % progress_increment == 0:
            print('progress: ', int(100*np.sum(checked)/natoms),'%')
        for j in range(nvoxels):
            if flags[j] == 0: continue
            if min_radius < radial_voxel[j] < max_radius:
                if atom_edgeflag == 1:
                    voxel = voxels[j]
                    for x1i, y1i, z1i in periodic_atoms[i]:
                        if abs(voxel[0] - x1i) > vdw_radii: continue
                        elif abs(voxel[1] - y1i) > vdw_radii: continue
                        elif abs(voxel[2] - z1i) > vdw_radii: continue
                        atompos = np.array([x1i, y1i, z1i])
                        distance_ppp = np.sqrt(np.sum(np.square(atompos-voxels[j])))
                        if distance_ppp < vdw_radii:
                            flags[j] = 0; voxel_atomIndexs[j] = atom[6]; pbc_count += 1;
                else:
                    voxel = voxels[j]
                    if abs(voxel[0] - x1) > vdw_radii: continue
                    elif abs(voxel[1] - y1) > vdw_radii: continue
                    elif abs(voxel[2] - z1) > vdw_radii: continue
                    distance_fff = np.sqrt(np.sum(np.square(atompos-voxels[j])))
                    if distance_fff < vdw_radii:
                        flags[j] = 0; voxel_atomIndexs[j] = atom[6];
    return pbc_count, npossible_pbc, voxel_atomIndexs


############################################################
# numba njit compilation for spatial distribution (serial) #
############################################################
@numba.njit(parallel=False)
def spat_dist_serial(xpos, ypos, zpos, voxels, flags, voxel_volume):
    volx = np.zeros(xpos.shape[0], dtype=float); 
    voly = np.zeros(ypos.shape[0], dtype=float);
    volz = np.zeros(zpos.shape[0], dtype=float);
    progress_increment = 10; count = 0; nvoxels = voxels.shape[:1][0];
    for i in range(voxels.shape[:1][0]):
        count += 1
        if 100*count/nvoxels % progress_increment == 0:
            print('progress: ', int(100*count/nvoxels),'%')
        if flags[i] == 0: continue
        voxel = voxels[i]
        xindex = np.where(xpos == voxel[0])[0][0]
        yindex = np.where(ypos == voxel[1])[0][0]
        zindex = np.where(zpos == voxel[2])[0][0]
        volx[xindex] += voxel_volume
        voly[yindex] += voxel_volume
        volz[zindex] += voxel_volume
    return volx, voly, volz


##############################################################
# numba njit compilation for spatial distribution (parallel) #
##############################################################
@numba.njit(parallel=True)
def spat_dist_parallel(xpos, ypos, zpos, voxels, flags, voxel_volume):
    volx = np.zeros(xpos.shape[0], dtype=float); 
    voly = np.zeros(ypos.shape[0], dtype=float);
    volz = np.zeros(zpos.shape[0], dtype=float);
    for i in numba.prange(voxels.shape[:1][0]):
        if flags[i] == 0: continue
        voxel = voxels[i]
        xindex = np.where(xpos == voxel[0])[0][0]
        yindex = np.where(ypos == voxel[1])[0][0]
        zindex = np.where(zpos == voxel[2])[0][0]
        volx[xindex] += voxel_volume
        voly[yindex] += voxel_volume
        volz[zindex] += voxel_volume
    return volx, voly, volz


################################################################
# numba njit compilation for free volume distribution (serial) #
################################################################
@numba.njit(parallel=False)
def find_zero_index(arr):
    index = 0;
    for index in range(arr.shape[0]):
        if arr[index] == 0: break
    return index

@numba.njit(parallel=False)
def free_dist_serial(voxels, flags, max_voxel_size, xlo, xhi, ylo, yhi, zlo, zhi, scaled_images, radial_voxel):
    # Find if voxels periodic postions
    nvoxels = voxels.shape[:1][0];
    voxel_periodic_flags = np.zeros(nvoxels);
    for i in range(nvoxels):
        if flags[i] == 0: continue
        x, y, z = voxels[i]; flag = 0;
        if abs(x - xlo) <= max_voxel_size or abs(x + xlo) <= max_voxel_size: flag = 1
        elif abs(x - xhi) <= max_voxel_size or abs(x + xhi) <= max_voxel_size: flag = 1 
        elif abs(y - ylo) <= max_voxel_size or abs(y + ylo) <= max_voxel_size: flag = 1
        elif abs(y - yhi) <= max_voxel_size or abs(y + yhi) <= max_voxel_size: flag = 1
        elif abs(z - zlo) <= max_voxel_size or abs(z + zlo) <= max_voxel_size: flag = 1
        elif abs(z - zhi) <= max_voxel_size or abs(z + zhi) <= max_voxel_size: flag = 1
        voxel_periodic_flags[i] = flag
    
    # Find voxel connectivty
    numpy_graph = np.zeros(30*(nvoxels+1), dtype=np.int64).reshape((nvoxels+1, 30)) # should be 27, but 30 allows overflow
    internal_flags = np.copy(flags); progress_increment = 10; count = 0; nimages = 27
    for i in range(nvoxels):
        count += 1
        if 100*count/nvoxels % progress_increment == 0:
            print('progress: ', int(100*count/nvoxels),'%')
        if internal_flags[i] == 0: continue
        x1, y1, z1 = voxels[i]
        r1 = radial_voxel[i]
        min_radius = r1 - 2*max_voxel_size
        max_radius = r1 + 2*max_voxel_size
        periodic_postions = np.zeros(nimages*3).reshape((nimages, 3))
        for a in range(nimages):
            ixlx, iyly, izlz = scaled_images[a]
            periodic_postions[a][0] = x1+ixlx
            periodic_postions[a][1] = y1+iyly
            periodic_postions[a][2] = z1+izlz
        for j in range(i, nvoxels):
            if internal_flags[j] == 0 or i == j: continue
            if min_radius < radial_voxel[j] < max_radius:
                x2, y2, z2 = voxels[j]
                # Find non-periodic connectivity
                if voxel_periodic_flags[i] == 0 and voxel_periodic_flags[j] == 0:
                    if abs(x1-x2) <= max_voxel_size and abs(y1-y2) <= max_voxel_size and abs(z1-z2) <= max_voxel_size:
                        distance = np.sqrt(np.sum(np.square(voxels[i]-voxels[j])))
                        if distance <= max_voxel_size:
                            id1 = i+1; id2 = j+1;
                            gindex1 = find_zero_index(numpy_graph[id1])
                            gindex2 = find_zero_index(numpy_graph[id2])
                            if gindex1 < 29 and gindex2 < 29:
                                if not np.any(numpy_graph[id1] == id2):
                                    numpy_graph[id1][gindex1] = id2
                                if not np.any(numpy_graph[id2] == id1):
                                    numpy_graph[id2][gindex2] = id1
                else: # else find periodic connectivty
                    for x1i, y1i, z1i in periodic_postions:
                        if abs(x1i-x2) <= max_voxel_size and abs(y1i-y2) <= max_voxel_size and abs(z1i-z2) <= max_voxel_size:
                            voxels_ppp = np.array([x1i, y1i, z1i])
                            distance = np.sqrt(np.sum(np.square(voxels_ppp-voxels[j])))
                            if distance <= max_voxel_size:
                                id1 = i+1; id2 = j+1;
                                gindex1 = find_zero_index(numpy_graph[id1])
                                gindex2 = find_zero_index(numpy_graph[id2])
                                if gindex1 < 29 and gindex2 < 29:
                                    if not np.any(numpy_graph[id1] == id2):
                                        numpy_graph[id1][gindex1] = id2
                                    if not np.any(numpy_graph[id2] == id1):
                                        numpy_graph[id2][gindex2] = id1
                                break
        internal_flags[i] = 0
    return numpy_graph


##################################################################
# numba njit compilation for free volume distribution (parallel) #
##################################################################
@numba.njit(parallel=True)
def free_dist_parallel(voxels, flags, max_voxel_size, xlo, xhi, ylo, yhi, zlo, zhi, scaled_images, radial_voxel):
    # Find if voxels periodic postions
    nvoxels = voxels.shape[:1][0];
    voxel_periodic_flags = np.zeros(nvoxels);
    for i in range(nvoxels):
        if flags[i] == 0: continue
        x, y, z = voxels[i]; flag = 0;
        if abs(x - xlo) <= max_voxel_size or abs(x + xlo) <= max_voxel_size: flag = 1
        elif abs(x - xhi) <= max_voxel_size or abs(x + xhi) <= max_voxel_size: flag = 1 
        elif abs(y - ylo) <= max_voxel_size or abs(y + ylo) <= max_voxel_size: flag = 1
        elif abs(y - yhi) <= max_voxel_size or abs(y + yhi) <= max_voxel_size: flag = 1
        elif abs(z - zlo) <= max_voxel_size or abs(z + zlo) <= max_voxel_size: flag = 1
        elif abs(z - zhi) <= max_voxel_size or abs(z + zhi) <= max_voxel_size: flag = 1
        voxel_periodic_flags[i] = flag
    
    # Find voxel connectivty
    numpy_graph = np.zeros(30*(nvoxels+1), dtype=np.int64).reshape((nvoxels+1, 30)) # should be 27, but 30 allows overflow
    internal_flags = np.copy(flags); nimages = 27
    for i in numba.prange(nvoxels):
        if internal_flags[i] == 0: continue
        x1, y1, z1 = voxels[i]
        r1 = radial_voxel[i]
        min_radius = r1 - 2*max_voxel_size
        max_radius = r1 + 2*max_voxel_size
        periodic_postions = np.zeros(nimages*3).reshape((nimages, 3))
        for a in range(nimages):
            ixlx, iyly, izlz = scaled_images[a]
            periodic_postions[a][0] = x1+ixlx
            periodic_postions[a][1] = y1+iyly
            periodic_postions[a][2] = z1+izlz
        for j in range(i, nvoxels):
            if internal_flags[j] == 0 or i == j: continue
            if min_radius < radial_voxel[j] < max_radius:
                x2, y2, z2 = voxels[j]
                # Find non-periodic connectivity
                if voxel_periodic_flags[i] == 0 and voxel_periodic_flags[j] == 0:
                    if abs(x1-x2) <= max_voxel_size and abs(y1-y2) <= max_voxel_size and abs(z1-z2) <= max_voxel_size:
                        distance = np.sqrt(np.sum(np.square(voxels[i]-voxels[j])))
                        if distance <= max_voxel_size:
                            id1 = i+1; id2 = j+1;
                            gindex1 = find_zero_index(numpy_graph[id1])
                            gindex2 = find_zero_index(numpy_graph[id2])
                            if gindex1 < 29 and gindex2 < 29:
                                if not np.any(numpy_graph[id1] == id2):
                                    numpy_graph[id1][gindex1] = id2
                                if not np.any(numpy_graph[id2] == id1):
                                    numpy_graph[id2][gindex2] = id1
                else: # else find periodic connectivty
                    for x1i, y1i, z1i in periodic_postions:
                        if abs(x1i-x2) <= max_voxel_size and abs(y1i-y2) <= max_voxel_size and abs(z1i-z2) <= max_voxel_size:
                            voxels_ppp = np.array([x1i, y1i, z1i])
                            distance = np.sqrt(np.sum(np.square(voxels_ppp-voxels[j])))
                            if distance <= max_voxel_size:
                                id1 = i+1; id2 = j+1;
                                gindex1 = find_zero_index(numpy_graph[id1])
                                gindex2 = find_zero_index(numpy_graph[id2])
                                if gindex1 < 29 and gindex2 < 29:
                                    if not np.any(numpy_graph[id1] == id2):
                                        numpy_graph[id1][gindex1] = id2
                                    if not np.any(numpy_graph[id2] == id1):
                                        numpy_graph[id2][gindex2] = id1
                                break
        internal_flags[i] = 0
    return numpy_graph

        
        
##############################
# Class to anaylze atom size #
##############################
class analyze:
    def __init__(self, m, v, vdw_radius, boundary, max_voxel_size, compute_free_volume_distributions, mass_map, run_mode, files2write, probe_diameter, vdw_method, log):  
        #---------------------------#
        # Base free volume analysis #
        #---------------------------#          
        # set atomTypeIDs
        self.atomTypeIDs = {} # { atom Type ID : (element, mass) }
        if vdw_method == 'dict':
            for ID, element in enumerate(sorted(m.elements), 1):
                self.atomTypeIDs[ID] = (element, mass_map[element][0])
        elif vdw_method in ['class1', 'class2']:
            for ID in m.masses:
                masses = m.masses[ID]
                self.atomTypeIDs[ID] = (masses.element, masses.coeffs[0])
        else: log.error(f'ERROR vdw_method {vdw_method} not supported')
            
        # Generate numpy arrays  (use for vectorizing euclidean distances calculations)
        log.out('\n\nBuilding numpy arrays to compile code')
        voxels = v.numpy_voxels; flags = v.numpy_flags;  radial_voxel = v.numpy_radial; periodic_atoms = [] # [x, y, z]
        atoms = [] # [x, y, z, r, vdw, pbcflag]; pbcflag = 0 (not near edge); pbcflag = 1 (near edge);   
        for i in m.atoms:
            atom = m.atoms[i]
            x = atom.x; y = atom.y; z = atom.z;
            radial_atom = misc_func.compute_distance(x, y, z, v.cx, v.cy, v.cz)
            
            # Set atom vdw radii based on method
            if vdw_method == 'dict':
                vdw_radii = vdw_radius[atom.element]
            elif vdw_method in ['class1', 'class2']:
                sigma = m.pair_coeffs[atom.type].coeffs[1]
                vdw_radii = sigma/2
            else: log.error(f'ERROR vdw_method {vdw_method} not supported')
                
            if v.pflags.count('f') == 3: atom_edgeflag = False
            else: atom_edgeflag = misc_func.check_near_edge(x, y, z, vdw_radii, v.xlo, v.xhi, v.ylo, v.yhi, v.zlo, v.zhi)
            pbcflag = 0;
            if atom_edgeflag: pbcflag = 1;
            periodic_atoms.append(misc_func.find_periodic_postions(v.scaled_images, x, y, z, v.cx, v.cy, v.cz, Npos=12))
            atoms.append([x, y, z, radial_atom, vdw_radii, pbcflag, atom.type])
        atoms = np.array(atoms); periodic_atoms = np.array(periodic_atoms)
        
        # Used numpy/numba implementation to find atom volumes and free volumes
        start_time = time.time()
        if run_mode == 'numba':
            self.pbc_count, self.npossible_pbc, voxel_atomIndexs = atom_vol_calc_serial(atoms, periodic_atoms, radial_voxel, voxels, flags, max_voxel_size, probe_diameter/2)
        elif run_mode == 'numba-p':
            self.pbc_count, self.npossible_pbc, voxel_atomIndexs = atom_vol_calc_parallel(atoms, periodic_atoms, radial_voxel, voxels, flags, max_voxel_size, probe_diameter/2)
        else: log.error(f'ERROR input run_mode = {run_mode} not supported')
        execution_time = (time.time() - start_time)
        log.out(f'Atom volume calculation execution time: {execution_time} (seconds)')
        if files2write['write_atoms_free'] or files2write['write_atoms_only'] or files2write['write_bonds_free']:
            self.voxel_atomIDs = {}; 
            for i in np.ndindex(voxel_atomIndexs.shape):
                if voxel_atomIndexs[i] > 0:
                    pos = voxels[i[0]]; TypeID = voxel_atomIndexs[i];
                    self.voxel_atomIDs[i[0]+1] = (pos[0], pos[1], pos[2], TypeID)
            nvoxel_atomIDs = len(self.voxel_atomIDs)
        else: nvoxel_atomIDs = len(voxel_atomIndexs[voxel_atomIndexs > 0])

        # Find the freeIDs after all atomIDs have been iterate through
        if files2write['write_atoms_free'] or files2write['write_free_only'] or files2write['write_bonds_free'] or compute_free_volume_distributions:
            self.voxel_freeIDs = {ID:v.voxelIDs[ID] for ID, flag in enumerate(flags, 1) if flag > 0}
            nvoxel_freeIDs = len(self.voxel_freeIDs)
        else: nvoxel_freeIDs = len(flags[flags > 0])
        
        # Find simulation cell volume
        self.simulation_volume = v.lx*v.ly*v.lz
        
        # Find free volume
        self.voxel_volume = v.dx*v.dy*v.dz
        self.atom_volume = self.voxel_volume*nvoxel_atomIDs
        self.free_volume = self.voxel_volume*nvoxel_freeIDs
        self.pfree_volume = 100*(self.free_volume/self.simulation_volume)
        
        
        #-------------------------------------------#
        # Compute free volume spatial distributions #
        #-------------------------------------------#
        if files2write['write_spat_dis-x'] or files2write['write_spat_dis-y'] or files2write['write_spat_dis-z']: 
            start_time = time.time()
            if run_mode == 'numba': 
                log.out('\n\nFinding free volume spatial distributions (compiled - serial) ...')
                volx, voly, volz = spat_dist_serial(v.numpyx, v.numpyy, v.numpyz, voxels, flags, self.voxel_volume)
            elif run_mode == 'numba-p': 
                log.out('\n\nFinding free volume spatial distributions (compiled - parallel) ...')
                volx, voly, volz = spat_dist_parallel(v.numpyx, v.numpyy, v.numpyz, voxels, flags, self.voxel_volume)
            self.spatial_distributions = {'x':{v.numpyx[i]:volx[i] for i in range(v.numpyx.shape[:1][0])},
                                          'y':{v.numpyy[i]:voly[i] for i in range(v.numpyy.shape[:1][0])},
                                          'z':{v.numpyz[i]:volz[i] for i in range(v.numpyz.shape[:1][0])}}
            execution_time = (time.time() - start_time)
            log.out(f'Free volume spatial distributions execution time: {execution_time} (seconds)')
            
            
        #-------------------------------------------------------------------#
        # Free volume connectivity and distribution analysis (EXPEIRMENTAL) #
        #-------------------------------------------------------------------# 
        if compute_free_volume_distributions:
            log.out('\n\n----------------------------------------------')
            log.out('| Starting free volume distribution analysis |')
            log.out('----------------------------------------------')
    
            # numba find free volume connectivity
            start_time = time.time()
            if run_mode == 'numba': 
                log.out('Finding free volume voxelID connectivity (compiled - serial) ...')
                numpy_graph = free_dist_serial(voxels, flags, max_voxel_size, v.xlo, v.xhi, v.ylo, v.yhi, v.zlo, v.zhi, np.array(v.scaled_images), radial_voxel)
            elif run_mode == 'numba-p': 
                log.out('Finding free volume voxelID connectivity (compiled - parallel) ...')
                numpy_graph = free_dist_parallel(voxels, flags, max_voxel_size, v.xlo, v.xhi, v.ylo, v.yhi, v.zlo, v.zhi, np.array(v.scaled_images), radial_voxel)
            free_volume_voxel_graph = {ID:set() for ID in self.voxel_freeIDs} # { voxelID : set(bonded voxels) }
            for i in range(numpy_graph.shape[:1][0]):
                if sum(numpy_graph[i]) > 0 and i > 0:
                    for j in numpy_graph[i]:
                        if j > 0:
                            id1 = int(i); id2 = int(j);
                            free_volume_voxel_graph[id1].add(id2)
                            free_volume_voxel_graph[id2].add(id1)
            execution_time = (time.time() - start_time)
            log.out(f'Free volume distribution analysis execution time: {execution_time} (seconds)')
        
            # Find free volume clusters
            free_volume_clusters = set([]) # { (tuple of atoms in cluster1), (nclusters) }
            checked = {ID:False for ID in self.voxel_freeIDs}
            for ID in free_volume_voxel_graph:
                if checked[ID]: continue
                visited=set([ID]); queue=[ID];
                while queue:
                    s = queue.pop(0) 
                    for neighbor in free_volume_voxel_graph[s]:
                        if checked[neighbor]: continue
                        visited.add(neighbor)
                        queue.append(neighbor)
                        checked[neighbor]=True
                free_volume_clusters.add( tuple(sorted(visited)) )
            free_volume_clusters = sorted(free_volume_clusters, key=lambda x: x[0]) # Sort all clusters based on 1st atomID in cluster
            free_volume_clusters = sorted(free_volume_clusters, key=len, reverse=True) # Sort all clusters by number of atoms
            
            # Analyze free volume clusters
            class Info: pass # .size .volume .pvolume
            voxel_volume = v.dx*v.dy*v.dz; self.free_volume_clusters = {} # { clusterID : Info object }
            self.voxelclusterID2molID = {i:len(free_volume_clusters)+1 for i in self.voxel_freeIDs} # { voxelID : clustered }
            for ID, cluster in enumerate(free_volume_clusters, 1):
                for i in cluster:
                    self.voxelclusterID2molID[i] = ID
                xspan, yspan, zspan = misc_func.find_free_volume_cluster_span(cluster, self.voxel_freeIDs)
                I = Info()
                I.size = len(cluster)
                I.psize = 100*len(cluster)/len(self.voxel_freeIDs)
                I.volume = len(cluster)*voxel_volume
                I.pvolume = 100*(len(cluster)*voxel_volume)/self.simulation_volume
                I.xspan = xspan
                I.yspan = yspan
                I.zspan = zspan
                self.free_volume_clusters[ID] = I
                
                
            #print('self.pbc_bonded_count', self.pbc_bonded_count)
            #bonded = {len(free_volume_voxel_graph[i]) for i in free_volume_voxel_graph}
            #print(bonded)
            #print(len(self.free_volume_clusters))
            #print('count TRUE freeID_edgeflags: ', list(freeID_edgeflags.values()).count(True))