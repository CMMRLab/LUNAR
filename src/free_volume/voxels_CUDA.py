# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 6th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.free_volume.misc_func as misc_func
from numba import cuda
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
        
        # Generatre voxels on CPU or GPU based on inputs
        start_time = time.time()
        if m.CUDA_threads_per_block_voxels == 0:
            self.numpy_voxels, self.numpy_radial, self.numpy_flags, self.numpyx, self.numpyy, self.numpyz = generate_serial(self.nx, self.ny, self.nz, self.dx, self.dy, self.dz, self.xlo, self.ylo, self.zlo, self.cx, self.cy, self.cz)
        else:
            # CUDA Kernel to Generate voxels
            @cuda.jit
            def generate_CUDA(numpy_voxels, n_xyz, d_xyz, xyz_add):
                i, j = cuda.grid(2)
                if i < n_xyz[0] and j < n_xyz[1]:
                    for k in range(n_xyz[2]):
                        cuda.syncthreads()
                        numpy_voxels[i*n_xyz[1]*n_xyz[2] + j*n_xyz[2] + k][0] = i*d_xyz[0] + xyz_add[0]
                        numpy_voxels[i*n_xyz[1]*n_xyz[2] + j*n_xyz[2] + k][1] = j*d_xyz[1] + xyz_add[1]
                        numpy_voxels[i*n_xyz[1]*n_xyz[2] + j*n_xyz[2] + k][2] = k*d_xyz[2] + xyz_add[2]
                        cuda.syncthreads()
                        
            # CUDA Kernel to find radius from center
            @cuda.jit
            def radial_CUDA(numpy_voxels, numpy_radial, c_xyz):
                i = cuda.grid(1)
                if i < numpy_voxels.shape[:1][0]:
                    cuda.syncthreads()
                    diffx = c_xyz[0] - numpy_voxels[i][0]
                    diffy = c_xyz[1] - numpy_voxels[i][1]
                    diffz = c_xyz[2] - numpy_voxels[i][2]
                    numpy_radial[i] = ( diffx*diffx + diffy*diffy + diffz*diffz )**0.5
                    cuda.syncthreads()
            
            # Initialize data structs for GPU
            nvoxels = self.nx*self.ny*self.nz;
            xadd = self.dx/2 + self.xlo; yadd =  self.dy/2 + self.ylo; zadd = self.dz/2 + self.zlo;
            xyz_add = cuda.to_device(np.array([xadd, yadd, zadd]))
            n_xyz = cuda.to_device(np.array([self.nx, self.ny, self.nz]))
            d_xyz = cuda.to_device(np.array([self.dx, self.dy, self.dz]))
            c_xyz = cuda.to_device(np.array([self.cx, self.cy, self.cz]))
            numpy_voxels = cuda.to_device(np.zeros(3*nvoxels).reshape((nvoxels, 3)))
            numpy_radial = cuda.to_device(np.zeros(nvoxels)); 
            
            # Function to round threads_per_block to nears 4th
            def round2nearest(value):
                opt = [4, 8, 16, 32, 64, 128, 512, 1024]
                dif = [abs(i-value) for i in opt]
                idx = dif.index(min(dif))
                return opt[idx]
            
            # Generate voxels on GPU
            threads_per_block_xy = round2nearest( m.CUDA_threads_per_block_voxels**(1/3) )
            threadsperblock = (threads_per_block_xy, threads_per_block_xy)
            blockspergrid_x = int(math.ceil(self.nx / threadsperblock[0]))
            blockspergrid_y = int(math.ceil(self.ny / threadsperblock[1]))
            blockspergrid = (blockspergrid_x, blockspergrid_y)
            log.out(f'\n\nGenerating voxels ({run_mode}-{threads_per_block_xy}x{threads_per_block_xy}x{threads_per_block_xy}) ...')
            generate_CUDA[blockspergrid, threadsperblock](numpy_voxels, n_xyz, d_xyz, xyz_add)
            self.numpy_voxels = numpy_voxels.copy_to_host();
            
            # Finding radial voxels
            threads_per_block_voxels = m.CUDA_threads_per_block_voxels
            log.out(f'Finding radius of voxels from center ({run_mode}-{threads_per_block_voxels}) ...')
            blocks_per_grid_voxels = math.ceil(nvoxels / threads_per_block_voxels) 
            radial_CUDA[blocks_per_grid_voxels, threads_per_block_voxels](cuda.to_device(self.numpy_voxels), numpy_radial, c_xyz)
            self.numpy_radial = numpy_radial.copy_to_host();
            self.numpy_flags = np.ones(nvoxels)
            self.numpyx, self.numpyy, self.numpyz = unique_xyz(self.numpy_voxels)

        # Generate required python voxelIDs dict        
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


##########################################################
# numba njit compilation for unique x, y and z locations #
##########################################################
@numba.njit(parallel=False)
def unique_xyz(numpy_voxels):
    numpyx = np.unique(numpy_voxels[:,0])
    numpyy = np.unique(numpy_voxels[:,1])
    numpyz = np.unique(numpy_voxels[:,2])
    return numpyx, numpyy, numpyz


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
        else: log.out(f'ERROR vdw_method {vdw_method} not supported')
            
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
        
        #-----------------------------------------#
        # CUDA Kernel to compute atom/free volume #
        #-----------------------------------------#
        start_time = time.time()
        @cuda.jit
        def atom_vol_calc_CUDA(atoms, periodic_atoms, radial_voxel, voxels, flags, voxel_atomIndexs, sizes, pbc):
            i = cuda.grid(1)
            if i < atoms.shape[:1][0]:
                atom = atoms[i]
                x1 = atom[0]; y1 = atom[1]; z1 = atom[2]; radial_atom = atom[3];
                vdw_radii = atom[4] + sizes[1]; atom_edgeflag = atom[5];
                min_radius = radial_atom - vdw_radii - sizes[0]
                max_radius = radial_atom + vdw_radii + sizes[0]
                if atom_edgeflag == 1: cuda.atomic.add(pbc, 1, 1)
                for j in range(voxels.shape[:1][0]):
                    if flags[j] == 0: continue
                    if min_radius < radial_voxel[j] < max_radius:
                        if atom_edgeflag == 1:
                            voxel = voxels[j]
                            for x1i, y1i, z1i in periodic_atoms[i]:
                                if abs(voxel[0] - x1i) > vdw_radii: continue
                                elif abs(voxel[1] - y1i) > vdw_radii: continue
                                elif abs(voxel[2] - z1i) > vdw_radii: continue
                                dx = x1i - voxel[0]; dy = y1i - voxel[1]; dz = z1i - voxel[2];
                                distance_ppp = math.sqrt(dx*dx + dy*dy + dz*dz)
                                if distance_ppp < vdw_radii:
                                    flags[j] = 0; voxel_atomIndexs[j] = int(atom[6]);
                                    cuda.atomic.add(pbc, 0, 1)
                        else:
                            voxel = voxels[j]
                            if abs(voxel[0] - x1) > vdw_radii: continue
                            elif abs(voxel[1] - y1) > vdw_radii: continue
                            elif abs(voxel[2] - z1) > vdw_radii: continue
                            dx = atom[0] - voxel[0]; dy = atom[1] - voxel[1]; dz = atom[2] - voxel[2];
                            distance_fff = math.sqrt(dx*dx + dy*dy + dz*dz)
                            if distance_fff < vdw_radii:
                                flags[j] = 0; voxel_atomIndexs[j] = int(atom[6]);
                                
        # Host data struct setup
        nvoxels = voxels.shape[:1][0]; natoms = atoms.shape[:1][0];
        voxel_atomIndexs = cuda.to_device(np.zeros(nvoxels))
        pbc = cuda.to_device(np.array([0, 0])) # [pbc_count, npossible_pbc]
        sizes = cuda.to_device(np.array([max_voxel_size, probe_diameter/2]))
        flags = cuda.to_device(flags)
        
        # Finding volumes using CUDA
        threads_per_block_atoms = m.CUDA_threads_per_block_atoms
        blocks_per_grid_atoms = math.ceil(natoms / threads_per_block_atoms)
        log.out(f'Finding atom volumes and free volumes ({run_mode}-{threads_per_block_atoms}) ...')
        atom_vol_calc_CUDA[blocks_per_grid_atoms, threads_per_block_atoms](cuda.to_device(atoms), cuda.to_device(periodic_atoms), cuda.to_device(radial_voxel), cuda.to_device(voxels), flags, voxel_atomIndexs, sizes, pbc)
        flags = flags.copy_to_host(); voxel_atomIndexs = voxel_atomIndexs.copy_to_host();
        self.pbc_count = pbc[0]; self.npossible_pbc = pbc[1]
        execution_time = (time.time() - start_time)
        log.out(f'Atom volume calculation execution time: {execution_time} (seconds)')
        
        # post process CUDA results
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
            if m.CUDA_threads_per_block_voxels == 0:
                log.out('Finding free volume voxelID connectivity (compiled - parallel) ...')
                numpy_graph = free_dist_parallel(voxels, flags, max_voxel_size, v.xlo, v.xhi, v.ylo, v.yhi, v.zlo, v.zhi, np.array(v.scaled_images), radial_voxel)
            else:
                #-------------------------------------------------#
                # CUDA Kernel to compute free volume connectivity #
                # is slower then numba parallel so will not use   #
                # until further optimization can happen ....      #
                #-------------------------------------------------#
                # Generate intial data structs
                @cuda.jit
                def free_dist_pre_setup_CUDA(voxels, voxel_periodic_flags, scaled_images, max_voxel_size, box):
                    i = cuda.grid(1)
                    if i < voxels.shape[:1][0]:
                        # Set periodic flags
                        x, y, z = voxels[i]; flag = 0;
                        if abs(x - box[0]) <= max_voxel_size or abs(x + box[0]) <= max_voxel_size: flag = 1
                        elif abs(x - box[1]) <= max_voxel_size or abs(x + box[1]) <= max_voxel_size: flag = 1 
                        elif abs(y - box[2]) <= max_voxel_size or abs(y + box[2]) <= max_voxel_size: flag = 1
                        elif abs(y - box[3]) <= max_voxel_size or abs(y + box[3]) <= max_voxel_size: flag = 1
                        elif abs(z - box[4]) <= max_voxel_size or abs(z + box[4]) <= max_voxel_size: flag = 1
                        elif abs(z - box[5]) <= max_voxel_size or abs(z + box[5]) <= max_voxel_size: flag = 1
                        voxel_periodic_flags[i] = flag
                
                # Host data struct setup
                start_time = time.time()
                voxel_periodic_flags = cuda.to_device(np.zeros(nvoxels))
                
                # Generate free dist intial arrays
                threads_per_block_voxels = m.CUDA_threads_per_block_voxels
                blocks_per_grid_voxels = math.ceil(nvoxels / threads_per_block_voxels)
                log.out(f'Setting up free volume voxelID connectivity calculation ({run_mode}-{threads_per_block_voxels}) ...')
                box = cuda.to_device(np.array([v.xlo, v.xhi, v.ylo, v.yhi, v.zlo, v.zhi]))
                simages = cuda.to_device(np.array(v.scaled_images))
                free_dist_pre_setup_CUDA[blocks_per_grid_voxels, threads_per_block_voxels](cuda.to_device(voxels), voxel_periodic_flags, simages, max_voxel_size, box)
                voxel_periodic_flags = voxel_periodic_flags.copy_to_host();
    
                # Perform free volume connectivity anaylsis         
                @cuda.jit
                def free_dist_CUDA(voxels, radial_voxel, scaled_images, internal_flags, voxel_periodic_flags, numpy_graph, max_voxel_size):
                    i = cuda.grid(1)
                    if i < voxels.shape[:1][0] and internal_flags[i] != 0:
                        x1, y1, z1 = voxels[i]
                        r1 = radial_voxel[i]
                        min_radius = r1 - 2*max_voxel_size
                        max_radius = r1 + 2*max_voxel_size
                        for j in range(nvoxels):
                            if internal_flags[j] == 0 or i == j: continue
                            if min_radius < radial_voxel[j] < max_radius:
                                x2, y2, z2 = voxels[j]
                                # Find non-periodic connectivity
                                if voxel_periodic_flags[i] == 0 and voxel_periodic_flags[j] == 0:
                                    if abs(x1-x2) <= max_voxel_size and abs(y1-y2) <= max_voxel_size and abs(z1-z2) <= max_voxel_size:
                                        dx = (x1 - x2)**2
                                        dy = (y1 - y2)**2
                                        dz = (z1 - z2)**2
                                        distance = math.sqrt(dx + dy + dz)
                                        if distance <= max_voxel_size:
                                            id1 = i+1; id2 = j+1;
                                            for gindex1 in range(numpy_graph[id1].shape[0]):
                                                if numpy_graph[id1][gindex1] == 0: break
                                            for gindex2 in range(numpy_graph[id2].shape[0]):
                                                if numpy_graph[id2][gindex2] == 0: break
                                            if gindex1 < 29 and gindex2 < 29:
                                                for k in range(numpy_graph[id1].shape[0]):
                                                    if numpy_graph[id1][k] == id2: break
                                                if k+1 == numpy_graph[id1].shape[0]:
                                                    numpy_graph[id1][gindex1] = id2
                                                for l in range(numpy_graph[id2].shape[0]):
                                                    if numpy_graph[id2][l] == id1: break
                                                if l+1 == numpy_graph[id2].shape[0]:
                                                    numpy_graph[id2][gindex1] = id1
                                else: # else find periodic connectivty
                                    for a in range(scaled_images.shape[:1][0]):
                                        ixlx, iyly, izlz = scaled_images[a]
                                        x1i = x1+ixlx; y1i = y1+iyly; z1i = z1+izlz;
                                        if abs(x1i-x2) <= max_voxel_size and abs(y1i-y2) <= max_voxel_size and abs(z1i-z2) <= max_voxel_size:
                                            dx = (x1i - x2)**2
                                            dy = (y1i - y2)**2
                                            dz = (z1i - z2)**2
                                            distance = math.sqrt(dx + dy + dz)
                                            if distance <= max_voxel_size:
                                                id1 = i+1; id2 = j+1;
                                                for gindex1 in range(numpy_graph[id1].shape[0]):
                                                    if numpy_graph[id1][gindex1] == 0: break
                                                for gindex2 in range(numpy_graph[id2].shape[0]):
                                                    if numpy_graph[id2][gindex2] == 0: break
                                                if gindex1 < 29 and gindex2 < 29:
                                                    for k in range(numpy_graph[id1].shape[0]):
                                                        if numpy_graph[id1][k] == id2: break
                                                    if k+1 == numpy_graph[id1].shape[0]:
                                                        numpy_graph[id1][gindex1] = id2   
                                                    for l in range(numpy_graph[id2].shape[0]):
                                                        if numpy_graph[id2][l] == id1: break
                                                    if l+1 == numpy_graph[id2].shape[0]:
                                                        numpy_graph[id2][gindex1] = id1
                                                break
                        internal_flags[i] = 0
                
                # Find voxel connectivty
                numpy_graph = cuda.to_device(np.zeros(30*(nvoxels+1), dtype=np.int64).reshape((nvoxels+1, 30))) # should be 27, but 30 allows overflow
                internal_flags = np.copy(flags);
                log.out(f'Finding free volume voxelID connectivity ({run_mode}-{threads_per_block_voxels}) ...')
                scaled_images = np.array(v.scaled_images)
                free_dist_CUDA[blocks_per_grid_voxels, threads_per_block_voxels](cuda.to_device(voxels), cuda.to_device(radial_voxel), cuda.to_device(scaled_images), cuda.to_device(internal_flags), cuda.to_device(voxel_periodic_flags), numpy_graph, max_voxel_size)
                log.out('Copying graph from GPU back to host ...')
                numpy_graph = numpy_graph.copy_to_host()
                
            # Generate graph for CPU and pure python
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