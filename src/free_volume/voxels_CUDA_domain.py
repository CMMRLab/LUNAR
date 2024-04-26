# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
January 1st, 2024
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
                
        #-----------------------#
        # Generate voxels step1 #
        #-----------------------#
        # Compute nx, ny, nz number of voxels required to build 3D voxel grid
        self.nx = math.ceil(self.lx/max_voxel_size)
        self.ny = math.ceil(self.ly/max_voxel_size)
        self.nz = math.ceil(self.lz/max_voxel_size)
        self.dx = self.lx/self.nx; self.dy = self.ly/self.ny; self.dz = self.lz/self.nz;
        self.halfdx = self.dx/2; self.halfdy = self.dy/2; self.halfdz = self.dz/2
        
        # Set probe_diameter if shortcut flag exists
        self.probe_diameter = probe_diameter
        if probe_diameter == 'min-voxel':
            log.out("probe_diameter = 'min-voxel', finding minimum dimension of voxel ...")
            self.probe_diameter = min([self.dx, self.dy, self.dz])
            log.out(f'   updated probe_diameter = {self.probe_diameter}')
        
        # Find domain decomposition grid size
        self.domain_size = round(2.1*max(system_vdw_radii), 2) + 2*max_voxel_size + self.probe_diameter
        self.min_cell = 3*self.domain_size
        
        # Exiting conditions (for domain decomposition)
        if self.lx < self.min_cell or self.ly < self.min_cell or self.lz < self.min_cell:
            log.out('ERROR simulation cell dimensions to small to use CUDA-dd run_mode. Minimum');
            log.out(f'simulation cell dimensions for numba-dd or numba-ddp is {self.min_cell}x{self.min_cell}x{self.min_cell} angstroms. Current')
            log.error(f'simulation cell dimensions {self.lx}x{self.ly}x{self.ly}. Please use other run_mode.')
        
        # Find images and scaled images 
        self.images, self.pflags = misc_func.generate_iflags(m, boundary, log)
        self.scaled_images = [(ix*self.lx, iy*self.ly, iz*self.lz) for (ix, iy, iz) in self.images]
        np_scaled_images = np.array(self.scaled_images)
        

        
        # Generatre voxels on CPU or GPU based on inputs
        start_time = time.time()
        if m.CUDA_threads_per_block_voxels == 0:# or m.CUDA_threads_per_block_voxels >= 0:
            self.numpy_voxels, self.numpy_radial, self.numpy_flags, self.numpyx, self.numpyy, self.numpyz, self.domain, self.domain_graph, self.nxx, self.nyy, self.nzz, ndomains, self.voxel_domains = generate_serial(self.nx, self.ny, self.nz, self.dx, self.dy, self.dz, self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi, self.cx, self.cy, self.cz, self.lx, self.ly, self.lz, self.domain_size, np_scaled_images, run_mode)
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
            def radial_CUDA(numpy_voxels, domain, numpy_radial, voxel_domains, c_xyz):
                i = cuda.grid(1)
                if i < numpy_voxels.shape[:1][0]:
                    cuda.syncthreads()
                    diffx = c_xyz[0] - numpy_voxels[i][0]
                    diffy = c_xyz[1] - numpy_voxels[i][1]
                    diffz = c_xyz[2] - numpy_voxels[i][2]
                    numpy_radial[i] = ( diffx*diffx + diffy*diffy + diffz*diffz )**0.5
                    for domainID in range(domain.shape[:1][0]):
                        if domain[domainID][0] <= numpy_voxels[i][0] <= domain[domainID][1] and domain[domainID][2] <= numpy_voxels[i][1] <= domain[domainID][3] and domain[domainID][4] <= numpy_voxels[i][2] <= domain[domainID][5]:
                            voxel_domains[i] = domainID; break
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
            voxel_domains = cuda.to_device(np.zeros(nvoxels, dtype=np.int64))
            
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
            
            # Generate domain and find its connectivty using the CPU
            self.domain, self.domain_graph, self.nxx, self.nyy, self.nzz, ndomains = generate_domain(self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi, self.cx, self.cy, self.cz, self.lx, self.ly, self.lz, self.domain_size, np_scaled_images, run_mode)
            
            # Finding radial voxels
            threads_per_block_voxels = m.CUDA_threads_per_block_voxels
            log.out(f'Finding radius of voxels from center and assigning voxel to domains ({run_mode}-{threads_per_block_voxels}) ...')
            blocks_per_grid_voxels = math.ceil(nvoxels / threads_per_block_voxels) 
            radial_CUDA[blocks_per_grid_voxels, threads_per_block_voxels](cuda.to_device(self.numpy_voxels), cuda.to_device(self.domain), numpy_radial, voxel_domains, c_xyz)
            self.numpy_radial = numpy_radial.copy_to_host();
            self.voxel_domains = voxel_domains.copy_to_host();
            self.numpy_flags = np.ones(nvoxels)
            self.numpyx, self.numpyy, self.numpyz = unique_xyz(self.numpy_voxels)
            
        # Generate required pytonh voxelIDs dict
        self.voxelIDs = dict(enumerate(self.numpy_voxels, 1))  # { voxelID : (x, y, z) }
        execution_time = (time.time() - start_time)
        log.out(f'Voxel generation execution time: {execution_time} (seconds)')
        
        # Generate linked list of domains and voxels in domain
        start_time = time.time()
        log.out('  Generating linked list (serial)')
        self.linked_lst = [set() for i in range(ndomains)]
        nvoxels_domains = self.voxel_domains.shape[:1][0]
        count = 0; progress_increment = 10;
        for voxelID in range(nvoxels_domains):
            domainID = self.voxel_domains[voxelID]
            self.linked_lst[domainID].add(voxelID)
            count += 1
            if 100*count/nvoxels_domains % progress_increment == 0:
                log.out(f'    progress: {int(100*count/nvoxels_domains)} %')
                
        # Making linked list homogeneous
        log.out('  Making linked list homogeneous (serial)')
        nmax = len(max(self.linked_lst, key=len))+1; # need an extra zero to modify in GPU run
        count = 0; progress_increment = 10;
        self.linked_lst = [list(i) for i in self.linked_lst]
        for n, i in enumerate(self.linked_lst):
            diff_zeros = abs(len(i)-nmax)*[0]
            self.linked_lst[n].extend(diff_zeros)
            count += 1
            if 100*count/ndomains % progress_increment == 0:
                print('    progress: ', int(100*count/ndomains),'%')
        self.linked_lst = np.array(self.linked_lst)
        execution_time = (time.time() - start_time)
        log.out(f'Linked list generation execution time: {execution_time} (seconds)')
        
        # Set atom domain
        for i in m.atoms:
            atom = m.atoms[i]
            x = atom.x; y = atom.y; z = atom.z;
            atom.domainID = assign_domain(x, y, z, self.domain)
            
        
        
#######################################################
# numba njit compilation for assigning atom to domain #
#######################################################
@numba.njit(parallel=False)
def assign_domain(xpos, ypos, zpos, domain):
    try:
        d = np.where( np.logical_and(xpos >= domain[:,0], xpos <= domain[:,1]) & np.logical_and(ypos >= domain[:,2], ypos <= domain[:,3])  & np.logical_and(zpos >= domain[:,4], zpos <= domain[:,5]) )[0]
        domainID = int(np.min(d))
    except:
        domainID = 0;
        for i in range(domain.shape[:1][0]):
            d = domain[i]
            if d[0] <= xpos <= d[1] and d[2] <= ypos <= d[3] and d[4] <= zpos <= d[5]:
                domainID = i; break;
    return domainID


##########################################################
# numba njit compilation for unique x, y and z locations #
##########################################################
@numba.njit(parallel=False)
def unique_xyz(numpy_voxels):
    numpyx = np.unique(numpy_voxels[:,0])
    numpyy = np.unique(numpy_voxels[:,1])
    numpyz = np.unique(numpy_voxels[:,2])
    return numpyx, numpyy, numpyz

##################################################
# numba njit compilitation for generating domain #
##################################################
@numba.njit(parallel=False)
def generate_domain(xlo, xhi, ylo, yhi, zlo, zhi, cx, cy, cz, lx, ly, lz, domain_size, scaled_images, run_mode):
    # Find subdomain regions
    nxx = math.ceil(lx/domain_size)
    nyy = math.ceil(ly/domain_size)
    nzz = math.ceil(lz/domain_size)
    dxx = lx/nxx; dyy = ly/nyy; dzz = lz/nzz;
    halfdxx = dxx/2; halfdyy = dyy/2; halfdzz = dzz/2;
    xxadd = halfdxx + xlo; yyadd =  halfdyy + ylo; zzadd = halfdzz + zlo;
    ndomains = nxx*nyy*nzz; count = 0; 
    domain = np.zeros(10*ndomains).reshape((ndomains, 10)) #  ( (xlo, xhi, ylo, yhi, zlo, zhi, xc, yc, zc, r) )
    print('  Generating Domains (serial)')
    for ii in range(nxx):
        for jj in range(nyy):
            for kk in range(nzz):
                xc = ii*dxx + xxadd
                yc = jj*dyy + yyadd
                zc = kk*dzz + zzadd
                diffx = cx - xc
                diffy = cy - yc
                diffz = cz - zc
                domain[count][0] = xc - halfdxx
                domain[count][1] = xc + halfdxx
                domain[count][2] = yc - halfdyy
                domain[count][3] = yc + halfdyy
                domain[count][4] = zc - halfdzz
                domain[count][5] = zc + halfdzz
                domain[count][6] = xc
                domain[count][7] = yc
                domain[count][8] = zc
                domain[count][9] = np.sqrt(diffx*diffx + diffy*diffy + diffz*diffz)
                count += 1
                
    # Generate domain graph
    domain_graph = domain_connectivty_serial(ndomains, domain, domain_size, scaled_images)
    return domain, domain_graph, nxx, nyy, nzz, ndomains


###################################################################################
# numba njit compilation for generating voxels and creating linked list in serial #
###################################################################################
@numba.njit(parallel=False)
def generate_serial(nx, ny, nz, dx, dy, dz, xlo, xhi, ylo, yhi, zlo, zhi, cx, cy, cz, lx, ly, lz, domain_size, scaled_images, run_mode):
    print('\n\nGenerating voxels (compiled - serial) ...')
    # Find subdomain regions
    nxx = math.ceil(lx/domain_size)
    nyy = math.ceil(ly/domain_size)
    nzz = math.ceil(lz/domain_size)
    dxx = lx/nxx; dyy = ly/nyy; dzz = lz/nzz;
    halfdxx = dxx/2; halfdyy = dyy/2; halfdzz = dzz/2;
    xxadd = halfdxx + xlo; yyadd =  halfdyy + ylo; zzadd = halfdzz + zlo;
    ndomains = nxx*nyy*nzz; count = 0; 
    domain = np.zeros(10*ndomains).reshape((ndomains, 10)) #  ( (xlo, xhi, ylo, yhi, zlo, zhi, xc, yc, zc, r) )
    print('  Generating Domains (serial)')
    for ii in range(nxx):
        for jj in range(nyy):
            for kk in range(nzz):
                xc = ii*dxx + xxadd
                yc = jj*dyy + yyadd
                zc = kk*dzz + zzadd
                diffx = cx - xc
                diffy = cy - yc
                diffz = cz - zc
                domain[count][0] = xc - halfdxx
                domain[count][1] = xc + halfdxx
                domain[count][2] = yc - halfdyy
                domain[count][3] = yc + halfdyy
                domain[count][4] = zc - halfdzz
                domain[count][5] = zc + halfdzz
                domain[count][6] = xc
                domain[count][7] = yc
                domain[count][8] = zc
                domain[count][9] = np.sqrt(diffx*diffx + diffy*diffy + diffz*diffz)
                count += 1

    # Generate voxels
    nvoxels = nx*ny*nz; progress_increment = 10; count = 0;
    xadd = dx/2 + xlo; yadd =  dy/2 + ylo; zadd = dz/2 + zlo;
    numpy_voxels = np.zeros(3*nvoxels).reshape((nvoxels, 3))
    numpy_radial = np.zeros(nvoxels); numpy_flags = np.ones(nvoxels)
    voxel_domains = np.zeros(nvoxels, dtype=np.int64);
    domain_guess = 0
    print('  Generating voxels and assigning to domain (serial)')
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                xpos = i*dx + xadd
                ypos = j*dy + yadd
                zpos = k*dz + zadd
                diffx = cx - xpos
                diffy = cy - ypos
                diffz = cz - zpos
                
                # Use previous voxel domain to guess current voxel domain (speed-up)
                d = domain[domain_guess]
                if d[0] <= xpos <= d[1] and d[2] <= ypos <= d[3] and d[4] <= zpos <= d[5]:
                    domainID = domain_guess
                else:
                    domainID = assign_domain(xpos, ypos, zpos, domain)
                    domain_guess = domainID
                numpy_radial[count] = np.sqrt(diffx*diffx + diffy*diffy + diffz*diffz)                 
                numpy_voxels[count][0] = xpos
                numpy_voxels[count][1] = ypos
                numpy_voxels[count][2] = zpos
                voxel_domains[count] = domainID
                count += 1
                if 100*count/nvoxels % progress_increment == 0:
                    print('    progress: ', int(100*count/nvoxels),'%')
                    
    # Find unique x, y, and z locations
    numpyx = np.unique(numpy_voxels[:,0])
    numpyy = np.unique(numpy_voxels[:,1])
    numpyz = np.unique(numpy_voxels[:,2])
    
    # Generate domain graph
    domain_graph = domain_connectivty_serial(ndomains, domain, domain_size, scaled_images)
    return numpy_voxels, numpy_radial, numpy_flags, numpyx, numpyy, numpyz, domain, domain_graph, nxx, nyy, nzz, ndomains, voxel_domains


#####################################################################################
# numba njit compilation for generating voxels and creating linked list in parallel #
#####################################################################################
@numba.njit(parallel=True)
def generate_parallel(nx, ny, nz, dx, dy, dz, xlo, xhi, ylo, yhi, zlo, zhi, cx, cy, cz, lx, ly, lz, domain_size, scaled_images, run_mode):
    print('\n\nGenerating voxels (compiled - parallel) ...')
    # Find subdomain regions
    nxx = math.ceil(lx/domain_size)
    nyy = math.ceil(ly/domain_size)
    nzz = math.ceil(lz/domain_size)
    dxx = lx/nxx; dyy = ly/nyy; dzz = lz/nzz;
    halfdxx = dxx/2; halfdyy = dyy/2; halfdzz = dzz/2;
    xxadd = halfdxx + xlo; yyadd =  halfdyy + ylo; zzadd = halfdzz + zlo;
    ndomains = nxx*nyy*nzz; count = 0; 
    domain = np.zeros(10*ndomains).reshape((ndomains, 10)) #  ( (xlo, xhi, ylo, yhi, zlo, zhi, xc, yc, zc, r) )
    print('  Generating Domains (serial)')
    for ii in range(nxx):
        for jj in range(nyy):
            for kk in range(nzz):
                xc = ii*dxx + xxadd
                yc = jj*dyy + yyadd
                zc = kk*dzz + zzadd
                diffx = cx - xc
                diffy = cy - yc
                diffz = cz - zc
                domain[count][0] = xc - halfdxx
                domain[count][1] = xc + halfdxx
                domain[count][2] = yc - halfdyy
                domain[count][3] = yc + halfdyy
                domain[count][4] = zc - halfdzz
                domain[count][5] = zc + halfdzz
                domain[count][6] = xc
                domain[count][7] = yc
                domain[count][8] = zc
                domain[count][9] = np.sqrt(diffx*diffx + diffy*diffy + diffz*diffz)
                count += 1

    # Generate voxels
    nvoxels = nx*ny*nz;
    xadd = dx/2 + xlo; yadd =  dy/2 + ylo; zadd = dz/2 + zlo;
    numpy_voxels = np.zeros(3*nvoxels).reshape((nvoxels, 3))
    numpy_radial = np.zeros(nvoxels); numpy_flags = np.ones(nvoxels)
    voxel_domains = np.zeros(nvoxels, dtype=np.int64);
    domain_guess = 0
    print('  Generating voxels and assigning to domain (parallel)')
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                xpos = i*dx + xadd
                ypos = j*dy + yadd
                zpos = k*dz + zadd
                diffx = cx - xpos
                diffy = cy - ypos
                diffz = cz - zpos
                
                # Use previous voxel domain to guess current voxel domain (speed-up)
                d = domain[domain_guess]
                if d[0] <= xpos <= d[1] and d[2] <= ypos <= d[3] and d[4] <= zpos <= d[5]:
                    domainID = domain_guess
                else:
                    domainID = assign_domain(xpos, ypos, zpos, domain)
                    domain_guess = domainID
                numpy_radial[i*ny*nz+j*nz+k] = np.sqrt(diffx*diffx + diffy*diffy + diffz*diffz)                 
                numpy_voxels[i*ny*nz+j*nz+k][0] = xpos
                numpy_voxels[i*ny*nz+j*nz+k][1] = ypos
                numpy_voxels[i*ny*nz+j*nz+k][2] = zpos
                voxel_domains[i*ny*nz+j*nz+k] = domainID
                    
    # Find unique x, y, and z locations
    numpyx = np.unique(numpy_voxels[:,0])
    numpyy = np.unique(numpy_voxels[:,1])
    numpyz = np.unique(numpy_voxels[:,2])
    
    # Generate domain graph
    domain_graph = domain_connectivty_serial(ndomains, domain, domain_size, scaled_images)
    return numpy_voxels, numpy_radial, numpy_flags, numpyx, numpyy, numpyz, domain, domain_graph, nxx, nyy, nzz, ndomains, voxel_domains


###################################################################
# numba njit compilation for finding domain connectivty in serial #
###################################################################
@numba.njit(parallel=False)
def domain_connectivty_serial(ndomains, domain, domain_size, scaled_images):
    domain_graph = np.zeros(30*ndomains, dtype=np.int64).reshape((ndomains, 30))
    internal_flags = np.ones(ndomains)
    print('  Finding Domain connectivity (serial)')
    for i in range(ndomains):
        x1 = domain[i][6]; y1 = domain[i][7];
        z1 = domain[i][8]; r1 = domain[i][9];
        min_radius = r1 - 2*domain_size
        max_radius = r1 + 2*domain_size
        periodic_postions = np.zeros(27*3).reshape((27, 3))
        for a in range(27):
            ixlx, iyly, izlz = scaled_images[a]
            periodic_postions[a][0] = x1+ixlx
            periodic_postions[a][1] = y1+iyly
            periodic_postions[a][2] = z1+izlz
        for j in range(i, ndomains):
            if internal_flags[j] == 0 or i == j: continue
            x2 = domain[j][6]; y2 = domain[j][7];
            z2 = domain[j][8]; r2 = domain[j][9];
            if min_radius < r2 < max_radius:
                for x1i, y1i, z1i in periodic_postions:
                    if abs(x1i-x2) > domain_size: continue
                    elif abs(y1i-y2) > domain_size: continue
                    elif abs(z1i-z2) > domain_size: continue
                    gindex1 = np.min(np.where(domain_graph[i] == 0)[0])
                    gindex2 = np.min(np.where(domain_graph[j] == 0)[0])
                    if gindex1 < 29 and gindex2 < 29:
                        if j not in domain_graph[i]:
                            domain_graph[i][gindex1] = j
                        if i not in domain_graph[j]:
                            domain_graph[j][gindex2] = i
                    break
        internal_flags[i] = 0
    return domain_graph
                    

######################################################################
# numba njit compilation for finding atom/free volume in serial mode #
######################################################################
@numba.njit(parallel=False)
def atom_vol_calc_serial(atoms, periodic_atoms, radial_voxel, voxels, flags, max_voxel_size, domain_graph, probe_radius, linked_lst):
    print('Finding atom volumes and free volumes (compiled - domain decomposition serial) ...')
    progress_increment = 10; checked = np.zeros(len(atoms)); natoms = atoms.shape[:1][0];
    pbc_count = 0; npossible_pbc = 0; voxel_atomIndexs = np.zeros(voxels.shape[:1])
    for i in range(natoms):
        atom = atoms[i]; checked[i] = 1;
        x1 = atom[0]; y1 = atom[1]; z1 = atom[2]; radial_atom = atom[3];
        atompos = np.array([x1, y1, z1])
        vdw_radii = atom[4] + probe_radius; atom_edgeflag = atom[5];
        if atom_edgeflag == 1: npossible_pbc += 1
        domainID = int(atom[6])
        tmp = np.unique(domain_graph[domainID])
        atom_domains = np.append(tmp, domainID)
        min_radius = radial_atom - vdw_radii - max_voxel_size
        max_radius = radial_atom + vdw_radii + max_voxel_size
        if 100*np.sum(checked)/natoms % progress_increment == 0:
            print('progress: ', int(100*np.sum(checked)/natoms),'%')
        for domainID in atom_domains:
            voxelIDs = linked_lst[domainID]
            for j in voxelIDs:
                if flags[j] == 0: continue
                if min_radius < radial_voxel[j] < max_radius:
                    if atom_edgeflag == 1:
                        voxel = voxels[j]
                        for x1i, y1i, z1i in periodic_atoms[i]:
                            if abs(voxel[0] - x1i) > vdw_radii: continue
                            elif abs(voxel[1] - y1i) > vdw_radii: continue
                            elif abs(voxel[2] - z1i) > vdw_radii: continue
                            atompos = np.array([x1i, y1i, z1i])
                            distance_ppp = np.linalg.norm(atompos-voxels[j])
                            if distance_ppp < vdw_radii:
                                flags[j] = 0; voxel_atomIndexs[j] = atom[7]; pbc_count += 1;
                    else:
                        voxel = voxels[j]
                        if abs(voxel[0] - x1) > vdw_radii: continue
                        elif abs(voxel[1] - y1) > vdw_radii: continue
                        elif abs(voxel[2] - z1) > vdw_radii: continue
                        distance_fff = np.linalg.norm(atompos-voxels[j])
                        if distance_fff < vdw_radii:
                            flags[j] = 0; voxel_atomIndexs[j] = atom[7];
    return pbc_count, npossible_pbc, voxel_atomIndexs


########################################################################
# numba njit compilation for finding atom/free volume in parallel mode #
########################################################################
@numba.njit(parallel=True)
def atom_vol_calc_parallel(atoms, periodic_atoms, radial_voxel, voxels, flags, max_voxel_size, domain_graph, probe_radius, linked_lst):
    print('Finding atom volumes and free volumes (compiled - domain decomposition parallel) ...')
    progress_increment = 10; checked = np.zeros(len(atoms)); natoms = atoms.shape[:1][0];
    pbc_count = 0; npossible_pbc = 0; voxel_atomIndexs = np.zeros(voxels.shape[:1])
    for i in numba.prange(natoms):
        atom = atoms[i]; checked[i] = 1;
        x1 = atom[0]; y1 = atom[1]; z1 = atom[2]; radial_atom = atom[3];
        atompos = np.array([x1, y1, z1])
        vdw_radii = atom[4] + probe_radius; atom_edgeflag = atom[5];
        if atom_edgeflag == 1: npossible_pbc += 1
        domainID = int(atom[6])
        tmp = np.unique(domain_graph[domainID])
        atom_domains = np.append(tmp, domainID)
        min_radius = radial_atom - vdw_radii - max_voxel_size
        max_radius = radial_atom + vdw_radii + max_voxel_size
        if 100*np.sum(checked)/natoms % progress_increment == 0:
            print('progress: ', int(100*np.sum(checked)/natoms),'%')
        for domainID in atom_domains:
            voxelIDs = linked_lst[domainID]
            for j in voxelIDs:
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
                                flags[j] = 0; voxel_atomIndexs[j] = atom[7]; pbc_count += 1;
                    else:
                        voxel = voxels[j]
                        if abs(voxel[0] - x1) > vdw_radii: continue
                        elif abs(voxel[1] - y1) > vdw_radii: continue
                        elif abs(voxel[2] - z1) > vdw_radii: continue
                        distance_fff = np.sqrt(np.sum(np.square(atompos-voxels[j])))
                        if distance_fff < vdw_radii:
                            flags[j] = 0; voxel_atomIndexs[j] = atom[7];
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
def free_dist_serial(voxels, flags, max_voxel_size, xlo, xhi, ylo, yhi, zlo, zhi, scaled_images, radial_voxel, voxel_domains, domain_graph, linked_lst):
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
        domainID = voxel_domains[i]
        tmp = np.unique(domain_graph[domainID])
        domains = np.append(tmp, domainID)
        periodic_postions = np.zeros(nimages*3).reshape((nimages, 3))
        for a in range(nimages):
            ixlx, iyly, izlz = scaled_images[a]
            periodic_postions[a][0] = x1+ixlx
            periodic_postions[a][1] = y1+iyly
            periodic_postions[a][2] = z1+izlz
        for domainID in domains:
            voxelIDs = linked_lst[domainID]
            for j in voxelIDs:
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
def free_dist_parallel(voxels, flags, max_voxel_size, xlo, xhi, ylo, yhi, zlo, zhi, scaled_images, radial_voxel, voxel_domains, domain_graph, linked_lst):
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
        domainID = voxel_domains[i]
        tmp = np.unique(domain_graph[domainID])
        domains = np.append(tmp, domainID)
        periodic_postions = np.zeros(nimages*3).reshape((nimages, 3))
        for a in range(nimages):
            ixlx, iyly, izlz = scaled_images[a]
            periodic_postions[a][0] = x1+ixlx
            periodic_postions[a][1] = y1+iyly
            periodic_postions[a][2] = z1+izlz
        for domainID in domains:
            voxelIDs = linked_lst[domainID]
            for j in voxelIDs:
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
                                    #if id2 not in numpy_graph[id1]:
                                        numpy_graph[id1][gindex1] = id2
                                    if not np.any(numpy_graph[id2] == id1):
                                    #if id1 not in numpy_graph[id2]:
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
        atoms = [] # [x, y, z, r, vdw, pbcflag, domainID]; pbcflag = 0 (not near edge); pbcflag = 1 (near edge);  
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
            atoms.append([x, y, z, radial_atom, vdw_radii, pbcflag, atom.domainID, atom.type])
        atoms = np.array(atoms); periodic_atoms = np.array(periodic_atoms)
        
        #-----------------------------------------#
        # CUDA Kernel to compute atom/free volume #
        #-----------------------------------------#
        start_time = time.time()
        @cuda.jit
        def atom_vol_calc_CUDA(atoms, periodic_atoms, radial_voxel, voxels, domain_graph, linked_lst, flags, voxel_atomIndexs, sizes, pbc):
            i = cuda.grid(1)
            if i < atoms.shape[:1][0]:
                atom = atoms[i]
                x1 = atom[0]; y1 = atom[1]; z1 = atom[2]; radial_atom = atom[3];
                vdw_radii = atom[4] + sizes[1]; atom_edgeflag = atom[5];
                min_radius = radial_atom - vdw_radii - sizes[0]
                max_radius = radial_atom + vdw_radii + sizes[0]
                domainID = int(atom[6])
                atom_domains = domain_graph[domainID]
                atom_domains[-1] = domainID # extra zero is available to do this
                if atom_edgeflag == 1: cuda.atomic.add(pbc, 1, 1)
                for domainID in atom_domains:
                    voxelIDs = linked_lst[domainID]
                    for j in voxelIDs:
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
                                        flags[j] = 0; voxel_atomIndexs[j] = int(atom[7]);
                                        cuda.atomic.add(pbc, 0, 1)
                            else:
                                voxel = voxels[j]
                                if abs(voxel[0] - x1) > vdw_radii: continue
                                elif abs(voxel[1] - y1) > vdw_radii: continue
                                elif abs(voxel[2] - z1) > vdw_radii: continue
                                dx = atom[0] - voxel[0]; dy = atom[1] - voxel[1]; dz = atom[2] - voxel[2];
                                distance_fff = math.sqrt(dx*dx + dy*dy + dz*dz)
                                if distance_fff < vdw_radii:
                                    flags[j] = 0; voxel_atomIndexs[j] = int(atom[7]);
                                
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
        atom_vol_calc_CUDA[blocks_per_grid_atoms, threads_per_block_atoms](cuda.to_device(atoms), cuda.to_device(periodic_atoms), cuda.to_device(radial_voxel), cuda.to_device(voxels), cuda.to_device(v.domain_graph), cuda.to_device(v.linked_lst), flags, voxel_atomIndexs, sizes, pbc)
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
                numpy_graph = free_dist_parallel(voxels, flags, max_voxel_size, v.xlo, v.xhi, v.ylo, v.yhi, v.zlo, v.zhi, np.array(v.scaled_images), radial_voxel, v.voxel_domains, v.domain_graph, v.linked_lst)
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
                def free_dist_CUDA(voxels, radial_voxel, scaled_images, internal_flags, voxel_periodic_flags, voxel_domains, domain_graph, linked_lst, numpy_graph, max_voxel_size):
                    i = cuda.grid(1)
                    if i < voxels.shape[:1][0] and internal_flags[i] != 0:
                        x1, y1, z1 = voxels[i]
                        r1 = radial_voxel[i]
                        min_radius = r1 - 2*max_voxel_size
                        max_radius = r1 + 2*max_voxel_size  
                        domain_graph[voxel_domains[i]][-1] = voxel_domains[i] # extra zero is available to do this
                        for domainID in domain_graph[voxel_domains[i]]:
                            for j in linked_lst[domainID]:
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
                free_dist_CUDA[blocks_per_grid_voxels, threads_per_block_voxels](cuda.to_device(voxels), cuda.to_device(radial_voxel), cuda.to_device(scaled_images), cuda.to_device(internal_flags), cuda.to_device(voxel_periodic_flags), cuda.to_device(v.voxel_domains), cuda.to_device(v.domain_graph), cuda.to_device(v.linked_lst), numpy_graph, max_voxel_size)
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