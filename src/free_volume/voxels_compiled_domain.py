# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
June 10, 2026
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


###############################################################
# Local sparse triclinic helpers compatible with numba njit   #
# h = [lx, ly, lz, yz, xz, xy]                                #
###############################################################
def get_box_parameters(m):    
    xline = m.xbox_line.split(); yline = m.ybox_line.split(); zline = m.zbox_line.split()
    xlo = float(xline[0]); xhi = float(xline[1])
    ylo = float(yline[0]); yhi = float(yline[1])
    zlo = float(zline[0]); zhi = float(zline[1])
    lx = xhi - xlo; ly = yhi - ylo; lz = zhi - zlo
    yz = m.yz; xz = m.xz; xy = m.xy

    if lx == 0: lx = 1.0
    if ly == 0: ly = 1.0
    if lz == 0: lz = 1.0

    h = np.array([lx, ly, lz, yz, xz, xy], dtype=np.float64)
    h_inv = np.zeros(6, dtype=np.float64)
    h_inv[0] = 1.0/h[0]
    h_inv[1] = 1.0/h[1]
    h_inv[2] = 1.0/h[2]
    h_inv[3] = -h[3] / (h[1]*h[2])
    h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2])
    h_inv[5] = -h[5] / (h[0]*h[1])
    boxlo = np.array([xlo, ylo, zlo], dtype=np.float64)
    boxhi = np.array([xhi, yhi, zhi], dtype=np.float64)
    tilts = np.array([xy, xz, yz], dtype=np.float64)
    return h, h_inv, boxlo, boxhi, tilts


@numba.njit(parallel=False)
def frac2pos_njit(fx, fy, fz, h, boxlo):
    x = h[0]*fx + h[5]*fy + h[4]*fz + boxlo[0]
    y = h[1]*fy + h[3]*fz + boxlo[1]
    z = h[2]*fz + boxlo[2]
    return x, y, z


@numba.njit(parallel=False)
def pos2frac_njit(x, y, z, h_inv, boxlo):
    dx = x - boxlo[0]
    dy = y - boxlo[1]
    dz = z - boxlo[2]
    fx = h_inv[0]*dx + h_inv[5]*dy + h_inv[4]*dz
    fy = h_inv[1]*dy + h_inv[3]*dz
    fz = h_inv[2]*dz
    return fx, fy, fz


@numba.njit(parallel=False)
def wrap_frac_scalar(f):
    return f - math.floor(f)

@numba.njit(parallel=False)
def frac_delta_to_cart_njit(dfx, dfy, dfz, h):
    # Convert a fractional displacement, not a position, to Cartesian.
    # Uses boxlo = 0 implicitly. h = [lx, ly, lz, yz, xz, xy].
    dx = h[0]*dfx + h[5]*dfy + h[4]*dfz
    dy = h[1]*dfy + h[3]*dfz
    dz = h[2]*dfz
    return dx, dy, dz

@numba.njit(parallel=False)
def mic_distance_frac_njit(fx1, fy1, fz1, fx2, fy2, fz2, h):
    # Minimum-image displacement in fractional space, then true Cartesian length.
    dfx = fx2 - fx1
    dfy = fy2 - fy1
    dfz = fz2 - fz1
    dfx -= round(dfx)
    dfy -= round(dfy)
    dfz -= round(dfz)
    dx, dy, dz = frac_delta_to_cart_njit(dfx, dfy, dfz, h)
    return math.sqrt(dx*dx + dy*dy + dz*dz)


@numba.njit(parallel=False)
def clamp_int(i, lo, hi):
    if i < lo:
        return lo
    if i > hi:
        return hi
    return i


@numba.njit(parallel=False)
def fractional_domain_id(fx, fy, fz, nxx, nyy, nzz):
    # wrap into [0, 1) for periodic-style binning
    fx = wrap_frac_scalar(fx)
    fy = wrap_frac_scalar(fy)
    fz = wrap_frac_scalar(fz)

    ii = int(math.floor(fx*nxx))
    jj = int(math.floor(fy*nyy))
    kk = int(math.floor(fz*nzz))
    ii = clamp_int(ii, 0, nxx-1)
    jj = clamp_int(jj, 0, nyy-1)
    kk = clamp_int(kk, 0, nzz-1)
    return ii*nyy*nzz + jj*nzz + kk


@numba.njit(parallel=False)
def find_zero_index(arr):
    for i in range(arr.shape[0]):
        if arr[i] == 0:
            return i
    return arr.shape[0] - 1


@numba.njit(parallel=False)
def add_graph_edge(graph, id1, id2):
    if id1 == id2:
        return
    if not np.any(graph[id1] == id2):
        g1 = find_zero_index(graph[id1])
        if g1 < graph.shape[1] - 1:
            graph[id1][g1] = id2
    if not np.any(graph[id2] == id1):
        g2 = find_zero_index(graph[id2])
        if g2 < graph.shape[1] - 1:
            graph[id2][g2] = id1


@numba.njit(parallel=False)
def build_fractional_domain_graph(nxx, nyy, nzz, pflags):
    # Same shape convention as old code: 0 means empty. Domain id 0 is valid,
    # but the old code also used zeros as empty slots. This preserves behavior,
    # but downstream code should eventually move to -1 sentinels.
    ndomains = nxx*nyy*nzz
    graph = np.zeros((ndomains, 30), dtype=np.int64)

    px = pflags[0] == 1
    py = pflags[1] == 1
    pz = pflags[2] == 1

    for ii in range(nxx):
        for jj in range(nyy):
            for kk in range(nzz):
                id1 = ii*nyy*nzz + jj*nzz + kk
                for di in range(-1, 2):
                    ni = ii + di
                    if ni < 0 or ni >= nxx:
                        if px:
                            ni = ni % nxx
                        else:
                            continue
                    for dj in range(-1, 2):
                        nj = jj + dj
                        if nj < 0 or nj >= nyy:
                            if py:
                                nj = nj % nyy
                            else:
                                continue
                        for dk in range(-1, 2):
                            nk = kk + dk
                            if nk < 0 or nk >= nzz:
                                if pz:
                                    nk = nk % nzz
                                else:
                                    continue
                            id2 = ni*nyy*nzz + nj*nzz + nk
                            add_graph_edge(graph, id1, id2)
    return graph


@numba.njit(parallel=False)
def generate_fractional_serial(nx, ny, nz, nxx, nyy, nzz, h, h_inv, boxlo, pflags):
    print('\n\nGenerating voxels (fractional/triclinic - serial) ...')
    nvoxels = nx*ny*nz
    ndomains = nxx*nyy*nzz

    numpy_voxels = np.zeros((nvoxels, 3), dtype=np.float64)
    numpy_frac_voxels = np.zeros((nvoxels, 3), dtype=np.float64)
    numpy_radial = np.zeros(nvoxels, dtype=np.float64)
    numpy_flags = np.ones(nvoxels, dtype=np.float64)
    voxel_domains = np.zeros(nvoxels, dtype=np.int64)

    # Domain array is retained mostly for compatibility/debugging.
    # Columns: flox, fhix, floy, fhiy, floz, fhiz, fcx, fcy, fcz, cart_r_from_cell_center
    domain = np.zeros((ndomains, 10), dtype=np.float64)
    count = 0
    for ii in range(nxx):
        for jj in range(nyy):
            for kk in range(nzz):
                flox = ii/nxx; fhix = (ii + 1)/nxx
                floy = jj/nyy; fhiy = (jj + 1)/nyy
                floz = kk/nzz; fhiz = (kk + 1)/nzz
                fcx = (ii + 0.5)/nxx
                fcy = (jj + 0.5)/nyy
                fcz = (kk + 0.5)/nzz
                x, y, z = frac2pos_njit(fcx, fcy, fcz, h, boxlo)
                cx, cy, cz = frac2pos_njit(0.5, 0.5, 0.5, h, boxlo)
                dx = x - cx; dy = y - cy; dz = z - cz
                domain[count][0] = flox; domain[count][1] = fhix
                domain[count][2] = floy; domain[count][3] = fhiy
                domain[count][4] = floz; domain[count][5] = fhiz
                domain[count][6] = fcx; domain[count][7] = fcy; domain[count][8] = fcz
                domain[count][9] = math.sqrt(dx*dx + dy*dy + dz*dz)
                count += 1

    print('  Generating voxel centers and assigning fractional domains')
    count = 0
    progress_increment = 10
    cx, cy, cz = frac2pos_njit(0.5, 0.5, 0.5, h, boxlo)
    for i in range(nx):
        fx = (i + 0.5)/nx
        for j in range(ny):
            fy = (j + 0.5)/ny
            for k in range(nz):
                fz = (k + 0.5)/nz
                x, y, z = frac2pos_njit(fx, fy, fz, h, boxlo)
                dx = x - cx; dy = y - cy; dz = z - cz
                numpy_voxels[count][0] = x
                numpy_voxels[count][1] = y
                numpy_voxels[count][2] = z
                numpy_frac_voxels[count][0] = fx
                numpy_frac_voxels[count][1] = fy
                numpy_frac_voxels[count][2] = fz
                numpy_radial[count] = math.sqrt(dx*dx + dy*dy + dz*dz)
                voxel_domains[count] = fractional_domain_id(fx, fy, fz, nxx, nyy, nzz)
                count += 1
                if 100*count/nvoxels % progress_increment == 0:
                    print('    progress: ', int(100*count/nvoxels), '%')

    domain_graph = build_fractional_domain_graph(nxx, nyy, nzz, pflags)
    return numpy_voxels, numpy_frac_voxels, numpy_radial, numpy_flags, domain, domain_graph, ndomains, voxel_domains


@numba.njit(parallel=False)
def assign_domain_fractional_from_cart(x, y, z, h_inv, boxlo, nxx, nyy, nzz):
    fx, fy, fz = pos2frac_njit(x, y, z, h_inv, boxlo)
    return fractional_domain_id(fx, fy, fz, nxx, nyy, nzz)


##########################################################
# numba njit compilation for unique x, y and z locations #
##########################################################
@numba.njit(parallel=False)
def unique_xyz(numpy_voxels):
    numpyx = np.unique(numpy_voxels[:,0])
    numpyy = np.unique(numpy_voxels[:,1])
    numpyz = np.unique(numpy_voxels[:,2])
    return numpyx, numpyy, numpyz


############################
# Class to generate voxels #
############################
class generate:
    def __init__(self, m, max_voxel_size, boundary, vdw_radius, probe_diameter, vdw_method, run_mode, log):        
        # Get box dimensions        
        self.h, self.h_inv, self.boxlo, self.boxhi, self.tilts = get_box_parameters(m)
        self.cx, self.cy, self.cz = frac2pos_njit(0.5, 0.5, 0.5, self.h, self.boxlo)
        self.xlo, self.ylo, self.zlo = self.boxlo
        self.xhi, self.yhi, self.zhi = self.boxhi
        self.lx = self.boxhi[0] - self.boxlo[0]
        self.ly = self.boxhi[1] - self.boxlo[1]
        self.lz = self.boxhi[2] - self.boxlo[2]
        
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
                if vdw_method == 'class1':
                    vdw_radii = ( (2**(1/6))*sigma )/2
                if vdw_method == 'class2':
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
        self.dx = self.lx/self.nx
        self.dy = self.ly/self.ny
        self.dz = self.lz/self.nz
        
        # Set probe_diameter if shortcut flag exists
        self.probe_diameter = probe_diameter
        if probe_diameter == 'min-voxel':
            log.out("probe_diameter = 'min-voxel', finding minimum dimension of voxel ...")
            self.probe_diameter = min([self.dx, self.dy, self.dz])
            log.out(f'   updated probe_diameter = {self.probe_diameter}')

    
        # Find images and scaled images 
        self.images, self.pflags = misc_func.generate_iflags(m, boundary, log)
        self.scaled_images = []
        for ix, iy, iz in self.images:
            sx, sy, sz = frac2pos_njit(ix, iy, iz, self.h, np.array([0.0, 0.0, 0.0]))
            self.scaled_images.append((sx, sy, sz))
        
        # Generate voxels (numba method)
        # pflags as numeric for numba: 1 periodic, 0 non-periodic.
        pflags_num = np.zeros(3, dtype=np.int64)
        for i, flag in enumerate(self.pflags):
            if flag == 'p': pflags_num[i] = 1
                
        # Find domain decomposition grid size
        self.domain_size = round(2.1*max(system_vdw_radii), 2) + 2*max_voxel_size + self.probe_diameter
        self.min_cell = 3*self.domain_size
        
        # Fractional domain counts.
        a_len = self.h[0]
        b_len = math.sqrt(self.h[5]**2 + self.h[1]**2)
        c_len = math.sqrt(self.h[4]**2 + self.h[3]**2 + self.h[2]**2)
        self.nxx = max(1, math.ceil(a_len/self.domain_size))
        self.nyy = max(1, math.ceil(b_len/self.domain_size))
        self.nzz = max(1, math.ceil(c_len/self.domain_size))
        
        # Exiting conditions (for domain decomposition)
        if a_len < self.min_cell or b_len < self.min_cell or c_len < self.min_cell:
            log.out('WARNING simulation cell is small relative to domain_size.')
            log.out('Using at least one fractional domain in each direction.')
        
        start_time = time.time()
        if run_mode in ['numba-dd', 'numba-ddp']:
            # Keep serial first; parallelizing this part is not the bottleneck compared with atom/voxel marking.
            out = generate_fractional_serial(self.nx, self.ny, self.nz, self.nxx, self.nyy, self.nzz, self.h, self.h_inv, self.boxlo, pflags_num)
            (self.numpy_voxels, self.numpy_frac_voxels, self.numpy_radial, self.numpy_flags, self.domain, self.domain_graph, ndomains, self.voxel_domains) = out
        else:
            log.error(f'ERROR input run_mode = {run_mode} not supported')
        
        self.numpyx, self.numpyy, self.numpyz = unique_xyz(self.numpy_voxels)
        self.voxelIDs = dict(enumerate(self.numpy_voxels, 1))
        log.out(f'Voxel generation execution time: {time.time() - start_time} (seconds)')
        
        # Same linked-list layout as old code.
        start_time = time.time()
        log.out('  Generating linked list (serial)')
        self.linked_lst = [set() for _ in range(ndomains)]
        for voxelID in range(self.voxel_domains.shape[0]):
            domainID = int(self.voxel_domains[voxelID])
            self.linked_lst[domainID].add(voxelID)
        
        log.out('  Making linked list homogeneous (serial)')
        nmax = len(max(self.linked_lst, key=len))
        self.linked_lst = [list(i) for i in self.linked_lst]
        for n, entries in enumerate(self.linked_lst):
            #entries.extend([0]*(nmax - len(entries)))
            entries.extend([-1]*(nmax - len(entries)))
        self.linked_lst = np.array(self.linked_lst, dtype=np.int64)
        log.out(f'Linked list generation execution time: {time.time() - start_time} (seconds)')
        
        # Set atom domain using fractional binning.
        for i in m.atoms:
            atom = m.atoms[i]
            atom.domainID = assign_domain_fractional_from_cart(
                atom.x, atom.y, atom.z, self.h_inv, self.boxlo, self.nxx, self.nyy, self.nzz)



################################################
# Function to convert minimum-image fractional #
# coordinate differences into a Cartesian      #
# displacement squared (r^2) for orthogonal    #
# or restricted triclinic simulation cells.    #
################################################
@numba.njit(parallel=False)
def frac_delta_to_cart_distance2(dfx, dfy, dfz, h):
    # Minimum-image fractional delta -> Cartesian displacement squared.
    # h = [lx, ly, lz, yz, xz, xy] for LAMMPS restricted triclinic.
    dfx = dfx - round(dfx)
    dfy = dfy - round(dfy)
    dfz = dfz - round(dfz)
    dx = h[0]*dfx + h[5]*dfy + h[4]*dfz
    dy = h[1]*dfy + h[3]*dfz
    dz = h[2]*dfz
    return dx*dx + dy*dy + dz*dz


######################################################################
# numba njit compilation for finding atom/free volume in serial mode #
######################################################################
@numba.njit(parallel=False)
def atom_vol_calc_serial_frac(atoms_frac, frac_voxels, flags, domain_graph, probe_radius, linked_lst, h):
    print('Finding atom volumes and free volumes (compiled - domain decomposition - fractional/triclinic serial) ...')
    natoms = atoms_frac.shape[0]
    nvoxels = frac_voxels.shape[0]
    voxel_atomIndexs = np.zeros(nvoxels)
    pbc_count = 0
    npossible_pbc = 0
    progress_increment = 10

    for i in range(natoms):
        if 100*(i + 1)/natoms % progress_increment == 0:
            print('progress: ', int(100*(i + 1)/natoms), '%')

        atom = atoms_frac[i]
        fax = atom[0]
        fay = atom[1]
        faz = atom[2]
        cutoff = atom[3] + probe_radius
        cutoff2 = cutoff*cutoff
        domainID = int(atom[4])
        atom_type = atom[5]

        tmp = np.unique(domain_graph[domainID])
        atom_domains = np.append(tmp, domainID)

        for didx in range(atom_domains.shape[0]):
            domainID2 = int(atom_domains[didx])
            voxelIDs = linked_lst[domainID2]
            for jj in range(voxelIDs.shape[0]):
                j = int(voxelIDs[jj])
                if j < 0:
                    continue
                if flags[j] == 0:
                    continue

                fv = frac_voxels[j]
                d2 = frac_delta_to_cart_distance2(fv[0] - fax, fv[1] - fay, fv[2] - faz, h)
                if d2 < cutoff2:
                    flags[j] = 0
                    voxel_atomIndexs[j] = atom_type
                    # Every periodic overlap is handled by MIC now, so this is not
                    # directly comparable to the old edge-image counter.
                    pbc_count += 1

    return pbc_count, npossible_pbc, voxel_atomIndexs


########################################################################
# numba njit compilation for finding atom/free volume in parallel mode #
########################################################################
@numba.njit(parallel=True)
def atom_vol_calc_parallel_frac(atoms_frac, frac_voxels, flags, domain_graph, probe_radius, linked_lst, h):
    print('Finding atom volumes and free volumes (compiled - domain decomposition - fractional/triclinic parallel) ...')
    natoms = atoms_frac.shape[0]
    nvoxels = frac_voxels.shape[0]
    voxel_atomIndexs = np.zeros(nvoxels)
    pbc_count = 0
    npossible_pbc = 0

    # NOTE: This preserves the same race-condition profile as the old parallel
    # routine: multiple atoms can try to claim the same voxel. Occupied/free
    # flags remain correct, but voxel_atomIndexs can be non-deterministic for
    # voxels close to overlapping atoms. Use the serial version if deterministic
    # atom ownership is required.
    for i in numba.prange(natoms):
        atom = atoms_frac[i]
        fax = atom[0]
        fay = atom[1]
        faz = atom[2]
        cutoff = atom[3] + probe_radius
        cutoff2 = cutoff*cutoff
        domainID = int(atom[4])
        atom_type = atom[5]

        tmp = np.unique(domain_graph[domainID])
        atom_domains = np.append(tmp, domainID)

        for didx in range(atom_domains.shape[0]):
            domainID2 = int(atom_domains[didx])
            voxelIDs = linked_lst[domainID2]
            for jj in range(voxelIDs.shape[0]):
                j = int(voxelIDs[jj])
                if j < 0:
                    continue
                if flags[j] == 0:
                    continue

                fv = frac_voxels[j]
                d2 = frac_delta_to_cart_distance2(fv[0] - fax, fv[1] - fay, fv[2] - faz, h)
                if d2 < cutoff2:
                    flags[j] = 0
                    voxel_atomIndexs[j] = atom_type
                    pbc_count += 1

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


##################################################################
# numba njit compilation for free volume distribution (parallel) #
##################################################################
@numba.njit(parallel=True)
def free_dist_parallel_frac_indexed(flags, nx, ny, nz, pflags):
    nvoxels = nx*ny*nz
    numpy_graph = np.zeros((nvoxels + 1, 30), dtype=np.int64)

    for idx in numba.prange(nvoxels):
        if flags[idx] == 0:
            continue

        i = idx // (ny*nz)
        rem = idx - i*ny*nz
        j = rem // nz
        k = rem - j*nz

        id1 = idx + 1

        for di in range(-1, 2):
            ni = i + di
            if ni < 0 or ni >= nx:
                if pflags[0] == 1:
                    ni = ni % nx
                else:
                    continue

            for dj in range(-1, 2):
                nj = j + dj
                if nj < 0 or nj >= ny:
                    if pflags[1] == 1:
                        nj = nj % ny
                    else:
                        continue

                for dk in range(-1, 2):
                    if di == 0 and dj == 0 and dk == 0:
                        continue

                    nk = k + dk
                    if nk < 0 or nk >= nz:
                        if pflags[2] == 1:
                            nk = nk % nz
                        else:
                            continue

                    idx2 = ni*ny*nz + nj*nz + nk

                    # write one direction only; Python graph build later symmetrizes
                    if idx2 <= idx:
                        continue

                    if flags[idx2] == 0:
                        continue

                    gindex = find_zero_index(numpy_graph[id1])
                    if gindex < 29:
                        numpy_graph[id1][gindex] = idx2 + 1

    return numpy_graph

################################################################
# numba njit compilation for free volume distribution (serial) #
################################################################
@numba.njit(parallel=False)
def free_dist_serial_frac_indexed(flags, nx, ny, nz, pflags):
    nvoxels = nx*ny*nz
    numpy_graph = np.zeros((nvoxels + 1, 30), dtype=np.int64)

    for idx in range(nvoxels):
        if flags[idx] == 0:
            continue

        i = idx // (ny*nz)
        rem = idx - i*ny*nz
        j = rem // nz
        k = rem - j*nz

        id1 = idx + 1

        for di in range(-1, 2):
            ni = i + di
            if ni < 0 or ni >= nx:
                if pflags[0] == 1:
                    ni = ni % nx
                else:
                    continue

            for dj in range(-1, 2):
                nj = j + dj
                if nj < 0 or nj >= ny:
                    if pflags[1] == 1:
                        nj = nj % ny
                    else:
                        continue

                for dk in range(-1, 2):
                    if di == 0 and dj == 0 and dk == 0:
                        continue

                    nk = k + dk
                    if nk < 0 or nk >= nz:
                        if pflags[2] == 1:
                            nk = nk % nz
                        else:
                            continue

                    idx2 = ni*ny*nz + nj*nz + nk

                    # write one direction only; Python graph build later symmetrizes
                    if idx2 <= idx:
                        continue

                    if flags[idx2] == 0:
                        continue

                    gindex = find_zero_index(numpy_graph[id1])
                    if gindex < 29:
                        numpy_graph[id1][gindex] = idx2 + 1

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
        voxels = v.numpy_voxels
        frac_voxels = v.numpy_frac_voxels
        flags = v.numpy_flags
        atoms_frac = [] # [frac_x, frac_y, frac_z, vdw, domainID, atom_type]
        for i in m.atoms:
            atom = m.atoms[i]
            x = atom.x; y = atom.y; z = atom.z
            
            # Set atom vdw radii based on method
            if vdw_method == 'dict':
                vdw_radii = vdw_radius[atom.element]
            elif vdw_method in ['class1', 'class2']:
                sigma = m.pair_coeffs[atom.type].coeffs[1]
                if vdw_method == 'class1':
                    vdw_radii = ( (2**(1/6))*sigma )/2
                if vdw_method == 'class2':
                    vdw_radii = sigma/2
            else: log.out(f'ERROR vdw_method {vdw_method} not supported')
            
            fx, fy, fz = pos2frac_njit(x, y, z, v.h_inv, v.boxlo)
            fx = fx - math.floor(fx)
            fy = fy - math.floor(fy)
            fz = fz - math.floor(fz)
            atoms_frac.append([fx, fy, fz, vdw_radii, atom.domainID, atom.type])
        atoms_frac = np.array(atoms_frac)
        
        
        # Used numpy/numba implementation to find atom volumes and free volumes
        start_time = time.time()
        if run_mode == 'numba-dd':
            self.pbc_count, self.npossible_pbc, voxel_atomIndexs = atom_vol_calc_serial_frac(atoms_frac, frac_voxels, flags, v.domain_graph, probe_diameter/2, v.linked_lst, v.h)
        elif run_mode == 'numba-ddp':
            self.pbc_count, self.npossible_pbc, voxel_atomIndexs = atom_vol_calc_parallel_frac(atoms_frac, frac_voxels, flags, v.domain_graph, probe_diameter/2, v.linked_lst, v.h)
        else: log.error(f'ERROR input run_mode = {run_mode} not supported')
        execution_time = (time.time() - start_time)
        log.out(f'Atom volume calculation execution time: {execution_time} (seconds)')
        
        # Write files
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
        self.voxel_volume = self.simulation_volume / (v.nx*v.ny*v.nz)
        self.atom_volume = self.voxel_volume*nvoxel_atomIDs
        self.free_volume = self.voxel_volume*nvoxel_freeIDs
        self.pfree_volume = 100*(self.free_volume/self.simulation_volume)
        
        
        #-------------------------------------------#
        # Compute free volume spatial distributions #
        #-------------------------------------------#
        if files2write['write_spat_dis-x'] or files2write['write_spat_dis-y'] or files2write['write_spat_dis-z']: 
            start_time = time.time()
            if run_mode == 'numba-dd': 
                log.out('\n\nFinding free volume spatial distributions (compiled - serial) ...')
                volx, voly, volz = spat_dist_serial(v.numpyx, v.numpyy, v.numpyz, voxels, flags, self.voxel_volume)
            elif run_mode == 'numba-ddp': 
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
            pflags_num = np.zeros(3, dtype=np.int64)
            for i, flag in enumerate(v.pflags):
                if flag == 'p':  pflags_num[i] = 1
    
            # numba find free volume connectivity
            start_time = time.time()
            if run_mode == 'numba-dd': 
                log.out('Finding free volume voxelID connectivity (compiled - serial) ...')
                numpy_graph = free_dist_serial_frac_indexed(flags, v.nx, v.ny, v.nz, pflags_num)
            elif run_mode == 'numba-ddp': 
                log.out('Finding free volume voxelID connectivity (compiled - parallel) ...')
                numpy_graph = free_dist_parallel_frac_indexed(flags, v.nx, v.ny, v.nz, pflags_num)
            
            # Convert numpy graph to python dict
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
            self.free_volume_clusters = {} # { clusterID : Info object }
            self.voxelclusterID2molID = {i:len(free_volume_clusters)+1 for i in self.voxel_freeIDs} # { voxelID : clustered }
            for ID, cluster in enumerate(free_volume_clusters, 1):
                for i in cluster:
                    self.voxelclusterID2molID[i] = ID
                xspan, yspan, zspan = misc_func.find_free_volume_cluster_span(cluster, self.voxel_freeIDs)
                I = Info()
                I.size = len(cluster)
                I.psize = 100*len(cluster)/len(self.voxel_freeIDs)
                I.volume = len(cluster)*self.voxel_volume
                I.pvolume = 100*(len(cluster)*self.voxel_volume)/self.simulation_volume
                I.xspan = xspan
                I.yspan = yspan
                I.zspan = zspan
                self.free_volume_clusters[ID] = I
                
                
            #print('self.pbc_bonded_count', self.pbc_bonded_count)
            #bonded = {len(free_volume_voxel_graph[i]) for i in free_volume_voxel_graph}
            #print(bonded)
            #print(len(self.free_volume_clusters))
            #print('count TRUE freeID_edgeflags: ', list(freeID_edgeflags.values()).count(True))