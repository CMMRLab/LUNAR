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
from tqdm import tqdm
import math
import sys

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
        
        # Set probe_diameter if shortcut flag exists
        self.probe_diameter = probe_diameter
        if probe_diameter == 'min-voxel':
            log.out("probe_diameter = 'min-voxel', finding minimum dimension of voxel ...")
            self.probe_diameter = min([self.dx, self.dy, self.dz])
            log.out(f'   updated probe_diameter = {self.probe_diameter}')
        
        # Generate voxels 
        ID = 0; self.voxelIDs = {} # { voxelID : (x, y, z, radial distance from center of simulation cell) }
        self.radial_voxel = {} # {voxelID : r}
        self.postions = {'x':set(), 'y':set(), 'z':set()}
        xadd = self.halfdx + self.xlo; yadd = self.halfdy + self.ylo; zadd = self.halfdz + self.zlo;
        log.out('\n\nGenerating voxels:')
        for nx in tqdm(range(self.nx)):
            for ny in range(self.ny):
                for nz in range(self.nz):
                    ID += 1
                    xpos = nx*self.dx + xadd
                    ypos = ny*self.dy + yadd
                    zpos = nz*self.dz + zadd
                    r = misc_func.compute_distance(xpos, ypos, zpos, self.cx, self.cy, self.cz)
                    self.voxelIDs[ID] = (xpos, ypos, zpos)
                    self.radial_voxel[ID] = r
                    self.postions['x'].add(xpos)
                    self.postions['y'].add(ypos)
                    self.postions['z'].add(zpos)
        
        
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
            element2TypeIDs = {} # { element : atom Type ID }
            for ID, element in enumerate(sorted(m.elements), 1):
                self.atomTypeIDs[ID] = (element, mass_map[element][0])
                element2TypeIDs[element] = ID
        elif vdw_method in ['class1', 'class2']:
            for ID in m.masses:
                masses = m.masses[ID]
                self.atomTypeIDs[ID] = (masses.element, masses.coeffs[0])
        else: log.error(f'ERROR vdw_method {vdw_method} not supported')
        
        # Iterate through atoms and voxelIDs to find where atom size is localized
        log.out('\n\nFinding atom volumes and free volumes (standard library) ...')
        freeIDs = {i for i in v.voxelIDs}; self.voxel_atomIDs = {}; self.pbc_count = 0; self.npossible_pbc = 0; 
        for i in tqdm(m.atoms):
            atom = m.atoms[i]; ownedIDs = set(); periodic_switch = 0;
            x1 = atom.x; y1 = atom.y; z1 = atom.z; 
            
            # Set atom vdw radii based on method
            if vdw_method == 'dict':
                vdw_radii = vdw_radius[atom.element]
            elif vdw_method in ['class1', 'class2']:
                sigma = m.pair_coeffs[atom.type].coeffs[1]
                vdw_radii = sigma/2
            else: log.error(f'ERROR vdw_method {vdw_method} not supported')
            vdw_radii += probe_diameter/2
            
            radial_atom = misc_func.compute_distance(x1, y1, z1, v.cx, v.cy, v.cz)
            min_radius = radial_atom - vdw_radii - max_voxel_size
            max_radius = radial_atom + vdw_radii + max_voxel_size
            if v.pflags.count('f') < 3: 
                atom_edgeflag = misc_func.check_near_edge(x1, y1, z1, vdw_radii, v.xlo, v.xhi, v.ylo, v.yhi, v.zlo, v.zhi)
                periodic_postions = [ (x1+ixlx, y1+iyly, z1+izlz) for ixlx, iyly, izlz in v.scaled_images ]
                if atom_edgeflag: self.npossible_pbc += 1
            else: periodic_postions = [ (x1, y1, z1) ]; atom_edgeflag = False;
            for ID in freeIDs: 
                x2, y2, z2 = v.voxelIDs[ID]
                r = v.radial_voxel[ID]
                if min_radius < r < max_radius:
                    if not atom_edgeflag:
                        if abs(x1 - x2) > vdw_radii: continue
                        elif abs(y1 - y2) > vdw_radii: continue
                        elif abs(z1 - z2) > vdw_radii: continue
                        else: 
                            distance_fff = misc_func.compute_distance(x1, y1, z1, x2, y2, z2)
                            if distance_fff < vdw_radii:
                                if vdw_method == 'dict': TypeID = element2TypeIDs[atom.element]
                                elif vdw_method in ['class1', 'class2']: TypeID = atom.type
                                else: log.error(f'ERROR vdw_method {vdw_method} not supported')
                                self.voxel_atomIDs[ID] = (x2, y2, z2, TypeID)
                                ownedIDs.add(ID); continue
                    else:
                        for x1i, y1i, z1i in periodic_postions:
                            if abs(x1i - x2) > vdw_radii: continue
                            elif abs(y1i - y2) > vdw_radii: continue
                            elif abs(z1i - z2) > vdw_radii: continue
                            else:
                                distance_ppp = misc_func.compute_distance(x1i, y1i, z1i, x2, y2, z2)
                                if distance_ppp < vdw_radii:
                                    if vdw_method == 'dict': TypeID = element2TypeIDs[atom.element]
                                    elif vdw_method in ['class1', 'class2']: TypeID = atom.type
                                    else: log.error(f'ERROR vdw_method {vdw_method} not supported')
                                    self.voxel_atomIDs[ID] = (x2, y2, z2, TypeID)
                                    ownedIDs.add(ID); self.pbc_count += 1;
                                    if periodic_switch < 5:
                                        periodic_postions.remove( (x1i, y1i, z1i) )
                                        periodic_postions.insert(0, (x1i, y1i, z1i))
                                        periodic_switch += 1
                                    break

            # Remove "owned" voxelIDs before next iteration of atomIDs
            freeIDs = freeIDs - ownedIDs

        # Find the freeIDs after all atomIDs have been iterate through
        self.voxel_freeIDs = {ID:v.voxelIDs[ID] for ID in freeIDs}
        
        # Find simulation cell volume
        self.simulation_volume = v.lx*v.ly*v.lz
        
        # Find free volume
        self.voxel_volume = v.dx*v.dy*v.dz
        self.atom_volume = self.voxel_volume*len(self.voxel_atomIDs)
        self.free_volume = self.voxel_volume*len(self.voxel_freeIDs)
        self.pfree_volume = 100*(self.free_volume/self.simulation_volume)
        
        #-------------------------------------------#
        # Compute free volume spatial distributions #
        #-------------------------------------------#
        if files2write['write_spat_dis-x'] or files2write['write_spat_dis-y'] or files2write['write_spat_dis-z']: 
            log.out('\n\nFinding free volume spatial distributions ...')
            self.spatial_distributions = {'x':{i:0 for i in v.postions['x']},
                                          'y':{i:0 for i in v.postions['y']},
                                          'z':{i:0 for i in v.postions['z']}}
            for ID in tqdm(self.voxel_freeIDs):
                xpos, ypos, zpos = self.voxel_freeIDs[ID]
                self.spatial_distributions['x'][xpos] += self.voxel_volume
                self.spatial_distributions['y'][ypos] += self.voxel_volume
                self.spatial_distributions['z'][zpos] += self.voxel_volume
        
        
        #----------------------------------------------------#
        # Free volume connectivity and distribution analysis #
        #----------------------------------------------------# 
        if compute_free_volume_distributions:
            log.out('\n\n----------------------------------------------')
            log.out('| Starting free volume distribution analysis |')
            log.out('----------------------------------------------')
            
            # Determine "nearness" to edge. Only search periodic image if atom is this close to an edge.
            freeID_edgeflags = {} # { VoxelID : edge flag }
            for i in self.voxel_freeIDs:
                x, y, z = self.voxel_freeIDs[i]
                if v.pflags.count('f') == 3: freeID_edgeflags[i] = False # If system is non-periodic set all to False (speed-up)
                else: freeID_edgeflags[i] = misc_func.check_near_edge(x, y, z, max_voxel_size, v.xlo, v.xhi, v.ylo, v.yhi, v.zlo, v.zhi)
            
            # Find free volume connectivity (7 sec -> 3 sec)
            free_volume_voxel_graph = {ID:set() for ID in self.voxel_freeIDs} # { voxelID : set(bonded voxels) }
            self.pbc_bonded_count = 0; freeIDs = {i for i in self.voxel_freeIDs}
            log.out('Finding free volume voxelID connectivity (standard library) ...')
            for id1 in tqdm(self.voxel_freeIDs):
                if len(free_volume_voxel_graph[id1]) >= 26: continue
                x1, y1, z1 = self.voxel_freeIDs[id1]
                r1 = v.radial_voxel[id1]
                freeIDs.remove(id1);
                min_radius = r1 - 2*max_voxel_size
                max_radius = r1 + 2*max_voxel_size
                if freeID_edgeflags[id1]:
                    periodic_postions = misc_func.find_periodic_postions(v.scaled_images, x1, y1, z1, v.cx, v.cy, v.cz, Npos=12)
                for id2 in freeIDs:
                    if id1 == id2: continue
                    #if len(free_volume_voxel_graph[id2]) >= 26: continue
                    if min_radius < v.radial_voxel[id2] < max_radius:
                        x2, y2, z2 = self.voxel_freeIDs[id2]
                    
                        # Find non periodically bonded free volume voxelIDs
                        if not freeID_edgeflags[id1] or not freeID_edgeflags[id2]:
                            if abs(x1-x2) > max_voxel_size: continue
                            elif abs(y1-y2) > max_voxel_size: continue
                            elif abs(z1-z2) > max_voxel_size: continue
                            dx = (x1 - x2)**2
                            dy = (y1 - y2)**2
                            dz = (z1 - z2)**2
                            distance = math.sqrt(dx + dy + dz)
                            if distance <= max_voxel_size:
                                free_volume_voxel_graph[id1].add(id2)
                                free_volume_voxel_graph[id2].add(id1)
                        
                        # Find periodically bonded free volume voxelIDs
                        elif freeID_edgeflags[id1] and freeID_edgeflags[id2]:
                            for x1i, y1i, z1i in periodic_postions:
                                if abs(x1i - x2) > max_voxel_size: continue
                                elif abs(y1i - y2) > max_voxel_size: continue
                                elif abs(z1i - z2) > max_voxel_size: continue
                                dx = (x1i - x2)**2
                                dy = (y1i - y2)**2
                                dz = (z1i - z2)**2
                                distance = math.sqrt(dx + dy + dz)
                                if distance <= max_voxel_size:
                                    free_volume_voxel_graph[id1].add(id2)
                                    free_volume_voxel_graph[id2].add(id1)
                                    self.pbc_bonded_count += 1; break;

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
            self.voxelclusterID2molID = {i:len(free_volume_clusters)+1 for i in self.voxel_freeIDs} # { voxelID : clusterID }
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
            