# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.6
March 20th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.cell_builder.misc_functions as misc_functions
import random
import math
import os


# Function to check if number is even
def is_even(number):
    if number % 2 == 0:
        return True
    else: return False
    
# Function to generate x, y, and z images from domain flags 
def generate_images_from_domain(ni):
    images = [0]; nspan = 0;
    if ni > 1:
        if is_even(ni):
            n = math.ceil(ni/2); nspan = 0.5;
            images = [i for i in range(-n, n)]
        else:
            n = int((ni-1)/2); nspan = 0;
            images = [i for i in range(-n, n+1)]
    return images, nspan

# Function to compute distance when math.dist is not available
def compute_distance(x1, y1, z1, x2, y2, z2):
    dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
    return math.sqrt(dx*dx + dy*dy + dz*dz)


# Class to build subcell object complete with atoms, bonds, angles, dihedrals, impropers, and coeffs
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz .comment
class Bond: pass  # .type .atomids = [atom1id, atom2id]
class Angle: pass  # .type .atomids = [atom1id, atom2id, atom3id]
class Dihedral:pass  # .type .atomids = [atom1id, atom2id, atom3id, atom4id]
class Improper: pass  # .type .atomids = [atom1,atom2,atom3,atom4]
class Coeff_class: pass  # .type .coeffs = []
class constructor:
    def __init__(self, files, qtys, offset_coeff_types, duplicate, max_rotations, distance_scale, moleculespans,
                 seed, reset_molids, domain, molecule_insertion, occurrences, subcells, log, pflag=True, grouping=False):
        # Structure info
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.angles = {}  # {angle number : angle object}
        self.dihedrals = {}  # {dihedral number : dihedral object}
        self.impropers = {}  # {improper number : improper object}
        
        # Quantities
        self.natoms = 0
        self.natomtypes = 0
        self.nbonds = 0
        self.nbondtypes = 0
        self.nangles = 0
        self.nangletypes = 0
        self.ndihedrals = 0
        self.ndihedraltypes = 0
        self.nimpropers = 0
        self.nimpropertypes = 0
        self.nbondbond = 0
        self.nbondangle = 0
        self.nangleangle = 0
        self.nangleangletorsion = 0
        self.nendbondtorsion = 0
        self.nmiddlebondtorsion = 0
        self.nbondbond13 = 0
        self.nangletorsion = 0
        
        #--------------------------------------#
        # Find max span of lattice grid points #
        #--------------------------------------#
        #self.maxspan = math.ceil(distance_scale*max(moleculespans))
        self.maxspan = distance_scale*max(moleculespans)
        if self.maxspan <= 0: self.maxspan = 5 # set a 5 angstrom default if particles
        if distance_scale < 1: log.warn('WARNING using distance scale smaller then 1 may result in overlapped atoms/molecules. USE AT YOUR OWN RISK.')
        
        # Find total number of subcells needed 
        self.nsubcells = 0;
        for fileID in files:
            self.nsubcells += qtys[fileID]
        self.nsubcells = int(duplicate*self.nsubcells)
        if pflag: 
            log.out('\n\n\n----------------------------------------------')
            log.out('              Subcell info')
            log.out('----------------------------------------------')
            log.out(f'Total number of subcells required:  {self.nsubcells}')
            log.out(f'Lattice spacing for subcells     : {self.maxspan} angstrom')
        
        #--------------------------------#
        # Set up molIDs based on fileIDs #
        #--------------------------------#
        self.molids_based_on_files = {files[fileID].filename:fileID for fileID in files} # { filename : molID }
        if reset_molids == 'files' and occurrences == 0 and pflag:
            log.out('Using index of files to set every atomIDs molID as the file index')
            for filename in self.molids_based_on_files:
                log.out('    all atoms in {} are being assigned to molID {}'.format(filename, self.molids_based_on_files[filename]))
                
        #--------------------------------#
        # Set up molIDs based on offsets #
        #--------------------------------#
        self.molids_based_on_offsets = {} # { filename : molID offset }
        offset = 0
        if reset_molids == 'offset' and occurrences == 0:
            for fileID in files:
                m = files[fileID]
                try: molids = list({m.atoms[i].molid for i in m.atoms})
                except: molids = [fileID]
                self.molids_based_on_offsets[m.filename] = offset
                offset += max(molids)
            if pflag:
                for filename in self.molids_based_on_offsets:
                    log.out('    all atoms in {} have their molIDs offset by {}'.format(filename, self.molids_based_on_offsets[filename]))
                  
        #--------------------------------------------------------------#
        # Setting up offsets (default will be to keep all offsets as   #
        # zeros and only apply offsets if offset_coeff_types is True). #
        # Addtionally offsets will only be applied on occurance = 0 to #
        # avoid doubling offsets if group_monomers_locally is True     #
        #--------------------------------------------------------------#
        self.coeff_offsets = {} # { filename : {'atom':offset, 'bond':offset, ... } }
        tally = {'atom':0, 'bond':0, 'angle':0, 'dihedral':0, 'improper':0}
        for fileID in files:
            m = files[fileID]
            tmp = {'atom': tally['atom'],
                   'bond': tally['bond'],
                   'angle': tally['angle'],
                   'dihedral': tally['dihedral'],
                   'improper': tally['improper']}
            self.coeff_offsets[m.filename] = tmp
            if offset_coeff_types and occurrences == 0:
                tally['atom'] += len(m.masses)
                tally['bond'] += len(m.bond_coeffs)
                tally['angle'] += len(m.angle_coeffs)
                tally['dihedral'] += len(m.dihedral_coeffs)
                tally['improper'] += len(m.improper_coeffs)


        #---------------------------------------------#
        # Check for any zeros in the qty's, if create #
        # system and add to system accordingly        #
        #---------------------------------------------#
        adding_to_system = False; lattice_shift = {'x':0, 'y':0, 'z':0}; overlapped_lattices = [] # { (ix, iy, ix), ...}
        if 0 in list(qtys.values()) and not grouping:
            adding_to_system = True
            nsystem = list(qtys.values()).count(0)
            log.out(f'\n\n{nsystem} LAMMPS datafiles were found to have a qty of ZERO. Using these files to define a system to add molecules to.')
            box = {'xlo':[], 'xhi':[], 'ylo':[], 'yhi':[],'zlo':[], 'zhi':[],}
            for fileID in files:
                m = files[fileID]
                if qtys[fileID] == 0:
                    filename = os.path.basename(m.filename)
                    log.out(f'  Added {filename} to system.')
                    
                    # Get simulation cell info
                    xline = m.xbox_line.split(); yline = m.ybox_line.split(); zline = m.zbox_line.split();
                    box['xlo'].append(float(xline[0])); box['xhi'].append(float(xline[1]));
                    box['ylo'].append(float(yline[0])); box['yhi'].append(float(yline[1]));
                    box['zlo'].append(float(zline[0])); box['zhi'].append(float(zline[1]));
                    
                    # Add "molecule" to the system applying no xshift, yshift, zshift and no rotation about x, y, z
                    self.add_molecule_to_system(m, fileID, reset_molids, occurrences, 0, 0, 0, 0, 0, 0, log)
            
            # Redefine simulation box
            self.xlo = min(box['xlo']); self.xhi = max(box['xhi']);
            self.ylo = min(box['ylo']); self.yhi = max(box['yhi']);
            self.zlo = min(box['zlo']); self.zhi = max(box['zhi']);
            self.lx = self.xhi - self.xlo
            self.ly = self.yhi - self.ylo
            self.lz = self.zhi - self.zlo
            self.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(self.xlo, self.xhi, 'xlo', 'xhi')
            self.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(self.ylo, self.yhi, 'ylo', 'yhi')
            self.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(self.zlo, self.zhi, 'zlo', 'zhi')
            
            # Compute nx, ny, nz number of voxels required to build 3D voxel grid
            nx = math.floor(self.lx/self.maxspan)
            ny = math.floor(self.ly/self.maxspan)
            nz = math.floor(self.lz/self.maxspan)
            lx = round(self.lx, 3); ly = round(self.ly, 3); lz = round(self.lz, 3);
            llx = round(nx*self.maxspan, 3); lly = round(ny*self.maxspan, 3); llz = round(nz*self.maxspan, 3);
            
            # If the diff bewteen li - lli - maxspan is less then the tolerance, add an extra lattic point
            tolerance = 1.0
            if abs(lx - llx - self.maxspan) < tolerance: 
                nx += 1
                llx = round(nx*self.maxspan, 3)
            if abs(ly - lly - self.maxspan) < tolerance:
                ny += 1
                lly = round(ny*self.maxspan, 3)
            if abs(lz - llz - self.maxspan) < tolerance:
                nz += 1
                llz = round(nz*self.maxspan, 3)
            log.out("  After adding files with qty's of ZERO's, the simulation domain has been established as:")
            log.out(f'    x-dimensions: {self.xbox_line}  (lx: {lx}, x-lattice-span: {llx})')
            log.out(f'    y-dimensions: {self.ybox_line}  (ly: {ly}, y-lattice-span: {lly})')
            log.out(f'    z-dimensions: {self.zbox_line}  (lz: {lz}, z-lattice-span: {llz})')
            domain = '{}x{}x{}'.format(nx, ny, nz)
            log.out(f'  Resulting in the creation of a domain of {domain}')
            
            # Find center of system simulation cell and set lattice shift values
            cx = (self.xhi + self.xlo)/2
            cy = (self.yhi + self.ylo)/2
            cz = (self.zhi + self.zlo)/2
            lattice_shift['x'] = cx
            lattice_shift['y'] = cy
            lattice_shift['z'] = cz

        #------------------------------------------------------------------------#
        # Find image information (by finding max nimages needed to hold problem) #
        # If grouping is True default to cubic, else decide based on user input  #
        #------------------------------------------------------------------------#    
        if grouping or domain == 'cubic':            
            # Generate sub cells
            nimages = math.ceil( (self.nsubcells**(1/3)-1)/2 ) # formula to find lowest nimages needed
            images = set([]) # set to hold unique images
            for ix in range(-nimages, nimages+1):
                for iy in range(-nimages, nimages+1):
                    for iz in range(-nimages, nimages+1):            
                        images.add( (ix, iy, iz) )
        else:
            n = domain.split('x')
            if len(n) != 3:
                log.error(f"ERROR domain = {domain} which does not supply 3 directions. Example: domain = '4x4x2'")
            else:
                nx = int(n[0]); ny = int(n[1]); nz = int(n[2])
                imagesx, nspanx = generate_images_from_domain(nx)
                imagesy, nspany = generate_images_from_domain(ny)
                imagesz, nspanz = generate_images_from_domain(nz)   
                lattice_shift['x'] += nspanx*self.maxspan
                lattice_shift['y'] += nspany*self.maxspan
                lattice_shift['z'] += nspanz*self.maxspan
                nlattices = len(imagesx)*len(imagesy)*len(imagesz)
                if nlattices < self.nsubcells and not adding_to_system: 
                    log.error(f'ERROR required lattice points = {self.nsubcells}, domain = {domain} with {nlattices} lattices. Adjust domain')
                images = set([]) # set to hold unique images
                for ix in imagesx:
                    for iy in imagesy:
                        for iz in imagesz:    
                            images.add( (ix, iy, iz) )
                    
        # Sort images in ascending order by sum of iflags {IF nimages=1, then
        # (0, 0, 0), (0, 1, 0), (0, 0, -1), ... (1, 1, -1), (1,1,1) }. 
        images = list(images)
        if grouping: images = sorted(images, key=lambda x: sum([abs(i) for i in x]))
        
        #--------------------------------------------#
        # Generate 3D grid and populate with gridIDs #
        #--------------------------------------------#       
        # Only iterate to as many subcells as needed
        lattice_iflags = {} # { gridID : (ix, iy, iz) }
        grid = {} # { gridID : (x, y, z) }
        gridID = 0
        for ix, iy, iz in images:
            x = ix*self.maxspan + lattice_shift['x']
            y = iy*self.maxspan + lattice_shift['y']
            z = iz*self.maxspan + lattice_shift['z']
            
            # if adding_to_system, we need to check for any overlaps between existing atoms and lattice point
            if adding_to_system:
                overlap = False
                for i in self.atoms:
                    atom = self.atoms[i]
                    distance = compute_distance(x, y, z, atom.x, atom.y, atom.z)
                    if distance < self.maxspan/2:
                        overlap = True; break
                if overlap: 
                    overlapped_lattices.append((ix, iy, iz))
                    continue
            gridID += 1
            grid[gridID] = (x, y, z)
            lattice_iflags[gridID] = (ix, iy, iz)
        if adding_to_system and overlapped_lattices:
            log.out(f'  While generating lattice points {len(overlapped_lattices)} from domain {domain} could not be generated as')
            log.out('  they would result in placing a molecule to be overlapped with the existing system\n\n')
        
        #----------------------------------------------------------#
        # Find new simulation cell bounds (if not adding to system #
        # otherwise use the simulation cell bounds from above).    #
        #----------------------------------------------------------#
        if not adding_to_system:
            dimensions = {'x':[], 'y':[], 'z':[]}
            for ID in grid:
                dimensions['x'].append(grid[ID][0])
                dimensions['y'].append(grid[ID][1])
                dimensions['z'].append(grid[ID][2])
            self.xlo = min(dimensions['x']) - 0.5*self.maxspan
            self.xhi = max(dimensions['x']) + 0.5*self.maxspan
            self.ylo = min(dimensions['y']) - 0.5*self.maxspan
            self.yhi = max(dimensions['y']) + 0.5*self.maxspan
            self.zlo = min(dimensions['z']) - 0.5*self.maxspan
            self.zhi = max(dimensions['z']) + 0.5*self.maxspan
            self.lx = self.xhi - self.xlo
            self.ly = self.yhi - self.ylo
            self.lz = self.zhi - self.zlo
            self.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(self.xlo, self.xhi, 'xlo', 'xhi')
            self.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(self.ylo, self.yhi, 'ylo', 'yhi')
            self.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(self.zlo, self.zhi, 'zlo', 'zhi')
        
        #--------------------------------------------------#
        # Assign molecules that have specified Region for  #
        # solvation first to guarntee lattice point exists #
        #--------------------------------------------------#
        solvation_flag = {} # { fileID : boolean }
        for fileID in files:
            if hasattr(files[fileID],'Regions'):
                if files[fileID].Regions:
                    solvation_flag[fileID] = 1
                else: solvation_flag[fileID] = 0
            else: solvation_flag[fileID] = 0
        nsolvation_molecules = list(solvation_flag.values()).count(1)
        if nsolvation_molecules > 1: 
            log.error('ERROR Region = {} specified in more then one read-in file. This option is only supported for a single file.')
        if nsolvation_molecules == 1:
            # Reset files and qtys to ensure the first molecules assigned to lattices are the ones with regions
            solvation_files = {} # { fileID : file-object }
            solvation_qtys  = {} # { fileID : n-quantity }
            solvation_flag  = dict(sorted(solvation_flag.items(), key=lambda item: item[1], reverse=True))
            for i in solvation_flag:
                solvation_files[i] = files[i]
                solvation_qtys[i]  = qtys[i]
            files = solvation_files
            qtys = solvation_qtys
        
        #--------------------------------------------------#
        # Function to find gridIDs in a specified region   #
        # based on lattice iflags (for solvation problems) #
        #--------------------------------------------------#
        def get_eligible_lattices(Regions, gridIDs, lattice_iflags):
            xlo, xhi = Regions['x']; ylo, yhi = Regions['y']; zlo, zhi = Regions['z']
            possible_gridIDs = []
            for gridID in gridIDs:
                ix, iy, iz = lattice_iflags[gridID]
                if xlo <= ix <= xhi and ylo <= iy <= yhi and zlo <= iz <= zhi:
                    possible_gridIDs.append(gridID)
            return possible_gridIDs
        
        #------------------------------------------------------#
        # Start random mapping of dupID set fileIDs to gridIDs #
        #------------------------------------------------------#
        if seed > 0: random.seed(seed)
        dupID = 0 # to keep track of duplicationIDs (starts from 1)
        dup2file = {} # { dupID : fileID used in location }
        dup2grid = {} # { dupID : gridID }
        gridIDs = [i for i in grid] # will reduce as we assign gridID->dupID
        for fileID in files:
            if qtys[fileID] == 0: continue
            for qty in range(qtys[fileID]):
                for i in range(duplicate):
                    dupID += 1 # increment dupID
                    
                    # Check for regions option
                    m = files[fileID]
                    if gridIDs:
                        if m.Regions:                         
                            # Find possible gridIDs to assign molecule to
                            possible_gridIDs = get_eligible_lattices(m.Regions, gridIDs, lattice_iflags)
                            if not possible_gridIDs:
                                log.error('ERROR ran out of lattice points to assign molecules in region. Increase region size or decrease qty or duplicate.')
                            else:
                                random_index_solvation = random.randint(0, len(possible_gridIDs)-1)
                                solvation_gridID = possible_gridIDs[random_index_solvation]
                                random_index = gridIDs.index(solvation_gridID) 
                                                    
                        # else get completely random integer in len(gridIDs) to get gridID to assign to dupID
                        else: random_index = random.randint(0, len(gridIDs)-1)
                        
                        # create maps and delete random_index from gridIDs
                        dup2file[dupID] = fileID
                        dup2grid[dupID] = gridIDs[random_index]
                        if not grouping:
                            if occurrences == 0: molecule_insertion[fileID][0] += 1
                            #else: molecule_insertion[fileID][0] += qtys[fileID]
                        del gridIDs[random_index]
        
        #-----------------------------------------------------------------#
        # Log how many files of each where able to be added to the system #
        #-----------------------------------------------------------------#
        if adding_to_system and occurrences == 0:
            for fileID in molecule_insertion:
                created, desired = molecule_insertion[fileID]
                if created == desired and qtys[fileID] != 0:
                    log.out(f'Created {created} of {desired} molecules from {files[fileID].filename}')
                if created != desired and qtys[fileID] != 0:
                    log.warn(f'WARNING created {created} of {desired} molecules from {files[fileID].filename}')
                    
        #--------------------------------------------#
        # Start moving atoms to randomized locations #
        #--------------------------------------------#
        def generate_incremented_lst(s):  
            return [n*s['increment'] for n in range( int((s['end']-s['start'])/s['increment'])+1) ]  
        rx = generate_incremented_lst({'start':0, 'end':max_rotations['x'], 'increment': 0.1}) # increment by 0.1 degrees
        ry = generate_incremented_lst({'start':0, 'end':max_rotations['y'], 'increment': 0.1}) # increment by 0.1 degrees
        rz = generate_incremented_lst({'start':0, 'end':max_rotations['z'], 'increment': 0.1}) # increment by 0.1 degrees
        if pflag: log.out('\n\nGenerating new simulation cell ....')
        for ID in dup2grid:
            xshift, yshift, zshift = grid[dup2grid[ID]] # get x, y, z shift
            m = files[dup2file[ID]] # get m-object
            fileid = dup2file[ID]
            if m.Regions and m.Rotations: 
                rx_local = generate_incremented_lst({'start':0, 'end':m.Rotations['x'], 'increment': 0.1}) # increment by 0.1 degrees
                ry_local = generate_incremented_lst({'start':0, 'end':m.Rotations['y'], 'increment': 0.1}) # increment by 0.1 degrees
                rz_local = generate_incremented_lst({'start':0, 'end':m.Rotations['z'], 'increment': 0.1}) # increment by 0.1 degrees
                phi = 0; theta = 0; psi = 0; # intialize and update if max_rotation > 0
                if m.Rotations['x'] > 0: phi = rx[random.randint(0, len(rx_local)-1)]
                if m.Rotations['y'] > 0: theta = ry[random.randint(0, len(ry_local)-1)]
                if m.Rotations['z'] > 0: psi = rz[random.randint(0, len(rz_local)-1)]
            else:
                # Find random phi (rx), theta (ry), psi (rz)
                phi = 0; theta = 0; psi = 0; # intialize and update if max_rotation > 0
                if max_rotations['x'] > 0: phi = rx[random.randint(0, len(rx)-1)]
                if max_rotations['y'] > 0: theta = ry[random.randint(0, len(ry)-1)]
                if max_rotations['z'] > 0: psi = rz[random.randint(0, len(rz)-1)]
                
            # Rotate molecule and add to system
            m = misc_functions.rotate_molecule(m, phi, theta, psi)
            self.add_molecule_to_system(m, fileid, reset_molids, occurrences, xshift, yshift, zshift, phi, theta, psi, log)
                
        #-------------------------------#
        # Update force field parameters #
        #-------------------------------#
        self.update_force_field(files, subcells, grouping, pflag, log)
        
        #---------------------------------------------------------------------#
        # If adding_to_system, check that no atoms were placed to be periodic #
        #---------------------------------------------------------------------#
        if adding_to_system: self. wrap_periodic_atoms()
        
    
    #---------------------#
    # Define some methods #
    #---------------------#
    # Method to wrap any atoms that may be periodic due to adding to system
    def wrap_periodic_atoms(self):
        for i in self.atoms:
            atom = self.atoms[i]
            if atom.x < self.xlo:
                atom.x += self.lx
                atom.ix += 1
            if atom.x > self.xhi:
                atom.x -= self.lx
                atom.ix -= 1
            if atom.y < self.ylo:
                atom.y += self.ly
                atom.iy += 1
            if atom.y > self.yhi:
                atom.y -= self.ly
                atom.iy -= 1
            if atom.z < self.zlo:
                atom.z += self.lz
                atom.iz += 1
            if atom.z > self.zhi:
                atom.z -= self.lz
                atom.iz -= 1
        return
    
    # Method to add molecule from m class to system
    def add_molecule_to_system(self, m, fileid, reset_molids, occurrences, xshift, yshift, zshift, phi, theta, psi, log):
        # Start moving atoms and reseting atomIDs
        atomID_map = {} # {orginal atomID : new atomID }
        for i in m.atoms:
            atom = m.atoms[i]
            self.natoms += 1 # increment atom count
            atomID_map[i] = self.natoms
            
            # Set molIDs
            if reset_molids == 'files':
                if occurrences == 0:
                    molid = self.molids_based_on_files[m.filename]
                else:
                    try: molid = atom.molid
                    except: molid = 1
            elif reset_molids == 'offset':
                if occurrences == 0:
                    offset = self.molids_based_on_offsets[m.filename]
                    try: molid = atom.molid
                    except: molid = 1
                    molid += offset
                else:
                    try: molid = atom.molid
                    except: molid = 1
            elif reset_molids == 'clusters':
                try: molid = atom.molid
                except: molid = 1
            elif reset_molids == 'skip':
                try: molid = atom.molid
                except: molid = 1
            else: log.error(f'ERROR request unsupported reset_molids {reset_molids} option.')
            
            # Save atom info
            a = Atom()
            a.type = atom.type + self.coeff_offsets[m.filename]['atom']
            a.charge = atom.charge
            a.x = atom.x + xshift
            a.y = atom.y + yshift
            a.z = atom.z + zshift
            a.comment = atom.comment
            try: a.fileid = atom.fileid # atom might already have a fileid if grouping
            except: a.fileid = fileid # for possible .mol2 resID setting to vis mols in VMD
            a.rotation = [phi, theta, psi]
            a.location =[xshift, yshift, zshift]
            a.filename = m.filename
            a.molid = molid
            if xshift == 0 and yshift == 0 and zshift == 0 and phi == 0 and theta == 0 and phi == 0: 
                try:
                    a.ix = atom.ix
                    a.iy = atom.iy
                    a.iz = atom.iz
                except:
                    a.ix = 0 # reset to zero
                    a.iy = 0 # reset to zero
                    a.iz = 0 # reset to zero        
            else:
                a.ix = 0 # reset to zero
                a.iy = 0 # reset to zero
                a.iz = 0 # reset to zero
            # if atom has attribute group, pass group info
            if hasattr(atom,'group'):
                a.group = atom.group
            self.atoms[self.natoms] = a
            
        # Setting new bonds
        for i in m.bonds:
            bond = m.bonds[i]
            self.nbonds += 1 # increment bond count
            id1, id2 = bond.atomids
            b = Bond()
            b.type = bond.type + self.coeff_offsets[m.filename]['bond']
            b.atomids = [atomID_map[id1], atomID_map[id2]]
            self.bonds[self.nbonds] = b
            
        # Setting new angles
        for i in m.angles:
            angle = m.angles[i]
            self.nangles += 1 # increment angle count
            id1, id2, id3 = angle.atomids
            a = Angle()
            a.type = angle.type + self.coeff_offsets[m.filename]['angle']
            a.atomids = [atomID_map[id1], atomID_map[id2], atomID_map[id3]]
            self.angles[self.nangles] = a
            
        # Setting new dihedrals
        for i in m.dihedrals:
            dihedral = m.dihedrals[i]
            self.ndihedrals += 1 # increment dihedral count
            id1, id2, id3, id4 = dihedral.atomids
            d = Dihedral()
            d.type = dihedral.type + self.coeff_offsets[m.filename]['dihedral']
            d.atomids = [atomID_map[id1], atomID_map[id2], atomID_map[id3], atomID_map[id4]]
            self.dihedrals[self.ndihedrals] = d
            
        # Setting new impropers
        for i in m.impropers:
            improper = m.impropers[i]
            self.nimpropers += 1 # increment improper count
            id1, id2, id3, id4 = improper.atomids
            i = Improper()
            i.type = improper.type + self.coeff_offsets[m.filename]['improper']
            i.atomids = [atomID_map[id1], atomID_map[id2], atomID_map[id3], atomID_map[id4]]
            self.impropers[self.nimpropers] = i   
        return
    
    # Method to handle force field
    def update_force_field(self, files, subcells, grouping, pflag, log):
        if grouping: self.filename = 'grp_local_subcell_{}'.format(len(subcells))
        else: self.filename = 'cell_builder_tmp_name'
        self.header = ''; 
        self.xy = 0
        self.xz = 0
        self.yz = 0
        self.masses = {}; self.mass_coeffs_style_hint = {};
        self.pair_coeffs = {}; self.pair_coeffs_style_hint = {};
        self.bond_coeffs = {}; self.bond_coeffs_style_hint = {};
        self.angle_coeffs = {}; self.angle_coeffs_style_hint = {};
        self.dihedral_coeffs = {}; self.dihedral_coeffs_style_hint = {};
        self.improper_coeffs = {}; self.improper_coeffs_style_hint = {};
        self.bondbond_coeffs = {}; self.bondbond_coeffs_style_hint = {};
        self.bondangle_coeffs = {}; self.bondangle_coeffs_style_hint = {};
        self.angleangletorsion_coeffs = {}; self.angleangletorsion_coeffs_style_hint = {};
        self.endbondtorsion_coeffs = {}; self.endbondtorsion_coeffs_style_hint = {};
        self.middlebondtorsion_coeffs = {}; self.middlebondtorsion_coeffs_style_hint = {};
        self.bondbond13_coeffs = {}; self.bondbond13_coeffs_style_hint = {};
        self.angletorsion_coeffs = {}; self.angletorsion_coeffs_style_hint = {};
        self.angleangle_coeffs = {}; self.angleangle_coeffs_style_hint = {};
        self.type_labels_flag = False
        self.atom_type_labels_forward = {}  # {atom type : atom type label}
        self.atom_type_labels_reverse = {}  # {atom type label: atom type}
        self.bond_type_labels_forward = {}  # {bond type: bond type label}
        self.bond_type_labels_reverse = {}  # {bond type label: bond type}
        self.angle_type_labels_forward = {}  # {angle type : angle type label}
        self.angle_type_labels_reverse = {}  # {angle type label : angle type}
        self.dihedral_type_labels_forward = {}  # {dihedral type : dihedral type label}
        self.dihedral_type_labels_reverse = {}  # {dihedral type label : dihedral type}
        self.improper_type_labels_forward = {}  # {improper type : improper type label}
        self.improper_type_labels_reverse = {}  # {improper type label: improper type}
        if pflag: log.out('\n\nFinding energy coeffs ....')
        for ID in files:
            m = files[ID]
            self.header = 'HEADER,' #m.header
            self.type_labels_flag = m.type_labels_flag
            
            # Find new masses and pair coeffs
            self.mass_coeffs_style_hint = m.mass_coeffs_style_hint
            self.pair_coeffs_style_hint = m.pair_coeffs_style_hint
            for i in m.masses:
                new_type = i + self.coeff_offsets[m.filename]['atom']
                self.masses[new_type] = m.masses[i]
                self.pair_coeffs[new_type] = m.pair_coeffs[i]
                if m.atom_type_labels_forward:
                    self.atom_type_labels_reverse[new_type] = m.atom_type_labels_reverse[i]
                    self.atom_type_labels_forward[m.atom_type_labels_reverse[i]] = new_type
            
            # Find new bond coeffs
            self.bond_coeffs_style_hint = m.bond_coeffs_style_hint
            for i in m.bond_coeffs:
                new_type = i + self.coeff_offsets[m.filename]['bond']
                self.bond_coeffs[new_type] = m.bond_coeffs[i]
                if m.bond_type_labels_forward:
                    self.bond_type_labels_reverse[new_type] = m.bond_type_labels_reverse[i]
                    self.bond_type_labels_forward[m.bond_type_labels_reverse[i]] = new_type
                
            # Find new angle, bondbond, and bondangle coeffs
            self.angle_coeffs_style_hint = m.angle_coeffs_style_hint
            self.bondbond_coeffs_style_hint = m.bondbond_coeffs_style_hint
            self.bondangle_coeffs_style_hint = m.bondangle_coeffs_style_hint
            for i in m.angle_coeffs:
                new_type = i + self.coeff_offsets[m.filename]['angle']
                self.angle_coeffs[new_type] = m.angle_coeffs[i]
                if i in m.bondbond_coeffs:
                    self.bondbond_coeffs[new_type] = m.bondbond_coeffs[i]
                if i in m.bondangle_coeffs:
                    self.bondangle_coeffs[new_type] = m.bondangle_coeffs[i]
                if m.angle_type_labels_forward:
                    self.angle_type_labels_reverse[new_type] = m.angle_type_labels_reverse[i]
                    self.angle_type_labels_forward[m.angle_type_labels_reverse[i]] = new_type
                
            # Find new dihedral, angleangletorsion, endbondtorsion,
            # middlebondtorsion, bonddbond13, and angletorsion coeffs
            self.dihedral_coeffs_style_hint = m.dihedral_coeffs_style_hint
            self.angleangletorsion_coeffs_style_hint = m.angleangletorsion_coeffs_style_hint
            self.endbondtorsion_coeffs_style_hint = m.endbondtorsion_coeffs_style_hint
            self.middlebondtorsion_coeffs_style_hint = m.middlebondtorsion_coeffs_style_hint
            self.bondbond13_coeffs_style_hint = m.bondbond13_coeffs_style_hint
            self.angletorsion_coeffs_style_hint = m.angletorsion_coeffs_style_hint
            for i in m.dihedral_coeffs:
                new_type = i + self.coeff_offsets[m.filename]['dihedral']
                self.dihedral_coeffs[new_type] = m.dihedral_coeffs[i]
                if i in m.angleangletorsion_coeffs:
                    self.angleangletorsion_coeffs[new_type] = m.angleangletorsion_coeffs[i]
                if i in m.endbondtorsion_coeffs:
                    self.endbondtorsion_coeffs[new_type] = m.endbondtorsion_coeffs[i]
                if i in m.middlebondtorsion_coeffs:
                    self.middlebondtorsion_coeffs[new_type] = m.middlebondtorsion_coeffs[i]
                if i in m.bondbond13_coeffs:
                    self.bondbond13_coeffs[new_type] = m.bondbond13_coeffs[i]
                if i in m.angletorsion_coeffs:
                    self.angletorsion_coeffs[new_type] = m.angletorsion_coeffs[i]
                if m.dihedral_type_labels_forward:
                    self.dihedral_type_labels_reverse[new_type] = m.dihedral_type_labels_reverse[i]
                    self.dihedral_type_labels_forward[m.dihedral_type_labels_reverse[i]] = new_type
                
            # Find new improper and angleangle coeffs
            self.improper_coeffs_style_hint = m.improper_coeffs_style_hint
            self.angleangle_coeffs_style_hint = m.angleangle_coeffs_style_hint
            for i in m.improper_coeffs:
                new_type = i + self.coeff_offsets[m.filename]['improper']
                self.improper_coeffs[new_type] = m.improper_coeffs[i]
                if i in m.angleangle_coeffs:
                    self.angleangle_coeffs[new_type] = m.angleangle_coeffs[i]
                if m.improper_type_labels_forward:
                    self.improper_type_labels_reverse[new_type] = m.improper_type_labels_reverse[i]
                    self.improper_type_labels_forward[m.improper_type_labels_reverse[i]] = new_type
                
        #---------------#
        # Update ntypes #
        #---------------#
        self.natomtypes = len(self.masses)
        self.nbondtypes = len(self.bond_coeffs)
        self.nangletypes = len(self.angle_coeffs)
        self.ndihedraltypes = len(self.dihedral_coeffs)
        self.nimpropertypes = len(self.improper_coeffs)
        return