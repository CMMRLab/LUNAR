# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.5
April 10th, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

##############################
# Import Necessary Libraries #
##############################
import src.io_functions as io_functions
import src.periodicity as periodicity
import src.read_lmp as read_lmp
import numpy as np
import math
import time
import sys
import os



###############################################
# class for finding and storing molecule info #
###############################################
class Info: pass # .size .mass .psize .pmass .atoms
class analysis:
    def __init__(self, topofile, N0, txtfile, fav, log=None): 
        start_time = time.time()       
        
        # Configure log (default is level='production', switch to 'debug' if debuging)
        if log is None:
            log = io_functions.LUNAR_logger()
        log.configure(level='production')
        #log.configure(level='debug')
        
        ########################################################
        # set version and print starting information to screen #
        ########################################################
        version = 'v1.5 / 10 April 2025'
        log.out(f'\n\nRunning cluster_analysis {version}')
        log.out(f'Using Python version {sys.version}')
        
        ##########################################################
        # Read LAMMPS datafile into memory as class "m" applying #
        # forward mapping of type labels (if applicable)         #
        ##########################################################
        if os.path.isfile(topofile):
            m = read_lmp.Molecule_File(topofile, method='forward', sections = ['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers'])
            log.out(f'Read in {m.filename}')
        else:
            log.error(f'ERROR lammps datafile: {topofile} does not exist')
        basename = topofile # Strip path from topofile
        if basename.endswith('data.gz'):
            basename = basename.replace('.gz', '') # remove .gz suffix
        basename = basename[:basename.rfind('.')] # Find file basename
        log.debug(f'basename = {basename}')
        
        
        #######################################
        # Unwrap atoms and get box dimensions #
        #######################################
        # Generate h and h_inv vector like LAMMPS does
        h, h_inv, box = periodicity.get_box_parameters(m)
        m = periodicity.unwrap_atoms(m)
        write_lmp = False
        if write_lmp:
            # Write LAMMPS datafile of unwrapped atoms to check how code did
            import src.write_lmp as write_lmp
            header = 'HEADER, Unwrap Testing'
            atom_style = 'full'
            include_type_labels = False
            write_lmp.file(m, basename+'_UNWRAPPED.data', header, atom_style, include_type_labels, log)
        
        
        ##############################################
        # Setup attributes of cluster_analysis class #
        ##############################################
        self.clusters = set([]) # { (tuple of atoms in cluster1), (nclusters) }
        self.info = {} # { molID : Info object }
        self.p = 0   # Extent of reaction (conversion)
        self.pg = 0  # Critical Extent of reaction (gel point)
        self.Xn = 0  # Degree of polymerization
        self.Mw = 0  # Weight-average molar mass
        self.Mn = 0  # Number-average molar mass
        self.Mz = 0  # Higher-average molar mass 
        self.Mz1 = 0 # Higher-average z+1 molar mass 
        self.RMW = 0 # Weight-averge reduced molecular weight
        self.filename = topofile
        nmol_initial_zeros = 10 # Set max molIDs to intialize w/zeros


        ##################################
        # Find clusters and analyze them #
        ##################################
        # Generate graph
        graph = {i:[] for i in m.atoms}
        for i in m.bonds:
            id1, id2 = m.bonds[i].atomids
            graph[id1].append(id2)
            graph[id2].append(id1)
        
        # Find clusters
        checked = {ID:False for ID in m.atoms}
        for ID in graph:
            if checked[ID]: continue
            visited=set([ID]); queue=[ID];
            while queue:
                s = queue.pop(0) 
                for neighbor in graph[s]:
                    if checked[neighbor]: continue
                    visited.add(neighbor)
                    queue.append(neighbor)
                    checked[neighbor]=True
            self.clusters.add( tuple(sorted(visited)) )
        self.clusters = sorted(self.clusters, key=lambda x: x[0]) # Sort all clusters based on 1st atomID in cluster
        self.clusters = sorted(self.clusters, key=len, reverse=True) # Sort all clusters by number of atoms

        
        # Find mapping of atomIDs to bondIDs
        atoms2bondIDs = {i:set() for i in m.atoms} # {atomIDs : set(bondIDs)}
        for i in m.bonds:
            id1, id2 = m.bonds[i].atomids
            atoms2bondIDs[id1].add(i)
            atoms2bondIDs[id2].add(i)
            
        # Map bonds to cluster index
        clusters_bonds = []
        for cluster in self.clusters:
            bondIDs = set()
            for i in cluster:
                bondIDs.update(atoms2bondIDs[i])
            clusters_bonds.append(bondIDs)
    
        # Analyze clusters
        natoms = [];  cmass = []; periodic_flags = []; rgs = []
        for cluster, bonds in zip(self.clusters, clusters_bonds):
            rgs.append(self.compute_radius_of_gyration(m, cluster))
            natoms.append(len(cluster));
            cmass.append(self.getmass(m, cluster))   
            periodic_flags.append(self.cluster_periodicity(m, cluster, bonds, h, h_inv))
                  
        # Summing mass and atoms                                                                                                                                                               
        self.mass_total = sum(cmass); self.size_total = sum(natoms);   

        # Initalize info with zeros for nmolIDs=nmol_initial_zeros
        for n in range(nmol_initial_zeros):  
            I = Info()
            I.rg = 7*[0]
            I.atoms = []
            I.size = 0
            I.mass = 0
            I.psize = 0
            I.pmass = 0
            I.periodic = [False, False, False]
            self.info[n+1] = I                                                                              
        
        # Log Info inside an oject acessible via molID in a dictionary  
        self.periodic_count = {'x':int(0), 'y':int(0), 'z':int(0)} # {'dir':count}   
        self.rg_log = {'Rg':[], 'Rxx':[], 'Ryy':[], 'Rzz':[], 'Rxy':[], 'Rxz':[], 'Ryz':[]}                         
        for n, (size, mass, periodic, rg) in enumerate(zip(natoms, cmass, periodic_flags, rgs)):
            # Compute pmass and psize if applicable
            psize = 0; pmass = 0; # Intialize as zeros and update
            if self.size_total > 0: psize = 100*size/self.size_total
            if self.mass_total > 0: pmass = 100*mass/self.mass_total
                
            # Log cluster info
            I = Info()
            I.rg = rg
            I.atoms = self.clusters[n]
            I.size = size
            I.mass = mass
            I.psize = psize
            I.pmass = pmass
            I.periodic = periodic
            self.info[n+1] = I
            
            # Count periodically spanning clusters
            if periodic[0]: self.periodic_count['x'] += int(1)
            if periodic[1]: self.periodic_count['y'] += int(1)
            if periodic[2]: self.periodic_count['z'] += int(1)
            
            # Log different Rgs
            self.rg_log['Rg'].append(rg[0])
            self.rg_log['Rxx'].append(rg[1])
            self.rg_log['Ryy'].append(rg[2])
            self.rg_log['Rzz'].append(rg[3])
            self.rg_log['Rxy'].append(rg[4])
            self.rg_log['Rxz'].append(rg[5])
            self.rg_log['Ryz'].append(rg[6])
            
        # Compute radius of gyration stats
        self.Rg_mean = self.compute_mean(self.rg_log['Rg'])
        self.Rg_stdev = self.compute_standard_deviation(self.rg_log['Rg'])
        
        self.Rxx_mean = self.compute_mean(self.rg_log['Rxx'])
        self.Rxx_stdev = self.compute_standard_deviation(self.rg_log['Rxx'])
        self.Ryy_mean = self.compute_mean(self.rg_log['Ryy'])
        self.Ryy_stdev = self.compute_standard_deviation(self.rg_log['Ryy'])
        self.Rzz_mean = self.compute_mean(self.rg_log['Rzz'])
        self.Rzz_stdev = self.compute_standard_deviation(self.rg_log['Rzz'])
        
        self.Rxy_mean = self.compute_mean(self.rg_log['Rxy'])
        self.Rxy_stdev = self.compute_standard_deviation(self.rg_log['Rxy'])
        self.Rxz_mean = self.compute_mean(self.rg_log['Rxz'])
        self.Rxz_stdev = self.compute_standard_deviation(self.rg_log['Rxz'])
        self.Ryz_mean = self.compute_mean(self.rg_log['Ryz'])
        self.Ryz_stdev = self.compute_standard_deviation(self.rg_log['Ryz'])
                
        # Find Mw, Mn, Mz, Mz1, and RMW from self.info[molID].mass. Equations from:
        #    Polymers: Chemistry and Physics of Modern Materials pg 8-9 &
        #    Multiscale Modeling for Virtual Manufacturing of Thermoset Composites
        ni = 0 # tally of number of clusters
        ni_mi1 = 0 # tally of ni_clusters*mi_cluster
        ni_mi2 = 0 # tally of ni_clusters*mi_cluster^2
        ni_mi3 = 0 # tally of ni_clusters*mi_cluster^3
        ni_mi4 = 0 # tally of ni_clusters*mi_cluster^4
        RMW_ni_mi = 0  # tally of ni_clusters*mi_cluster for RMW calculation
        RMW_ni_mi2 = 0 # tally of ni_clusters*mi_cluster^2 for RMW calculation
        for n in self.info:
            mass = self.info[n].mass
            ni += 1
            ni_mi1 += 1*mass
            ni_mi2 += 1*mass**2
            ni_mi3 += 1*mass**3
            ni_mi4 += 1*mass**4
            if n > 1: # skip over 1st largest cluster from tally
                RMW_ni_mi += 1*mass
                RMW_ni_mi2 += 1*mass**2
        # Only compute if div by zero is not possible else leave as default zero's
        if ni > 0: self.Mn = ni_mi1/ni
        if ni_mi1 > 0: self.Mw = ni_mi2/ni_mi1
        if ni_mi2> 0: self.Mz = ni_mi3/ni_mi2
        if ni_mi3 > 0: self.Mz1 = ni_mi4/ni_mi3
        if RMW_ni_mi > 0: self.RMW = RMW_ni_mi2/RMW_ni_mi 

        
        # Compute Xn and p from N for poly condensation reaction. Equations from:
        # Polymers: Chemistry and Physics of Modern Materials pg 32
        N = len([n for n in self.info if self.info[n].size > 0]) # current number of molecules
        if N0 > 0 and fav > 0:
            self.p = (2*(N0-N))/(N0*fav) # extent of reaction
            self.pg = 2/fav
        if self.p > 0 and fav > 0:
            self.Xn = 2/(2-self.p*fav) # degree of polymerization
        
        # Print out table of molIDs
        log.out('\n\n----------------------------------------------------------Cluster Analysis----------------------------------------------------------')
        log.out('{:^6} {:^10} {:^8} {:^15} {:^8} {:^32} {:28} {:^18}'.format('molID', 'size', '%size', 'mass', '%mass', 'GyRadius (Rg, Rxx, Ryy, Rzz)', 'GyRadius^2 (Rxy, Rxz, Ryz)', 'periodic (x, y, z)'))
        log.out('------------------------------------------------------------------------------------------------------------------------------------')  
        for n in self.info:
            molID = self.info[n]
            if molID.size > 0:
                pbc = ', '.join([str(i)[0] for i in molID.periodic])
                rg_ii = ', '.join(['{:.2f}'.format(i) for i in molID.rg[:4]])
                rg_ij = ', '.join(['{:.2f}'.format(i) for i in molID.rg[4:]])
                log.out('{:^6} {:^10} {:^8.2f} {:^15.2f} {:^8.2f} {:^32} {:^28} {:^18}'.format(n, molID.size, molID.psize, molID.mass, molID.pmass, rg_ii, rg_ij, pbc))
            
        # Print out Mn, Mw, Mz, Mz+1 Xn, and p
        log.out('\n\n------------------------Distributions------------------------')
        log.out('{} {:.4f}'.format('Extent of reaction aka converison      (p)   : ', self.p))
        log.out('{} {:.4f}'.format('Critical extent of reaction gel point  (pg)  : ', self.pg))
        log.out('{} {:.4f}'.format('Degree of Polymerization               (Xn)  : ', self.Xn))
        log.out('{} {:.4f}'.format('Weight-average molar mass              (Mw)  : ', self.Mw))
        log.out('{} {:.4f}'.format('Number-average molar mass              (Mn)  : ', self.Mn))
        log.out('{} {:.4f}'.format('Higher-average molar mass              (Mz)  : ', self.Mz))
        log.out('{} {:.4f}'.format('Higher-average molar mass              (Mz+1): ', self.Mz1))
        log.out('{} {:.4f}'.format('Weight-averge reduced molecular weight (RMW) : ', self.RMW))        
        log.out('{} {}'.format('Count of infinitely spanned clusters     (x) : ', self.periodic_count['x']))
        log.out('{} {}'.format('Count of infinitely spanned clusters     (y) : ', self.periodic_count['y']))
        log.out('{} {}'.format('Count of infinitely spanned clusters     (z) : ', self.periodic_count['z']))
        
        log.out('\n\n----------------------Radius-of-Gyration---------------------') 
        log.out('{} {:.4f}'.format('Radius of Gyration                  mean(Rg) : ', self.Rg_mean))
        log.out('{} {:.4f}'.format('Radius of Gyration                 stdev(Rg) : ', self.Rg_stdev))
        
        log.out('{} {:.4f}'.format('Radius of Gyration                 mean(Rxx) : ', self.Rxx_mean))
        log.out('{} {:.4f}'.format('Radius of Gyration                stdev(Rxx) : ', self.Rxx_stdev))
        log.out('{} {:.4f}'.format('Radius of Gyration                 mean(Ryy) : ', self.Ryy_mean))
        log.out('{} {:.4f}'.format('Radius of Gyration                stdev(Ryy) : ', self.Ryy_stdev))
        log.out('{} {:.4f}'.format('Radius of Gyration                 mean(Rzz) : ', self.Rzz_mean))
        log.out('{} {:.4f}'.format('Radius of Gyration                stdev(Rzz) : ', self.Rzz_stdev))
        log.out('{} {:.4f}'.format('Radius of Gyration^2               mean(Rxy) : ', self.Rxy_mean))
        log.out('{} {:.4f}'.format('Radius of Gyration^2              stdev(Rxy) : ', self.Rxy_stdev))
        log.out('{} {:.4f}'.format('Radius of Gyration^2               mean(Rxz) : ', self.Rxz_mean))
        log.out('{} {:.4f}'.format('Radius of Gyration^2              stdev(Rxz) : ', self.Rxz_stdev))
        log.out('{} {:.4f}'.format('Radius of Gyration^2               mean(Ryz) : ', self.Ryz_mean))
        log.out('{} {:.4f}'.format('Radius of Gyration^2              stdev(Ryz) : ', self.Ryz_stdev))

            
        # Write .txt file if desired
        if txtfile:
            log.write_logged(basename+'_clusters.txt')
            path = os.path.dirname(os.path.abspath(topofile))
            log.out(f'\n\nAll outputs can be found in {path} directory')
                
        # Print completion of code
        log.out('\n\nNormal program termination\n\n')
        
        # Script run time
        execution_time = (time.time() - start_time)
        log.out('Execution time in seconds: ' + str(execution_time))
        
        # Show number of warnings and errors
        log.out_warnings_and_errors()
            
        
    # Method to compute cluster mass                                  
    def getmass(self, m, cluster):                                                                                                                  
        return sum([m.masses[m.atoms[i].type].coeffs[0] for i in cluster])
    
    # Method to find cluster periodicity
    def cluster_periodicity(self, m, cluster, bonds, h, h_inv):
        # Initialize Flags
        periodic_flags = [False, False, False]
        
        # Check bond lengths
        half_box_lambda_space = 0.5
        for i in bonds:
            id1, id2 = m.bonds[i].atomids
            atom1 = m.atoms[id1]
            atom2 = m.atoms[id2]
            pos1 = np.array([atom1.x, atom1.y, atom1.z])
            pos2 = np.array([atom2.x, atom2.y, atom2.z])
            delta = pos2 - pos1
            frac = delta @ h_inv
            dx, dy, dz = frac
            if abs(dx) >= half_box_lambda_space: periodic_flags[0] = True
            if abs(dy) >= half_box_lambda_space: periodic_flags[1] = True
            if abs(dz) >= half_box_lambda_space: periodic_flags[2] = True    
        return periodic_flags
    
    # Method to compute the radius of gyration of a cluster
    # https://en.wikipedia.org/wiki/Radius_of_gyration
    # https://docs.lammps.org/compute_gyration.html
    def compute_radius_of_gyration(self, m, cluster):
        # Find center of mass
        M, xm, ym, zm = 0, 0, 0, 0
        for i in cluster:
            atom = m.atoms[i]
            mass = m.masses[atom.type].coeffs[0]
            M += mass
            xm += mass*atom.x
            ym += mass*atom.y
            xm += mass*atom.z
        com_x = xm/M
        com_y = ym/M
        com_z = zm/M
        
        # Compute radius of gyration
        Rg2 = 0
        Rg2_xx, Rg2_yy, Rg2_zz = 0, 0, 0
        Rg2_xy, Rg2_xz, Rg2_yz = 0, 0, 0
        for i in cluster:
            atom = m.atoms[i]
            mass = m.masses[atom.type].coeffs[0]
            
            rx = atom.x - com_x
            ry = atom.y - com_y
            rz = atom.z - com_z
            Rg2 += mass*(rx*rx + ry*ry + rz*rz)
            
            Rg2_xx += mass*(rx*rx)
            Rg2_yy += mass*(ry*ry)
            Rg2_zz += mass*(rz*rz)
            
            Rg2_xy += mass*(rx*ry)
            Rg2_xz += mass*(rx*rz)
            Rg2_yz += mass*(ry*rz)
        
        Rg = math.sqrt(Rg2/M)
        Rg_xx = math.sqrt(Rg2_xx/M)
        Rg_yy = math.sqrt(Rg2_yy/M)
        Rg_zz = math.sqrt(Rg2_zz/M)
        Rg2_xy = Rg2_xy/M
        Rg2_xz = Rg2_xz/M
        Rg2_yz = Rg2_yz/M
        Rfull = [Rg, Rg_xx, Rg_yy, Rg_zz, Rg2_xy, Rg2_xz, Rg2_yz]
        return Rfull
    
    #######################################################################
    # Methods for computing statistic values to avoid numpy as dependancy #
    #######################################################################
    def compute_mean(self, data):
        if data:
            return sum(data)/len(data)
        else:
            return 0
     
    def compute_variance(self, data):
        if data:
          mean = self.compute_mean(data)
          deviations = [(x - mean)**2 for x in data]
          variance = sum(deviations)/len(data)
          return variance
        else:
            return 0
     
    def compute_standard_deviation(self, data):
      variance = self.compute_variance(data)
      return math.sqrt(variance)  
    
