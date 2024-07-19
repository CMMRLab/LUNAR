# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 2.1 (2.N's no longer uses Jake's ring code to find cycles')
June 6th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import math
import time


#########################################################
# Function for finding rings/cycles in molecular system #
#########################################################
def find_cycles(m, find_rings, log):
    #------------------------------------------------------------------------------#
    # Generate graph which will be used to find all cycles in the molecular system #
    #------------------------------------------------------------------------------#
    graph = {i:[] for i in m.atoms}
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        graph[id1].append(id2); graph[id2].append(id1);
    
    #----------------------------------------------------------------------------------------#
    # DFS search to find start and end point to use for finding cycles. Implementation:      # 
    # https://stackoverflow.com/questions/40833612/find-all-cycles-in-a-graph-implementation #
    # Meant to be called from within a loop like: for path in dfs(graph, start, end):        #
    # but could be called as paths = list(dfs(graph, start, end)) to get a list of all paths #
    #----------------------------------------------------------------------------------------#
    def dfs(graph, start, end):
        stack = [(start, [])]
        while stack:
            vertex, path = stack.pop()
            if path and vertex == end: yield path; continue;
            for neighbor in graph[vertex]:
                if neighbor in path: continue;
                stack.append((neighbor, path+[neighbor]))
        
    #--------------------------------------------------------------------------------#
    # Functions to reduce the size of a graph to speed up execution of dfs algorithm #
    #--------------------------------------------------------------------------------#
    # Function to find neighs away
    def find_Nneighs_away(atomID, maxdepth, graph):
        neighbors=[[] for i in range(maxdepth)]; 
        neighbors[0]=graph[atomID]; visited=graph[atomID]+[atomID];
        for n, i in enumerate(neighbors):
            for j in neighbors[n]:
                if n+1 < maxdepth:
                    new = []
                    for k in graph[j]:
                        if k not in visited:
                            visited.append(k); new.append(k);
                    neighbors[n+1].extend(new)
        return neighbors
    
    # Function to generate reduced graph
    def reduced_graph(atomID, maxdepth, graph, m, elements2walk):
        rgraph = {}; rgraph[atomID] = graph[atomID];
        neighs = find_Nneighs_away(atomID, maxdepth, graph)
        for n, depth in enumerate(neighs, 1): 
            for ID in depth:
                if n < len(neighs) and m.atoms[ID].element in elements2walk: rgraph[ID] = graph[ID]
                else: rgraph[ID] = []
        return rgraph
    
    #-------------------------------------------------------------------------------------------------#
    # Find rings/cycles in the molecular system based on user inputs (will produce tuples of atomIDs) #
    #-------------------------------------------------------------------------------------------------#
    elements2walk = find_rings['elements2walk']; rings2check = find_rings['rings2check'];
    maxringsize = math.ceil( 0.5*(max(rings2check)+1) ); cycles = set([]) # hold all unqiue cycles
    log.out('Finding rings ...')
    start_time = time.time()
    time_per_atom = 0 # DETDA is about 0.0001 (if over 0.25 generate output log)
    progress_increment = 5; count = 0; natoms = len(graph);
    nevery = math.ceil(natoms*(progress_increment/100))
    checked = {i:False for i in graph}
    use_checked = False; log_progress = False
    for n, atomID in enumerate(graph, 1):
        # log progress
        count += 1
        if  time_per_atom > 0.25:
            log_progress = True
            use_checked = True
        if log_progress and count % nevery == 0 or log_progress and count in [1, natoms]:
            current_percent = int(100*count/natoms)
            elapsed_time = '{:.4f}'.format(time.time() - start_time)
            message = '    Taking longer then expected ... Completed: {} %. Elapsed time: {} seconds'.format(current_percent, elapsed_time)
            log.out(message)
        
        # Find rings
        if use_checked and checked[atomID]: continue
        if len(graph[atomID]) == 1 or m.atoms[atomID].element not in elements2walk: continue 
        rgraph = reduced_graph(atomID, maxringsize, graph, m, elements2walk) # Find reduced graph for speedup
        for path in dfs(rgraph, start=atomID, end=atomID): # Find all cycles in the reduced graph (rgraph)
            if len(path) in rings2check and len(path) >= 2: cycles.add(tuple(sorted(path)))    
            if use_checked:
                for i in path:
                    checked[i] = True
        time_per_atom = (time.time() - start_time)/n
    if use_checked: log.out('Ring analysis was slow. Sped-up by checking if atom was walked prior from previous paths.')
    return list(sorted(cycles)), graph


##########################################
# Class to perform a fused ring analysis #
##########################################
class Data: pass # .size .mass .pmass .psize .formula
class fused_rings:
    def __init__(self, m, find_rings, ringIDs, graph, log):
        self.atom2fusedringIDs = {i:0 for i in m.atoms} # { atomID : [lst of ringIDs atom belongs too ]}
        self.clusters = [] # [ lst of sets to find fused clusters ]
        self.bonded = [] # [ (ringID1, ringID2), ... bonded fused rings ]
        self.data = {} # { cluster-id : Data object }
        self.atoms = set([])  # set of atoms in fused structures
        self.bonds = set([])  # set of bonds in fused structures
        ninitial_zeros = 10 # Set max FusedIDs to intialize w/zeros
        
        #------------------------------------------------#
        # Find fused rings and assign "bond" via ringIDs #
        #------------------------------------------------#
        fused2check = find_rings['fused2check']
        id2_ringIDs = [i for i in ringIDs]
        ringed_graph = {ID:[] for ID in ringIDs}
        checked = {} # for finding clusters in next step
        log.out('Finding fused-ring clusters ...')
        for id1 in ringIDs:
            checked[id1] = False
            ring1 = ringIDs[id1]
            id2_ringIDs.remove(id1)
            if len(ring1) not in fused2check: continue
            for id2 in id2_ringIDs:
                ring2 = ringIDs[id2]
                if id1 == id2 or len(ring2) not in fused2check: continue
                shared_atomIDs = set(ring1).intersection(ring2)
                if len(shared_atomIDs) >= 2: # If ring has two or more shared atomIDs, it must be fused
                    self.bonded.append( tuple(sorted([id1, id2])) )
                    ringed_graph[id1].append(id2); ringed_graph[id2].append(id1);

        #-----------------------------------------#
        # Perform cluster analysis on fused rings #
        #-----------------------------------------#
        clusters = set([])
        for ID in ringIDs:
            if checked[ID]: continue
            visited=set([ID]); queue=[ID];
            while queue:
                s = queue.pop(0) 
                for neighbor in ringed_graph[s]:
                    if checked[neighbor]: continue
                    visited.add(neighbor)
                    queue.append(neighbor)
                    checked[neighbor]=True
            clusters.add( tuple(sorted(visited)) )
        clusters = sorted(clusters, key=lambda x: x[0]) # Sort all clusters based on 1st ringID in cluster
        clusters = sorted(clusters, key=len, reverse=True) # Sort all clusters by number of rings
        
        # Sort clusters based on number of atoms in each cluster
        index2atoms = {n:set() for n in range(len(clusters))} # { index : atoms in cluster }
        for n, cluster in enumerate(clusters):
            for ID in cluster:
                for atomID in ringIDs[ID]:
                    index2atoms[n].add(atomID)
        index2natoms = {i:len(index2atoms[i]) for i in index2atoms} # { index : natoms in cluster }
        index2natoms = dict(sorted(index2natoms.items(), key=lambda x:x[1], reverse=True)) # sort dict by value
        for index in index2natoms:
            self.clusters.append(clusters[index])
        
        #----------------------------------#
        # Functions to analyze fused rings #
        #----------------------------------#
        # Function to compute fused ring mass
        def cluster_size_mass(cluster):
            size = 0; mass = 0; atoms = set();
            for ID in cluster:
                for atomID in ringIDs[ID]:
                    atoms.add(atomID)
            for atomID in atoms:
                mass += m.masses[m.atoms[atomID].type].coeffs[0]; size += 1;
            return size, round(mass, 2)
        
        # Funtion to find fused ring formula
        def fused_ring_formula(cluster):
            formulas = []; formula = ''
            for ringID in cluster:
                atoms = ringIDs[ringID]
                elements = [m.atoms[atomid].element for atomid in atoms];
                base_elements = list(sorted({i for i in elements})); tmp = ''
                for element in base_elements:
                    tmp += '{}{}-'.format(element, elements.count(element))
                formulas.append(tmp[:-1])
            base_formulas = list(sorted({i for i in formulas}))
            for ring in base_formulas:
                formula += '{}:{} '.format(formulas.count(ring), ring)
            return formula
        
        # Function to find bonds in fused cluster
        def find_bonds_in_fused_cluster(cluster):
            atoms = set([]); bonds = set([]);
            for ringID in cluster:
                atoms.update(ringIDs[ringID])
                for id1 in ringIDs[ringID]:
                    for id2 in graph[id1]:
                        bonds.add( tuple(sorted([id1, id2])) )
            return atoms, bonds
        
        # Initalize info with zeros for nFusedIDs=ninitial_zeros
        for n in range(ninitial_zeros):  
            d = Data()
            d.ringIDs = []; d.nrings = 0;
            d.prings = 0; d.size = 0;
            d.mass = 0; d.formula = 'Initial';
            d.pmass = 0; d.psize = 0;
            self.data[n+1] = d   
            
        #------------------------#
        # Analyze fused clusters #
        #------------------------#
        total_rings = len(ringIDs)
        for n, cluster in enumerate(self.clusters, 1):
            size, mass = cluster_size_mass(cluster)
            formula = fused_ring_formula(cluster)
            atoms, bonds = find_bonds_in_fused_cluster(cluster)
            nrings = len(cluster)
            
            # Compute pmass, psize, and prings
            pmass = 0; psize = 0; prings = 0;
            if m.total_system_mass > 0: pmass = round(100*mass/m.total_system_mass, 2)
            if m.total_system_mass > 0: psize = round(100*size/m.total_system_size, 2)
            if total_rings > 0: prings = round(100*nrings/total_rings, 2)
            
            # Save info for for log file
            d = Data()
            d.ringIDs = cluster
            d.nrings = nrings
            d.prings = prings
            d.size = size
            d.mass = mass
            d.formula = formula
            d.pmass = pmass
            d.psize = psize
            self.data[n] = d
            self.atoms.update(atoms)
            self.bonds.update(bonds)
            
            # Add fusedringID "n" to each atomID in each fused ringed cluster
            for ID in cluster:
                for atomID in ringIDs[ID]:
                    self.atom2fusedringIDs[atomID] = n
        

##################################################################################
# Class to call find_cycles and fused_rings and combined results into one object #
##################################################################################
class Data1: pass #  .count .pcount
class Partitioned: pass #  .size .mass .pmass .psize
class DataRingIDs: pass # .size .mass .pmass .psize .formula
class ring_analysis:
    def __init__(self, m, find_rings, log):
        self.atoms_partitioned = {i:0 for i in m.atoms} # { atom number : partitioned ring size }
        self.partitioned_count = 0  # Count of atoms that have been partitioned
        self.rings_set = {} # { atom number : {set of rings atom is in (NO duplicates)} }
        self.rings_lst = {} # { atom number : [list of rings atom is in (duplicates)] }
        self.cycles = [] # {(id1, id2, id3, ...ringsize), ... nwalked rings}
        self.count = {}  # { ring size : count of rings }
        self.total = 0   # number of rings in the molecular system
        if 'all' in find_rings['rings2check']: find_rings['rings2check'] = [i+3 for i in range(len(m.atoms))]

        #----------------------------------------------------#
        # Find rings/cycles info and set per atom quantities #
        #----------------------------------------------------#
        # Find all cycles in system based in inputs found in find_rings
        self.cycles, self.graph = find_cycles(m, find_rings, log)
        self.total = len(self.cycles)
        self.ringIDs = {n:list(cycle) for n, cycle in enumerate(self.cycles, 1)} # {cycle Index/ringID : [lst of atomIDs]}
        self.ringformulas = {i:[] for i in m.atoms}  # { atomID : [ringformulas] }
        self.atom2ringIDs = {} # { atomID : [lst of ringIDs atom belongs too ]}
        self.clusters = {}     # { ringID : Data object }

        # Create data set with atom and ring size                                        
        for ID in m.atoms:
            self.atom2ringIDs[ID] = []
            self.rings_set[ID] = set([])
            self.rings_lst[ID] = []

        # Loop through cycles to set ringsize atom belongs to and tally ringsizes
        self.count = {i:0 for i in find_rings['rings2check']}
        for cycle in self.cycles:
            self.count[len(cycle)] += 1
            for ID in cycle:
                self.rings_lst[ID].append(len(cycle))
                self.rings_set[ID].add(len(cycle))
                
        # Add ringID to self.atom2ringIDs dict
        for ringID in self.ringIDs:
            ring = self.ringIDs[ringID]
            for atomID in ring:
                self.atom2ringIDs[atomID].append(ringID)
                
        # Partition atoms that belong to multiple rings or keep as zero if not in ring
        for i in self.rings_set:
            rings = self.rings_set[i]
            if len(rings) > 1: self.partitioned_count += 1
            if int(6) in rings: self.atoms_partitioned[i] = 6
            elif int(5) in rings: self.atoms_partitioned[i] = 5
            elif int(7) in rings: self.atoms_partitioned[i] = 7
            elif int(4) in rings: self.atoms_partitioned[i] = 4
            elif int(3) in rings: self.atoms_partitioned[i] = 3
            elif int(8) in rings: self.atoms_partitioned[i] = 8
            elif len(rings) > 0: self.atoms_partitioned[i] = min(rings)
                
        # Find fused rings if user wants
        if find_rings['fused-rings']: self.fused = fused_rings(m, find_rings, self.ringIDs, self.graph, log)

        #--------------------------------------------------------------------#
        # Analyze all aspects of the "ringed" nature of the molecular system #
        #--------------------------------------------------------------------#
        # New data structure to hold all of the system ring data
        self.data = {} # { ring size : Data Object -> partioned object dict with 'element' keys}
        elements_all = m.elements + ['all'] # Add all option to elements
        for i in self.count:  
            if self.total == 0: pcount = 0
            else: pcount = round(100*self.count[i]/self.total, 2)
            d = Data1()
            d.count = self.count[i]; d.pcount = pcount;
            tmp_dict = {} # { 'element symbol' : Partitioned object }
            for j in elements_all:
                p = Partitioned()
                p.size = 0; p.mass = 0;
                p.pmass = 0; p.psize = 0;
                tmp_dict[j] = p
            d.partitioned = tmp_dict
            self.data[i] = d

        # Go through and tally partioned rings data
        for i in self.atoms_partitioned:
            ring = self.atoms_partitioned[i]
            element = m.atoms[i].element
            mass = m.masses[m.atoms[i].type].coeffs[0]
            if ring not in find_rings['rings2check'] or ring == 0: continue # skip rings not of interest or zeroes
            
            # Tally mass and natoms (size) per element and 'all' key
            self.data[ring].partitioned[element].size += 1
            self.data[ring].partitioned[element].mass += mass
            self.data[ring].partitioned['all'].size += 1
            self.data[ring].partitioned['all'].mass += mass
            
        # Find percents and update
        for i in self.data:
            for j in self.data[i].partitioned:
                # Find mass and size and then compute pmass and psize
                pmass = 0; psize = 0;
                mass = self.data[i].partitioned[j].mass
                size = self.data[i].partitioned[j].size
                if m.total_system_mass > 0: pmass = round(100*mass/m.total_system_mass, 2)
                if m.total_system_size > 0: psize = round(100*size/m.total_system_size, 2)
                
                # Update self.data[ring size].partitioned['element']
                self.data[i].partitioned[j].pmass = pmass
                self.data[i].partitioned[j].psize = psize
                
        #---------------------------#
        # Find ringed clusters info #
        #---------------------------#
        # Function to find cluster formula
        def find_cluster_formula(cluster, m):
            elements = [m.atoms[atomid].element for atomid in cluster];
            base_elements = list(sorted({i for i in elements})); formula = '';
            for element in base_elements:
                formula += '{}{}-'.format(element, elements.count(element))
            return formula[:-1]

        # Find ringed clusters info from self.cycles
        for n, ring in enumerate(self.cycles, 1):
            mass = sum([m.masses[m.atoms[i].type].coeffs[0] for i in ring])
            formula = find_cluster_formula(ring, m)
            size = len(ring); pmass = 0; psize = 0;
            if m.total_system_mass > 0: pmass = round(100*mass/m.total_system_mass, 2)
            if m.total_system_size > 0: psize = round(100*size/m.total_system_size, 2)
            dr = DataRingIDs()
            dr.size = size
            dr.mass = mass
            dr.pmass = pmass
            dr.psize = psize
            dr.atoms = ring
            dr.formula = formula
            self.clusters[n] = dr
            for i in ring:
                self.ringformulas[i].append(formula)
                
                
#####################################################
# Function to call all ring codes and add data to m #
#####################################################
def add_ring_data2m(m, find_rings, log):
    
    # Call rings class from above that will call other auxiliary ring class
    rings_data = ring_analysis(m, find_rings, log)
    
    # add rings and find_rings attributes to m
    m.rings = rings_data; m.find_rings = find_rings;
    
    # add .rings .ring .ringID .ringformula to m.atoms[ID]
    # .rings = [lst or ringsizes atom belongs too]
    # .ring = integer value of partioned ring
    # .ringID = ringID (zero if not logical)
    # .ringformula = ring formula (blank if not logical)
    for i in m.atoms:
        atom = m.atoms[i]
        atom.rings = sorted(rings_data.rings_set[i])
        atom.ring = rings_data.atoms_partitioned[i]
        

        # Get ringID for atom if it is logical. It is not logical if the atom
        # does not belong to a ring or if the atom belongs to two rings in a
        # fused structure. So ringID is only for atoms that belong to 1-ring.
        if len(rings_data.atom2ringIDs[i]) == 1:
            ringID = rings_data.atom2ringIDs[i][0] # only use the 1st ringID
            ringformula = rings_data.clusters[ringID].formula # only use 1st ring forumla
            atom.ringID =  ringID 
            atom.ringformula = ringformula
        else: atom.ringID = 0; atom.ringformula = '';
        
        # If fused rings are found add fusedringID to atom; else set as zero
        if find_rings['fused-rings']: atom.fusedringID = rings_data.fused.atom2fusedringIDs[i]
        else: atom.fusedringID = 0
    return m