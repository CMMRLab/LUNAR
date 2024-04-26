# -*- coding: utf-8 -*-
"""
Revision 1.4
January 4th, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
Questions? jdkemppa@mtu.edu
"""
##############################
# Import Necessary Libraries #
##############################          
from collections import OrderedDict                                                                    
                                                                                                               

#############################
### cluster analysis code ###
#############################
class Data:
    pass # .size .mass .pmass .psize .formula

class cluster_analysis:                                                                                                           
    def __init__(self, m, log):
        self.data = {} # { cluster-id : Data object }
        self.atoms = {} # { cluster-id : {set of atoms in cluster id }}
        self.formula = {} # { atom-id : formula atom belongs too }
        self.molids = {} # {atom-id : molid}
        self.clusters = set([]) # { (tuple of atoms in cluster1), (nclusters) }
        
        # Function to find cluster formula
        def find_cluster_formula(cluster, m):
            elements = [m.atoms[atomid].element for atomid in cluster]
            base_elements = list(sorted({i for i in elements}))
            formula = ''
            for element in base_elements:
                formula += '{}{}-'.format(element, elements.count(element))
            return formula[:-1]
            
        # Generate graph
        graph = {i:[] for i in m.atoms}
        for i in m.bonds:
            id1, id2 = m.bonds[i].atomids
            graph[id1].append(id2)
            graph[id2].append(id1)
            
        # Find clusters
        log.out('Finding molecules ....')
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
        self.clusters = sorted(self.clusters, key=lambda x: x[0])    # Sort all clusters based on 1st atomID in cluster
        self.clusters = sorted(self.clusters, key=len, reverse=True) # Sort all clusters by number of atoms
        
        # Find molids after sorting clusters. molid will be count + 1 where the largest
        # cluster will be molid = 1 and smallet cluster will be molid = nclusters in system
        self.natoms = [];  self.cmass = []; 
        for molID, cluster in enumerate(self.clusters, 1):
            self.cmass.append(sum([m.masses[m.atoms[i].type].coeffs[0] for i in cluster])); self.natoms.append(len(cluster)); 
            for ID in cluster:
                self.molids[ID] = molID

        # Summing mass and atoms                                                                                                                                                               
        self.mass_total = sum(self.cmass); self.size_total = sum(self.natoms); 

        # Save cluster info                               
        for n, (size, mass) in enumerate(zip(self.natoms, self.cmass)):  
            # Compute pmass and psize if applicable and find cluster formula
            psize = 0; pmass = 0; # Intialize as zeros and update
            if self.size_total > 0: psize = round(100*size/self.size_total, 2)
            if self.mass_total > 0: pmass = round(100*mass/self.mass_total, 2)
            formula = find_cluster_formula(self.clusters[n], m)
            
            # Save info for for log file
            d = Data()
            d.size = size
            d.mass = mass
            d.pmass = pmass
            d.psize = psize
            d.formula = formula
            self.data[n + 1] = d
            
            # atom atoms info to atoms object
            self.atoms[ n + 1 ] = self.clusters[n]
            
            # Loop through atoms in cluster[n] and assign formula
            for i in self.clusters[n]:
                self.formula[i] = formula



##############################################################
# Function to call cluster_analysis class and update m class #
##############################################################
def add_molecule_data2m(m, log):
    
    # Class cluster_analysis class
    molecules = cluster_analysis(m, log)
    
    # Add molecules attribute to m
    m.molecules = molecules
    
    # Loop through m.atoms and add molecule instance and update molid
    class Molecule: pass # .mass .size .formula
    for i in m.atoms:
        atom = m.atoms[i]
        
        # Find molid and update
        molid = molecules.molids[i]
        atom.molid = molid
        
        # Create a molecule instance for each atom
        M = Molecule()
        M.formula = molecules.formula[i]
        M.mass = molecules.data[molid].mass
        M.size = molecules.data[molid].size
        atom.molecule = M

    return m