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
        log.out('Finding molecules ...')
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


##############################################################################
# Class to recreate m class with desired info, but leave out small molecules #
##############################################################################
class Atom:
    pass # .type .molid .charge .x .y .z .ix. iy .iz .comment .element .molecule

class Bond:
    pass  # .type .atomids = [atom1id, atom2id]
    
class cluster_removal:
    def __init__(self, m, log, delete_atoms, method):
        # Provide attributes of original m class that are not atoms or bonds related
        # ReaxFF Specific info
        self.reaxff_flag = m.reaxff_flag
        if self.reaxff_flag: self.reaxff = m.reaxff # Only add reaxff attribute if reaxff_flag
        
        # Bond dists Specific info
        self.bonddist_flag = m.bonddist_flag
        if self.bonddist_flag: self.bonds_via_dist = m.bonds_via_dist # Only add bond_stats attribute if bonddist_flag
        
        # General info
        self.molecules = m.molecules
        self.elements = m.elements
        self.filename = m.filename
        self.header = m.header
        self.xbox_line = m.xbox_line
        self.ybox_line = m.ybox_line
        self.zbox_line = m.zbox_line
        self.xy = m.xy
        self.xz = m.xz
        self.yz = m.yz
        self.masses = m.masses
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.bond_coeffs = {} # { bondtype-id : coeffs class}
        self.velocities = {}  # {atom number : tuple of velocities}
        self.natoms = 0 # Initialize and update later
        self.nbonds = 0 # Initialize and update later
        self.nbondtypes = 1 # No unique bond-typing will be performed
        self.natomtypes = m.natomtypes
        self.kept_molecules = {} # {molecule-formula : count of molecule formula}
        self.delete_atoms = delete_atoms
        self.total_system_mass = 0 # Initialize and update later
        self.total_system_size = 0 # Initialize and update later
        
        
        ################################
        # Find atoms the meet criteria #
        ################################
        kept_atoms = []; idcounter = 0; # kept atoms log and idcounter to set new ids
        kept_flags = {i:False for i in m.atoms}
        atomid_map = {} # {old-atomid:new-atomid}
        m.atoms = dict(OrderedDict(sorted(m.atoms.items()))) # Re-order m.atoms to keep ids as close as possible
        for i in m.atoms:
            atom = m.atoms[i]
            size = atom.molecule.size
            mass = atom.molecule.mass
            
            # If method is 'extract' use >= for selection
            if method == 'extract':
                # If delete method is mass search using mass
                if delete_atoms['method'] == 'mass':
                    if mass >= delete_atoms['criteria']:
                        if delete_atoms['criteria'] == 0:
                            idcounter = i
                        else: idcounter += 1
                        atomid_map[i] = idcounter
                        kept_atoms.append(i)
                        kept_flags[i]=True
                        
                # If delete method is mass search using mass
                elif delete_atoms['method'] == 'size':
                    if size >= delete_atoms['criteria']:
                        if delete_atoms['criteria'] == 0:
                            idcounter = i
                        else: idcounter += 1
                        atomid_map[i] = idcounter
                        kept_atoms.append(i)
                        kept_flags[i]=True
                        
                # else raise exception
                else: log.error(f"ERROR delete_atoms['method'] method not supported:  {delete_atoms['method']}")
                
            # If method is 'by-product' use < for selection
            elif method == 'by-product':
                # If delete method is mass search using mass
                if delete_atoms['method'] == 'mass':
                    if mass < delete_atoms['criteria']:
                        if delete_atoms['criteria'] == 0:
                            idcounter = i
                        else: idcounter += 1
                        atomid_map[i] = idcounter
                        kept_atoms.append(i)
                        kept_flags[i]=True
                        
                # If delete method is mass search using mass
                elif delete_atoms['method'] == 'size':
                    if size < delete_atoms['criteria']:
                        if delete_atoms['criteria'] == 0:
                            idcounter = i
                        else: idcounter += 1
                        atomid_map[i] = idcounter
                        kept_atoms.append(i)
                        kept_flags[i]=True
                        
                # else raise exception
                else: log.error(f"ERROR delete_atoms['method'] method not supported:  {delete_atoms['method']}")
                
            # else cluster_removal method not supported
            else: log.error(f"ERROR cluster_removal method not supported:  {method}")
            
            
        #############################################
        # Build self.atoms instance with kept atoms #
        #############################################
        # sort kept atoms to keep atomids as close as possible when reseting atomids
        kept_molids = set(); kept_atoms = sorted(kept_atoms) # If no volatiles are removed the ids will stay the same
        for i in kept_atoms:
            self.atoms[atomid_map[i]] = m.atoms[i]
            kept_molids.add(self.atoms[atomid_map[i]].molid)
            try: velocity = m.velocities[i]
            except: velocity = (0, 0, 0)
            self.velocities[atomid_map[i]] = velocity
            
            
        #############################################
        # Build self.bonds instance with kept atoms #
        #############################################
        idcounter = 0;
        for i in m.bonds:
            id1, id2 = m.bonds[i].atomids
            if kept_flags[id1] and kept_flags[id2]:
                idcounter += 1
                id1_mapped = atomid_map[id1]
                id2_mapped = atomid_map[id2]
                
                # Recreate self.bonds from m.bonds (reseting all bondtypes to 1)
                B = Bond()
                B.type = 1 # set all bondtypes as 1
                B.atomids = sorted([id1_mapped, id2_mapped])
                self.bonds[idcounter] = B
                
                
        ###########################
        # Create self.bond_coeffs #
        ###########################
        class Coeff_class: pass  # .type .coeffs = []
        # Set only 1 bond coeff type
        C = Coeff_class()
        C.type = 1
        C.coeffs = []
        self.bond_coeffs[1] = C
        
        
        ######################################
        # Update self.natoms and self.nbonds #
        ######################################
        self.natoms = len(self.atoms); self.nbonds = len(self.bonds);
        
        
        #######################
        # Find kept_molecules #
        #######################
        kept_molids = sorted(kept_molids)
        kept_formulas = [self.molecules.data[i].formula for i in kept_molids]
        unqiue_formulas = sorted(set(kept_formulas), key=len) # Find all unqiue formulas
        for i in unqiue_formulas:
            self.kept_molecules[i] = kept_formulas.count(i)
            
            
        ############################################################
        # Update self.total_system_mass and self.total_system_size #
        ############################################################
        for i in self.atoms:
            self.total_system_mass += self.masses[self.atoms[i].type].coeffs[0]
            self.total_system_size += 1