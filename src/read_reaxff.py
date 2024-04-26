# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 3.0
June 2nd, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import math


#######################################
# Class to read ReaxFF bondorder file #
#######################################
class Atoms: pass # .type .nb .abo .nlp .q
class Steps: pass # .atoms .bonds
class read_BO_file:
    def __init__(self, bondfile, log):
        self.steps = {} # { step : Step OBJECT }
        self.pairs = set([]) # { sorted(ID1, ID2), nbonds } # GLOBAL ACROSS ALL TIMESTEPS
        self.atoms = {} # { atomID : Atom OBJECT } # GLOBAL ACROSS ALL TIMESTEPS
        
        # Open and read file ReaxFF bondfile
        with open(bondfile, 'r') as f:
            
            # Itializatizing flags
            data_flag = False
            
            # Loop through file
            log.out('Reading ReaxFF Bond order file and finding bonds ....')
            for line in f:
                
                # split and strip line
                string = line
                line = line.strip()
                line = line.split()
                
                # Setting flags
                if '#' in string:
                    data_flag = False
                    if 'Timestep' in line:
                        step = int(line[2]) # timestep of data
                        atoms = {} # { atomID    : atoms object }
                        bonds = {} # { (ID1,ID2) : [lst of BOs] }
                    continue
                elif '#' not in string:
                    data_flag = True
                    
                # Find data
                if data_flag and len(line) >= 3:                        
                    # Find values that are always constant
                    atomID = int(line[0])
                    atomtypeID = int(line[1])
                    nb = int(line[2])
                    
                    # log atom info
                    a = Atoms()
                    a.type = atomtypeID
                    a.nb = nb
                    a.abo = float(line[-3])
                    a.nlp = float(line[-2])
                    a.q = float(line[-1])
                    atoms[atomID] = a
                    self.atoms[atomID] = a
                    
                    # Find pair(s)/bo(s) and then find bond(s) and log bo(s)
                    pairs = line[3:nb+3]; bos = line[nb+4:2*nb+4];
                    for bo, ID in zip(bos, pairs):
                        ID = int(ID); bo = float(bo);
                        if ID < atomID: bond = (ID, atomID)
                        else: bond = (atomID, ID)
                        if bond in bonds: bonds[bond].append(bo)
                        else: bonds[bond] = [bo]
                        self.pairs.add(bond)
                
                # Log atoms and bonds dict of current step
                s = Steps()
                s.atoms = atoms
                s.bonds = bonds
                self.steps[step] = s


################################################################
# Statistic's functions to use for analyzing bond stats. *NOTE #
# not using numpy to make this code have zero dependancies*    #
################################################################
def compute_mean(data):
  return sum(data)/len(data)
 
def compute_variance(data):
  mean = compute_mean(data)
  deviations = [(x - mean)**2 for x in data]
  variance = sum(deviations)/len(data)
  return variance
 
def compute_standard_deviation(data):
  variance = compute_variance(data)
  return math.sqrt(variance)


#####################################################
# Function to create bonds by averaging bond orders #
# over timesteps and applying bonding constraints   #
#####################################################
class Stats: pass # .count .avg .min .max .std .cutoff
class avgs: pass # .nlps .abos
class create_bonds:
    def __init__(self, bondfile, bond_info, bondorder, log, maxstep='all'):
        self.timesteps = []      # list of timesteps found
        self.bonds = []          # list of found bonds
        self.flagged_bonds = []  # list of flagged bonds
        self.statistics = {}     # { tuple bond type : stats object }
        self.abo_stats = {}      # { element symbol : stats object }
        
        # Read file
        BO = read_BO_file(bondfile, log)
        timesteps = sorted(list(BO.steps.keys()))
        
        # Update maxstep if maxstep defaulted to 'all'
        if maxstep == 'all': maxstep = max(timesteps)
        self.timesteps = [i for i in timesteps if i <= maxstep]
        
        # Create sorted bondorder and statistic dictionaries
        sorted_bo = {} # { tuple(sorted(elem1, elem2)) : BO cut-off }
        statistic = {} # { tuple(sorted(elem1, elem2)) : [lst of avgs] }
        for bond in bondorder:
            sortedbond = bond
            if bond != 'unknown': sortedbond = tuple(sorted(bond))
            sorted_bo[sortedbond] = bondorder[bond]; statistic[sortedbond] = [];
            
        # Set division factor as 2*ntimesteps, b/c bo's will appear twice for same pair
        div = 2*len(self.timesteps)
        
        # Iterate through steps to find all bonds data
        bonds = {pair:[] for pair in BO.pairs}
        abos = {}; graph = {};
        for atom in BO.atoms:
            abos[atom] = []; graph[atom] = [];
        for step in BO.steps:
            stepatoms = BO.steps[step].atoms
            stepbonds = BO.steps[step].bonds
            if step > maxstep: continue
            for atom in stepatoms:
                abos[atom].append(stepatoms[atom].abo)
            for bond in stepbonds:
                bonds[bond].extend( stepbonds[bond] )
        
        # Average together Bond orders from all read in timesteps of Bond order data
        bond_BO_avg = {i:sum(bonds[i])/div for i in bonds} # { tuple(id1, id2): bo time-average }
        
        # Generate graph
        for ID1, ID2 in bond_BO_avg:
            graph[ID1].append(ID2)
            graph[ID2].append(ID1)
        
        # Start finding bonds that meet criteria
        bonds = set([]); flagged_bonds = set([]); atoms = set(); elements = set();
        for ID1 in graph:  
            # Find average bond orders for bonded atoms and sort dictionary in descending order of BOs
            avg_bo = {} # { (ID1,ID2) : avgBO }
            for ID2 in graph[ID1]:
                if ID1 < ID2: bond = (ID1, ID2)
                else: bond = (ID2, ID1)
                avg_bo[bond] = bond_BO_avg[bond]
            avg_bo = dict(sorted(avg_bo.items(), key=lambda x:x[1], reverse=True ))
            
            # Find element and max_nb based on atom type
            element, maxnb = bond_info[BO.atoms[ID1].type]
            atoms.add(ID1); elements.add(element);
            
            # Create bonds based on user-defined cut-offs
            countnb = 0
            for bond in avg_bo:
                avgBO = avg_bo[bond]; countnb += 1;
                
                # Find elements 1 and 2 to determine bond type
                element1, maxnb1 = bond_info[BO.atoms[bond[0]].type]
                element2, maxnb2 = bond_info[BO.atoms[bond[1]].type]
                if element1 < element2: bondtype = (element1, element2)
                else: bondtype = (element2, element1)
                
                # Find minBO (Intialized as 0.3 and updated if found)
                minBO = 0.3 # Default if not found in user input
                if bondtype in sorted_bo: minBO = sorted_bo[bondtype]
                elif 'unknown' in sorted_bo: minBO = sorted_bo['unknown']
                
                # If number of bonds is less than specified and BO is higher then minimum create the bond
                if countnb <= maxnb and avgBO >= minBO:
                    bonds.add(bond); statistic[bondtype].append(avgBO)
                else: flagged_bonds.add(bond)
                
        # Removing any flagged bonds if they were created      
        bonds = list(set(bonds) - set(flagged_bonds)) 
        
        # Sorting bonds and flagged bonds                
        bonds = sorted(bonds); flagged_bonds = sorted(flagged_bonds) 
        self.bonds = bonds; self.flagged_bonds = flagged_bonds;
        
        # Find bond order statistics
        for i in statistic:
            lst = statistic[i]
    
            # Find bondtype name elem1-elem2 or 'unknown'
            if i == 'unknown': bondtype = i
            else: bondtype = '{}-{}'.format(i[0], i[1])
            
            # Find cut-off used
            cutoff = 0.3 # Default
            if i in sorted_bo: cutoff = sorted_bo[i]

            # Add to stats if bondtype exists in system
            if len(lst)/2 > 0:
                s = Stats()
                s.count = int(len(lst)/2)
                s.avg = '{:.4f}'.format( compute_mean(lst) )
                s.min = '{:.4f}'.format( min(lst) )
                s.max = '{:.4f}'.format( max(lst) )
                s.std = '{:.4f}'.format( compute_standard_deviation(lst) )
                s.cutoff = cutoff
                self.statistics[bondtype] = s
                
        # Find all abos per element type and update self.abo_stats
        cutoffs_abo = {} # {element symbol : cut-off }
        element_count = {i:0 for i in elements}
        elements_abo = {i:[] for i in elements}
        for ID in atoms:
            element, maxnb = bond_info[BO.atoms[ID].type]
            element_count[element] += 1
        for ID1, ID2 in self.bonds:
            element1, maxnb1 = bond_info[BO.atoms[ID1].type]
            element2, maxnb2 = bond_info[BO.atoms[ID2].type]
            elements_abo[element1].extend(abos[ID1])
            elements_abo[element2].extend(abos[ID2])
            cutoffs_abo[element1] = maxnb1
            cutoffs_abo[element2] = maxnb2
        for element in elements_abo:
            lst = elements_abo[element]
            s = Stats()
            s.count = element_count[element]
            s.avg = '{:.4f}'.format( compute_mean(lst) )
            s.min = '{:.4f}'.format( min(lst) )
            s.max = '{:.4f}'.format( max(lst) )
            s.std = '{:.4f}'.format( compute_standard_deviation(lst) )
            s.cutoff = cutoffs_abo[element]
            self.abo_stats[element] = s