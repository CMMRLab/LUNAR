# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
August 6th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

##############################
# Import Necessary Libraries #
##############################
import src.bond_react_merge.auto_gen_map_file as agmf
import os


################################################################
# Class to find all info needed to write a bond/react map-file #
################################################################
class find:
    def __init__(self, pre, post, pre2post, fileN_tagged, log, ndepth_update=0):
        #################################
        # Find pre molecules and molids #
        #################################
        pre_molecules, pre_molids, pre_formulas, pre_cluster_info = agmf.clusters(pre, log, pflag=False, rxntype='Pre')
        pre_atoms, pre_types, pre_elements, pre_graph_full, pre_neigh_types, pre_neigh_elements = agmf.connections(pre)
        
        
        #################################################################################################
        # Try topolgoy matching between different files and different molecular formulas from each file #
        #################################################################################################
        for premolID in pre_cluster_info: 
            # Get basic cluster info
            pre_cluster = pre_cluster_info[premolID]; pre_atomIDs = pre_cluster.atoms;
            pre_graph = {i:pre_graph_full[i] for i in pre_atomIDs}; found_edgeIDs = [];
            pre_edgeIDs = [i for i in pre2post.pre_edge if i in list(pre_atomIDs)]; 
            match_info = {} # { tuple(filename, formula) : {'matched':[matched edgeIDs], 'map':{}, 'class':m, 'found': T or F} }
            for edgeID in pre_edgeIDs:
                pre_edgeType = pre_types[edgeID]; pre_BFS_atoms = agmf.bfs_neighs_depth(edgeID, pre_graph)
                pre_BFS_types = pathIDs2atomtypes(pre_BFS_atoms, pre_types);
                for file in fileN_tagged:
                    m = fileN_tagged[file]; filename = os.path.basename(m.filename); formula2atomIDs = get_formula2molatomIDs(m, log);
                    for formula in formula2atomIDs:
                        match_key = tuple([filename, formula])
                        if match_key not in match_info: match_info[match_key] = { 'matched':[], 'map':{}, 'class':m, 'found': False, 'match_atoms': []}
                        atomIDs2match = formula2atomIDs[formula][0] # use the 1st tuple of atomIDs from each formula to try for match
                        match_atoms, match_types, match_elements, match_graph, match_neigh_types, match_neigh_elements = agmf.connections(m)
                        for matchID in atomIDs2match:
                            if match_types[matchID] == pre_edgeType and matchID not in match_info[match_key]['matched']: # 1st check if type is the same
                                match_BFS_atoms = agmf.bfs_neighs_depth(matchID, match_graph)
                                match_BFS_types = pathIDs2atomtypes(match_BFS_atoms, match_types)#[:len(pre_BFS_types)]
                                match = check_preBFS2matchBFS_typelst(pre_BFS_types, match_BFS_types, pre_BFS_atoms, match_BFS_atoms, pre_types, match_types, pre_graph, match_graph, ringsym=False)
                                if match:
                                    match_info[match_key]['matched'].append(matchID); match_info[match_key]['found'] = True;
                                    match_info[match_key]['class'] = m; match_info[match_key]['map'][edgeID] = matchID;
                                    match_info[match_key]['match_atoms'] = atomIDs2match; found_edgeIDs.append(edgeID); break;
                    
                    
            ################################################
            # Loop through match_info to find equivalences #
            ################################################
            count = 0;
            for pair in match_info:
                info = match_info[pair]
                if info['found'] and count == 0:
                    # Update any edgeID that equivs were found (control of depth was control at that stage, so update all that are in equivs)
                    equivs, depths = get_equivs(pre, info['class'], info['map'], info['match_atoms'], ndepth_update); count += 1;
                    file_molecule = '{}/{}'.format(pair[0], pair[1])
                    log.out('    {} {} {} {}'.format('Some rxn edgeIDs in molID: ', premolID, ' matched with file ', file_molecule))
                    for i in equivs:
                        l = equivs[i]; depth = depths[i];
                        mapped_charge = info['class'].atoms[l].charge
                        mapped_comment = info['class'].atoms[l].comment
                        
                        # Update pre-rxn template
                        atom = pre.atoms[i]
                        oldcharge = atom.charge
                        atom.charge = mapped_charge # Update this instance of charge as well
                        atom.mapped_charge = mapped_charge
                        atom.mapped_comment = 'Charge  mapping: {:>10.6f}  -> {:>10.6f} mapped from atomID ({}): {:^3} in {} file at a Depth: {:^3} from edgeIDs.'.format(
                                               oldcharge, mapped_charge, mapped_comment, l, file_molecule, depth)
                        
                        # Assume pre-2-post rxn equivalences are correct and update
                        j = pre2post.map[i]
                        atom = post.atoms[j]
                        oldcharge = atom.charge
                        atom.charge = mapped_charge # Update this instance of charge as well
                        atom.mapped_charge = mapped_charge
                        atom.mapped_comment = 'Charge  mapping: {:>10.6f}  -> {:>10.6f} mapped from atomID ({}): {:^3} in {} file at a Depth: {:^3} from edgeIDs.'.format(
                                               oldcharge, mapped_charge, mapped_comment, l, file_molecule, depth)
                if count > 0: break;
                    
            #############################################
            # Update any edgeID that could not be found #
            #############################################
            for i in pre_edgeIDs:
                if i not in found_edgeIDs:
                    log.warn(f'    WARNING edgeID : {i} equivalent not found and charge not updated')

                    # Update pre-rxn template
                    atom = pre.atoms[i]
                    atom.mapped_charge = atom.charge
                    atom.mapped_comment = 'Charge  mapping: {:>10.6f}  -> {:>10.6f} FAILED, CHARGE KEPT THE SAME.'.format(atom.charge, atom.charge)
                    
                    # Assume pre-2-post rxn equivalences are correct and update
                    j = pre2post.map[i]
                    atom = post.atoms[j]
                    atom.mapped_charge = atom.charge
                    atom.mapped_comment = 'Charge  mapping: {:>10.6f}  -> {:>10.6f} FAILED, CHARGE KEPT THE SAME.'.format(atom.charge, atom.charge)
            

######################################################
# Function to help find information and equivalences #
######################################################
# Function to get equivs from pre-rxn to matched topology
def get_equivs(pre, matched, mapped_edges, match_atoms, ndepth_update):
    # Getpre and matched info
    pre_atoms, pre_types, pre_elements, pre_graph, pre_neigh_types, pre_neigh_elements = agmf.connections(pre)
    match_atoms, match_types, match_elements, match_graph, match_neigh_types, match_neigh_elements = agmf.connections(matched)
    
    # Initialze equivs at edge atoms
    equivs = {} # { pre-atomID : matched-atomID equivalet } 
    depths = {} # { pre-atomID : depth from edge atom }
    for i in mapped_edges:
        j = mapped_edges[i]; equivs[i] = j; depths[i] = 0;
    
    # If ndepth_update > 0 create reduced graph and start matching at each depth from an edge atom
    if ndepth_update > 0:
        bonds2remove = []; # Find bonds2remove to split molecule at edge atoms
        for i in mapped_edges:
            j = mapped_edges[i]; intial_1st_neigh = '';
            pre_atom_first_neigh = pre_graph[i][0] # Assume only 1-bonded neigh on edge atom from pre-rxn
            pre_type_first_neigh = pre_types[pre_atom_first_neigh] # Get atom type
            pre_neighneigh_types = sorted([pre_types[s] for s in pre_graph[pre_atom_first_neigh]])
            for k in match_graph[j]:
                match_neighneigh_types = sorted([match_types[s] for s in match_graph[k]])
                if match_types[k] == pre_type_first_neigh: intial_1st_neigh = k
                if pre_neighneigh_types == match_neighneigh_types: intial_1st_neigh = k; break;
            if intial_1st_neigh != '': # Find bonds2remove if 1st neighs of match do not equal pre-1st-neighs
                for l in match_graph[j]:
                    if l != intial_1st_neigh: 
                        bonds2remove.append( sorted([j, l]) )
    
        # Find reduce data structs
        reducedbonds = []; reducegraph = {i:[] for i in match_atoms};
        for i in matched.bonds:
            bond = sorted(matched.bonds[i].atomids)
            if bond not in bonds2remove and bond[0] in match_atoms and bond[1] in match_atoms:
                reducedbonds.append(bond)
        for id1, id2 in reducedbonds:
            reducegraph[id1].append(id2); reducegraph[id2].append(id1);
            
        # Iterate through each mapped edges to find equivs emaniting from that location
        pre_edges = list(mapped_edges.keys()); matched_edges = list(mapped_edges.values());
        for edgeID in mapped_edges:
            pre_BFS_atoms = agmf.bfs_neighs_depth(edgeID, pre_graph)[:ndepth_update]
            match_BFS_atoms = agmf.bfs_neighs_depth(mapped_edges[edgeID], reducegraph)[:ndepth_update]
            
            # Start matching at each depth
            for n in range(min([len(pre_BFS_atoms), len(match_BFS_atoms)])):
                pre_depth_atoms = [i for i in pre_BFS_atoms[n] if i not in pre_edges]
                match_depth_atoms = [i for i in match_BFS_atoms[n] if i not in matched_edges]
                
                # Match atoms at each depth
                used_preIDs = []; used_matchIDs = []; nequivs = min([len(pre_depth_atoms), len(match_depth_atoms)]);
                for preN in range(nequivs):
                    preID = pre_depth_atoms[preN]
                    pre_type = pre_types[preID]; pre_neightypes = sorted([pre_types[i] for i in pre_graph[preID]])
                    for matchN in range(nequivs):
                        matchID = match_depth_atoms[matchN]
                        match_type = match_types[matchID]; match_neightypes = sorted([match_types[i] for i in match_graph[matchID]])
                        if pre_type == match_type and pre_neightypes == match_neightypes and preID not in used_preIDs and matchID not in used_matchIDs:
                            used_preIDs.append(preID); used_matchIDs.append(matchID); equivs[preID] = matchID; depths[preID] = n+1; break;
    return equivs, depths

# Function to check preBFS to matchBFS lst
def check_preBFS2matchBFS_typelst(pre_BFS_types, match_BFS_types, pre_BFS_atoms, match_BFS_atoms, pre_types, match_types, pre_graph, match_graph, ringsym=False):
    try:
        match = False; found = {i:False for i in range(len(pre_BFS_types))}
        for n, depth_types in enumerate(pre_BFS_types): 
            pre_counts = {i:depth_types.count(i) for i in set(depth_types)} # { atomType : count of atomType }
            for pre_type in pre_counts:
                pre_count = pre_counts[pre_type]
                match_count = match_BFS_types[n].count(pre_type)
                
                # If count is accetable look deeper into if it is a topological match
                if pre_count <= match_count:
                    pre_atoms = pre_BFS_atoms[n]; match_atoms = match_BFS_atoms[n];
                    for l in pre_atoms:
                        pre_neigh_types = sorted([pre_types[s] for s in pre_graph[l]])
                        for t in match_atoms:
                            match_neigh_types = sorted([match_types[s] for s in match_graph[t]])
                            if pre_neigh_types == match_neigh_types: found[n] = True
        
        # Check found for only True values in the dictionary
        found_values = set(list(found.values()))
        if len(found_values) == 1 and True in found_values:  match = True
    except: match = False
    return match 

# Function to get path atom types
def pathIDs2atomtypes(path, types):
    path_types = []
    for depth in path:
        path_types.append(sorted([types[i] for i in depth]))
    return path_types    
                         
# Function to get general formula-to-atomIDs in molecule (to only search a certain molecular formula once for speed up)
def get_formula2molatomIDs(m, log):
    molecules, molids, formulas, cluster_info = agmf.clusters(m, log, pflag=False, rxntype='Pre', advancedformula=True)
    unique_formulas = set(formulas); formula2atomIDs = {i:[] for i in unique_formulas} # { formula : [lst of atomIDs] }
    for molID in cluster_info: 
        # Get basic cluster info
        cluster = cluster_info[molID]
        atomIDs = cluster.atoms
        formula = cluster.formula    
        formula2atomIDs[formula].append(atomIDs)
    return formula2atomIDs