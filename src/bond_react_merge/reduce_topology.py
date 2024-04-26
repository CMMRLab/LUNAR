# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
July 31st, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

##############################
# Import Necessary Libraries #
##############################
import src.bond_react_merge.auto_gen_map_file as agmf
import src.atom_removal.rm_atoms as rm_atoms


# Function to reduce template down
def template(pre, post, BondingIDs, CreateIDs, Reduce, Remove, Keep, log):

    # Generate graphs and find IDs to search outward from
    pre_atoms, pre_types, pre_elements, pre_graph, pre_neigh_types, pre_neigh_elements = agmf.connections(pre)
    post_atoms, post_types, post_elements, post_graph, post_neigh_types, post_neigh_elements = agmf.connections(post)
    
    # Find reduce method based on length of lists
    if len(Reduce) == 5:
        preID1, preID2, postID1, postID2, depth1 = Reduce; 
        depth2 = depth1; BFS_reduce = True;
    elif len(Reduce) == 2 and len(BondingIDs) == 4:
        preID1, preID2, postID1, postID2 = BondingIDs
        depth1, depth2 = Reduce; BFS_reduce = True;
    else: BFS_reduce = False;
    
    # Find pre/post Keep and Remove atoms
    pre_keep = []; post_keep = [];
    if Keep:
        if len(Keep) % 2 == 0:
            middle = int(len(Keep)/2)
            pre_keep = Keep[:middle]
            post_keep = Keep[middle:]
        else: log.out('Keep option used but odd number of atomIDs provided')
    pre_remove = []; post_remove = []; remove_flag = False
    if Remove:
        if len(Remove) % 2 == 0:
            middle = int(len(Remove)/2)
            pre_remove = Remove[:middle]
            post_remove = Remove[middle:]
            remove_flag = True
        else: log.out('Remove option used but odd number of atomIDs provided')


    # If BFS_reduce finding pre/post atoms to keep and delete
    pre_reduced = ''; post_reduced = ''; atoms_deleted = False;
    if BFS_reduce:
        atoms_deleted = True
        # Find preIDs to keep/remove and regenerage m class as *_reduced
        preIDs2keep1 = find_Nneighs_away(preID1, pre_graph, depth1)
        preIDs2keep2 = find_Nneighs_away(preID2, pre_graph, depth2)
        preIDs2keep = sorted(set(preIDs2keep1 + preIDs2keep2 + [preID1] + [preID2] + pre_keep))
        preIDs2remove = [i for i in pre.atoms if i not in preIDs2keep] + pre_remove
        log.out('  pre-rxn atom removal:')
        pre_reduced = rm_atoms.constructor(pre, preIDs2remove, log, method='atomIDs', pflag=True)
        
        # Find postIDs to keep/remove and regenerage m class as *_reduced
        if len(Reduce) == 2 and len(BondingIDs) == 4:
            try:
                post_graph[postID1].remove(postID2)
                post_graph[postID2].remove(postID1)
            except:
                log.error(f'ERROR BondingIDs or ReduceIDs post-bond {postID1}-{postID2} does not exist. Update BondingIDs or ReduceIDs')
        postIDs2keep1 = find_Nneighs_away(postID1, post_graph, depth1)
        postIDs2keep2 = find_Nneighs_away(postID2, post_graph, depth2)
        postIDs2keep = sorted(set(postIDs2keep1 + postIDs2keep2 + [postID1] + [postID2] + post_keep + CreateIDs))
        postIDs2remove = [i for i in pre.atoms if i not in postIDs2keep] + post_remove
        log.out('  post-rxn atom removal:')
        post_reduced = rm_atoms.constructor(post, postIDs2remove, log, pflag=True)
        
    # elif remove_flag is independant of Reduce
    elif remove_flag:
        atoms_deleted = True
        pre_reduced = rm_atoms.constructor(pre, pre_remove, log, pflag=True)
        post_reduced = rm_atoms.constructor(post, post_remove, log, pflag=True)
    
    # Update header info after removing atoms (if applicable)
    if atoms_deleted:
        new_BondingIDs = []; new_CreateIDs = []; new_Reduce = []; strings = [];
        if BondingIDs:
            for n, i in enumerate(BondingIDs):
                if n <= 1: new_BondingIDs.append(pre_reduced.atomID_map[i])
                else: new_BondingIDs.append(post_reduced.atomID_map[i])
            strings.append('BondingIDs = {};'.format(str(new_BondingIDs)))
        if len(Reduce) == 5:
            for n, i in enumerate(Reduce):
                if n <= 1: new_Reduce.append(pre_reduced.atomID_map[i])
                elif n <= 3: new_Reduce.append(post_reduced.atomID_map[i])
                else: new_Reduce.append(i)
            strings.append('Reduce = {};'.format(str(new_Reduce)))
        if len(Reduce) == 2: strings.append('Reduce = {};'.format(str(Reduce)))
        if CreateIDs: 
            new_CreateIDs = [post_reduced.atomID_map[i] for i in CreateIDs]
            strings.append('CreateIDs = {};'.format(str(new_CreateIDs)))
        if strings:
            strings.insert(0, 'bond_react_merge.py reduce or remove option used. Updating Header with new atomIDs >')
            pre_reduced.header = ' '.join(strings)

        # Print warning if for some reason pre and post have differnt number of atoms in the template        
        if pre_reduced.natoms != post_reduced.natoms:
            log.warn('WARNING template reduce failed due to different number of pre-atoms to post atoms. Try adjusting depth variable.')
    return pre_reduced, post_reduced


# Function to find neighs away
def find_Nneighs_away(atomid, graph, N):
    neighbors = [[] for i in range(N)]; 
    neighbors[0] = graph[atomid]
    visited = graph[atomid] + [atomid]
    for n, i in enumerate(neighbors):
        for j in neighbors[n]:
            if n+1 < N:
                new = []
                for k in graph[j]:
                    if k not in visited:
                        visited.append(k)
                        new.append(k)
                neighbors[n+1].extend(new)
    return [ID for depth in neighbors for ID in depth]


# Function to find rxn pairs
def find_rxn_pairs(merge):
    template_pairs = {} # {1:[pre1, post1], 2:[pre2, post2], ... N:[] }
    for file in merge:    
        # strip last number of filename. Assumes tags either:
        tmpname = ''.join([i for i in file if i.isalpha()])
        tmpid = ''.join([i for i in file if i.isdigit()])
        
        # Log pairs if tmpid exists, else create key/valued list if tmpname is not data
        if tmpname != 'data':
            if tmpid in template_pairs:
                template_pairs[tmpid].append(file)
            else: template_pairs[tmpid] = [file]
            
    # Loop through template_pairs and try finding map
    return template_pairs