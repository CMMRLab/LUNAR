# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 16th, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""


# Function to find neighbor_ids and update mm.atoms[ID] instance
def neighbors(mm, Nth_neigh_depth):
    
    # generate graph of 1st neighbors
    mm.graph = {i:[] for i in mm.atoms} # { atomid : [lst-of-1st-neighs]}
    
    # add 1st-neighs to graph
    for i in mm.bonds:
        id1, id2 = mm.bonds[i].atomids
        mm.graph[id1].append(id2)
        mm.graph[id2].append(id1)
    
    # Function to find neighs away
    def find_Nneighs_away(atomid, N, graph):
        
        # Initalize as many neighbors sub-lists as desired
        # [ [1st-neigh-lst], [2nd-neigh-lst], [N-neigh-lst], ...]
        neighbors = [[] for i in range(N)]; 
        
        # Initialize 1st neighs
        neighbors[0] = graph[atomid]
    
        # Initilaize visted with source node and 1st neighs
        visited = graph[atomid] + [atomid]
        
        # Loop through and find as many neighbors as desired
        for n, i in enumerate(neighbors):
            for j in neighbors[n]:
                if n+1 < N:
                    new = []
                    for k in graph[j]:
                        if k not in visited:
                            visited.append(k)
                            new.append(k)
                    neighbors[n+1].extend(new)
                    
        return {n+1:i for n, i in enumerate(neighbors)} # regenerate as dict with key as depth
    
    # Loop through mm.atoms and find neighbor dict:
    # {N-neigh-depth: [lst-of-neighs]}
    # {1st-neighs,   2nd-neighs,  Nth-neighs,}
    # {1: [1, 2, 3], 2: [4, 5, 6], ....}
    for i in mm.atoms:
        atom = mm.atoms[i]
        atom.neighbor_ids = find_Nneighs_away(i, Nth_neigh_depth, mm.graph)
        
    return mm
