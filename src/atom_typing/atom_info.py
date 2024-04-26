# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 16th, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Funtion to add needed infor for final atom-typing
def add(mm):
    
    # Function to build info lst = [element, ringsize, nb]
    def build_neighinfo_lst(mm, atomid):
        dictionary = {} # {neigh-depth : [sorted lst of neigh info]}
        for depth in mm.atoms[atomid].neighbor_ids:
            tmp_lst = []
            for i in mm.atoms[atomid].neighbor_ids[depth]:
                tmp_lst.append(mm.atoms[i].info)
                
            # Sort tmp_lst by nb [2], then ringsize [1] and then by element [0]
            tmp_lst = sorted(tmp_lst, key=lambda x: x[2])  
            tmp_lst = sorted(tmp_lst, key=lambda x: x[1])  
            tmp_lst = sorted(tmp_lst, key=lambda x: x[0])  
            
            # Add to dictionary
            dictionary[depth] = tmp_lst
        return dictionary
    
    
    # Loop through mm.atoms and add new attributes of .nb and .info
    for i in mm.atoms:
        atom = mm.atoms[i]
        neighIDs = atom.neighbor_ids
        
        # Add new instances of nb and info
        atom.nb = len(neighIDs[1])
        atom.info = [atom.element, atom.ring, atom.nb]
        
    # Loop through mm.atoms and add neighbor_info attributes
    for i in mm.atoms:
        atom = mm.atoms[i]
        atom.neighbor_info = build_neighinfo_lst(mm, i)

    return mm