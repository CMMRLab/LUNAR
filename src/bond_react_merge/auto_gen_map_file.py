# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
April 14, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

##############################
# Import Necessary Libraries #
##############################
import src.bond_react_merge.graph_theory as gt
from collections import Counter
import os


################################################################
# Class to find all info needed to write a bond/react map-file #
################################################################
class find:
    def __init__(self, pre, post, filename, pairname, dataN_types2nb_map, log):
        self.map = {} # { pre-atomid : postatomid }
        self.pre_edge = {} # { pre-atomid's determined as edge atoms : comment }
        self.post_edge = {} # { post-atomid's determined as edge atoms : comment }
        self.cost = {} # { tuple(pre-atomid, postatomid): COST DICT}
        self.InitiatorIDs = {'ID1':'FAILED', 'ID2':'FAILED'} # Initialize and update later if two molecules exist in pre rxn-template
        self.deleteIDs = {} # { pre-rxn atom-id : comment }
        self.iterations = 0 # convergence for equivalences
        self.max_cost = 0 # max accetable cost iterated to
        pre_filename = os.path.basename(pre.filename)
        post_filename = os.path.basename(post.filename)
        
        
        ##################################################
        # Get pre headers info: BondingIDs and CreateIDs #
        ##################################################
        BondingIDs, CreateIDs, Reduce, Remove, Keep = get_lmp_header_info(pre, log)
        # Updating BondingIDs from Reduce
        if len(Reduce) == 5 and len(BondingIDs) == 0:
            log.out('Updating BondingIDs from Reduce for Equivalence matching')
            for n, i in enumerate(Reduce):
                if n <= 3: BondingIDs.append(i)
        self.BondingIDs = BondingIDs; self.CreateIDs = CreateIDs; self.Reduce = Reduce;
        
        
        ##############################################################################
        # Check for errors before proceding. If an error is found the code will exit #
        ##############################################################################
        error_checks(pre, post, CreateIDs, pre_filename, post_filename, log)
        
        
        ##########################################################################################################################################################################
        # find pre/post atoms=[lst of atomids]; types={atomid:atom-type}; elements={atomid:element}; neighs={atomid: [lst of neighs]}; neigh_types={atomid: [lst of atomtypes]}; #
        # neigh_elements={atomid: [lst of elements]}; 2nd_types={1: {2: ['cp', 'oc']}}; 2nd_elementss={1: {2: ['C', 'O']}};                                                      #
        ##########################################################################################################################################################################
        pre_atoms, self.pre_types, self.pre_elements, pre_graph, pre_neigh_types, pre_neigh_elements = connections(pre)
        post_atoms, self.post_types, self.post_elements, post_graph, post_neigh_types, post_neigh_elements = connections(post)
        if len(BondingIDs) == 4: # If remove the "bond" from the post_graph if using BondingIDs option
            try:
                post_graph[BondingIDs[2]].remove(BondingIDs[3])
                post_graph[BondingIDs[3]].remove(BondingIDs[2])
            except: log.error(f'ERROR BondingIDs or ReduceIDs post-bond {BondingIDs[2]}-{BondingIDs[3]} does not exist. Update BondingIDs or ReduceIDs')
        
        
        ######################################
        # find pre/post molecules and molids #
        ######################################
        pre_molecules, self.pre_molids, pre_formulas, pre_cluster_info = clusters(pre, log, pflag=True, rxntype='Pre', advancedformula=False)
        post_molecules, self.post_molids, post_formulas, post_cluster_info = clusters(post, log, pflag=True, rxntype='Post', advancedformula=False)
        
        
        ##################################################################################################################################
        # Find edge atoms for pre/post (This will be based on nb for atom-types in read dataN tags vs preN/postN tagged files. If no     #
        # info is in dataN_types2nb_map it will be based on if the element is not in elements2skip and if it only has 1 bonded neighbor) #
        ##################################################################################################################################
        elements2skip = ['Br', 'Ca', 'Cl', 'F', 'H'] # elements known to have one bonded atom that are usualy terminal
        # Find pre-rxn edge atoms
        for i in pre_atoms:
            nb = len(pre_graph[i]); element = self.pre_elements[i]; atom_type = self.pre_types[i];
            if atom_type in dataN_types2nb_map:
                
                # If nb not in nb_dataN it must be an edge atom or the user messed something up either way, look it as an edge atom.
                if nb not in dataN_types2nb_map[atom_type] and nb >= 1:
                    if nb == 1:
                        comment = '{} (   nb==1 and nb does not match most common nb for the atom-type from the dataN tagged files     ) '.format(atom_type)
                        self.pre_edge[i] = comment
                    if nb > 1: # Warn user if possible edgeID has more then 1 bonded neighbor, since its likely an issue on their end
                        log.warn(f'    WARNING atomID {i} was found to be a possible edge atom since the number of bonded neighbors differ from')
                        log.out(f'    other {atom_type} types, but has more then 1-bonded neighbor. Rejecting atomID {i} from edge atom consideration.\n')
                
            # elif only try guessing based on known terminating elements if len(dataN_types2nb_map) == 0
            elif element not in elements2skip and nb == 1 and len(dataN_types2nb_map) == 0:
                comment = '{} (   nb==1 and not a known terminating element of: {}   ) '.format(atom_type, ', '.join(elements2skip))
                self.pre_edge[i] = comment
        
        # Find post-rxn edge atoms
        for i in post_atoms:
            nb = len(post_graph[i]); element = self.post_elements[i]; atom_type = self.post_types[i]
            if atom_type in dataN_types2nb_map:
                
                # If nb not in nb_dataN it must be an edge atom or the user messed something up either way, look it as an edge atom.
                if nb not in dataN_types2nb_map[atom_type] and nb == 1:
                    comment = '{} (   nb>=1 and nb does not match most common nb for the atom-type from the dataN tagged files     ) '.format(atom_type)
                    self.post_edge[i] = comment
                
            # elif only try guessing based on known terminating elements if len(dataN_types2nb_map) == 0
            elif element not in elements2skip and nb == 1 and len(dataN_types2nb_map) == 0:
                comment = '{} (   nb==1 and not a known terminating element of: {}   ) '.format(atom_type, ', '.join(elements2skip))
                self.post_edge[i] = comment

        
        ############################################################################
        # Loop through pre_atoms and check against post_atoms and assign costs ... #
        ############################################################################
        cost_tally = {} # { tuple(pre-atomid : postatomid): COST TALLY}
        cost_dict =  {} # { tuple(pre-atomid : postatomid): COST DICT}
        pretypes = [self.pre_types[i] for i in self.pre_types]
        for pre_atom in pre_atoms:
            for post_atom in post_atoms:
                if post_atom in CreateIDs: continue
                # Find pair, pre/post comment
                pair = tuple([pre_atom, post_atom])
                pre_comment = pre.atoms[pre_atom].comment
                post_comment = post.atoms[post_atom].comment
                
                # Find pre/post element
                pre_element = pre_comment[pre_comment.rfind('/')+1:].capitalize()
                post_element = post_comment[post_comment.rfind('/')+1:].capitalize()
                
                # Find pre/post atom-type
                pre_atomtype = pre_comment[:pre_comment.rfind('/')]
                post_atomtype = post_comment[:post_comment.rfind('/')]
                
                # Find pre/post molID to assign any possible deleteID costs
                pre_molid = self.pre_molids[pre_atom]
                post_molid = self.post_molids[post_atom]
                
                # Find cummulative neighs
                pre_cumulative_neighs = bfs_neighs_depth(pre_atom, pre_graph)
                post_cumulative_neighs = bfs_neighs_depth(post_atom, post_graph)
                
                # Only procede if pre/post elements are the same (enforce chemistry to reduce errors and minimize equivalence search space)
                if pre_element == post_element and pair not in cost_tally:
                
                    # Intialize tmpcost dict with zeros and update later on
                    tmpcost = {'atom-type':0, 'nb':0, 'neigh_1st_atom-types':0, 'neigh_1st_elements':0,
                               'edge':0, 'new-type':0, 'depth_from_edge':0, 'deleteID':0, 'ID':0, 'molID':0,
                               'cumulative_neigh_types':0, 'cumulative_neigh_elements':0, 'cumulative_neigh_nb':0}
                    
                    #########################
                    # Start assigning costs #
                    #########################                    
                    # type_cost = difference value based on atom-types between pre atomid and post atomid (Look for difference in
                    # atomtypes. This will come with a moderate penality cost of 1*natoms since this will help match NON-REACTIVE atoms)
                    if pre_atomtype == post_atomtype: type_cost = 0
                    else: type_cost = 10
                    
                    # nb_cost = number of bonded costs and is difference of bonded atoms (Look for same relative connectivity, assuming
                    # connectivity should change minimally)
                    if post_atom not in CreateIDs:
                        nb_cost = int(abs(len(pre_graph[pre_atom]) - len(post_graph[post_atom])) )
                    else: nb_cost = 0
                    
                    # neigh_1st_elements_cost = is how many 1st neigh atom types match between pre to post
                    neigh_1st_elements_cost = diff_lst(pre_neigh_elements[pre_atom], post_neigh_elements[post_atom])
                     
                    # neigh_1st_atomtypes_cost = is how many 1st neigh atom types match between pre to post
                    neigh_1st_atomtypes_cost = diff_lst(pre_neigh_types[pre_atom], post_neigh_types[post_atom])
                    
                    # cummulative neigh costs: atom-types, elements, and nbs
                    max_neigh_depth = min( [len(pre_cumulative_neighs), len(post_cumulative_neighs)] ) - 1
                    cumulative_neigh_types_cost = 0; cumulative_neigh_elements_cost = 0; cumulative_neigh_nb_cost = 0;
                    for n in range(max_neigh_depth):
                        # Difference in cummulative atom types
                        cumulative_pre_types = [self.pre_types[i] for i in pre_cumulative_neighs[n]]
                        cumulative_post_types = [self.post_types[i] for i in post_cumulative_neighs[n] if i not in CreateIDs]
                        cumulative_neigh_types_cost += diff_lst(cumulative_pre_types, cumulative_post_types)

                        # Difference in cummulative elements
                        cumulative_pre_elements = [self.pre_elements[i] for i in pre_cumulative_neighs[n]]
                        cumulative_post_elements = [self.post_elements[i] for i in post_cumulative_neighs[n] if i not in CreateIDs]
                        cumulative_neigh_elements_cost += diff_lst(cumulative_pre_elements, cumulative_post_elements)
                        
                        # Difference in cummulative nbs
                        cumulative_pre_nb = [len(pre_graph[i]) for i in pre_cumulative_neighs[n]]
                        cumulative_post_nb = [len(post_graph[i]) for i in post_cumulative_neighs[n] if i not in CreateIDs]
                        cumulative_neigh_types_cost += diff_lst(cumulative_pre_nb, cumulative_post_nb)
                    if len(post_cumulative_neighs) != len(pre_cumulative_neighs):
                        ndiff = 10*abs(len(pre_cumulative_neighs) - len(post_cumulative_neighs))
                        cumulative_neigh_types_cost += ndiff
                        cumulative_neigh_elements_cost += ndiff
                        cumulative_neigh_nb_cost += ndiff
                    
                    # deleteID cost if post_molid > 1 assign a large cost to help deleteIDs be assigned lass
                    if post_molid > 1: deleteID_cost = 20
                    else: deleteID_cost = 0
                    
                    # molID cost to assign all atoms from smaller molecules last
                    if pre_molid > 1: molID_cost = 1
                    else: molID_cost = 0
                    
                    # new_type_cost = cost if post_atomtype is not in pretypes lst (this should be a pretty heavy cost
                    # to ensure that any new atom types get assigned last since this will be the atomtypes that dont
                    # have a clear match and making the equivalence assignment last on these types will reduce the number
                    # of possibilities to increase the probability of getting it correct).
                    if post_atomtype in pretypes: new_type_cost = 0
                    else: new_type_cost = 10
                    
                    # Depth from edge cost which assigns a penality based on the difference of each atoms (pre/post)
                    # shortest path depth(s) to all edges atoms found prior to assigning cost fitting
                    pre_depth_from_edge = []; post_depth_from_edge = [];
                    for edge_atom in self.pre_edge:
                        shortest_path = gt.find_shortest_path(pre_graph, pre_atom, edge_atom)
                        if shortest_path: pre_depth_from_edge.append(len(shortest_path))
                    for edge_atom in self.post_edge:
                        shortest_path = gt.find_shortest_path(post_graph, post_atom, edge_atom)
                        if shortest_path: post_depth_from_edge.append(len(shortest_path))
                    pre_depth_from_edge = sorted(pre_depth_from_edge); post_depth_from_edge = sorted(post_depth_from_edge);
                    depth_from_edge = 2*diff_lst(pre_depth_from_edge, post_depth_from_edge) # Weight it by 2 since this is important
                    
                    # edge_cost = a penality if both atoms are not edge atoms
                    if pre_atom in self.pre_edge and post_atom not in self.post_edge: 
                        if nb_cost <= 1 and type_cost == 0: edge_cost = 10
                        else: edge_cost = 10
                    elif pre_atom not in self.pre_edge and post_atom in self.post_edge:
                        if nb_cost <= 1 and type_cost == 0: edge_cost = 10
                        else: edge_cost = 10
                    else: edge_cost = 0
                    
                    # ID_cost is a penality if atomID in pre/post differ. This should not be an extreme cost
                    # since it may cause issues for some templates. The purpose of this cost is to keep atomIDs
                    # as close as possible to one another, especially in very symmertric templates.
                    if pre_atom == post_atom: ID_cost = 0
                    else: ID_cost = round(abs(pre_atom-post_atom)/2, 2) 
                    
                    # Update tmpcost with new found costs
                    tmpcost['neigh_1st_atom-types'] = neigh_1st_atomtypes_cost
                    tmpcost['neigh_1st_elements'] = neigh_1st_elements_cost
                    tmpcost['cumulative_neigh_types'] = cumulative_neigh_types_cost
                    tmpcost['cumulative_neigh_elements'] = cumulative_neigh_elements_cost
                    tmpcost['cumulative_neigh_nb'] = cumulative_neigh_nb_cost
                    tmpcost['atom-type'] = type_cost
                    tmpcost['new-type'] = new_type_cost
                    tmpcost['edge'] = edge_cost
                    tmpcost['nb'] = nb_cost
                    tmpcost['depth_from_edge'] = depth_from_edge
                    tmpcost['deleteID'] = deleteID_cost
                    tmpcost['ID'] = ID_cost
                    tmpcost['molID'] = molID_cost

                    # If pair not already in cost_tally add to cost_tally and cost_dict
                    # Also make sure edge_cost is zero so edges are mapped properly
                    if pair not in cost_tally and edge_cost == 0:
                        cost_tally[pair] = sum(list(tmpcost.values()))
                        cost_dict[pair] = tmpcost

                        
        ##################################################################################################
        # Find equivalences from costs, initialize self.map and self.cost with zeros and update later on #
        ##################################################################################################
        # Intialize self.map and self.cost with zeroes
        for i in pre_atoms:
            self.map[i] = 0; 
            self.cost[(i, 0)] = {'atom-type':0, 'nb': 0, 'neigh_1st_atom-types':0, 'neigh_1st_elements':0,
                                 'edge':0, 'new-type':0, 'depth_from_edge':0, 'deleteID':0, 'ID':0, 'molID':0,
                                 'cumulative_neigh_types':0, 'cumulative_neigh_elements':0, 'cumulative_neigh_nb':0}
                        
        # Dynamically check all possible equivalences based on cost approach with two Methods: 1) if user supplies BondingIDs 2) general cost minimizer
        used_postids = set(CreateIDs); used_preids = set(); paired = set([]); # Data-structs to keep track of what has been used
        while_count = 0; max_while = 250*len(pre_atoms);  # while count and max loop break condition
        acceptable_cost = 0; cost_increment = 1; # acceptable cost (tally); increment 
        update = 2; # update condtion to increment acceptable_cost by (>=2 seems to work best, gives N-chances to match per acceptable_cost)
        #--------------------------------------------------------------------------------------------------------------------#
        # Method 1 (if user provideds BondingIDs get intial maps from BFS radial searching and let Method 2 handle the rest) #
        #--------------------------------------------------------------------------------------------------------------------#
        if len(BondingIDs) == 4:
            # Get pre-rxn/post mol1 and mol2 and set BondingIDs equivs and log associated info
            log.out(f'  BondingIDs {str(BondingIDs)} will be used as much as possible in a BFS mapping pattern')
            pre_mol1_BFS = bfs_neighs_depth(BondingIDs[0], pre_graph); pre_mol2_BFS = bfs_neighs_depth(BondingIDs[1], pre_graph);
            post_mol1_BFS = bfs_neighs_depth(BondingIDs[2], post_graph); post_mol2_BFS = bfs_neighs_depth(BondingIDs[3], post_graph);
            self.map[BondingIDs[0]] = BondingIDs[2]; self.cost[(BondingIDs[0], BondingIDs[2])] = {'BondingIDs-BFS-pattern-Depth-0-mol1':0}
            self.map[BondingIDs[1]] = BondingIDs[3]; self.cost[(BondingIDs[1], BondingIDs[3])] = {'BondingIDs-BFS-pattern-Depth-0-mol2':0}
            used_preids.add(BondingIDs[0]); used_preids.add(BondingIDs[1]); used_postids.add(BondingIDs[2]); used_postids.add(BondingIDs[3]);
            keys = ['atom-type', 'nb', 'neigh_1st_elements', 'neigh_1st_atom-types', 'cumulative_neigh_types'] 
            
            # Map pre_mol1 to post_mol1
            min_depth = min([len(pre_mol1_BFS), len(post_mol1_BFS)]); tmp_count = 0; tmp_cost = 0; 
            for depth in range(min_depth):
                pre_neighs = pre_mol1_BFS[depth]; post_neighs = post_mol1_BFS[depth]; max_map = min([len(pre_neighs), len(post_neighs)]); tmp_used = [];
                while len(tmp_used) < max_map:
                    tmp_count += 1;
                    if tmp_count % update: tmp_cost += cost_increment
                    for i in pre_neighs:
                        tmptally, tmpdicts = reduce_and_sort_costdict(i, cost_dict, keys, paired, post_graph)
                        for pair in tmptally:
                            if pair[1] not in post_neighs: continue
                            preid, postid = pair; cost = tmptally[pair];
                            if preid not in used_preids and postid not in used_postids and cost <= acceptable_cost:
                                used_preids.add(preid);  used_postids.add(postid); tmp_used.append(preid); paired.add(pair); self.map[preid] = postid;
                                self.cost[pair] = tmpdicts[pair]; self.cost[pair] = {'BondingIDs-BFS-pattern-Depth-'+str(depth+1)+'-mol1':0}; break;
                    if tmp_count >= max_while/10: break
                
            # Map pre_mol2 to post_mol2
            min_depth = min([len(pre_mol2_BFS), len(post_mol2_BFS)]); tmp_count = 0; tmp_cost = 0; 
            for depth in range(min_depth):
                pre_neighs = pre_mol2_BFS[depth]; post_neighs = post_mol2_BFS[depth]; max_map = min([len(pre_neighs), len(post_neighs)]); tmp_used = [];
                while len(tmp_used) < max_map:
                    tmp_count += 1;
                    if tmp_count % update: tmp_cost += cost_increment
                    for i in pre_neighs:
                        tmptally, tmpdicts = reduce_and_sort_costdict(i, cost_dict, keys, paired, post_graph)
                        for pair in tmptally:
                            if pair[1] not in post_neighs: continue
                            preid, postid = pair; cost = tmptally[pair];
                            if preid not in used_preids and postid not in used_postids and cost <= acceptable_cost:
                                used_preids.add(preid);  used_postids.add(postid); tmp_used.append(preid); paired.add(pair); self.map[preid] = postid;
                                self.cost[pair] = tmpdicts[pair]; self.cost[pair] = {'BondingIDs-BFS-pattern-Depth-'+str(depth+1)+'-mol2':0}; break;
                    if tmp_count >= max_while/10: break
        
        #--------------------------------------------------------#
        # Method 2 (Default if user does not provide BondingIDs) #
        #--------------------------------------------------------#
        #sorted_cost_tally = dict(sorted(cost_tally.items(), key=lambda x:x[1], reverse=False))  # [0=keys;1=values] sort by cost
        #pre_atoms = [i[0] for i in sorted_cost_tally] # iterate through atomIDs with increasing costs
        #pre_atoms = sorted(pre_atoms, reverse=True)
        # sorted_cost_tally = dict(sorted(cost_tally.items(), key=lambda x:x[1], reverse=False))  # [0=keys;1=values] sort by cost
        # for i in sorted_cost_tally:
        #     keys = ['atom-type', 'nb', 'neigh_1st_elements', 'edge', 'depth_from_edge', 'deleteID', 'cumulative_neigh_types',
        #             'cumulative_neigh_elements', 'molID', 'neigh_1st_atom-types']
        #     tmptally, tmpdicts = reduce_and_sort_costdict(i[0], cost_dict, keys, paired, post_graph)
        #     if i[0] == 1:
        #         print()
        #         print(i, tmptally[i])
        #         print(tmpdicts[i])
        while len(used_preids) < len(pre_atoms):
            while_count += 1
            if while_count % update:
                acceptable_cost += cost_increment
            for i in pre_atoms:
                keys = ['atom-type', 'nb', 'neigh_1st_elements', 'edge', 'depth_from_edge', 'deleteID', 'cumulative_neigh_types',
                        'cumulative_neigh_elements', 'molID', 'neigh_1st_atom-types']                 
                if len(set(pre_formulas)) == 1: keys.append('ID') # If pre-rxn template has the same formula, assume the symmetry and assign 'ID' to constrain the domain
                tmptally, tmpdicts = reduce_and_sort_costdict(i, cost_dict, keys, paired, post_graph)
                for pair in tmptally:
                    preid, postid = pair; cost = tmptally[pair];
                    if preid not in used_preids and postid not in used_postids and cost <= acceptable_cost:
                        used_preids.add(preid); used_postids.add(postid); paired.add(pair);
                        self.map[preid] = postid; self.cost[pair] = tmpdicts[pair]; break
            if while_count >= max_while: 
                log.warn(' WARNING Equivalence map timed out. NOT ALL Equivalences could be found for {} -> {}'.format(pairname[0], pairname[1]));  break

        # Update convergence criteria: self.iterations with while_count and self.max cost with most current acceptable_cost
        self.iterations = while_count; self.max_cost = acceptable_cost;
        

        ###################################################################################################################
        # Find initiator ids from pre_molecules. First try for a two molecule search else warn about manual update status #
        ###################################################################################################################
        # If two molecules are found in pre rxn molecules assign intiatorIDs else print and use defaults from above
        if len(pre_molecules) == 2:
            #--------------------------------------------------------------------------------------------------------------------------------------------#
            # Method 1 for setting InitiatorIDs, by searching 4-6 neighs deep from any edge atom and using more unqiue type (Default, but kind of basic) #
            #--------------------------------------------------------------------------------------------------------------------------------------------#
            # Check that edgeIDs were found, if some edgeIDs were found, try method 1
            if self.pre_edge:
                self.InitiatorIDs = {} # Reset InitiatorIDs 
                for molid, mol in enumerate(pre_molecules):
                    # Find all bonds in molecule, generate molecule graph, find edge atoms, and find N-neighs
                    # intiator ids. Lastly, find the furtherst atomid from and of the edge atoms and used as an IntiatorID
                    mol_bonds = [pre.bonds[i].atomids for i in pre.bonds if check_bonds(pre.bonds[i].atomids, mol)];
                    mol_graph = gen_graph(list(mol), mol_bonds); mol_edges = [i for i in self.pre_edge if i in mol];
                    mol_neigh = {i:{n:j for n, j in enumerate(bfs_neighs_depth(i, mol_graph)[:6], 1)} for i in mol if i in mol_edges} # Find neighs 6 deep from edge atoms
                    
                    # Find all neighbors from all edges that are at least 4-neighs deep from any edge
                    neighs_3_deep_edge = []
                    for i in mol_neigh:
                        depth = mol_neigh[i]
                        for j in depth:
                            if j >= 3:
                                neighs_3_deep_edge.extend(depth[j])
                    neighs_3_deep_edge = sorted(set(neighs_3_deep_edge)) # Remove duplicates and sort
                    if neighs_3_deep_edge:
                        molecule_types = [self.pre_types[i] for i in neighs_3_deep_edge if i not in self.pre_edge]
                        if not molecule_types: molecule_types = [self.pre_types[i] for i in neighs_3_deep_edge]
                        reverse_map = {self.pre_types[i]:i for i in neighs_3_deep_edge}
                        types_sorted = [item for items, c in Counter(molecule_types).most_common() for item in [items] * c]

                        most_unqiue_type = types_sorted[-1]
                        most_unqiue_id = reverse_map[most_unqiue_type]
                        comment = '{:^5}    pre-molid: {:^5}    Count of type on pre-molid (that are at least 4-deep from any edge atom): {:^5}'.format(most_unqiue_type, molid+1,  molecule_types.count(most_unqiue_type))
                        self.InitiatorIDs[most_unqiue_id] = comment
                        
                #---------------------------------------------------------------------------------------------------------------------------------------------#
                # Method 2 for setting InitiatorIDs, by comparing bonds to find which post bond join two pre molecules (will only override Method 1 if found) #
                #---------------------------------------------------------------------------------------------------------------------------------------------#
                # Check if any map is still zero, if none are zero, try method2 
                matches = [self.map[i] for i in self.map]
                if 0 not in matches:
                    # Get pre bonds and post bonds (with post bonds mapped onto pre bonds)
                    pre_bonds = [tuple(sorted(pre.bonds[i].atomids)) for i in pre.bonds]
                    reverse_map = {self.map[i]:i for i in self.map}; post_bonds = [];
                    for i in post.bonds:
                        id1, id2 = post.bonds[i].atomids
                        post_bonds.append( tuple(sorted( [reverse_map[id1], reverse_map[id2]] )) )
                    
                    # Find differences in bonds and loop through to find if any link the two molecules
                    possible_InitiatorIDs = []; diff_bonds = set(post_bonds) - set(pre_bonds);
                    for id1, id2 in diff_bonds:
                        if id1 in pre_molecules[0] and id2 in pre_molecules[1]: possible_InitiatorIDs.append( sorted([id1, id2]) )
                        if id1 in pre_molecules[1] and id2 in pre_molecules[0]: possible_InitiatorIDs.append( sorted([id1, id2]) )
                    
                    # If len(possible_InitiatorIDs) >= 1 Reset InitiatorIDs
                    if len(possible_InitiatorIDs) >= 1:
                        self.InitiatorIDs = {} # Reset InitiatorIDs
                        index2use = 0 # default will be to use index zero
                        for n, IDs in enumerate(possible_InitiatorIDs):
                            Elements = [self.pre_elements[i] for i in IDs]
                            if 'H' not in Elements: index2use = n
                        
                        # Use index2use to decide which InitiatorID set to use
                        IDs = sorted(possible_InitiatorIDs[index2use])
                        otherIDs = [i for n, i in enumerate(possible_InitiatorIDs) if n != index2use]
                        if otherIDs: additional_comment = ' other possible linking atomIDs {}'.format(str(otherIDs))
                        else: additional_comment = ''
                        self.InitiatorIDs[IDs[0]] =  'links pre-rxn molecule molid 1 -> 2' + additional_comment
                        self.InitiatorIDs[IDs[1]] =  'links pre-rxn molecule molid 1 -> 2' + additional_comment
            
            # Reset InitiatorIDs based on supplied BondingIDs (if present and not from Reduce)
            if len(BondingIDs) == 4 and len(Reduce) in [0, 2]:
                self.InitiatorIDs = {} # Reset InitiatorIDs 
                self.InitiatorIDs[BondingIDs[0]] =  'BondingID set in pre-rxn LAMMPS datafile HEADER for mol1'
                self.InitiatorIDs[BondingIDs[1]] =  'BondingID set in pre-rxn LAMMPS datafile HEADER for mol2'
                
          
        #---------------------------------------------------------------------------------------------------#
        # else warn that pre template does not have two molecules and user must update IntiatarIDs manually #
        #---------------------------------------------------------------------------------------------------#
        else: log.warn(f'   WARNING {pre_filename} does not contain two molecules. {filename} map file will need manual updating of InitiatorIDs')


        #################################################################################################################
        # Find deleteIDs if post-rxn has two molecules, using the smallest molecule and assume equivalences are correct #
        #################################################################################################################
        # Check that all equivs have been found (all post-rxn-ids are none zeros)
        post_rxn_ids = [self.map[i] for i in self.map]
        if post_rxn_ids.count(0) == 0: # If count of 0 in post_rxn_ids equals zero all equiv's were assigned and are assumed to be correct
            for i in self.map:
                postid = self.map[i]; molid = self.post_molids[postid];
                
                # IF molid is greater then 1 assume it is a by-product and add to deleteIDs dict (molids are sorted by number of atoms,
                # so molid > 1 means that molid 1 is the largest molecule fragment and is assumed to be the molecule fragment to keep)
                if molid > 1:
                    comment = 'Pre-rxn-atomID: {:^5} has a Post-rxn-atomID: {:^5} that belongs to a Post-molID: {:^5} and has been tagged for deletion'.format(i, postid, molid)
                    self.deleteIDs[i] = comment
                    
        # Warn user if the auto-detetection of deleteIDs found deleteIDs
        if self.deleteIDs:
            log.warn(f'    WARNING It was detected that {post_filename} post-rxn has more then one molecule. This engaged the')
            log.out(f'    auto-deleteIDs detection and the written {filename} contains deleteIDs for the following pre-rxn-ids:')
            for i in self.deleteIDs:
                log.out('    {}'.format(self.deleteIDs[i]))
                
                
###########################################################
# Function to "order" atoms about an atomID. For example: #
#   atomID = 3                                            #
#   atoms = [1,2,3,4,5,6]                                 #
# The returned list of atoms will be:                     #
#  [3, 2, 4, 1, 5, 6]                                     #
###########################################################
def order_atomIDs_around_atomID(atomID, atoms):
    atomID_diffs = {i:abs(i-atomID) for i in sorted(atoms)} # {atomID : differents between atomIDs}
    atomID_diffs = dict(sorted(atomID_diffs.items(), key=lambda x:x[1], reverse=False))  # [0=keys;1=values]
    ordered_atoms = list(atomID_diffs.keys())
    return ordered_atoms
                
                
################################################################################
# Function to reduce inner dict and sort based atomids in pair and lowest cost #
################################################################################
def reduce_and_sort_costdict(atomid, cost_dict, keys, paired, post_graph):            
    tmptally = {} # { tuple(pre-atomid : postatomid): COST TALLY}
    tmpdicts = {} # { tuple(pre-atomid : postatomid): COST DICT}
    for pair in cost_dict:
        if pair[0] == atomid and pair not in list(paired):
            tmp = {i:cost_dict[pair][i] for i in cost_dict[pair] if i in keys}
            tally = sum([cost_dict[pair][i] for i in cost_dict[pair] if i in keys])
            tmpdicts[pair] = tmp;  tmptally[pair] = tally;
            
    # Check for duplicate tally and if duplicates found add the atomID cost
    tallies = list(tmptally.values())
    duplicates = set([x for x in tallies if tallies.count(x) > 1])
    if duplicates:
        for pair in tmptally:
            if tmptally[pair] in duplicates:
                tmptally[pair] += cost_dict[pair]['ID']
                tmpdicts[pair]['ID'] = cost_dict[pair]['ID']
                
    # Sort tmptally then create tmpdicts from tmptally since tmptally has already been sorted
    # and tmpdicts will inherit tmptally sequence (ASSUMES PYTHON 3.7 OR GREATER)
    tmptally = dict(sorted(tmptally.items(), key=lambda x:x[0][0], reverse=False)) # [0=keys;1=values][0=preID] 
    tmptally = dict(sorted(tmptally.items(), key=lambda x:x[0][1], reverse=False)) # [0=keys;1=values][1=postID] 
    tmptally = dict(sorted(tmptally.items(), key=lambda x:x[1], reverse=False))    # [0=keys;1=values] sort by cost 
    tmpdicts = {j:tmpdicts[j] for j in tmptally} # ReBuild tmpdicts from order in tmpy tally
    return tmptally, tmpdicts
                

####################################################
# Function to get info from LAMMPS datafile header #
####################################################
def get_lmp_header_info(m, log):
    BondingIDs = []; CreateIDs = []; Reduce = []; Remove = []; Keep = []; filename = os.path.basename(m.filename);
    
    # Function to parse header and get list
    def get_list(line):
        lst = []; string = ''; bracket = False; equal = False;
        count1 = 0; count2 = 0;
        for i in line[-1]:
            if i == '=': equal = True; continue
            if i == '[': bracket = True; count1 += 1; continue;
            if i == ']': bracket = False; count2 += 1; continue;
            if count1 > 1 or count2 > 1: bracket = False
            if bracket and equal: string += i
        string = string.split(',') # split by commas
        for ID in string:
            try: lst.append(int(ID))
            except: pass
        return lst
    
    # Try getting BondingIDs
    if 'BondingIDs' in m.header:
        line = m.header.split('BondingIDs') # Find BondingIDs string 
        BondingIDs = get_list(line) # Get list of info
        if BondingIDs: log.out(f'{filename} had the following BondingIDs specified in the header: {str(BondingIDs)}')
        
    # Try getting reduce
    if 'Reduce' in m.header:
        line = m.header.split('Reduce') # Find Reduce string
        Reduce = get_list(line) # Get list of info
        if Reduce: log.out(f'{filename} had the following Reduce specified in the header: {str(Reduce)}')
        
    # Try getting Remove
    if 'Remove' in m.header:
        line = m.header.split('Remove') # Find Remove string
        Remove = get_list(line) # Get list of info
        if Remove: log.out(f'{filename} had the following Remove specified in the header: {str(Remove)}')
        
    # Try getting Keep
    if 'Keep' in m.header:
        line = m.header.split('Keep') # Find Keep string
        Keep = get_list(line) # Get list of info
        if Keep: log.out(f'{filename} had the following Keep specified in the header: {str(Keep)}')
    
    # Try getting CreateIDs 
    if 'CreateIDs' in m.header:
        line = m.header.split('CreateIDs') # Find CreateIDs string 
        CreateIDs = get_list(line) # Get list of info
        if CreateIDs: log.out(f'{filename} had the following CreateIDs specified in the header: {str(CreateIDs)}')
    return BondingIDs, CreateIDs, Reduce, Remove, Keep
    

########################################################################
# Function to check for errors before attempting to find map file data #
########################################################################
def error_checks(pre, post, CreateIDs, pre_filename, post_filename, log):    
    # Check that same number of atoms exist in each file pair
    if pre.natoms != post.natoms-len(CreateIDs):
        log.error(f'ERROR attempting to generate map files, but {pre_filename} and {post_filename} have different number of atoms')
                
    # Find elements in each file and exit for each file that does not contain atom-type/element comment style
    pre_elements = []; post_elements = []; unique_elements = set();
    for i in pre.atoms:
        comment = pre.atoms[i].comment
        split = comment.split('/')
        if '/' not in comment or not split[-1].isalpha():
            log.error(f'ERROR generate_map_file requires atom-type/element comment style. {pre_filename} has this style: {comment} for atomID {i}')
        elif split[-1].islower():
            log.out(f'ERROR generate_map_file requires atom-type/element comment style. {pre_filename} has this style: {comment} for atomID {i}');
            log.out('The element is assigned as an atom-type (all2lmp default if it can not determine the element type. Indicates incorrect atom typing).')
            log.error('Please check for missing coeffs and if none are missing assign an element sybmol in place of the atom-type.')
        else: pre_elements.append(comment[comment.rfind('/')+1:].capitalize()); unique_elements.add(comment[comment.rfind('/')+1:].capitalize());
    for i in post.atoms:
        comment = post.atoms[i].comment
        split = comment.split('/')
        if i in CreateIDs: continue
        if '/' not in comment or not split[-1].isalpha():
            log.error(f'ERROR generate_map_file requires atom-type/element comment style. {post_filename} has this style: {comment} for atomID {i}')
        elif split[-1].islower():
            log.out(f'ERROR generate_map_file requires atom-type/element comment style. {post_filename} has this style: {comment} for atomID {i}');
            log.out('The element is assigned as an atom-type (all2lmp default if it can not determine the element type. Indicates incorrect atom typing).')
            log.error('Please check for missing coeffs and if none are missing assign an element sybmol in place of the atom-type.')
        else: post_elements.append(comment[comment.rfind('/')+1:].capitalize()); unique_elements.add(comment[comment.rfind('/')+1:].capitalize());
        
    # sort unique_elements and loop through to make sure count is consistant between files
    unique_elements = sorted(unique_elements)
    for i in unique_elements:
        pre_elem_count = pre_elements.count(i)
        post_elem_count = post_elements.count(i)
        if pre_elem_count != post_elem_count:
            log.out(f'ERROR generate_map_file requires same number of elements in pre/post files. {pre_filename} and {post_filename} have different count of {i} elements');
            log.out(f'  {pre_filename} has {set(pre_elements)} elements')
            log.error(f'  {post_filename} has {set(post_elements)} elements')
    return
            
            
###############################
# Graph theory type functions #
###############################
# Function to check if bonds are in molecule atoms
def check_bonds(bond, mol):
    id1, id2 = bond; return_boolean = False
    if id1 in mol and id2 in mol:
        return_boolean = True
    return return_boolean

# Function to generate graph
def gen_graph(mol, bonds):
    # Intialize and add to graph
    graph = {i:[] for i in mol}
    for id1, id2 in bonds:
        graph[id1].append(id2)
        graph[id2].append(id1)
    return graph

# Function to find all neighs
def bfs_neighs_depth(atomID, graph):
    # Use BFS and iterate through all atoms to find neighbors list of lists
    neighbors = [graph[atomID]]; visited = set(graph[atomID] + [atomID]);
    for n, neighs in enumerate(neighbors):
        for neigh in neighs:
            neighbors.append([]); tmp = [];
            for adjacent in graph[neigh]:
                if adjacent in visited: continue
                tmp.append(adjacent)
                visited.add(adjacent)
            neighbors[n+1].extend(tmp)
            
    # Reduce neighbors to only neighs without empty lists.
    reduced = [neighs for neighs in neighbors if neighs != []]
    return reduced

            
#####################################
# Function for finding bonded atoms #
#####################################
def connections(data):    
    # Intialize atoms dictionary
    graph = {}; atoms = []; elements = {}; atomtypes = {};
    for i in data.atoms:
        # Find atomtype and element
        comment = data.atoms[i].comment
        atomtype = comment.split('/')[0]
        element = comment[comment.rfind('/')+1:].capitalize()
        
        # Create data-structs
        atoms.append(i); graph[i] = []
        atomtypes[i] = atomtype
        elements[i] = element
        
    # add in connect atoms to bonded dictionary
    for i in data.bonds:
        id1, id2 = data.bonds[i].atomids        
        graph[id1].append(id2)
        graph[id2].append(id1)
    
    # Find neigh_types and neigh_elements
    neigh_types = {}; neigh_elements = {}; # {atomid : atomtype or atom element}
    for i in graph:
        tmp_types = []; tmp_elements = [];
        for j in graph[i]:
            tmp_types.append(atomtypes[j])
            tmp_elements.append(elements[j])
        neigh_types[i] = sorted(tmp_types); neigh_elements[i] = sorted(tmp_elements);
    return sorted(atoms), atomtypes, elements, graph, neigh_types, neigh_elements


############################################
# Function for finding molecules in system #
############################################
def clusters(m, log, pflag, rxntype='none', advancedformula=False):
    # Function to find cluster formula
    def find_cluster_formula(cluster, graph, m, advancedformula):
        elements = [m.atoms[i].comment[m.atoms[i].comment.rfind('/')+1:].capitalize() for i in cluster]
        base_elements = list(sorted({i for i in elements})); formula = '';
        for element in base_elements:
            formula += '{}{}-'.format(element, elements.count(element))
        formula = formula[:-1]
        
            
        # Advanced nb ending
        if advancedformula:
            nb = {i:{} for i in base_elements} # { element:{nb:count} }
            formula += '/'
            for i in cluster:
                element = m.atoms[i].comment[m.atoms[i].comment.rfind('/')+1:].capitalize()
                bonded = len(graph[i])
                if bonded in nb[element]: nb[element][bonded] += 1
                else: nb[element][bonded] = 0
            for element in nb:
                formula += '{'
                for bonded in nb[element]:
                    formula += '{}:{}-'.format(bonded, nb[element][bonded])
                formula = formula[:-1]
                formula += '}'
        return formula

    # Generate graph
    graph = {i:[] for i in m.atoms}
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        graph[id1].append(id2)
        graph[id2].append(id1)
        
    # Initialize clusters and molids
    checked = {ID:False for ID in m.atoms}; clusters = set([]); molids = {};
    for ID in graph:
        molids[ID] = 1
        if checked[ID]: continue
        visited=set([ID]); queue=[ID];
        while queue:
            s = queue.pop(0) 
            for neighbor in graph[s]:
                if checked[neighbor]: continue
                visited.add(neighbor)
                queue.append(neighbor)
                checked[neighbor]=True
        clusters.add( tuple(sorted(visited)) )

    # Sort clusters in a very unique fashion
    clusters = list(clusters)
    clusters = sorted(clusters, key=lambda x: x[0]) # Sort all clusters based on 1st atomID in cluster
    clusters = sorted(clusters, key=len, reverse=True) # Sort all clusters by number of atoms
    
    # Find molids after sorting clusters. molid will be count + 1 where the largest
    # cluster will be molid = 1 and smallet cluster will be molid = nclusters in system
    natoms = []; formulas = []; info = {} # {molID : info object}
    class Info: pass # .size .mass .psize .pmass .atoms
    for n, cluster in enumerate(clusters):
        # Update molIDs
        for atomID in cluster:
            molids[atomID] = n+1
    
        # if advancedformula and number of clusters is less then 50
        # add molID to string as well to help near_edge_charges.py
        formula = find_cluster_formula(cluster, graph, m, advancedformula)
        if advancedformula and len(clusters) <= 50:
            formula += '/molID:{}'.format(n+1)
        
        # Analyze cluster
        size = len(cluster)
        natoms.append(size)
        formulas.append(formula)
        

        
        # Log cluster info
        I = Info()
        I.atoms = clusters[n]
        I.formula = formula
        I.size = size
        info[n+1] = I
        
    # Print if pflag
    if pflag:
        log.out('\n  ----------------------{:^4} rxn Cluster Analysis-----------------------'.format(rxntype))
        log.out('  {:^10} {:^15} {:^20} {:^20}'.format('molID', 'Molecule Size', '%Size of Molecule', 'formula'))
        log.out('  ----------------------------------------------------------------------'); size_total = sum(natoms);                              
        for n in range(len(natoms)):
            natom = '{: >6}'.format(natoms[n]); psize = '{:.2f}'.format(round(100*natoms[n]/size_total, 2))  
            log.out('  {:^10} {:^15} {:^20} {:^20}'.format(n+1, natom, psize, formulas[n]))   
        log.out('')                                                                       
    return clusters, molids, formulas, info


#################################################################
# Function to count differing elements from lst of pre/post     #
# will be used for counting different 1st neigh types/elements  #
#################################################################
def diff_lst(pre_lst, post_lst):
    count = 0; pre_lst = sorted(pre_lst); post_lst = sorted(post_lst)
    for n, pre in enumerate(pre_lst):
        try: 
            post = post_lst[n]
            if pre != post:
                count += 1
        except: count += 1
    return int(count)


######################################################
# Function to find nb map for every atom-type in the # 
# corresponding dataN tagged files (will be used to  #
# help make sure edge atoms are correct)             #
######################################################
def dataN_type2nb_map(datafile, dataN_types2nb_lst):    
    atoms, atomtypes, elements, bonded, neigh_types, neigh_elements = connections(datafile)
    for ID in atomtypes:
        atomtype = atomtypes[ID]; nb = len(bonded[ID]);
        if atomtype in dataN_types2nb_lst:
            dataN_types2nb_lst[atomtype].append(nb)
        else: dataN_types2nb_lst[atomtype] = [nb]
    return dataN_types2nb_lst


##################################
# Function to write new map file #
##################################
def write_map(filename, pre2post, title, version, comment_flag):
    
    # function to print_dict
    def string_dict(dictionary):
        string = ''
        for key in dictionary:
            string += '{:<2} {}: {:^5} {:>2}'.format('', key, dictionary[key], '')
        return string
    
    # Open and write file
    with open(filename,'w') as f:
        f.write(f'{title[0:220]}\n\n') # Make max header length of 220 characters 

        # Write Nequivalences and NedgeIDS
        if pre2post.pre_edge: f.write(f'{len(pre2post.pre_edge)} edgeIDs\n')
        f.write(f'{len(pre2post.map)} equivalences\n')
        
        # Write NdeleteIDs if they exists
        if pre2post.deleteIDs: f.write(f'{len(pre2post.deleteIDs)} deleteIDs\n')
        if pre2post.CreateIDs: f.write(f'{len(pre2post.CreateIDs)} createIDs\n')
        f.write('\n')
        
        # Write InitiatorIDs
        f.write('InitiatorIDs\n\n')
        for i in sorted(pre2post.InitiatorIDs):
            comment = '{:^5} {:^5}'.format('#', pre2post.InitiatorIDs[i])
            
            # Write to file
            if comment_flag:
                f.write('{:>2} {:^5}\n'.format(i, comment))
            else: f.write('{:>2}\n'.format(i))
        
        # Write EdgeIDs
        if pre2post.pre_edge:
            f.write('\nEdgeIDs\n\n')
            for i in sorted(pre2post.pre_edge):
                comment = '{:^5} {:^5}'.format('#', pre2post.pre_edge[i])
                
                # Write to file
                if comment_flag:
                    f.write('{:>2} {:^5}\n'.format(i, comment))
                else: f.write('{:>2}\n'.format(i))
            
        # Write CreateIDs
        if pre2post.CreateIDs:
            f.write('\nCreateIDs\n\n')
            for i in sorted(pre2post.CreateIDs):
                comment = '{:^5} {:^5} {}'.format('#', 'user input ', pre2post.post_types[i])
                f.write('{:>2} {:^5}\n'.format(i, comment))
                
        # Write deleteIDs if they exists
        if pre2post.deleteIDs:
            f.write('\nDeleteIDs\n\n')
            for i in pre2post.deleteIDs:
                comment = '{:^5} {:^5}'.format('#', pre2post.deleteIDs[i])
                
                # Write to file
                if comment_flag:
                    f.write('{:>2} {:^5}\n'.format(i, comment))
                else: f.write('{:>2}\n'.format(i))
            
        # Write Equivalences
        if comment_flag: f.write('\nEquivalences # ( pre -> post )\n\n')
        else: f.write('\nEquivalences\n\n')
        for i in pre2post.map:
            j = pre2post.map[i]
            
            # Build comment for data logging
            costdict = pre2post.cost[(i, j)];
            totalcost = sum(list(costdict.values()));
            breakdowncost = string_dict(costdict);
            
            # Find pre/post types and molids
            pretype = pre2post.pre_types[i];
            premolid = pre2post.pre_molids[i]
            try: posttype = pre2post.post_types[j];
            except: posttype = 'N/A';
            try: postmolid = pre2post.post_molids[j]
            except: postmolid = 0
            types = '{:^5} {:^5} -> {:^5}'.format('types: ', pretype, posttype);
            molid = '{:^5} {:^5} -> {:^5}'.format('molids: ', premolid, postmolid)
            cost_total = '{:^5} {:>5.3f}'.format('Total cost: ', totalcost)
            cost_breakdown = '( {:^5} )'.format(breakdowncost)
            comment = '# {:<20} {:>15} {:>20} {:>100}'.format(types, molid, cost_total, cost_breakdown)
            
            # Write to file
            if comment_flag:
                f.write('{:>2} {:^15} {}\n'.format(i, j, comment))
            else: f.write('{:>2} {:^15}\n'.format(i, j))      
    return