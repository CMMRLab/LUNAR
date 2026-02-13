# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 2nd, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

##############################
# Import Necessary Libraries #
##############################
import src.bond_react_merge.auto_gen_map_file as agmf
import src.bond_react_merge.graph_theory as gt
import src.io_functions as io_functions
import os

# Function to add info the pre/post molecule files
def add2molecules(molecule_file_options, pre, post, pre2post, newfile, log):
    
    # Loop through options and apply what ever options are desired
    for option in molecule_file_options:
        option = option.split('/') # split options
        
        # if option[0] == 'fragment' apply fragment options
        if option[0] == 'fragment':
            fragmentID = option[1] # set fragmentID
            
            # fragment 'custom_charges' option
            if option[2] == 'custom_charges':
                # if depth_from_edge_ in option[3] apply depth from edge option
                if 'depth_from_edge_' in option[3]:
                    tmp = option[3].split('_');  depth = int(tmp[-1]);
                    pre = depth_from_edge(option, fragmentID, depth, pre, pre2post, newfile, log, rxntype='pre')
                    post = depth_from_edge(option, fragmentID, depth, post, pre2post, newfile, log, rxntype='post')
                    
                # elif equiv_cost_0 in option[3] apply depth from edge option
                elif 'equiv_cost_' in option[3]:
                    tmp = option[3].split('_');  cost = int(tmp[-1]);
                    pre, post = equiv_cost(option, fragmentID, depth, pre, post, pre2post, cost, newfile, log)

        # elif option[0] == 'molecule' apply molecule options
        elif option[0] == 'molecule':
            pre =  molecule(pre, option, newfile, log)
            post =  molecule(post, option, newfile, log)
    
    return pre, post


################################################
# Function to perfrom depth_from_edge_N option #
################################################
def depth_from_edge(option, fragmentID, depth, rxn, pre2post, newfile, log, rxntype):
    # Initialize fragment_atoms
    fragment_atoms = []
    
    #################################
    # find pre molecules and molids #
    #################################
    rxn_molecules, rxn_molids, rxn_formulas, rxn_info = agmf.clusters(rxn, log, pflag=False, rxntype='Pre', advancedformula=False)
    
    #############################################################################
    # Find pre molecule graphs from pre_molecules and start matching topologies #
    #############################################################################
    for molid, mol in enumerate(rxn_molecules):            
        # Find molecule bonds and graph
        mol_bonds = [rxn.bonds[i].atomids for i in rxn.bonds if agmf.check_bonds(rxn.bonds[i].atomids, list(mol))]
        mol_graph = agmf.gen_graph(list(mol), mol_bonds)
        
        # Find molecule edge ids for this molecule and create dict of atoms in mol to keep track of shortest distances to edge atoms
        if rxntype == 'pre': edge = pre2post.pre_edge
        elif rxntype == 'post': edge = pre2post.post_edge
        edgeIDs = [i for i in edge if i in list(mol)]
        mol_atoms = {i:[] for i in mol}
        
        # Loop through mol_atoms and edgeIDs to find shortest distances of atom from edgeIDs and add to valued-lst
        for i in mol_atoms:
            for j in edgeIDs:
                if i != j:
                    shortest_path = gt.find_shortest_path(mol_graph, i, j, path=[])
                    mol_atoms[i].append(len(shortest_path))
                else:
                    mol_atoms[i].append(0)
                
        # Loop through mol_atoms and see if any path is to short to be append atom to fragment_atoms
        for i in mol_atoms:
            # Check if all paths are greater then the cut-off depth
            path_lengths = mol_atoms[i]
            if path_lengths:
                if min(path_lengths) > depth+1:
                    fragment_atoms.append(i)
                
    # Check if rxn has dict fragments and add new info to rxn accordingly
    class Fragments:
        pass # .comment .atoms
        
    if hasattr(rxn,'fragments'):
        f = Fragments()
        f.comment = '/'.join(option)
        f.atoms = sorted(fragment_atoms)
        rxn.fragments[fragmentID] = f
    else:
        rxn.fragments = {} # {fragmentID : fragment class}
        f = Fragments()
        f.comment = '/'.join(option)
        f.atoms = sorted(fragment_atoms)
        rxn.fragments[fragmentID] = f
        
    # Print outcome
    basename = io_functions.get_basename(rxn.filename, newfile=newfile, character=':', pflag=True)
    newname = '{}.lmpmol'.format(basename)
    log.out(f'  Adding molecule FragmentID {fragmentID} based on {"/".join(option)} to')
    log.out(f'  {newname} with a depth-cutoff: {depth}')
    log.out('')
    return rxn


################################################
# Function to perfrom depth_from_edge_N option #
################################################
def equiv_cost(option, fragmentID, depth, pre, post, pre2post, cost, newfile, log):
    # Initialize pre/post fragment_atoms
    pre_fragment_atoms = []
    post_fragment_atoms = []
    
    # Loop through map and find equiv and cost
    for i in pre2post.map:
        j = pre2post.map[i]

        # Find total cost
        costdict = pre2post.cost[(i, j)];
        totalcost = sum(list(costdict.values()));
        
        # If totalcost > cost log accordingly
        if totalcost > cost:
            pre_fragment_atoms.append(i)
            post_fragment_atoms.append(j)
            
    # Check if rxn has dict fragments and add new info to rxn accordingly
    class Fragments:
        pass # .comment .atoms
    
    # Add to pre class
    if hasattr(pre,'fragments'):
        f = Fragments()
        f.comment = '/'.join(option)
        f.atoms = sorted(pre_fragment_atoms)
        pre.fragments[fragmentID] = f
    else:
        pre.fragments = {} # {fragmentID : fragment class}
        f = Fragments()
        f.comment = '/'.join(option)
        f.atoms = sorted(pre_fragment_atoms)
        pre.fragments[fragmentID] = f
        
    # Add to post class
    if hasattr(post,'fragments'):
        f = Fragments()
        f.comment = '/'.join(option)
        f.atoms = sorted(post_fragment_atoms)
        post.fragments[fragmentID] = f
    else:
        post.fragments = {} # {fragmentID : fragment class}
        f = Fragments()
        f.comment = '/'.join(option)
        f.atoms = sorted(post_fragment_atoms)
        post.fragments[fragmentID] = f
        
    # Print outcome
    pre_basename = io_functions.get_basename(pre.filename, newfile=newfile, character=':', pflag=True)
    pre_newname = '{}.lmpmol'.format(pre_basename)
    post_basename = io_functions.get_basename(post.filename, newfile=newfile, character=':', pflag=True)
    post_newname = '{}.lmpmol'.format(post_basename)
    log.out(f'  Adding molecule FragmentID {fragmentID} based on {"/".join(option)} to')
    log.out(f'  {pre_newname} and {post_newname} with a cost-cutoff: {cost}')
    log.out('')
    return pre, post    

#########################################
# Function to add molecule info to file #
#########################################
def molecule(rxn, option, newfile, log):
    
    #################################
    # find pre molecules and molids #
    #################################
    rxn_molecules, rxn_molids, rxn_formulas, rxn_info = agmf.clusters(rxn, log, pflag=False, rxntype='Pre', advancedformula=False)
    
    # Check if rxn has dict fragments and add new info to pre accordingly
    class Molecules:
        pass # .comment .type
    
    # Generate atoms dict
    atoms = {} # {atomID : Molecules Object}
    atomids = sorted(list(rxn_molids.keys()))
    for i in atomids:
        M = Molecules()
        comment = rxn.atoms[i].comment
        atomtype = comment.split('/')[0]
        M.type = atomtype
        M.comment = '/'.join(option)
        M.molid = rxn_molids[i]
        atoms[i] = M
        
    if hasattr(rxn,'molecules'):
        M = Molecules()
        M.comment = '/'.join(option)
        M.atoms = atoms
        rxn.molecules = M
    else:
        M = Molecules()
        M.comment = '/'.join(option)
        M.atoms = atoms
        rxn.molecules = M
        
    # Print outcome
    basename = io_functions.get_basename(rxn.filename, newfile=newfile, character=':', pflag=True)
    newname = '{}.lmpmol'.format(basename)
    log.out(f'  Adding molecules section based on {"/".join(option)} to {newname}')
    log.out('')
    return rxn