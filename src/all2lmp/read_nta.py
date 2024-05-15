# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
October 11th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
# Funtion to read nta file
def atoms(nta_file, m, log):
    
    ####################################
    # Opening and reading the frc file #
    ####################################
    with open(nta_file, 'r') as f:
        nta = {} # { atomid : new atom type } (For style id)
        name = {} # {atomid : nta:NAME}
        types = {} # { atomtype : new atom type } (For style type)
        edge = {} # { atomid : [lst of edge atoms] }
        charges = {} # { atomid : charge } (For charge OPT were OPT=id or type or nta)
        tmpq = {'id':{}, 'type':{}, 'nta':{}} # tmp charges to assign to atomid key in charges after read
        neutralize = {'all': False, 'bond-inc': False, 'user-defined': False, 'zero': False} # Intialize all as False and update if found in file
        type_equivs = {} # { set typed : equiv to map set type to }
        remove = {'angle-nta':[], 'dihedral-nta':[], 'improper-nta':[], 
                  'angle-ID':[], 'dihedral-ID':[], 'improper-ID':[], 
                  'zero':{'angle':False, 'dihedral':False, 'improper':False}}
        
        style = '' # set as empty string and update once found
        charge_style = '' # set as empty string and update once found
        
        # Initializing flags    
        atomtypes_flag = False
        edgeatoms_flag = False
        charges_flag = False
        neutralize_flag = False
        equivs_flag = False
        remove_flag = False
        
        # Looping through each line of the nta file
        linenumber = 0;
        for line in f:
            linenumber += 1
            if linenumber == 1: continue
        
            # Strip comment's
            line = line.split('#')[0]
            line = line.rstrip()
            
            # Split line
            line = line.strip()
            line_split = line.split()

            # Setting flags
            if line == '':# or not any(chr.isdigit() for chr in line):
                atomtypes_flag = False
                edgeatoms_flag = False
                charges_flag = False
                neutralize_flag = False
                equivs_flag = False
                remove_flag = False
            elif 'style' in line:
                style = line_split[1] # should be in second index location
                atomtypes_flag = True
                edgeatoms_flag = False
                charges_flag = False
                neutralize_flag = False
                equivs_flag = False
                remove_flag = False
                continue
            elif 'edge' in line and 'id' in line:
                edgeatoms_flag = True
                atomtypes_flag = False
                charges_flag = False
                neutralize_flag = False
                equivs_flag = False
                remove_flag = False
                continue
            elif 'charge' in line and 'neutralize' not in line:
                charges_flag = True
                atomtypes_flag = False
                edgeatoms_flag = False
                neutralize_flag = False
                equivs_flag = False
                remove_flag = False
                charge_style = line_split[1] # should be in second index location
                continue
            elif 'remove' in line:
                equivs_flag = False
                neutralize_flag = False
                charges_flag = False
                atomtypes_flag = False
                edgeatoms_flag = False
                remove_flag = True
                angle_id = False; dihedral_id = False; improper_id = False;
                angle_nta = False; dihedral_nta = False; improper_nta = False;
            if 'neutralize' in line:
                neutralize_flag = True
                charges_flag = False
                atomtypes_flag = False
                edgeatoms_flag = False
                equivs_flag = False
                remove_flag = False
            if 'equivs' in line:
                equivs_flag = True
                neutralize_flag = False
                charges_flag = False
                atomtypes_flag = False
                edgeatoms_flag = False
                remove_flag = False


            
            # Finding were atom types exist in the file
            if atomtypes_flag and len(line_split) >= 2:
                # Set nta based on style id
                if style == 'id':
                    try: atomid = int(line_split[0])
                    except: atomid = line_split[0]
                    atom_type = line_split[1]
                    if ':' in atom_type:
                        nta[atomid] = atom_type[:atom_type.rfind(':')] # strip any ':NAME ending'
                    else: nta[atomid] = atom_type
                    name[atomid] = atom_type
                
                # Set nta based on style type
                elif style == 'type':
                    # Only valid file for style type is a LAMMPS .data or .dat file
                    if m.filename.endswith('data') or m.filename.endswith('dat') or m.filename.endswith('pdb'):
                        try: typeid = int(line_split[0])
                        except: typeid = line_split[0]
                        atom_type = line_split[1]
                        types[typeid] = atom_type
                    else: log.error('{} {} {}'.format('ERROR style type used in', nta_file, ' and read in topofile is not a LAMMPS datafile or manually edit pdb file'))
                # Exit if style is not type or id
                else: log.error('{} {} {}'.format('ERROR requesting style of unsupported type in:', nta_file, ' only supported style are id or type'))
                    
                
            # Finding where edge id maybe in file
            elif edgeatoms_flag and len(line_split) >= 2:
                atomid = int(line_split[0])
                edgetypes = line_split[1:] # Take everything after 1st index
                
                # If edgetypes is not empty add to edge dict
                if edgetypes:
                    edge[atomid] = edgetypes
                    
            # Finding where charge maybe in file and add info to tmpq based on specific charge style
            elif charges_flag and len(line_split) >= 2:
                try: key = int(line_split[0]) # try getting integer for id and type style
                except: key = str(line_split[0]) # execept its fate as a string for nta style
                try: charge = float(line_split[1])
                except: charge = line_split[1]
                tmpq[charge_style].update({key:charge})
                #print(line_split)
                
            # Finding neutralize information
            elif neutralize_flag and len(line_split) >= 4:
                if line_split[3] == 'user-defined': neutralize['user-defined'] = True
                if line_split[3] == 'bond-inc': neutralize['bond-inc'] = True
                if line_split[3] == 'all': neutralize['all'] = True
                if line_split[3] == 'zero': neutralize['zero'] = True
                
            # Finding equivs
            elif equivs_flag and len(line_split) >= 2:
                type_equivs[line_split[0]] = line_split[1]
                
            # Finding remove flags
            elif remove_flag:
                # Flag information
                if 'angle-remove' in line:
                    if 'zero' in line:
                        remove['zero']['angle'] = True
                    if 'nta' in line:
                        angle_id = False; dihedral_id = False; improper_id = False;
                        angle_nta = True; dihedral_nta = False; improper_nta = False;
                        continue
                    if 'id' in line:
                        angle_id = True; dihedral_id = False; improper_id = False;
                        angle_nta = False; dihedral_nta = False; improper_nta = False;
                        continue
                elif 'dihedral-remove' in line:
                    if 'zero' in line:
                        remove['zero']['dihedral'] = True
                    if 'nta' in line:
                        angle_id = False; dihedral_id = False; improper_id = False;
                        angle_nta = False; dihedral_nta = True; improper_nta = False;
                        continue
                    if 'id' in line:
                        angle_id = False; dihedral_id = True; improper_id = False;
                        angle_nta = False; dihedral_nta = False; improper_nta = False;
                        continue
                elif 'improper-remove' in line:
                    if 'zero' in line:
                        remove['zero']['improper'] = True
                    if 'nta' in line:
                        angle_id = False; dihedral_id = False; improper_id = False;
                        angle_nta = False; dihedral_nta = False; improper_nta = True;
                        continue
                    if 'id' in line:
                        angle_id = False; dihedral_id = False; improper_id = True;
                        angle_nta = False; dihedral_nta = False; improper_nta = False;
                        continue
                    
                # Find information
                if angle_id and len(line_split) == 3:
                    IDs = []
                    for i in line_split:
                        try: ID = int(i)
                        except: ID = i
                        IDs.append(ID)
                    remove['angle-ID'].append(tuple(IDs))
                if angle_nta and len(line_split) == 3:
                    remove['angle-nta'].append(tuple(line_split))
                if dihedral_id and len(line_split) == 4:
                    IDs = []
                    for i in line_split:
                        try: ID = int(i)
                        except: ID = i
                        IDs.append(ID)
                    remove['dihedral-ID'].append(tuple(IDs))
                if dihedral_nta and len(line_split) == 4:
                    remove['dihedral-nta'].append(tuple(line_split))
                if improper_id and len(line_split) == 5:
                    IDs = []
                    for i in line_split:
                        try: ID = int(i)
                        except: ID = i
                        IDs.append(ID)
                    remove['improper-ID'].append(tuple(IDs))
                if improper_nta and len(line_split) == 5:
                    remove['improper-nta'].append(tuple(line_split))
                
        ###############################################################
        # If style was type build nta dict by looping through m.atoms #
        ###############################################################
        if style == 'type':
            # Check that number of atom types in datafile are consistent with number of atom types in nta_file
            natomtypes_datafile = m.natomtypes; natomtypes_nta_file = len(types);
            if natomtypes_datafile != natomtypes_nta_file:
                log.error('{} {} {} {}'.format('ERROR number of atom types in: ', m.filename, ' are inconsistent with number of atom types in: ', nta_file))
                
            # If check above passes build nta dict
            for atomid in m.atoms:
                # Find type based on type in m.atoms and type in types
                atom = m.atoms[atomid]
                atom_type = types[atom.type]
                if ':' in atom_type:
                    nta[atomid] = atom_type[:atom_type.rfind(':')] # strip any ':NAME ending'
                else: nta[atomid] = atom_type
                name[atomid] = atom_type
                
        #####################################
        # use type_equivs to map atom types #
        #####################################
        if type_equivs:
            equivs = ',  '.join(list(type_equivs.keys()))
            log.out(f'{nta_file} had equivs flag. Applying equivalent mapping to [ {equivs} ] atom types')
            equivs_count = {i:0 for i in type_equivs}
            for i in nta:
                if nta[i] in type_equivs:
                    equivs_count[nta[i]] += 1
                    equiv = type_equivs[nta[i]]
                    nta[i] = equiv
                    name[i] = equiv
            for i in equivs_count:
                log.out(' - {} -> {} mapping occured on {} atoms'.format(i, type_equivs[i], equivs_count[i]))
                
                
        #######################
        # Perform some Checks #
        #######################
        # Check that length of nta matchs length of m.atoms. If not exit
        if len(nta) != len(m.atoms):
            log.error('{} {} {} {}'.format('ERROR number of atom ids in: ', m.filename, ' are inconsistent with number of atom ids in: ', nta_file))
            
        # Check that every atom id in m.atoms can be mapped onto an atom id from nta dict
        for i in nta:
            try: m.atoms[i]
            except: log.error('{} {} {} {} {} {} {}'.format('ERROR atom id: ', i, ' in ', nta_file, ' cannot be mapped onto an atom type in ', m.filename, ' please make ids consistant'))
                
        # Check that every atom id in m.bonds[i].atomids can be mapped onto an atom id from nta dict
        for i in m.bonds:
            id1, id2 = m.bonds[i].atomids
            
            # Check that atom id1 can be mapped onto a atom type
            try: nta[id1]
            except: log.error('{} {} {} {} {} {} {} {} {} {}'.format('ERROR atom id: ', id1, ' in bond: ', id1, id2, ' in ',  m.filename, ' cannot be mapped onto an atom type in ', nta_file, ' please make ids consistant'))

            # Check that atom id2 can be mapped onto a atom type
            try: nta[id2]
            except: log.error('{} {} {} {} {} {} {} {} {} {}'.format('ERROR atom id: ', id2, ' in bond: ', id1, id2, ' in ',  m.filename, ' cannot be mapped onto an atom type in ', nta_file, ' something wrong with atomIDs in bond'))
                
                
        #############################
        # Add charge info to charge #
        #############################
        # Add charge to charges based on "id" style
        if tmpq['id']:
            for ID in tmpq['id']:
                charges[ID] = tmpq['id'][ID] 
                
        # Add charge to charges based on "type" style
        if tmpq['type']:
            # Only valid file for style type is a LAMMPS .data or .dat file
            if m.filename.endswith('data') or m.filename.endswith('dat'):
                # If check above passes update charge based in type
                for atomid in m.atoms:
                    atom = m.atoms[atomid]
                    # IF atomid has a typeID in tmpq update charges
                    if atom.type in tmpq['type']:
                        charges[atomid] = tmpq['type'][atom.type]
            # Else warn user and exit
            else: log.error('{} {} {}'.format('ERROR charge type used in', nta_file, ' and read in topofile is not a LAMMPS datafile'))
                
        # Add charge to charges based on "nta" style
        if tmpq['nta']:
            for ID in nta:
                if nta[ID] in tmpq['nta']:
                    charges[ID] = tmpq['nta'][nta[ID]]   
                    
                    
        ###############################################################################
        # Check neutralize dictionary that only one neutralization method was defined #
        ###############################################################################
        neutralize_booleans = [neutralize[i] for i in neutralize]
        if neutralize_booleans.count(True) > 1:
            log.error(f'ERROR can only defined one "neutralize system charge" method in the {nta_file}. Currently {neutralize_booleans.count(True)} are defined.')
    return nta, name, edge, charges, neutralize, remove