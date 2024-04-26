# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
March 1st, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
def atoms(cta_file, m, log):
    
    # Opening and reading the frc file
    with open(cta_file, 'r') as f:
        cta = {} # { atomid : new atom type } (For style id)
        types = {} # { atomtype : new atom type } (For style type)
        style = '' # set as empty string and update once found
        
        # Initializing flags    
        atomtypes_flag = False
        
        # Looping through each line of the frc file
        for line in f:
            
            # Strip comment's
            line = line.split('#')[0]
            line = line.rstrip()
            
            # Split line
            line = line.strip()
            line_split = line.split()

            
            # Setting flags
            if line == '':
                atomtypes_flag = False
            elif 'style' in line:
                style = line_split[1] # should be in second index location
                atomtypes_flag = True
                continue

            
            # Finding were atom types exist in the file
            if atomtypes_flag:
                # Set nta based on style id
                if style == 'id':
                    atomid = int(line_split[0])
                    atom_type = line_split[1]
                    cta[atomid] = atom_type
                
                # Set nta based on style type
                elif style == 'type':
                    typeid = int(line_split[0])
                    atom_type = line_split[1]
                    types[typeid] = atom_type
                
                # Exit if style is not type or id
                else:
                    log.error('ERROR requesting style of unsupported type in: {} only supported style are id or type'.format(cta_file))
                    
                    
        # If style was type build cta dict by looping through m.atoms
        if style == 'type':
            # Check that number of atom types in datafile are consistent with number of atom types in nta_file
            natomtypes_datafile = m.natomtypes; natomtypes_nta_file = len(types);
            if natomtypes_datafile != natomtypes_nta_file:
                log.error('ERROR number of atom types in: {} are inconsistent with number of atom types in: {}'.format(m.filename, cta_file))
                
            # If check above passes build nta dict
            for atomid in m.atoms:
                # Find type based on type in m.atoms and type in types
                atom = m.atoms[atomid]
                atom_type = types[atom.type]
                cta[atomid] = atom_type
            
        # Check that length of nta matchs length of m.atoms. If not exit
        if len(cta) != len(m.atoms):
            log.error('{} {} {} {}'.format('ERROR number of atom ids in: ', m.filename, ' are inconsistent with number of atom ids in: ', cta_file))
            
        # Check that every atom id in m.atoms can be mapped onto an atom id from nta dict
        for i in cta:
            try: m.atoms[i]
            except:
                log.error('{} {} {} {} {} {} {}'.format('ERROR atom id: ', i, ' in ', cta_file, ' cannot be mapped onto an atom type in ', m.filename, ' please make ids consistant'))
                
        # Check that every atom id in m.bonds[i].atomids can be mapped onto an atom id from nta dict
        for i in m.bonds:
            id1, id2 = m.bonds[i].atomids
            
            # Check that atom id1 can be mapped onto a atom type
            try: cta[id1]
            #if id1 not in nta_atomids:
            except:
                log.error('{} {} {} {} {} {} {} {} {} {}'.format('ERROR atom id: ', id1, ' in bond: ', id1, id2, ' in ',  m.filename, ' cannot be mapped onto an atom type in ', cta_file, ' please make ids consistant'))

            # Check that atom id2 can be mapped onto a atom type
            try: cta[id2]
            #if id2 not in nta_atomids:
            except:
                log.error('{} {} {} {} {} {} {} {} {} {}'.format('ERROR atom id: ', id2, ' in bond: ', id1, id2, ' in ',  m.filename, ' cannot be mapped onto an atom type in ', cta_file, ' something wrong with atomIDs in bond'))
    
    return cta



    
        
