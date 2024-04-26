# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
Decemeber 15th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
#############################
# Function to create chunks #
#############################
def divide_chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

        
###########################################
# Function to add whitepsace where needed #
###########################################
def manage_whitespace(record, nchars, justification, log):
    string = str(record)
    if len(string) > nchars: 
        log.error(f'ERROR write_pdb.py failed to write .pdb file. Character overflow, max chars = {nchars}; record chars = {len(string)}')
    
    # Manage white characters
    nempty = nchars - len(string)
    if nempty > 0:
        empty = ' '*nempty
        if justification == 'left':
            string = '{}{}'.format(string, empty)
        elif justification == 'right':
            string = '{}{}'.format(empty, string)
        else: log.error(f'ERROR justification method {justification} not supported')
    return string

############################################
# Function to manage whitespace for floats #
############################################
def manage_floats(position, nchars, log):
    string = '{:.3f}'.format(position)
    if len(string) > nchars: 
        log.error(f'ERROR write_pdb.py failed to write .pdb file. Position character overflow, max chars = {nchars}; record chars = {len(string)}')
        
    # Manage white characters
    nempty = nchars - len(string)
    if nempty > 0:
        empty = ' '*nempty
        string = '{}{}'.format(empty, string)
    return string

#########################################################
# Function to insert string into string. NOTE lo starts #
# at 1 to make it easier to align w/ pdb file rules     #
#########################################################
def insert_into_string(lo, hi, string, chars):
    lst = [i for i in string]
    if lo == hi:
        lst[lo-1] = chars
    else:
        for n, i in enumerate(range(lo-1, hi)):
            try: lst[i] = chars[n]
            except: pass
    string = ''.join(lst)
    return string


################################
# Function to write .pdb files #
################################
def file(m, doc_title, version, ff_name, pdb_file, log):
    # Generate graph
    graph = {i:set() for i in m.atoms}
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        graph[id1].add(id2)
        graph[id2].add(id1)
        
    atom_names = {} # { atomID : atom name }
    if pdb_file == 'types':
        log.out('Writing optional .pdb file in types format')
        for i in m.atoms:
            atom_names[i] = m.atoms[i].nta_type
    elif pdb_file == 'typeIDs':
        log.out('Writing optional .pdb and .nta file in typeIDs format')
        atomtypes_lst = sorted({m.atoms[i].nta_type for i in m.atoms})
        nta_types_forward = {n:i for n, i in enumerate(atomtypes_lst, 1)}
        nta_types_reverse = {i:n for n, i in enumerate(atomtypes_lst, 1)}
        for i in m.atoms:
            atom_names[i] = nta_types_reverse[m.atoms[i].nta_type]
        
        # write ntafile
        with open(doc_title+'.nta','w') as f: 
            f.write(f'New type assignment built with atom_typing.py {version} w/{ff_name} atom-types\n')
            f.write('style type\n') # Format will be type
            for i in nta_types_forward:
                f.write('{:<6} {:<15}\n'.format(i, nta_types_forward[i]))

    
    # Writing pdb file
    with open(doc_title+'.pdb','w') as f: 
        # Write intro
        f.write(f'TITLE    Built with atom_typing.py {version} w/{ff_name} atom-types\n') 
        f.write('REMARK   Atom Name columns: 13-16, will be used for atomType or atomTypeID\n') 
        f.write('REMARK   1         2         3         4         5         6         7         8\n') 
        f.write('REMARK78901234567890123456789012345678901234567890123456789012345678901234567890\n') 
        
        # Write atoms
        for i in m.atoms:
            atom = m.atoms[i]; string = 80*' '
            record_type = manage_whitespace('HETATM', 6, 'left', log) # 1-6
            string = insert_into_string(1, 6, string, record_type) # insert into string
            
            atom_serial_number = manage_whitespace(i, 5, 'right', log) # 7-11
            string = insert_into_string(7, 11, string, atom_serial_number) # insert into string
            
            atom_name = manage_whitespace(atom_names[i], 4, 'left', log) # 13-16
            string = insert_into_string(13, 16, string, atom_name) # insert into string
            
            location_indicator = manage_whitespace('', 1, 'right', log) # 17 (leaving empty for now...)
            string = insert_into_string(17, 17, string, location_indicator) # insert into string
            
            residue_name = manage_whitespace('', 3, 'right', log) # 18-20 (leaving empty for now...)
            string = insert_into_string(18, 20, string, residue_name) # insert into string
            
            chain_identifier = manage_whitespace('', 1, 'right', log) # 22 (leaving empty for now...)
            string = insert_into_string(22, 22, string, chain_identifier) # insert into string
            
            residue_sequence = manage_whitespace(atom.molid, 4, 'right', log) # 23-26
            string = insert_into_string(23, 26, string, residue_sequence) # insert into string
            
            residue_insertion_code = manage_whitespace('', 1, 'right', log) # 27 (leaving empty for now...)
            string = insert_into_string(27, 27, string, residue_insertion_code) # insert into string
            
            x = manage_floats(atom.x, 8, log) # 31-38
            string = insert_into_string(31, 38, string, x) # insert into string
            
            y = manage_floats(atom.y, 8, log) # 39-46
            string = insert_into_string(39, 46, string, y) # insert into string
            
            z = manage_floats(atom.z, 8, log) # 47-54
            string = insert_into_string(47, 54, string, z) # insert into string
            
            occupancy = manage_whitespace('', 6, 'right', log) # 55-60 (leaving empty for now...)
            string = insert_into_string(55, 60, string, occupancy) # insert into string
            
            temp_factor = manage_whitespace('', 6, 'right', log) # 61-66 (leaving empty for now...)
            string = insert_into_string(61, 66, string, temp_factor) # insert into string
            
            segment_identifier = manage_whitespace('', 6, 'right', log) # 73-76 (leaving empty for now...)
            string = insert_into_string(73, 76, string, segment_identifier) # insert into string
            
            element = manage_whitespace(atom.element, 2, 'left', log) # 77-78
            string = insert_into_string(77, 78, string, element) # insert into string
            f.write(f'{string}\n')
            
        # Write connects
        graphIDs = sorted(list(graph.keys()))
        for i in graphIDs:
            if graph[i]:
                bonded = sorted(graph[i])
                chunks = list(divide_chunks(bonded, 4))
                for chunk in chunks:
                    string = 31*' '; bondIDs = 25*' '; start_index = 0; char_span = 5;
                    for ID_local, ID in enumerate(chunk):
                        lo_index = int(ID_local*char_span + start_index + 1)
                        hi_index = int(ID_local*char_span + start_index + char_span)
                        ID = manage_whitespace(ID, 5, 'right', log)
                        bondIDs = insert_into_string(lo_index, hi_index, bondIDs, ID) # insert into string
                    
                    record_type = manage_whitespace('CONECT', 6, 'left', log) # 1-6
                    string = insert_into_string(1, 6, string, record_type) # insert into string
                
                    atom_serial_number = manage_whitespace(i, 5, 'right', log) # 7-11
                    string = insert_into_string(7, 11, string, atom_serial_number) # insert into string
                    string = insert_into_string(12, 31, string, bondIDs) # insert into string
                    f.write(f'{string}\n')
        f.write('END\n')
    return 