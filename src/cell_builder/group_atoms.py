# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
April 11th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import os


# Function to convert groupIDs to ints
def lstofstrs_2_setofints(linestr, lstofstrs, log):
    setofints = set()
    for i in lstofstrs:
        if ':' in i:
            try:
                sequence = i.split(':')
                if len(sequence) == 2:
                    start = int(sequence[0])
                    end = int(sequence[1])
                    inc = 1
                elif len(sequence) == 3:
                    start = int(sequence[0])
                    end = int(sequence[1])
                    inc = int(sequence[2])
                else: log.error(f'ERROR invalid entry {i} in line: {linestr}. Use either A:B or A:B:C sequence')
                for j in range(start, end+1, inc):
                    setofints.add(int(j))
            except: log.error(f'ERROR invalid entry {i} in line: {linestr}. Use either A:B or A:B:C sequence')  
        else:
            try: setofints.add(int(i))
            except: log.error(f'ERROR invalid entry {i} in line: {linestr}')
    return setofints
    
# Function to read groups file
def read_groups_file(filename, log):
    groups = {} # {filename : { (groupID, style) : set(IDs)} }
    datafile = 'NOT-DEFINED' # Default filename to assign groups to
    strings = '' # strings from looping through file
    valid_styles = ['id', 'type']
    
    # Open file and get strings (use: print(repr(strings)) to see '\n' chars)
    with open(filename, 'r') as f:
        for line in f:
            if line == '' or line.startswith('#'): continue
            strings += line
    
    # Find '&' and '\n' (and and newline chars). There could be a whitespace
    # between them like '& \n' or '&   \n' or they could be next to each other
    # like '&\n'. Thus this needs to be parsed properly.
    and_newline_characters = set(); tmp_chars = ''
    and_char = False; newline_char = False;
    for i in strings:
        if i == '&': and_char = True
        if i == '\n': newline_char = True
        if and_char:
            tmp_chars += i
        if newline_char:
            and_char = False; newline_char = False;
            if tmp_chars != '':
                and_newline_characters.add(tmp_chars)
            tmp_chars = ''
    
    # Loop through and remove any '&' and '\n' grouping
    # to effectively join the lines together.
    for and_newline in and_newline_characters:
        strings = strings.replace(and_newline, '')
    
    # Finally split strings by '\n' character and parse
    lines = [i for i in strings.split('\n') if i != '']
    for linestr in lines:
        # Strip comment's and split by whitespace
        line = linestr.split('#')[0]
        line = line.rstrip()
        line = line.split()

        # Start parsing groups
        if 'read_data' in line and len(line) >= 2:
            datafile = line[1]
            groups[datafile] = {}
        if 'group' in line and len(line) >= 3:
            groupID = line[1]
            style   = line[2]
            if style in valid_styles:
                if datafile == 'NOT-DEFINED':
                    log.error(f'ERROR undefied datafile for groupID {groupID} in {filename}. Use read_data FILENAME.data to assign filename above group.')
                ids = lstofstrs_2_setofints(linestr, line[3:], log)
                if (groupID, style) in groups[datafile]:
                    groups[datafile][(groupID, style)].update(ids)
                else:
                    groups[datafile][(groupID, style)] = ids
            else: log.error(f'ERROR unsupported group style {style}')
    return groups

# Function to add "group" attribute to atoms in the read-in files
def add_group_to_atoms(files, grouping_file, log):
    # Read groups from file
    groups = read_groups_file(grouping_file, log)
    log.out(f'Read-in {grouping_file}')
    for i in groups:
        log.debug('')
        log.debug(i)
        for j in groups[i]:
            log.debug(f'{j} {groups[i][j]}')
            
    # Apply groups to atoms
    for fileID in files:
        m = files[fileID]
        filename = os.path.basename(m.filename)
        for i in m.atoms:
            m.atoms[i].group = 'unassigned'
        for file in groups:
            if file == filename:
                for groupID, style in groups[file]:
                    values = groups[file][(groupID, style)]
                    for i in m.atoms:
                        atom = m.atoms[i]
                        if style == 'type' and atom.type in values:
                            atom.group = groupID
                        if style == 'id' and i in values:
                            atom.group = groupID
    
    # Get all defined groups and perform a debug check
    log.debug('\n\nChecking groups have been assigned')
    defined_groups = set()
    for fileID in files:
        m = files[fileID]
        log.debug('')
        log.debug(m.filename)
        for i in m.atoms:
            atom = m.atoms[i]
            defined_groups.add(atom.group)
            log.debug(f'{i} {atom.group}')
    log.debug('\n\nend of groups\n\n')
    log.out(f'Grouping file {grouping_file} generated {len(defined_groups)} groups:')
    
    # log groups that where generated
    defined_groups = sorted(defined_groups)
    for grp in defined_groups:
        log.out(f' - {grp}')
    return files, defined_groups

# Function to create chunks and create print string
def divide_chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# Write groups file
def write(newfile, defined_groups, m, log, write_group_lmp_file=True):
    lmp_script = 'include.{}_groups.script'.format(newfile)
    log.out(f'\n\nWriting {lmp_script} script')
    with open(lmp_script, 'w') as f:
        f.write('# LAMMPS script to use w/ the "include" command to define groups of atoms.\n')
        f.write('# The following group names have been defined:\n')
        for name in defined_groups:
            f.write('#    - {}\n'.format(name))
            
        # Start writing groups
        f.write('\n\n')
        groups = {i:set() for i in defined_groups}
        for i in m.atoms:
            group = m.atoms[i].group
            groups[group].add(i)
        IDs_per_row = 20
        for group in groups:
            atomIDs = list(sorted(groups[group]))
            chunks = list(divide_chunks(atomIDs, IDs_per_row))
            buffer_str = 'group {}'.format(group)
            buffer_len = len(buffer_str) + 3 # 3 spaces required
            buffer = ''.join([' ' for i in range(buffer_len)])
            for n, chunk in enumerate(chunks):
                chunk = [str(i) for i in chunk]
                if n == 0: lst = ['group'] + [group] + ['id'] + chunk 
                else: lst = [buffer] + chunk
                if n+1 < len(chunks): lst = lst + ['&']
                f.write('{}\n'.format(' '.join(lst)))
            f.write('\n\n')
            
    # Write a simple LAMMPS datafile to visualize the groups of atoms
    if write_group_lmp_file:
        filename = '{}_groups.data'.format(newfile)
        with open(filename,'w') as f: 
            # Write header
            f.write('HEADER, cell_builder.py simple datafile to visualize grouped atoms\n\n') # Make max header length of 220 characters 
            
            f.write(f'{m.natoms} atoms\n')
            if m.nbonds > 0: f.write(f'{m.nbonds} bonds\n')
            f.write('\n')
            
            # Write structure quantity types
            f.write(f'{len(groups)} atom types\n')
            if m.nbonds > 0: f.write('1 bond types\n')
            
            # Write simulation cell
            f.write('\n{0}\n{1}\n{2}\n'.format(m.xbox_line, m.ybox_line, m.zbox_line))
            if m.xy != 0 or m.xz != 0 or m.yz != 0:
                f.write('{:>12.9f} {:^9.9f} {:^9.9f} {} {} {}\n'.format(m.xy, m.xz, m.yz, 'xy', 'xz', 'yz'))
            
            # Write massses
            f.write('\nMasses\n\n')
            for n, group in enumerate(groups, 1): 
                f.write('{:^3} {:^10.5f} # {}\n'.format(n, 1, group))
                    
            # Write bond coeffs
            if m.nbonds > 0:
                # Write bond coeffs section
                f.write('\nBond Coeffs\n\n')
                comment = '{:^2} {:5}'.format('#', 'generic bond type')
                f.write('{:^3} {:^2}\n'.format(1, comment))
                
                
            # Write atoms            
            f.write('\nAtoms # full\n\n') 
            for n, group in enumerate(groups, 1):
                for i in groups[group]:
                    atom = m.atoms[i]
                    f.write('{:^6} {:^4} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2}\n'.format(i, 1, n, atom.charge, atom.x, atom.y, atom.z, atom.ix, atom.iy, atom.iz))
                
            # Write bonds
            if m.nbonds > 0:
                f.write("\nBonds\n\n")
                for i in m.bonds:
                    bond = m.bonds[i]
                    id1, id2 = bond.atomids
                    f.write('{:^2} {:^2} {:^2} {:^2}\n'.format(i, 1, id1, id2))
    return