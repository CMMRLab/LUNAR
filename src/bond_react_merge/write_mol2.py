# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 8th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Import Necessary Libraries
from collections import OrderedDict


def file(name, m, file, version):
    
    # Find box dimensions to remove periodic boundary conditions
    xline = m.xbox_line.split()
    yline = m.ybox_line.split()
    zline = m.zbox_line.split()
    lx = float(xline[1])-float(xline[0])
    ly = float(yline[1])-float(yline[0])
    lz = float(zline[1])-float(zline[0])
    
    # set max_x, max_y, max_z w/ minimum image convention
    max_x = lx/2
    max_y = ly/2
    max_z = lz/2
    
    # Find PBC bonds
    flagged_bondids = []; flagged_bond_atoms = []; # to set molids as zeros for PBC bonding atoms
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        pbc_flag = False
        
        # Find x1, y1, z1, x2, y2, and z2
        x1 = m.atoms[id1].x
        y1 = m.atoms[id1].y
        z1 = m.atoms[id1].z
        x2 = m.atoms[id2].x
        y2 = m.atoms[id2].y
        z2 = m.atoms[id2].z
        
        # Find differences
        diffx = abs(x2 - x1)
        diffy = abs(y2 - y1)
        diffz = abs(z2 - z1)
        
        # if bond fails minimum image convention it is periodic
        if diffx > max_x:
            pbc_flag = True
    
        # if bond fails minimum image convention it is periodic
        if diffy > max_y:
            pbc_flag = True
    
        # if bond fails minimum image convention it is periodic
        if diffz > max_z:
            pbc_flag = True
            
        # If pbc_flag the bond was found to be periodic and flag its id
        if pbc_flag:
            flagged_bondids.append(i)
            flagged_bond_atoms.append(id1)
            flagged_bond_atoms.append(id2)
            
    
    # Find molids and bonds
    molids = sorted({m.atoms[i].molid for i in m.atoms})
    bonding_atoms = []; nbonds = 0;
    for i in m.bonds:
        if i not in flagged_bondids:
            nbonds += 1
            id1, id2 = m.bonds[i].atomids
            bonding_atoms.append(id1)
            bonding_atoms.append(id2)
    # Find new number of atoms
    natoms = 0;
    for i in m.atoms:
        natoms += 1
            
            

    # Writing new file with bonds information
    # https://chemicbook.com/2021/02/20/mol2-file-format-explained-for-beginners-part-2.html
    with open(name+'.mol2','w') as f: 
        
        # Write molecule section
        f.write('@<TRIPOS>MOLECULE\n')
        header = '{} > bond_react_merge: {} .mol2 file (filetag: {})'.format(m.header, version, file)
        f.write(f'{header[-220:len(header)]}\n') # Make max header length of 220 characters 
        f.write(f'  {natoms} {nbonds}    {len(molids)}    0    0\n')
        f.write('SMALL\n')
        f.write('NO_CHARGES\n')
        f.write('****\n')
        f.write('Energy = 0\n')
        
        # Write Atoms info
        f.write('\n@<TRIPOS>ATOM\n')
        m.atoms = dict(OrderedDict(sorted(m.atoms.items())))
        for i in m.atoms:
            atom = m.atoms[i]
            comment = atom.comment
            
            ###################
            # Find atoms info #
            ###################
            # if / in comment get element, else use atom type 1st character
            if '/' in comment:
                elem = comment[comment.rfind('/')+1:].capitalize()
            else:
                elem = comment[0].capitalize()
            
            # Find atom positions
            x = '{:>17.4f}'.format(atom.x) # float point
            y = '{:>10.4f}'.format(atom.y) # float point
            z = '{:>10.4f}'.format(atom.z) # float point
            
            # Find sub structure names and nomenclature (otherwise known as residues)
            try: subst_id = atom.molid
            except: subst_id = 0
            subst_name = '****' # The name of the substructure containing the atom (string)
            subst_name = '{:>4}'.format('m'+str(subst_id)) # Will give access through [VMD ResName Coloring]
            charge = '{:>10.4f}'.format(atom.charge)
            
            # Write atoms info
            f.write('{:>7} {} {} {} {} {} {:>7} {:>7} {}\n'.format(i, elem, x, y, z, elem, subst_id, subst_name, charge))
            
        # Write Bonds info
        f.write('@<TRIPOS>BOND\n')
        newid = 0 # For skipping periodic bonds
        m.bonds = dict(OrderedDict(sorted(m.bonds.items())))
        for i in m.bonds:
            bond = m.bonds[i]
            id1, id2 = bond.atomids

                        
            # Set bond_type as 1 for now...
            # Possible options:
            #    1 = single
            #    2 = double
            #    3 = triple
            #    am = amide
            #    ar = aromatic
            #    du = dummy
            #    un = unknown (cannot be determined from the parameter tables)
            #    nc = not connected
            bond_type = '1'
            
            # Write bonds info
            if i not in flagged_bondids:
                newid += 1
                f.write('{:>6} {:>6} {:>6} {:>6}\n'.format(newid, id1, id2, bond_type))