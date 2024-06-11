# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
June 11th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
# Function to write out subcells into a LAMMPS dump file
def subcells2dump(basename, subcells):
    with open(basename+'_grp_local_subcells.dump','w') as f:
        for n, m in enumerate(subcells):
            f.write('ITEM: TIMESTEP\n')
            f.write('{}\n'.format(n))
            f.write('ITEM: NUMBER OF ATOMS\n')
            f.write('{}\n'.format(len(m.atoms)))
            
            xline = m.xbox_line.split()
            yline = m.ybox_line.split()
            zline = m.zbox_line.split()
            f.write('ITEM: BOX BOUNDS pp pp pp\n')
            f.write('{} {}\n'.format(xline[0], xline[1]))
            f.write('{} {}\n'.format(yline[0], yline[1]))
            f.write('{} {}\n'.format(zline[0], zline[1]))
            
            f.write('ITEM: ATOMS x y z type id q\n')
            for i in m.atoms:
                atom = m.atoms[i]
                f.write('{:^15.9f} {:^15.9f} {:^15.9f} {:^3} {:^3} {:^10.6f}\n'.format(atom.x, atom.y, atom.z, atom.type, i, atom.charge))
    return