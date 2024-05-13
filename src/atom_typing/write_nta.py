# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
May 13th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to write log file
def file(mm, name, version, ff_name, include_comments_nta):
    
    # write file
    with open(name, 'w') as f:
        filename = name[:name.rfind('.')]
        f.write(f'New type assignment file for {filename}.data > atom_typing {version} w/{ff_name} atom-types\n')
        if 'general' in ff_name:
            equivs = sorted({mm.atoms[i].nta_type for i in mm.atoms})
            if include_comments_nta:
                f.write('equivs # please manually set equivs. format is two columns with column-1 = set-atom-type-in-"style id" and column-2 = atom-type-to-map-to\n') # If general is used assume user wants to set equivs
            else:
                f.write('equivs')
            for i in equivs:
                if include_comments_nta:
                    comment = '# set-atom-type-in-"style id"  atom-type-to-map-to'
                else: comment = ''
                if ff_name == 'general:0':
                    f.write('{:<20} {}\n'.format(i, comment))
                elif ff_name == 'general:1':
                    f.write('{:<40} {}\n'.format(i, comment))
                elif ff_name == 'general:2':
                    f.write('{:<60} {}\n'.format(i, comment))
                elif ff_name == 'general:3':
                    f.write('{:<80} {}\n'.format(i, comment))
                elif ff_name == 'general:4':
                    f.write('{:<100}  {}\n'.format(i, comment))
                else:
                    f.write('{:<20} {}\n'.format(i, comment))
            f.write('\n')
        f.write('style id\n') # Format will be id since every atomid has been given a unqiue atom-type
        for i in mm.atoms:
            atom = mm.atoms[i]
            if include_comments_nta:
                comment = '{:^2} {:5}'.format('#', atom.nta_info)
            else: comment = ''
            if ff_name == 'general:0':
                f.write('{:<6} {:<20} {}\n'.format(i, atom.nta_type, comment))
            elif ff_name == 'general:1':
                f.write('{:<6} {:<40} {}\n'.format(i, atom.nta_type, comment))
            elif ff_name == 'general:2':
                f.write('{:<6} {:<60} {}\n'.format(i, atom.nta_type, comment))
            elif ff_name == 'general:3':
                f.write('{:<6} {:<80} {}\n'.format(i, atom.nta_type, comment))
            elif ff_name == 'general:4':
                f.write('{:<6} {:<100} {}\n'.format(i, atom.nta_type, comment))
            else:
                f.write('{:<6} {:<15} {}\n'.format(i, atom.nta_type, comment))
    return
