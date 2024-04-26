# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
September 19th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to write log file
def file(mm, name, version, ff_name):
    
    # write file
    with open(name, 'w') as f:
        filename = name[:name.rfind('.')]
        f.write(f'New type assignment file for {filename}.data > atom_typing {version} w/{ff_name} atom-types\n')
        if 'general' in ff_name:
            equivs = sorted({mm.atoms[i].nta_type for i in mm.atoms})
            f.write('equivs # please manually set equivs. format is two columns with column-1 = set-atom-type-in-"style id" and column-2 = atom-type-to-map-to\n') # If general is used assume user wants to set equivs
            for i in equivs:
                if ff_name == 'general:0':
                    f.write('{:<20} # set-atom-type-in-"style id"  atom-type-to-map-to\n'.format(i))
                elif ff_name == 'general:1':
                    f.write('{:<40} # set-atom-type-in-"style id"  atom-type-to-map-to\n'.format(i))
                elif ff_name == 'general:2':
                    f.write('{:<60} # set-atom-type-in-"style id"  atom-type-to-map-to\n'.format(i))
                elif ff_name == 'general:3':
                    f.write('{:<80} # set-atom-type-in-"style id"  atom-type-to-map-to\n'.format(i))
                elif ff_name == 'general:4':
                    f.write('{:<100} # set-atom-type-in-"style id"  atom-type-to-map-to\n'.format(i))
                else:
                    f.write('{:<20} # set-atom-type-in-"style id"  atom-type-to-map-to\n'.format(i))
            f.write('\n')
        f.write('style id\n') # Format will be id since every atomid has been given a unqiue atom-type
        for i in mm.atoms:
            atom = mm.atoms[i]
            comment = '{:^2} {:5}'.format('#', atom.nta_info)
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
