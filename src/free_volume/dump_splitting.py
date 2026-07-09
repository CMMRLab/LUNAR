# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
July 9, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.read_dump as read_dump
import math

#-------------------------------------------------------------#
# Function for stringing together float values for parameters #
#-------------------------------------------------------------#
def string_parameters(coeff):
    string = ''
    for i in coeff:
        if isinstance(i, float):
            string += '{:>16.8f}'.format(i)
        elif isinstance(i, int):
            string += '{:>16}'.format(i)
        else:
            string += '{:>16}'.format(i)
    return string


#---------------------------------------------------------------------#
# Function to take read topofile and dump file to write new data file #
#---------------------------------------------------------------------#
def topo_dump_to_data(m, box, atoms, filename):
    # Writing m file with m atoms
    with open(filename,'w') as f: 
        # Write header
        header = str(m.header) +  ' > LUNAR/free_volume.py dump splitting option'
        f.write(f'{header[-220:len(header)]}\n\n') # Make max header length of 220 characters 
        
        # Write structure quantities
        f.write(f'{len(atoms)} atoms\n')
        
        # Write structure quantity types
        f.write(f'{m.natomtypes} atom types\n')
        
        # write box size
        f.write('\n')
        f.write(f"{box['xlo']: .16e} {box['xhi']: .16e} xlo xhi\n")
        f.write(f"{box['ylo']: .16e} {box['yhi']: .16e} ylo yhi\n")
        f.write(f"{box['zlo']: .16e} {box['zhi']: .16e} zlo zhi\n")
        if box['xy'] != 0.0 or box['xz'] != 0.0 or box['yz'] != 0.0:
            f.write(f"{box['xy']: .16e} {box['xz']: .16e} {box['yz']: .16e} xy xz yz\n")
            
        # Write massses
        #f.write(f'Masses  # {m.mass_coeffs_style_hint}\n\n')
        f.write('\nMasses\n\n')
        type_ids = set()
        for i in m.masses: 
            type_ids.add( i )
            coeff = m.masses[i]
            parms = coeff.coeffs
            if coeff.type != 'N/A':
                f.write('{:^3} {:^10.5f} # {}\n'.format(i, parms[0], coeff.type))
            else:
                f.write('{:^3} {:^10.5f}\n'.format(i, parms[0]))
            
        # Write pair coeffs
        if m.pair_coeffs:
            if m.pair_coeffs_style_hint != 'N/A':
                f.write(f'\nPair Coeffs  # {m.pair_coeffs_style_hint}\n\n')
            else:
                f.write('\nPair Coeffs\n\n')
            for i in m.pair_coeffs: 
                pair = m.pair_coeffs[i]
                if pair.type != 'N/A':
                    f.write('{:^3} {} # {}\n'.format(i, string_parameters(pair.coeffs), pair.type))
                else:
                    f.write('{:^3} {}\n'.format(i, string_parameters(pair.coeffs)))
                    
        # Write atoms in full charge style
        f.write('\nAtoms # charge\n\n')   
        for n, atom in enumerate(atoms):

            type_id = atom.get('type', 0)
            if type_id not in type_ids:
                raise Exception(f'ERROR trying to write a dump file from a datafile FF and atomTypeID: {type_id} from dump file is not in the corresponding datafile.')
                
            i  = atom.get('id', n)
            q  = atom.get('q', 0.0)
            x  = atom.get('x', 0.0)
            y  = atom.get('y', 0.0)
            z  = atom.get('z', 0.0)
            ix = atom.get('ix', 0.0)
            iy = atom.get('iy', 0.0)
            iz = atom.get('iz', 0.0)
            f.write('{:^6} {:^2} {:^10.6f} {:^15.9f} {:^15.9f} {:^15.9f} {:^2} {:^2} {:^2}\n'.format(i, type_id, q, x, y, z, ix, iy, iz)) 
    return


#-----------------------------------------------#
# Dump Settings parsing that are a string like: #
#  "start=1; end=10; nevery=1; style=step"      #
#  "start=1; end=10; nevery=1; style=section"   #
#-----------------------------------------------#
# Function to parse misc string into dict of settings 
def get_misc_setting(setting_string):
    # Setup the globals namespace to limit scope of what eval() can do
    allowed_builtins = ['min','max','sum','abs','len','map','range','reversed']
    copied_builtins = globals()['__builtins__'].copy()
    globals_dict = {'e':math.e, 'pi':math.pi}
    globals_dict['__builtins__'] = {key:copied_builtins[key] for key in allowed_builtins}
    
    # Parse misc string
    setting = {} # {keyword:float or int or Boolean}
    tmp1 = setting_string.split(';')
    for tmp2 in tmp1:
        tmp3 = tmp2.split('=')
        if len(tmp3) >= 2:
            i = tmp3[0].strip()
            try: j = eval(tmp3[1], globals_dict)
            except: j = str(tmp3[1])
            setting[i] = j
    return setting
    
    
# Test reading dumpfile
if __name__ == "__main__": 
    dumpfile = '../../EXAMPLES/free_volume/wildcard_DETDA/detda_rep_1_free_volume_test.dump'
    dump = read_dump.read_dumpfile(dumpfile)
    
    box, atoms = dump.parse_atoms_and_box(0)
    print(f"{box['xlo']: .16e} {box['xhi']: .16e} xlo xhi")
    print(f"{box['ylo']: .16e} {box['yhi']: .16e} ylo yhi")
    print(f"{box['zlo']: .16e} {box['zhi']: .16e} zlo zhi")
    if box['xy'] != 0.0 or box['xz'] != 0.0 or box['yz'] != 0.0:
        print(f"{box['xy']: .16e} {box['xz']: .16e} {box['yz']: .16e} xy xz yz")
        
    for atom in atoms:
        print(atom)

    
    
    