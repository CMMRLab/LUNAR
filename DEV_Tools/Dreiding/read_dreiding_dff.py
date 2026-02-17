# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
February 17, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

This code will read the info in all2lmp_dreiding.dff file
into dictionaries and classes. This is a completely unqiue
code to read basic info in the all2lmp_dreiding.dff file
format developed by Josh and then a secondary code will
use this to read basic dreiding info and make a formal
.frc file that can be read into all2lmp. This .frc file
should theortically be compatable with material studio
as well.
"""


##############################################################
### Class to Read and Parse all2lmp specific Dreiding file ###
##############################################################

class Atom_types:
    pass # .ver .ref .type .mass, .element, .connection, .comment
    
class Bond:
    pass # .ver .ref .ri .n .ke (set by global Ke from above section)
    
class Angle:
    pass # .ver .ref .theta  .ke (set by global Ke from above section)
    
class Torsion:
    pass # .ver .ref .k .phi0
    
class Nonbond:
    pass # .ver .ref .r0  .d0  .xi
    
class Inversion:
    pass # .ver .ref .r0  .d0
    
class Hbond:
    pass # .ver .ref .dhb .rhb .theta
    
# Function to check if variable is a float
def check_float(variable):
    try:
        float(variable)
        return_boolean = True
    except:
        return_boolean = False
    return return_boolean


class read:
    def __init__(self, filename):
        # set frc file type
        self.type = 'fix-bond'
        
        # atom types info/bondics
        self.atom_types = {}      # {atom type : atom object}
        self.bonds = {}           # {atom type : bond object}
        self.angles = {}          # {atom type : angle object}
        self.torsions = {}        # {tuple(atom types) : torsion object}
        self.inversions = {}      # {atom type : inversion object}
        self.nonbonds = {}        # {atom type : nonbond object}
        self.hbond_gasteiger = {} # {tuple(atom types) : hbond object}
        self.hbond_nocharge = {}  # {tuple(atom types) : hbond object}
        
        
        # Open and read file
        with open(filename, 'r') as f:
            
            # Initializing flags    
            atomtypes_flag = False
            bonds_flag = False
            angles_flag = False
            torsion_flag = False
            inversion_flag = False
            nonbond_flag = False
            hbond_type = ''
            
            # Loop through each line
            for line in f:
                
                # Strip comment's
                line = line.split('!')[0]
                line = line.rstrip()
                
                # Strip line
                line = line.strip()
                string = line # to access string info
                line = line.split()
                
                # Setting flags if '#' character is in string of 1st element in list
                # as False to break each section flag from previous iteration
                if len(line) > 0 and line[0][0] == '#':
                    atomtypes_flag = False
                    bonds_flag = False
                    angles_flag = False
                    torsion_flag = False
                    inversion_flag = False
                    nonbond_flag = False
                    hbond_flag = False

                    
                # atom types info
                if '#atom_types' in line:
                    atomtypes_flag = True
                    
                # bond info
                elif '#harmonic_bond' in string or '#morse_bond' in string and 'DREIDING' in string:
                    bonds_flag = True
                    
                # angle info
                elif '#harmonic_angle' in string and 'DREIDING' in string:
                    angles_flag = True
                    
                # torsion info
                elif '#torsion' in string and 'DREIDING' in string:
                    torsion_flag = True
                    
                # inversion info
                elif '#out_of_plane' in string and 'DREIDING' in string:
                    inversion_flag = True
                    
                # nonbond info
                elif '#nonbond' in string and 'DREIDING' in string:
                    nonbond_flag = True
                    
                # nonbond info
                elif '#nonbond' in string and 'hbonding' in string:
                    hbond_flag = True
                    
                    
                # Find force field info based on flags
                if atomtypes_flag:
                    # Check if len(line) >= 5 and if line[0] is a float
                    if len(line) >= 5 and check_float(line[0]):
                        a = Atom_types()
                        ver = float(line[0])
                        ref = int(line[1])
                        Type = line[2]
                        a.type = Type
                        a.mass = float(line[3])
                        a.element = line[4]
                        a.connection = int(line[5])
                        comments = ''
                        for i in range(6, len(line), 1):
                            comments += ' ' + line[i]
                        a.comment = comments
                        a.ver = ver
                        a.ref = ref
                        self.atom_types[Type] = a
                        
                # Find bond info
                elif bonds_flag:
                    # Find delta and k values from harmonic bond if line is not empty
                    if '#harmonic_bond' in string:
                        harmonic_section = True; morse_section = False;
                    elif '#morse_bond' in string:
                        harmonic_section = False; morse_section = True;
                    
                    # Find delta_harmonic and ke_harmonic
                    if harmonic_section:
                        if line and line[1] == 'delta':
                            delta_harmonic = float(line[-1])
                        if line and line[1] == 'Ke':
                            Ke_harmonic = float(line[-1])
                            
                    # Find delta_morse and ke_morse and de_morse
                    if morse_section:
                        if line and line[1] == 'delta':
                            delta_morse = float(line[-1])
                        if line and line[1] == 'Ke':
                            Ke_morse = float(line[-1])
                        if line and line[1] == 'De':
                            De_morse = float(line[-1])

                    
                    # Check if len(line) >= 5 and if line[0] is a float
                    if len(line) >= 5 and check_float(line[0]):
                        b = Bond()
                        ver = float(line[0])
                        ref = int(line[1])
                        Type = line[2]
                        ri = float(line[3])
                        n = float(line[4])
                        b.type = Type
                        b.ri = ri
                        b.n = n
                        b.ke_harmonic = Ke_harmonic
                        b.delta_harmonic = delta_harmonic
                        b.ke_morse = Ke_morse
                        b.de_morse = De_morse
                        b.delta_morse = delta_morse
                        b.ver = ver
                        b.ref = ref
                        self.bonds[Type] = b
                        
                # Find angle info
                elif angles_flag:                    
                    # Find Kijk
                    if line and line[1] == 'Kijk':
                        Kijk = float(line[-1])
                        #print('Kijk', Kijk)
                        
                    # Check if len(line) == 4 and if line[0] is a float
                    if len(line) == 4 and check_float(line[0]) and check_float(line[3]):
                        #print(line)
                        ag = Angle()
                        ver = float(line[0])
                        ref = int(line[1])
                        Type = line[2]
                        theta = float(line[3])
                        ag.type = Type
                        ag.theta = theta
                        ag.ke = Kijk
                        ag.ver = ver
                        ag.ref = ref
                        self.angles[Type] = ag
                        
                # Find torsion info
                elif torsion_flag:
                    # Check if len(line) == 7 and if line[6] is a float
                    if len(line) == 7 and check_float(line[4]) and check_float(line[6]):
                        #print(line)
                        t = Torsion()
                        ver = float(line[0])
                        ref = int(line[1])
                        Typei = line[2]
                        Typej = line[3]
                        v = float(line[4])
                        n = int(line[5])
                        phi0 = float(line[6])
                        t.v = v
                        t.n = n
                        t.phi0 = phi0
                        t.ver = ver
                        t.ref = ref
                        self.torsions[(Typei, Typej)] = t
                        
                        
                # Find inversion info
                elif inversion_flag:
                    # Check if len(line) == 5 and if line[3] and line[4] is a float
                    if len(line) == 5 and check_float(line[3]) and check_float(line[4]):
                        #print(line)
                        i = Inversion()
                        ver = float(line[0])
                        ref = int(line[1])
                        Type = line[2]
                        k = float(line[3])
                        phi0 = float(line[4])
                        i.k = k
                        i.phi0 = phi0
                        i.ver = ver
                        i.ref = ref
                        self.inversions[Type] = i
                        
                        
                # Find nonbond info
                elif nonbond_flag:
                    # Check if len(line) == 5 and if line[3] and line[4] is a float
                    if len(line) == 6 and check_float(line[3]) and check_float(line[4]):
                        #print(line)
                        nb = Nonbond()
                        ver = float(line[0])
                        ref = int(line[1])
                        Type = line[2]
                        r0 = float(line[3])
                        d0 = float(line[4])
                        xi = float(line[5])
                        nb.r0 = r0
                        nb.d0 = d0
                        nb.xi = xi
                        nb.ver = ver
                        nb.ref = ref
                        self.nonbonds[Type] = nb
                        
                elif hbond_flag:
                    if 'hbonding-' in string:
                        if 'Gasteiger' in string:
                            hbond_type = 'Gasteiger'
                            hbonds_dict = self.hbond_gasteiger
                        if 'no-charge' in string:
                            hbond_type = 'no-charge'
                            hbonds_dict = self.hbond_nocharge
                        continue
                    
                    if hbond_type in ['Gasteiger', 'no-charge'] and len(line) >= 7 and check_float(line[4]) and check_float(line[5]) and check_float(line[6]):
                        hb = Hbond()
                        ver = float(line[0])
                        ref = int(line[1])
                        type_i = line[2]
                        type_j = line[3]
                        dhb    = float(line[4])
                        rhb    = float(line[5])
                        theta  = float(line[6])
                        hb.dhb = dhb
                        hb.rhb = rhb
                        hb.theta = theta
                        hb.ver = ver
                        hb.ref = ref
                        hbonds_dict[(type_i, type_j)] = hb
                        
                        
                        
                    
                    
                    
                    
                
                
                
##########################    
### Testing reading ff ###    
##########################
if __name__ == '__main__':
    
    # frc FF file to read
    file = 'all2lmp_dreiding.dff'    
    ff = read(file)
    
    
    # set pflag
    pflag = True
    pflag = False
    
    # Print if pflag:
    if pflag:
    
        print('\n\n')
        ##############
        # atom types #
        ##############        
        i = 'H_'
        atom = ff.atom_types[i]
        print('----------------------atom types testing---------------------')
        print('type    mass       element     connection    comment')
        print('-------------------------------------------------------------')
        print(f'{atom.type}\t\t{atom.mass}\t\t{atom.element}\t\t\t{atom.connection}\t\t\t{atom.comment}\n\n')
        
        
        ########################
        # harmonic/morse bonds #
        ########################        
        i = 'H_'
        bond = ff.bonds[i]
        print('----------------------------------------------bond testing-----------------------------------------------')
        print('type    Ri          n       Ke-harmonic    delta-harmonic       Ke-morse    delta-morse     De-morse')
        print('---------------------------------------------------------------------------------------------------------')
        print(f'{bond.type}\t\t{bond.ri}\t\t{bond.n}\t\t{bond.ke_harmonic}\t\t\t{bond.delta_harmonic}\t\t\t\t{bond.ke_morse}\t\t{bond.delta_morse}\t\t\t{bond.de_morse}\n\n')
        
        
        ###############
        # angle types #
        ###############        
        i = 'H_'
        angle = ff.angles[i]
        print('------angle testing------')
        print('type    Ke          Theta')
        print('-------------------------')
        print(f'{angle.type}\t\t{angle.ke}\t\t{angle.theta}\n\n')
        
        
        #################
        # torsion types #
        #################        
        i = 'Sp3'; j = 'Sp3'
        tor = ff.torsions[(i,j)]
        print('---------torsion testing---------')
        print('typei   typej    v      n    phi0')
        print('---------------------------------')
        print(f'{i}\t\t{j}\t\t{tor.v}\t\t{tor.n}\t{tor.phi0}\n\n')
        
        
        ###################
        # inversion types #
        ###################        
        i = 'Planar'
        inv = ff.inversions[i]
        print('------inversion testing-----')
        print('type        k          phi0')
        print('----------------------------')
        print(f'{i}\t\t{inv.k}\t\t{inv.phi0}\n\n')
        
        
        #################
        # nonbond types #
        #################        
        i = 'H_'
        nb = ff.nonbonds[i]
        print('--nonbond testing--')
        print('type    r0      d0')
        print('-------------------')
        print(f'{i}\t\t{nb.r0}\t\t{nb.d0}\n\n')