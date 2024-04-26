# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 12:19:49 2023

@author: jdkem


Atom Type Definitions:
    parm-type   typeID  classID atomType      comment      atomic-number    mass     nb
    atom          108     48      CA      "Phenol C-OH"    6                12.011    3
    
Van der Waals Parameters:
    parm-type     type or classID        sigma      epsilon
    vdw                  1               2.9400     0.0610
    
Bond Stretching Parameters:
    parm-type     typeID1  typeID2        k2         r0
    bond            1        2          367.00     1.3800
    
Angle Bending Parameters:
    parm-type     typeID1  typeID2  typeID3         k       theta0
    angle           25        1        25         33.00     109.47
    
Torsional Parameters:
    parm-type     typeID1  typeID2  typeID3   typeID4    k1     ?1 1     k2      ?1  2    k3   ?1  3
    torsion       0           2        2        2       -2.500 0.0 1     1.250 180.0 2   3.100 0.0 3
    
Improper Torsional Parameters:
    parm-type     typeID1  typeID2  typeID3   typeID4    kchi    Chi0   n 
    imptors       0          0        3          4       21.000  180.0  2
    
Atomic Partial Charge Parameters:
    parm-type   classID          partial charge
    charge        1                -0.2200
"""
class Atom_types:
    pass # .type .mass, .atomic_number .connection .comment .typeID .classID
    
class Quadratic_bonds:
    pass # .r0, .k2
    
class Quadratic_angles:
    pass # .theta0, .k2
    
class Torsion_1_opls:
    pass # .k1 .k2 .k3 .k4
    
class OOP:
    pass # .kphi, .n, .phi0
    
class Pair_coeff:
    pass # .r, .eps  


###########################################################################
# Open and read .prm file to store parameters in dictionaries and classes #
###########################################################################
class prm:
    def __init__(self, filename):
        self.atom_types = {}          # { tuple(atom-type, typeID) : atom object}
        self.charges = {}             # { typeID : partial charge}
        self.quadratic_bonds = {}     # { tuple(classID1, classID2) : quadratic bond object}
        self.quadratic_angles = {}    # { tuple(classID1, classID2, classID3) : quadratic angle object}
        self.torsion_1_opls = {}      # { tuple(classID1, classID2, classID3, classID4) : torsion1 object}
        self.out_of_plane = {}        # { tuple(classID1, classID2, classID3, classID4) : oop object}
        self.pair_coeffs_12_6 = {}    # { ID (Type or Class) : pair coeff object}
        
        with open(filename, 'r') as f:
            for wholeline in f:
                wholeline = wholeline.rstrip()
                line = wholeline.split()
                if len(line) < 1: continue # skip short lines
                if len(line) >= 1 and '#' in line[0]: continue # skip comment lines
                
                # Find atom types
                if line[0] == 'atom' and len(line) >= 8:
                    typeID = int(line[1]); i = line[3];
                    a = Atom_types()
                    a.TypeID = int(line[1])
                    a.classID = int(line[2])
                    a.atomic_number = int(line[-3])
                    a.mass = float(line[-2])
                    a.connection = int(line[-1])
                    a.comment = parse_atom_type_comments(wholeline)
                    a.type = i
                    self.atom_types[(i, typeID)] = a
                    
                # Find partial charges
                if line[0] == 'charge' and len(line) >= 3:
                    self.charges[int(line[1])] = float(line[2])
                    
                # Find bonds
                if line[0] == 'bond' and len(line) >= 5:
                    qb = Quadratic_bonds()
                    classID1 = int(line[1]); classID2 = int(line[2]);
                    qb.k2 = float(line[3])
                    qb.r0 = float(line[4])
                    self.quadratic_bonds[(classID1,classID2)] = qb
                    
                # Find angles
                if line[0] == 'angle' and len(line) >= 6:
                    qa = Quadratic_angles()
                    classID1 = int(line[1]); classID2 = int(line[2]); classID3 = int(line[3]);
                    qa.k2 = float(line[4])
                    qa.theta0 = float(line[5])
                    self.quadratic_angles[(classID1,classID2,classID3)] = qa
                    
                # Find torsion
                if line[0] == 'torsion' and len(line) >= 14:
                    t1o = Torsion_1_opls()
                    classID1 = int(line[1]); classID2 = int(line[2]); classID3 = int(line[3]);  classID4 = int(line[4]);
                    t1o.k1 = float(line[5])
                    t1o.k2 = float(line[8])
                    t1o.k3 = float(line[11])
                    t1o.k4 = float(0.0)
                    self.torsion_1_opls[(classID1,classID2,classID3,classID4)] = t1o
                    
                # Find oops
                if line[0] == 'imptors' and len(line) >= 8:
                    o = OOP()
                    classID1 = int(line[1]); classID2 = int(line[2]); classID3 = int(line[3]);  classID4 = int(line[4]);
                    o.kchi = float(line[5])
                    o.chi0 = float(line[6])
                    o.n = int(line[7])
                    self.out_of_plane[(classID1,classID2,classID3,classID4)] = o
                    
                # Find vdw
                if line[0] == 'vdw' and len(line) >= 4:
                    p = Pair_coeff()
                    ID = int(line[1]);
                    p.r = float(line[2])
                    p.eps = float(line[3])
                    self.pair_coeffs_12_6[ID] = p
                    
                    
# Parser for reading atom type comments
def parse_atom_type_comments(line):
    string = ''; flag = False; count = 0;
    bounding_characters = ['"', "'"] # " or ' bounding chars
    for i in line:
        if i in bounding_characters: flag = True; count += 1; continue;
        if count >= 2: flag = False
        if flag: string += i
    return string.rstrip()


###########################    
### Testing reading prm ###    
###########################
if __name__ == '__main__':
    filename = 'oplsaa.prm'    
    parms = prm(filename)
    pflag = True # print flag
    if pflag:
        # Atom types test for oplsaa.prm
        print('\n\n-------------------------------------atom types testing-------------------------------------')
        print('{:^10} {:^10} {:^10} {:^10} {:^15} {:^15} {:^15}'.format('type', 'typeID', 'classID', 'mass', 'atomic_number', 'connection', 'comment'))
        print('--------------------------------------------------------------------------------------------')
        atom = parms.atom_types[('OH', 96)]
        print('{:^10} {:^10} {:^10} {:^10} {:^15} {:^15} {:^15}'.format(atom.type, atom.TypeID, atom.classID, atom.mass, atom.atomic_number, atom.connection, atom.comment))
        atom = parms.atom_types[('CT', 100)]
        print('{:^10} {:^10} {:^10} {:^10} {:^15} {:^15} {:^15}'.format(atom.type, atom.TypeID, atom.classID, atom.mass, atom.atomic_number, atom.connection, atom.comment))
        atom = parms.atom_types[('CA', 163)]
        print('{:^10} {:^10} {:^10} {:^10} {:^15} {:^15} {:^15}'.format(atom.type, atom.TypeID, atom.classID, atom.mass, atom.atomic_number, atom.connection, atom.comment))
        
        # vdw test for oplsaa.prm
        print('\n\n----------------vdw testing----------------')
        print('{:^10} {:^15} {:^15}'.format('typeID', 'r', 'eps'))
        print('-------------------------------------------')
        vdw = parms.pair_coeffs_12_6[1]
        print('{:^10} {:^15} {:^15}'.format(1, vdw.r, vdw.eps))
        vdw = parms.pair_coeffs_12_6[4]
        print('{:^10} {:^15} {:^15}'.format(4, vdw.r, vdw.eps))
        vdw = parms.pair_coeffs_12_6[12]
        print('{:^10} {:^15} {:^15}'.format(12, vdw.r, vdw.eps))
        
        # bond test for oplsaa.prm
        print('\n\n----------------bond testing---------------')
        print('{:^10} {:^15} {:^15}'.format('classID', 'k2', 'r0'))
        print('-------------------------------------------')
        bond = parms.quadratic_bonds[(1,2)]
        print('{:^10} {:^15} {:^15}'.format('(1,2)', bond.k2, bond.r0))
        bond = parms.quadratic_bonds[(6,20)]
        print('{:^10} {:^15} {:^15}'.format('(6,20)', bond.k2, bond.r0))
        bond = parms.quadratic_bonds[(13,50)]
        print('{:^10} {:^15} {:^15}'.format('(13,50)', bond.k2, bond.r0))
        
        # angle test for oplsaa.prm
        print('\n\n---------------angle testing---------------')
        print('{:^10} {:^15} {:^15}'.format('classIDs', 'k2', 'theta0'))
        print('-------------------------------------------')
        angle = parms.quadratic_angles[(25,1,25)]
        print('{:^10} {:^15} {:^15}'.format('(25,1,25)', angle.k2, angle.theta0))
        angle = parms.quadratic_angles[(12,3,12)]
        print('{:^10} {:^15} {:^15}'.format('(12,3,12)', angle.k2, angle.theta0))
        angle = parms.quadratic_angles[(4,3,57)]
        print('{:^10} {:^15} {:^15}'.format('(4,3,57)', angle.k2, angle.theta0))
        
        # torsion test for oplsaa.prm
        print('\n\n------------------------------torsion testing------------------------------')
        print('{:^15} {:^15} {:^15} {:^15} {:^15}'.format('classIDs', 'k1', 'k2', 'k3', 'k4'))
        print('---------------------------------------------------------------------------')
        tor = parms.torsion_1_opls[(0,2,2,2)]
        print('{:^15} {:^15} {:^15} {:^15} {:^15}'.format('(0,2,2,2)', tor.k1, tor.k2, tor.k3, tor.k4))
        tor = parms.torsion_1_opls[(24,3,13,21)]
        print('{:^15} {:^15} {:^15} {:^15} {:^15}'.format('(24,3,13,21)', tor.k1, tor.k2, tor.k3, tor.k4))
        tor = parms.torsion_1_opls[(47,3,24,45)]
        print('{:^15} {:^15} {:^15} {:^15} {:^15}'.format('(47,3,24,45)', tor.k1, tor.k2, tor.k3, tor.k4))
        
        # oop test for oplsaa.prm
        print('\n\n------------------------oop testing------------------------')
        print('{:^15} {:^15} {:^15} {:^15}'.format('classIDs', 'kchi', 'chi0', 'n'))
        print('-----------------------------------------------------------')
        oop = parms.out_of_plane[(0,0,3,4)]
        print('{:^15} {:^15} {:^15} {:^15}'.format('(0,0,3,4)', oop.kchi, oop.chi0, oop.n))
        oop = parms.out_of_plane[(0,0,3,52)]
        print('{:^15} {:^15} {:^15} {:^15}'.format('(0,0,3,52)', oop.kchi, oop.chi0, oop.n))
        oop = parms.out_of_plane[(0,0,24,0)]
        print('{:^15} {:^15} {:^15} {:^15}'.format('(0,0,24,0)', oop.kchi, oop.chi0, oop.n))