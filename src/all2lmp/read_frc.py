# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.5
May 3rd, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

This script will read an entire .frc file and save information into classes for easy access.
Each class is accesible via a key to a dictionary value. The key is a tuple of the atom types
that make up the coeff. The key will only store the exact ordering of atom types in each coeff
so, forwards and reverse checking will need to be preformed when it comes time to fill in each
coeff parameters.

If there is empty spaces due to symmetry in any of the parameter sections the empty space will
be used as a symmetric coeff and will be read in accordingly. If there are duplicates of atom
types for each coeff, the script will dynamically check to make sure the tuple of atom types
belongs to the highest possible version of the read in parameters. This operation is performed
on the for the following coeffs:
    - Equivalences
    - Auto-equivalences
    - Bond-inc
    - Pair coeffs
    - Bond coeffs (quadratic and quartic)
    - Bond coeffs (morse and auto morse)
    - Angle coeffs (quadratic and quartic)
    - Dihedral coeffs (Torsion 1 and Torsion 3)
    - Improper coeffs (quadratic and quartic)
    - Bondbond
    - Bondbond13
    - Bondangle
    - Angleangle
    - Endbondtorsion
    - Middlebondtorsion
    - Angletorsion
    - Angleangletorsion

ForceField specific sections that need to be kept track of during reading in and execution of 
the rest of the code:
    
    Class2 FF's:
        ----------------------
        | IFF/PCFF sections: |
        ---------------------
            #atom_types            cff91
            #equivalence           cff91
            #auto_equivalence      cff91_auto
            
            #nonbond(9-6)          cff91
            #bond_increments       cff91_auto
            
            #quartic_bond          cff91
            #quadratic_bond        cff91_auto
            
            #quartic_angle         cff91
            #quadratic_angle       cff91_auto
                                
            #torsion_3             cff91
            #torsion_1             cff91_auto
            
            #wilson_out_of_plane   cff91
            #wilson_out_of_plane   cff91_auto
            
            #bond-bond             cff91
            #bond-bond_1_3         cff91
            #bond-angle            cff91
            #angle-angle           cff91
            #end_bond-torsion_3    cff91
            #middle_bond-torsion_3 cff91
            #angle-torsion_3       cff91
            #angle-angle-torsion_1 cff91
         
        ---------------------
        | compass sections: |
        ---------------------
            #atom_types            compass
            #equivalence           compass
            
            #nonbond(9-6)          compass
            #bond_increments       compass
            
            #quartic_bond          compass
            
            #quartic_angle         compass
            
            #torsion_3             compass
            
            #wilson_out_of_plane   compass
            
            
            #bond-bond             compass
            #bond-bond_1_3         compass
            #bond-angle            compass
            #end_bond-torsion_3    compass
            #middle_bond-torsion_3 compass
            #angle-torsion_3       compass
            #angle-angle           compass
            #angle-angle-torsion_1 compass
                
    Class1 FF's:
        -------------------------
        | cvff/clayff sections: |
        -------------------------
            #atom_types	           cvff
            #equivalence	       cvff 
            #auto_equivalence	   cvff_auto
            
            #nonbond(12-6)	       cvff 
            #bond_increments       cvff
            
            #morse_bond	           cvff
            #morse_bond            cvff_auto
            
            #quadratic_bond	       cvff
            #quadratic_bond        cvff_auto
            
            #quadratic_angle       cvff
            #quadratic_angle       cvff_auto
            
            #torsion_1	           cvff 
            #torsion_1             cvff_auto
            
            #out_of_plane	       cvff 
            #out_of_plane	       cvff_auto
            
    Class0 FF's:
        --------------------
        | oplsaa sections: |
        --------------------
            #atom_types                 cvff
            #equivalence                cvff
            
            #nonbond(12-6)              cvff
            #bond_increments            cvff
            
            #morse_bond                 cvff_auto
            
            #quadratic_bond             cvff
            #quadratic_bond             cvff_auto
            
            #quadratic_angle            cvff
            #quadratic_angle            cvff_auto
            
            #torsion_1                  opls
            #torsion_1                  cvff_auto
            
            #out_of_plane               cvff
            #out_of_plane               cvff_auto
            #out_of_plane-out_of_plane  cvff
            
            #bond-bond                  cvff
            #bond-angle                 cvff
            #angle-angle-torsion_1      cvff
            #angle-angle                cvff
            
    ClassD FF:
        ----------------------
        | DREIDING sections: |
        ----------------------
            #atom_types	           cvff
            #equivalence	       cvff 
            
            #nonbond(12-6)	       cvff 
            
            #morse_bond	           cvff
            
            #quadratic_bond	       cvff
            
            #quadratic_angle       cvff
            
            #torsion_1	           cvff 
            
            #out_of_plane          DREIDING
"""
class Atom_types:
    pass # .type .mass, .element, .connection, .comment
     
class Equivalences:
    pass # .type, .nonb, .bond, .angle, .torsion, .oop
    
class Auto_equivalences:
    pass # .type, .nonb, .bond_inct, .bond .angle_end, .angle_apex, .torsion_end, .torsion_center, .oop_end, .oop_center
    
class Bond_incs:
    pass # .ij .ji .ver .ref
    
class Morse_bonds:
    pass # .r0, .d .alpha .ver .ref
        
class Quadratic_bonds:
    pass # .r0, .k2 .ver .ref
    
class Quartic_bonds:
    pass # .r0, .k2, .k3, .k4 .ver .ref
    
class Quadratic_angles:
    pass # .theta0, .k2 .ver .ref

class Quartic_angles:
    pass # .theta0, .k2, .k3, .k4 .ver .ref

class Torsion1:
    pass # .kphi, .n, .phi0 .ver .ref

class Torsion3:
    pass # .v1, .phi1, .v2, .phi2, .v3, .phi3 .ver .ref
    
class Wilson_out_of_plane:
    pass # .kchi, .chi0 .ver .ref
    
class DREIDING_OOP:
    pass # .phi0 .kl .ver .ref

class Pair_coeff:
    pass # .r, .eps  .ver .ref
    
class Bondbond:
    pass # .kb_bp

class Bondbond13:
    pass # .kb_bp
    
class Bondangle:
    pass # .kb_theta .kbp_theta  
    
class Angleangle:
    pass # .k_theta_thetap
    
class Endbondtorsion:
    pass # .l_f1, .l_f2, .l_f3, .r_f1, .r_f2, .r_f3
    
class Middlebondtorsion:
    pass # .f1 .f2 .f3  
    
class Angletorsion:
    pass # .l_f1, .l_f2, .l_f3, .r_f1, .r_f2, .r_f3 
    
class Angleangletorsion:
    pass # .k_ang_ang_tor 
    
# Function to check if variable is a float
def check_float(variable):
    try:
        float(variable)
        return_boolean = True
    except:
        return_boolean = False
    return return_boolean

###########################################################################
# Open and read .frc file to store parameters in dictionaries and classes #
###########################################################################
class forcefield_file:
    def __init__(self, frc_file, log):
        print_higher_versions = True # option to print higher versions (True or False)
        
        # set frc file type and log filename
        self.type = 'fix-bond'
        self.filename = frc_file
        
        # atom types and equivs/auto-equivs parameters
        self.atom_types = {}               # class1/class2     {atom type : atom object}
        self.equivalences = {}             # class1/class2     {equivalence : equivalence object}
        self.auto_equivalences = {}        # class1/class2     {auto equivalence : auto equivalence object}
        
        # bond types parameters
        self.bond_increments = {}          # class1/class2     { (type_i, type_j) : delta_ij or delta_ji (order dependant on input)}
        self.quadratic_bonds_auto = {}     # class1 auto-equiv { (type_i, type_j) : quadratic bond object}
        self.morse_bonds = {}              # class1 equivalent { (type_i, type_j) : morse bond object}
        self.morse_bonds_auto = {}         # class1 auto-equiv { (type_i, type_j) : morse bond object}
        self.quadratic_bonds = {}          # class2 auto-equiv { (type_i, type_j) : quadratic bond object}
        self.quartic_bonds = {}            # class2 equivalent { (type_i, type_j) : quartic bond object}

        # angle types parameters
        self.quadratic_angles_auto = {}    # class1 auto-equiv { (type_i, type_j, type_k) : quadratic angle object}
        self.quadratic_angles = {}         # class2 auto-equiv/class1 equivalent{ (type_i, type_j, type_k) : quadratic angle object}
        self.quartic_angles = {}           # class2 equivalent { (type_i, type_j, type_k) : quadratic angle object}
        
        # dihedral types parameters
        self.torsion_1_auto = {}           # class1 auto-equiv { (type_i, type_j, type_k, type_l) : torsion1 object}
        self.torsion_1 = {}                # class2 auto-equiv/class1 equivalent { (type_i, type_j, type_k, type_l) : torsion1 object}
        self.torsion_1_opls = {}           # class0 equivalent { (type_i, type_j, type_k, type_l) : torsion1 object}
        self.torsion_3 = {}                # class2 equivalent { (type_i, type_j, type_k, type_l) : torsion3 object}
        
        # improper types parameters
        self.out_of_plane = {}             # class1 equivalent { (type_i, type_j, type_k, type_l) : oop object}
        self.out_of_plane_auto = {}        # class1 auto-equiv { (type_i, type_j, type_k, type_l) : oop object}
        self.wilson_out_of_plane = {}      # class2 equivalent { (type_i, type_j, type_k, type_l) : oop object}
        self.wilson_out_of_plane_auto = {} # class2 auto-equiv { (type_i, type_j, type_k, type_l) : oop object}
        self.out_of_plane_DREIDING = {}    # classd equivalent { (type_i, type_j, type_k, type_l) : oop object}
        
        # atom types parameters
        self.pair_coeffs_9_6 = {}          # class1 equivalent {atom type : pair coeff object}
        self.pair_coeffs_12_6 = {}         # class2 equivalent {atom type : pair coeff object}
        
        # crossterm types parameters
        self.bondbond = {}                 # class2 equivalent { (type_i, type_j, type_k) : bondbond object}
        self.bondbond13 = {}               # class2 equivalent { (type_i, type_j, type_k, type_l) : bondbond13 object}
        self.bondangle = {}                # class2 equivalent { (type_i, type_j, type_k, type_l) : bondangle object}
        self.angleangle = {}               # class2 equivalent { (type_i, type_j, type_k, type_l) : angleangle object}
        self.endbondtorsion = {}           # class2 equivalent { (type_i, type_j, type_k, type_l) : endbondtorsion object}
        self.middlebondtorsion = {}        # class2 equivalent { (type_i, type_j, type_k, type_l) : middlebondtorsion object}
        self.angletorsion = {}             # class2 equivalent { (type_i, type_j, type_k, type_l) : angletorsion object}
        self.angleangletorsion = {}        # class2 equivalent { (type_i, type_j, type_k, type_l) : angleangletorsion object}
        self.torsiontorsion = {}           # class2 equivalent { (type_i, type_j, type_k, type_l, type_m) : torsiontorsion object}
        
        # Dictionary to log higher versions
        higher_versions_log = {'equivalences':             set([]), 'auto_equivalences':        set([]),
                               'bond_increments':          set([]), 'quadratic_bonds':          set([]),
                               'quadratic_bonds auto':     set([]), 'morse_bonds':              set([]),
                               'morse_bonds auto':         set([]), 'quartic_bonds':            set([]),
                               'quadratic_angles':         set([]), 'quadratic_angles auto':    set([]),
                               'quartic_angles':           set([]),  'torsion_1':                set([]),
                               'torsion_1_auto':           set([]),  'torsion_1_opls':           set([]),
                               'torsion_3':                set([]),  'out_of_plane':             set([]),
                               'out_of_plane auto':        set([]),  'wilson_out_of_plane':      set([]),
                               'wilson_out_of_plane_auto': set([]),  'pair_coeffs_9_6':          set([]),
                               'pair_coeffs_12_6':         set([]),  'bondbond':                 set([]),
                               'bondbond13':               set([]),  'bondangle':                set([]),
                               'angleangle':               set([]),  'endbondtorsion':           set([]),
                               'middlebondtorsion':        set([]),  'angletorsion':             set([]),
                               'angleangletorsion':        set([]),  'torsiontorsion':           set([])}

        ####################################
        # Opening and reading the frc file #
        ####################################
        with open(frc_file, 'r') as f:
            # atom types and equivs/auto-equivs and bondics parameters
            atomtypes_flag = False
            equivalence_flag = False
            autoequivalence_flag = False
            bond_inc_flag = False
            
            # bond types parameters
            quadraticbond_flag = False
            quadraticbond_auto_flag = False
            quarticbond_flag = False
            morsebond_flag = False
            morsebond_auto_flag = False
            
            # angle types parameters
            quadraticangle_flag = False
            quadraticangle_auto_flag = False
            quarticangle_flag = False
            
            # dihedral types parameters
            torsion1_flag = False
            torsion1_auto_flag = False
            torsion1_opls_flag = False
            torsion3_flag = False
            
            # improper types parameters
            wilson1_flag = False
            wilson1_auto_flag = False
            class1_oop_flag = False
            class1_oop_flag_auto = False
            DREINDING_OOP_flag = False
            
            # atom types parameters
            paircoeff_flag = False
            pair_coeff_type_9_6 = False
            pair_coeff_type_12_6 = False
            
            # crossterm types parameters
            bondbond_flag = False
            bondbond13_flag = False
            bondangle_flag = False
            angleangle_flag = False
            endbondtorsion_flag = False
            middlebondtorsion_flag = False
            angletorsion_flag = False
            angleangletorsion_flag = False           
            
            # Looping through each line of the frc file
            for line in f:

                # Split line and strip line
                line = line.strip()
                string = line # For accessing string of each line
                line = line.split()
                
                
                # Setting flags if '#' character is in string of 1st element in list
                # as False to break each section flag from previous iteration
                if len(line) > 0 and line[0][0] == '#':
                    # atom types and equivs/auto-equivs and bondics parameters
                    atomtypes_flag = False
                    equivalence_flag = False
                    autoequivalence_flag = False
                    bond_inc_flag = False
                    
                    # bond types parameters
                    quadraticbond_flag = False
                    quadraticbond_auto_flag = False
                    quarticbond_flag = False
                    morsebond_flag = False
                    morsebond_auto_flag = False
                    
                    # angle types parameters
                    quadraticangle_flag = False
                    quadraticangle_auto_flag = False
                    quarticangle_flag = False
                    
                    # dihedral types parameters
                    torsion1_flag = False
                    torsion1_auto_flag = False
                    torsion1_opls_flag = False
                    torsion3_flag = False
                    
                    # improper types parameters
                    wilson1_flag = False
                    wilson1_auto_flag = False
                    class1_oop_flag = False
                    class1_oop_flag_auto = False
                    DREINDING_OOP_flag = False
                    
                    # atom types parameters
                    paircoeff_flag = False
                    pair_coeff_type_9_6 = False
                    pair_coeff_type_12_6 = False
                    
                    # crossterm types parameters
                    bondbond_flag = False
                    bondbond13_flag = False
                    bondangle_flag = False
                    angleangle_flag = False
                    endbondtorsion_flag = False
                    middlebondtorsion_flag = False
                    angletorsion_flag = False
                    angleangletorsion_flag = False
                
                # atom types and equivs/auto-equivs and bondics parameters
                if '#atom_types' in line:
                    atomtypes_flag = True
                elif '#equivalence' in line:
                    equivalence_flag = True
                elif '#auto_equivalence' in line:
                    autoequivalence_flag = True
                elif '#bond_increments' in line:
                    bond_inc_flag = True
                
                # bond coeffs parameters
                elif '#quadratic_bond' in line:
                    # class2 auto-equiv quadraticbond_flag
                    if '#quadratic_bond' in line and 'cff91_auto' in string:
                        quadraticbond_flag = True
                    # class1 equivalent quadraticangle_flag
                    elif '#quadratic_bond' in line and 'cvff' in string and 'auto' not in string:
                        quadraticbond_flag = True
                    # class1 auto-equiv quadraticangle_auto_flag
                    elif '#quadratic_bond' in line and 'cvff' in string and 'auto' in string:
                        quadraticbond_auto_flag = True
                # class1 morse bond
                elif '#morse_bond' in line:
                    # class1 morse bond
                    if '#morse_bond' in line and 'cvff' in string and 'auto' not in string:
                        morsebond_flag = True
                    # class1 morse bond auto
                    elif '#morse_bond' in line and 'cvff' in string and 'auto' in string:
                        morsebond_auto_flag = True
                # class2 quartic bond
                elif '#quartic_bond' in line:
                    quarticbond_flag = True
                
                # angle coeffs parameters
                elif '#quadratic_angle' in line:
                    # class2 auto-equiv quadraticangle_flag
                    if '#quadratic_angle' in line and 'cff91_auto' in string:
                        quadraticangle_flag = True
                    # class1 equivalent quadraticangle_flag
                    elif '#quadratic_angle' in line and 'cvff' in string and 'auto' not in string:
                        quadraticangle_flag = True
                    # class1 auto-equiv quadraticangle_auto_flag
                    elif '#quadratic_angle' in line and 'cvff' in string and 'auto' in string:
                        quadraticangle_auto_flag = True
                # class2 quartic angle
                elif '#quartic_angle' in line:
                    quarticangle_flag = True                   
                
                # dihedral coeffs parameters
                elif '#torsion_1' in line:
                    # class2 auto-equiv torsion1_flag
                    if '#torsion_1' in line and 'cff91_auto' in string:
                        torsion1_flag = True
                    # class1 equivalent torsion1_flag
                    elif '#torsion_1' in line and 'cvff' in string and 'auto' not in string:
                        torsion1_flag = True
                    # class1 auto-equiv torsion1_flag
                    elif '#torsion_1' in line and 'cvff' in string and 'auto' in string:
                        torsion1_auto_flag = True 
                    # class0 equivalent for oplsa
                    elif '#torsion_1' in line and 'opls' in string:
                        torsion1_opls_flag = True
                # class2 quartic dihedral
                elif '#torsion_3' in line:
                    torsion3_flag = True
                
                # improper coeffs parameters
                elif '#wilson_out_of_plane' in line and 'auto' not in string:
                    wilson1_flag = True
                elif '#wilson_out_of_plane' in line and 'auto' in string:
                    wilson1_auto_flag = True
                elif '#out_of_plane' in line and 'auto' not in string and 'DREIDING' not in string:
                    class1_oop_flag = True
                elif '#out_of_plane' in line  and 'auto' in string and 'DREIDING' not in string:
                    class1_oop_flag_auto = True
                elif '#out_of_plane' in string and 'DREIDING':
                    DREINDING_OOP_flag = True
                
                # pair coeffs parameters
                elif '#nonbond(9-6)' in line or '#nonbond(12-6)' in line:
                    paircoeff_flag = True
                    if '#nonbond(9-6)' in line:
                        pair_coeff_type_9_6 = True
                    if '#nonbond(12-6)' in line:
                        pair_coeff_type_12_6 = True
                
                # crossterm coeffs parameters
                elif '#bond-bond' in line and '#bond-bond_1_3' not in line:
                    bondbond_flag = True
                elif '#bond-bond_1_3' in line:# and '#bond-bond' not in line:
                    bondbond13_flag = True
                elif '#bond-angle' in line: 
                    bondangle_flag = True
                elif '#angle-angle' in line and '#angle-angle-torsion_1' not in line: 
                    angleangle_flag = True
                elif '#end_bond-torsion_3' in line: 
                    endbondtorsion_flag = True
                elif '#middle_bond-torsion_3' in line: 
                    middlebondtorsion_flag = True
                elif '#angle-torsion_3' in line: 
                    angletorsion_flag = True
                elif '#angle-angle-torsion_1' in line: 
                    angleangletorsion_flag = True
                    
                # Finding atom types parameters    
                if atomtypes_flag and len(line) >= 5 and check_float(line[0]):
                    a = Atom_types()
                    ver = float(line[0]); ref = float(line[1]); i = line[2];
                    a.type = i
                    a.mass = float(line[3])
                    a.element = line[4]
                    try: # Some .frc files have connections column and others do not
                        a.connection = int(line[5])
                        comments = ' '.join([line[i] for i in range(6, len(line), 1)])
                    except:
                        a.connection = 'n.u.' # n.u. stands for not used
                        comments = comments = ' '.join([line[i] for i in range(5, len(line), 1)])
                    a.comment = comments
                    a.ver = ver
                    a.ref = ref
                    self.atom_types[i] = a
                
                # Find equivalences parameters
                elif equivalence_flag and len(line) >= 8 and check_float(line[0]):
                    e = Equivalences()
                    ver = float(line[0]);ref = float(line[1]); i = line[2];
                    e.type = i
                    e.nonb = line[3]
                    e.bond = line[4]
                    e.angle = line[5]
                    e.torsion = line[6]
                    e.oop = line[7]
                    e.ver = ver
                    e.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if i not in self.equivalences: self.equivalences[i] = e
                    elif i in self.equivalences:
                        logged_version = self.equivalences[i].ver
                        higher_versions_log['equivalences'].add(i)
                        if ver >= logged_version: self.equivalences[i] = e
                
                # Find autoequivalences parameters  
                elif autoequivalence_flag and len(line) >= 12 and check_float(line[0]):
                    ae = Auto_equivalences()
                    ver = float(line[0]); ref = float(line[1]); i = line[2];
                    ae.type = i
                    ae.nonb = line[3]
                    ae.bond_inct = line[4]
                    ae.bond = line[5]
                    ae.angle_end = line[6]
                    ae.angle_apex = line[7]
                    ae.torsion_end = line[8]
                    ae.torsion_center = line[9]
                    ae.oop_end = line[10]
                    ae.oop_center = line[11]
                    ae.ver = ver
                    ae.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if i not in self.auto_equivalences: self.auto_equivalences[i] = ae
                    elif i in self.auto_equivalences:
                        logged_version = self.auto_equivalences[i].ver
                        higher_versions_log['auto_equivalences'].add(i)
                        if ver >= logged_version: self.auto_equivalences[i] = ae
                
                # Find bond increments parameters   
                elif bond_inc_flag and len(line) >= 6 and check_float(line[0]):
                    bi = Bond_incs()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3];
                    bi.ij = float(line[4])
                    bi.ji = float(line[5])
                    bi.ver = ver
                    bi.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j) not in self.bond_increments: self.bond_increments[(i,j)] = bi
                    elif (i,j) in self.bond_increments:
                        logged_version = self.bond_increments[(i,j)].ver
                        higher_versions_log['bond_increments'].add(tuple([i,j]))
                        if ver >= logged_version: self.bond_increments[(i,j)] = bi
                    
                # Find quadratic bond parameters
                elif quadraticbond_flag and len(line) >= 6 and check_float(line[0]):
                    qb = Quadratic_bonds()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3];
                    qb.r0 = float(line[4])
                    qb.k2 = float(line[5])
                    qb.ver = ver
                    qb.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j) not in self.quadratic_bonds: self.quadratic_bonds[(i,j)] = qb
                    elif (i,j) in self.quadratic_bonds:
                        logged_version = self.quadratic_bonds[(i,j)].ver
                        higher_versions_log['quadratic_bonds'].add(tuple([i,j]))
                        if ver >= logged_version: self.quadratic_bonds[(i,j)] = qb
                               
                # Find quadratic bond auto parameters
                elif quadraticbond_auto_flag and len(line) >= 6 and check_float(line[0]):
                    qb = Quadratic_bonds()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3];
                    qb.r0 = float(line[4])
                    qb.k2 = float(line[5])
                    qb.ver = ver
                    qb.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j) not in self.quadratic_bonds_auto: self.quadratic_bonds_auto[(i,j)] = qb
                    elif (i,j) in self.quadratic_bonds_auto:
                        logged_version = self.quadratic_bonds_auto[(i,j)].ver
                        higher_versions_log['quadratic_bonds auto'].add(tuple([i,j]))
                        if ver >= logged_version: self.quadratic_bonds[(i,j)] = qb
                            
                # Find morse bond parameters
                elif morsebond_flag and len(line) >= 7 and check_float(line[0]):
                    mb = Morse_bonds()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3];
                    mb.r0 = float(line[4])
                    mb.d = float(line[5])
                    mb.alpha = float(line[6])
                    mb.ver = ver
                    mb.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j) not in self.morse_bonds: self.morse_bonds[(i,j)] = mb
                    elif (i,j) in self.morse_bonds:
                        logged_version = self.morse_bonds[(i,j)].ver
                        higher_versions_log['morse_bonds'].add(tuple([i,j]))
                        if ver >= logged_version: self.morse_bonds[(i,j)] = mb
                            
                # Find morse bond parameters
                elif morsebond_auto_flag and len(line) >= 7 and check_float(line[0]):
                    mb = Morse_bonds()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3];
                    mb.r0 = float(line[4])
                    mb.d = float(line[5])
                    mb.alpha = float(line[6])
                    mb.ver = ver
                    mb.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j) not in self.morse_bonds_auto: self.morse_bonds_auto[(i,j)] = mb
                    elif (i,j) in self.morse_bonds_auto:
                        logged_version = self.morse_bonds_auto[(i,j)].ver
                        higher_versions_log['morse_bonds auto'].add(tuple([i,j]))
                        if ver >= logged_version: self.morse_bonds_auto[(i,j)] = mb
                            
                # Find quartic bond parameters
                elif quarticbond_flag and len(line) >= 8 and check_float(line[0]):
                    qb = Quartic_bonds()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3];
                    qb.r0 = float(line[4])
                    qb.k2 = float(line[5])
                    qb.k3 = float(line[6])
                    qb.k4 = float(line[7])
                    qb.ver = ver
                    qb.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j) not in self.quartic_bonds: self.quartic_bonds[(i,j)] = qb
                    elif (i,j) in self.quartic_bonds:
                        logged_version = self.quartic_bonds[(i,j)].ver
                        higher_versions_log['quartic_bonds'].add(tuple([i,j]))
                        if ver >= logged_version: self.quartic_bonds[(i,j)] = qb                        
                            
                # Find quadratic angle parameters
                elif quadraticangle_flag and len(line) >= 7 and check_float(line[0]):
                    qa = Quadratic_angles()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4];
                    qa.theta0 = float(line[5])
                    qa.k2 = float(line[6])
                    qa.ver = ver
                    qa.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k) not in self.quadratic_angles: self.quadratic_angles[(i,j,k)] = qa
                    elif (i,j,k) in self.quadratic_angles:
                        logged_version = self.quadratic_angles[(i,j,k)].ver
                        higher_versions_log['quadratic_angles'].add(tuple([i,j,k]))
                        if ver >= logged_version: self.quadratic_angles[(i,j,k)] = qa
                            
                # Find quadratic angle parameters
                elif quadraticangle_auto_flag and len(line) >= 7 and check_float(line[0]):
                    qa = Quadratic_angles()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4];
                    qa.theta0 = float(line[5])
                    qa.k2 = float(line[6])
                    qa.ver = ver
                    qa.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k) not in self.quadratic_angles_auto: self.quadratic_angles_auto[(i,j,k)] = qa
                    elif (i,j,k) in self.quadratic_angles_auto:
                        logged_version = self.quadratic_angles_auto[(i,j,k)].ver
                        higher_versions_log['quadratic_angles auto'].add(tuple([i,j,k]))
                        if ver >= logged_version: self.quadratic_angles_auto[(i,j,k)] = qa
                    
                # Find quartic angle parameters   
                elif quarticangle_flag and len(line) >= 9 and check_float(line[0]):
                    qa = Quartic_angles()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4];
                    qa.theta0 = float(line[5])
                    qa.k2 = float(line[6])
                    qa.k3 = float(line[7])
                    qa.k4 = float(line[8])
                    qa.ver = ver
                    qa.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k) not in self.quartic_angles: self.quartic_angles[(i,j,k)] = qa
                    elif (i,j,k) in self.quartic_angles:
                        logged_version = self.quartic_angles[(i,j,k)].ver
                        higher_versions_log['quartic_angles'].add(tuple([i,j,k]))
                        if ver >= logged_version: self.quartic_angles[(i,j,k)] = qa
                    
                # Find torsion 1 parameters   
                elif torsion1_flag and len(line) >= 9 and check_float(line[0]):
                    t1 = Torsion1()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    t1.kphi = float(line[6])
                    t1.n = int(line[7])
                    t1.phi0 = float(line[8])
                    t1.ver = ver
                    t1.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.torsion_1: self.torsion_1[(i,j,k,l)] = t1
                    elif (i,j,k,l) in self.torsion_1:
                        logged_version = self.torsion_1[(i,j,k,l)].ver
                        higher_versions_log['torsion_1'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.torsion_1[(i,j,k,l)] = t1
                            
                # Find torsion 1 auto parameters
                elif torsion1_auto_flag and len(line) >= 9 and check_float(line[0]):
                    t1 = Torsion1()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    t1.kphi = float(line[6])
                    t1.n = int(line[7])
                    t1.phi0 = float(line[8])
                    t1.ver = ver
                    t1.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.torsion_1_auto: self.torsion_1_auto[(i,j,k,l)] = t1
                    elif (i,j,k,l) in self.torsion_1_auto:
                        logged_version = self.torsion_1_auto[(i,j,k,l)].ver
                        higher_versions_log['torsion_1_auto'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.torsion_1_auto[(i,j,k,l)] = t1
                                
                # Find torsion 1 opls parameters
                elif torsion1_opls_flag and len(line) >= 9 and check_float(line[0]):
                    t1 = Torsion1()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    t1.k1 = float(line[6])
                    t1.k2 = float(line[7])
                    t1.k3 = float(line[8])
                    t1.k4 = float(line[9])
                    t1.ver = ver
                    t1.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.torsion_1_opls: self.torsion_1_opls[(i,j,k,l)] = t1
                    elif (i,j,k,l) in self.torsion_1_opls:
                        logged_version = self.torsion_1_opls[(i,j,k,l)].ver
                        higher_versions_log['torsion_1_opls'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.torsion_1_opls[(i,j,k,l)] = t1
                    
                # Find torsion 3 parameters   
                elif torsion3_flag and len(line) >= 12 and check_float(line[0]):
                    t3 = Torsion3()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    t3.v1 = float(line[6])
                    t3.phi1 = float(line[7])
                    t3.v2 = float(line[8])
                    t3.phi2 = float(line[9])
                    t3.v3 = float(line[10])
                    t3.phi3 = float(line[11])
                    t3.ver = ver
                    t3.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.torsion_3: self.torsion_3[(i,j,k,l)] = t3
                    elif (i,j,k,l) in self.torsion_3:
                        logged_version = self.torsion_3[(i,j,k,l)].ver
                        higher_versions_log['torsion_3'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.torsion_3[(i,j,k,l)] = t3
                    
                # Find wilson out of plane parameters    
                elif wilson1_flag and len(line) >= 8 and check_float(line[0]):
                    w = Wilson_out_of_plane()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    w.kchi = float(line[6])
                    w.chi0 = float(line[7])
                    w.ver = ver
                    w.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.wilson_out_of_plane: self.wilson_out_of_plane[(i,j,k,l)] = w
                    elif (i,j,k,l) in self.wilson_out_of_plane:
                        logged_version = self.wilson_out_of_plane[(i,j,k,l)].ver
                        higher_versions_log['wilson_out_of_plane'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.wilson_out_of_plane[(i,j,k,l)] = w
                
                # Find wilson out of plane auto parameters   
                elif wilson1_auto_flag and len(line) >= 8 and check_float(line[0]):
                    w = Wilson_out_of_plane()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    w.kchi = float(line[6])
                    w.chi0 = float(line[7])
                    w.ver = ver
                    w.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.wilson_out_of_plane_auto: self.wilson_out_of_plane_auto[(i,j,k,l)] = w
                    elif (i,j,k,l) in self.wilson_out_of_plane_auto:
                        logged_version = self.wilson_out_of_plane_auto[(i,j,k,l)].ver
                        higher_versions_log['wilson_out_of_plane_auto'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.wilson_out_of_plane_auto[(i,j,k,l)] = w
                            
                # Find class1 oop parameters
                elif class1_oop_flag and len(line) >= 9 and check_float(line[0]):
                    w = Wilson_out_of_plane()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    w.kchi = float(line[6])
                    w.n = int(line[7])
                    w.chi0 = float(line[8])
                    w.ver = ver
                    w.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.out_of_plane: self.out_of_plane[(i,j,k,l)] = w
                    elif (i,j,k,l) in self.out_of_plane:
                        logged_version = self.out_of_plane[(i,j,k,l)].ver
                        higher_versions_log['out_of_plane'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.out_of_plane[(i,j,k,l)] = w
                                
                # Find class1 oop auto parameters
                elif class1_oop_flag_auto and len(line) >= 9 and check_float(line[0]):
                    w = Wilson_out_of_plane()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    w.kchi = float(line[6])
                    w.n = int(line[7])
                    w.chi0 = float(line[8])
                    w.ver = ver
                    w.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.out_of_plane_auto: self.out_of_plane_auto[(i,j,k,l)] = w
                    elif (i,j,k,l) in self.out_of_plane_auto:
                        logged_version = self.out_of_plane_auto[(i,j,k,l)].ver
                        higher_versions_log['out_of_plane auto'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.out_of_plane_auto[(i,j,k,l)] = w
                                
                # Find DREIDING oop parameters
                elif DREINDING_OOP_flag and len(line) >= 8 and check_float(line[6]):
                    DOOP = DREIDING_OOP()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    DOOP.kl = float(line[6])
                    DOOP.phi0 = float(line[7])
                    DOOP.ver = ver
                    DOOP.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.out_of_plane_DREIDING: self.out_of_plane_DREIDING[(i,j,k,l)] = DOOP
                    elif (i,j,k,l) in self.out_of_plane_DREIDING:
                        logged_version = self.out_of_plane_DREIDING[(i,j,k,l)].ver
                        #higher_versions_log['out_of_plane'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.out_of_plane_DREIDING[(i,j,k,l)] = DOOP
                    
                # Find pair coeff parameters   
                elif paircoeff_flag and len(line) >= 5 and check_float(line[0]):
                    # Find 9-6
                    if pair_coeff_type_9_6:
                        p = Pair_coeff()
                        ver = float(line[0]); ref = float(line[1]); i = line[2];
                        p.r = float(line[3])
                        p.eps = float(line[4])
                        p.ver = ver
                        p.ref = ref
                        
                        # If parameter not in dictionary add parameter to dictionary; elif log highest version
                        if i not in self.pair_coeffs_9_6: self.pair_coeffs_9_6[i] = p
                        elif i in self.pair_coeffs_9_6:
                            logged_version = self.pair_coeffs_9_6[i].ver
                            higher_versions_log['pair_coeffs_9_6'].add(i)
                            if ver >= logged_version: self.pair_coeffs_9_6[i] = p
                                
                    # Find 12-6
                    if pair_coeff_type_12_6:
                        p = Pair_coeff()
                        ver = float(line[0]); ref = float(line[1]); i = line[2];
                        p.A = float(line[3])
                        p.B = float(line[4])
                        p.ver = ver
                        p.ref = ref
                        
                        # If parameter not in dictionary add parameter to dictionary; elif log highest version
                        if i not in self.pair_coeffs_12_6: self.pair_coeffs_12_6[i] = p
                        elif i in self.pair_coeffs_12_6:
                            logged_version = self.pair_coeffs_12_6[i].ver
                            higher_versions_log['pair_coeffs_12_6'].add(i)
                            if ver >= logged_version: self.pair_coeffs_12_6[i] = p
                    
                # Find bondbond parameters    
                elif bondbond_flag and len(line) >= 6 and check_float(line[0]):
                    b = Bondbond()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4];
                    b.kb_bp = float(line[5])
                    b.ver = ver
                    b.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k) not in self.bondbond: self.bondbond[(i,j,k)] = b
                    elif (i,j,k) in self.bondbond:
                        logged_version = self.bondbond[(i,j,k)].ver
                        higher_versions_log['bondbond'].add(tuple([i,j,k]))
                        if ver >= logged_version: self.bondbond[(i,j,k)] = b
                    
                # Find bondbond13 parameters    
                elif bondbond13_flag and len(line) >= 7 and check_float(line[0]):
                    b = Bondbond13()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    b.kb_bp = float(line[6])
                    b.ver = ver
                    b.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.bondbond13: self.bondbond13[(i,j,k,l)] = b
                    elif (i,j,k,l) in self.bondbond13:
                        logged_version = self.bondbond13[(i,j,k,l)].ver
                        higher_versions_log['bondbond13'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.bondbond13[(i,j,k,l)] = b
                    
                # Find bondangle parameters   
                elif bondangle_flag and len(line) >= 6 and check_float(line[0]):
                    ba = Bondangle()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4];
                    ba.kb_theta = float(line[5])
                    ba.ver = ver
                    ba.ref = ref
                    try: ba.kbp_theta = float(line[6])
                    except: ba.kbp_theta = float(line[5])
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k) not in self.bondangle: self.bondangle[(i,j,k)] = ba
                    elif (i,j,k,l) in self.bondangle:
                        logged_version = self.bondangle[(i,j,k)].ver
                        higher_versions_log['bondangle'].add(tuple([i,j,k]))
                        if ver >= logged_version: self.bondangle[(i,j,k)] = ba
                    
                # Find angleangle parameters  
                elif angleangle_flag and len(line) >= 7 and check_float(line[0]):
                    aa = Angleangle()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    aa.k_theta_thetap = float(line[6])
                    aa.ver = ver
                    aa.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.angleangle: self.angleangle[(i,j,k,l)] = aa
                    elif (i,j,k,l) in self.angleangle:
                        logged_version = self.angleangle[(i,j,k,l)].ver
                        higher_versions_log['angleangle'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.angleangle[(i,j,k,l)] = aa
                    
                # Find end bond torsion parameters   
                elif endbondtorsion_flag and len(line) >= 9 and check_float(line[0]):
                    ebt = Endbondtorsion()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    ebt.l_f1 = float(line[6])
                    ebt.l_f2 = float(line[7])
                    ebt.l_f3 = float(line[8])
                    ebt.ver = ver
                    ebt.ref = ref
                    try:
                        ebt.r_f1 = float(line[9])
                        ebt.r_f2 = float(line[10])
                        ebt.r_f3 = float(line[11])
                    except:
                        ebt.r_f1 = float(line[6])
                        ebt.r_f2 = float(line[7])
                        ebt.r_f3 = float(line[8])
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.endbondtorsion: self.endbondtorsion[(i,j,k,l)] = ebt
                    elif (i,j,k,l) in self.endbondtorsion:
                        logged_version = self.endbondtorsion[(i,j,k,l)].ver
                        higher_versions_log['endbondtorsion'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.endbondtorsion[(i,j,k,l)] = ebt
                    
                # Find middle torsion parameters   
                elif middlebondtorsion_flag and len(line) >= 9 and check_float(line[0]):
                    mbt = Middlebondtorsion()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    mbt.f1 = float(line[6])
                    mbt.f2 = float(line[7])
                    mbt.f3 = float(line[8])
                    mbt.ver = ver
                    mbt.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.middlebondtorsion: self.middlebondtorsion[(i,j,k,l)] = mbt
                    elif (i,j,k,l) in self.middlebondtorsion:
                        logged_version = self.middlebondtorsion[(i,j,k,l)].ver
                        higher_versions_log['middlebondtorsion'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.middlebondtorsion[(i,j,k,l)] = mbt
                    
                # Find angle torsion parameters   
                elif angletorsion_flag and len(line) >= 9 and check_float(line[0]):
                    at = Angletorsion()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    at.l_f1 = float(line[6])
                    at.l_f2 = float(line[7])
                    at.l_f3 = float(line[8])
                    at.ver = ver
                    at.ref = ref
                    try:
                        at.r_f1 = float(line[9])
                        at.r_f2 = float(line[10])
                        at.r_f3 = float(line[11])
                    except:
                        at.r_f1 = float(line[6])
                        at.r_f2 = float(line[7])
                        at.r_f3 = float(line[8])
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.angletorsion: self.angletorsion[(i,j,k,l)] = at
                    elif (i,j,k,l) in self.angletorsion:
                        logged_version = self.angletorsion[(i,j,k,l)].ver
                        higher_versions_log['angletorsion'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.angletorsion[(i,j,k,l)] = at
                
                # Find angle angle torsion parameters  
                elif angleangletorsion_flag and len(line) >= 7 and check_float(line[0]):
                    aat = Angleangletorsion()
                    ver = float(line[0]); ref = float(line[1]);
                    i = line[2]; j = line[3]; k = line[4]; l = line[5];
                    aat.k_ang_ang_tor = float(line[6])
                    aat.ver = ver
                    aat.ref = ref
                    
                    # If parameter not in dictionary add parameter to dictionary; elif log highest version
                    if (i,j,k,l) not in self.angleangletorsion: self.angleangletorsion[(i,j,k,l)] = aat
                    elif (i,j,k,l) in self.angleangletorsion:
                        logged_version = self.angleangletorsion[(i,j,k,l)].ver
                        higher_versions_log['angleangletorsion'].add(tuple([i,j,k,l]))
                        if ver >= logged_version: self.angleangletorsion[(i,j,k,l)] = aat

        # print out higher logged versions
        if print_higher_versions:
            # Map to get dictionary that is associated with each key in higher_versions_log
            coeff_map = {'pair_coeffs_9_6':self.pair_coeffs_9_6, 'pair_coeffs_12_6':self.pair_coeffs_12_6, 'equivalences':self.equivalences,
                         'auto_equivalences':self.auto_equivalences, 'bond_increments':self.bond_increments, 'quadratic_bonds':self.quadratic_bonds,
                         'quadratic_bonds auto':self.quadratic_bonds_auto, 'morse_bonds':self.morse_bonds, 'morse_bonds auto':self.morse_bonds_auto,
                         'quartic_bonds':self.quartic_bonds, 'quadratic_angles':self.quadratic_angles, 'quadratic_angles auto':self.quadratic_angles_auto,
                         'quartic_angles':self.quartic_angles, 'torsion_1':self.torsion_1, 'torsion_1_auto':self.torsion_1_auto, 'torsion_3':self.torsion_3,
                         'wilson_out_of_plane':self.wilson_out_of_plane, 'wilson_out_of_plane_auto':self.wilson_out_of_plane_auto,
                         'bondbond':self.bondbond, 'bondbond13':self.bondbond13, 'bondangle':self.bondangle, 'angleangle':self.angleangle,
                         'endbondtorsion':self.endbondtorsion, 'middlebondtorsion':self.middlebondtorsion, 'angletorsion':self.angletorsion,
                         'angleangletorsion':self.angleangletorsion, 'torsiontorsion':self.torsiontorsion, 'torsion_1_opls':self.torsion_1_opls,
                         'out_of_plane': self.out_of_plane, 'out_of_plane auto': self.out_of_plane_auto}
            
            # Find any higher versions that have been logged and print
            for i in higher_versions_log:
                for j in higher_versions_log[i]:
                    log.out('Using higher version of parameters for #{} {} version: {} '.format(i, ' '.join(j), '{:.2f}'.format(coeff_map[i][j].ver))) 
            log.out('')