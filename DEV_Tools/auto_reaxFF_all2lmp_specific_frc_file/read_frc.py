# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.4
December 19th, 2022
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
it logs belongs to the higher possible version of the read in parameters. This operation is
done for the following coeffs:
    - Equivalences
    - Auto-equivalences
    - Bond-inc
    - Pair coeffs
    - Bond coeffs (quadratic and quartic)
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

Version 1.0 would just take the last of the duplicates and store the parameter data for the
coeffs associated with it. It was determined by comparison to msi2lmp scripts that logging
the higher version was the better approach when determining which parameters to use. Thus 1.1
addresses this issue.

Version 1.2 adds more flexible header spacing and white space between different .frc files,
making it more flexible to nuances from different .frc file spacing.



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
    
class Torsiontorsion:
    pass # Current IFF .frc file does not have this section filled out
    
    
# Function to check if variable is a float
def check_float(variable):
    try:
        float(variable)
        return_boolean = True
    except:
        return_boolean = False
    return return_boolean

# Open and read .frc file and store info in dictionaries and classes
class forcefield_file:
    def __init__(self, frc_file):
        # atom types info/bondics
        self.atom_types = {}                # class1/class2     {atom type : atom object}
        self.equivalences = {}              # class1/class2     {equivalence : equivalence object}
        self.auto_equivalences = {}         # class1/class2     {auto equivalence : auto equivalence object}
        
        # bond types info
        self.bond_increments = {}           # class1/class2     { (type_i, type_j) : delta_ij or delta_ji (order dependant on input)}
        self.quadratic_bonds_auto = {}      # class1 auto-equiv { (type_i, type_j) : quadratic bond object}
        self.morse_bonds = {}                # class1 equivalent { (type_i, type_j) : morse bond object}
        self.morse_bonds_auto = {}           # class1 auto-equiv { (type_i, type_j) : morse bond object}
        self.quadratic_bonds = {}           # class2 auto-equiv { (type_i, type_j) : quadratic bond object}
        self.quartic_bonds = {}             # class2 equivalent { (type_i, type_j) : quartic bond object}

        # angle types info
        self.quadratic_angles_auto = {}     # class1 auto-equiv { (type_i, type_j, type_k) : quadratic angle object}
        self.quadratic_angles = {}          # class2 auto-equiv/class1 equivalent{ (type_i, type_j, type_k) : quadratic angle object}
        self.quartic_angles = {}            # class2 equivalent { (type_i, type_j, type_k) : quadratic angle object}
        
        # dihedral types info
        self.torsion_1_auto = {}            # class1 auto-equiv { (type_i, type_j, type_k, type_l) : torsion1 object}
        self.torsion_1 = {}                 # class2 auto-equiv/class1 equivalent { (type_i, type_j, type_k, type_l) : torsion1 object}
        self.torsion_3 = {}                 # class2 equivalent { (type_i, type_j, type_k, type_l) : torsion3 object}
        
        # improper types info
        self.out_of_plane = {}              # class1 equivalent { (type_i, type_j, type_k, type_l) : oop object}
        self.out_of_plane_auto = {}         # class1 auto-equiv { (type_i, type_j, type_k, type_l) : oop object}
        self.wilson_out_of_plane = {}       # class2 equivalent { (type_i, type_j, type_k, type_l) : oop object}
        self.wilson_out_of_plane_auto = {}  # class2 auto-equiv { (type_i, type_j, type_k, type_l) : oop object}
        
        # atom types info
        self.pair_coeffs_9_6 = {}           # class1 equivalent {atom type : pair coeff object}
        self.pair_coeffs_12_6 = {}          # class2 equivalent {atom type : pair coeff object}
        
        # crossterm types info
        self.bondbond = {}                  # class2 equivalent { (type_i, type_j, type_k) : bondbond object}
        self.bondbond13 = {}                # class2 equivalent { (type_i, type_j, type_k, type_l) : bondbond13 object}
        self.bondangle = {}                 # class2 equivalent { (type_i, type_j, type_k, type_l) : bondangle object}
        self.angleangle = {}                # class2 equivalent { (type_i, type_j, type_k, type_l) : angleangle object}
        self.endbondtorsion = {}            # class2 equivalent { (type_i, type_j, type_k, type_l) : endbondtorsion object}
        self.middlebondtorsion = {}         # class2 equivalent { (type_i, type_j, type_k, type_l) : middlebondtorsion object}
        self.angletorsion = {}              # class2 equivalent { (type_i, type_j, type_k, type_l) : angletorsion object}
        self.angleangletorsion = {}         # class2 equivalent { (type_i, type_j, type_k, type_l) : angleangletorsion object}
        self.torsiontorsion = {}            # class2 equivalent { (type_i, type_j, type_k, type_l, type_m) : torsiontorsion object}
        
        # Dictionary to log higher versions
        higher_versions_log = {'equivalences':             set([]),
                               'auto_equivalences':        set([]),
                               'bond_increments':          set([]),
                               'quadratic_bonds':          set([]),
                               'quadratic_bonds auto':     set([]),
                               'morse_bonds':              set([]),
                               'morse_bonds auto':         set([]),
                               'quartic_bonds':            set([]),
                               'quadratic_angles':         set([]),
                               'quadratic_angles auto':    set([]),
                               'quartic_angles':           set([]),
                               'torsion_1':                set([]),
                               'torsion_1_auto':           set([]),
                               'torsion_3':                set([]),
                               'out_of_plane':             set([]),
                               'out_of_plane auto':        set([]),
                               'wilson_out_of_plane':      set([]),
                               'wilson_out_of_plane_auto': set([]),
                               'pair_coeffs_9_6':          set([]),
                               'pair_coeffs_12_6':         set([]),
                               'bondbond':                 set([]),
                               'bondbond13':               set([]),
                               'bondangle':                set([]),
                               'angleangle':               set([]),
                               'endbondtorsion':           set([]),
                               'middlebondtorsion':        set([]),
                               'angletorsion':             set([]),
                               'angleangletorsion':        set([]),
                               'torsiontorsion':           set([])}
  

        # Opening and reading the frc file
        with open(frc_file, 'r') as f:
            
            # Initializing flags    
            # atom types info/bondics
            atomtypes_flag = False
            equivalence_flag = False
            autoequivalence_flag = False
            bond_inc_flag = False
            
            # bond types info
            quadraticbond_flag = False
            quadraticbond_auto_flag = False
            quarticbond_flag = False
            morsebond_flag = False
            morsebond_auto_flag = False
            
            # angle types info
            quadraticangle_flag = False
            quadraticangle_auto_flag = False
            quarticangle_flag = False
            
            # dihedral types info
            torsion1_flag = False
            torsion1_auto_flag = False
            torsion3_flag = False
            
            # improper types info
            wilson1_flag = False
            wilson1_auto_flag = False
            class1_oop_flag = False
            class1_oop_flag_auto = False
            
            # atom types info
            paircoeff_flag = False
            pair_coeff_type_9_6 = False
            pair_coeff_type_12_6 = False
            
            # crossterm types info
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
                    # atom types info/bondics
                    atomtypes_flag = False
                    equivalence_flag = False
                    autoequivalence_flag = False
                    bond_inc_flag = False
                    
                    # bond types info
                    quadraticbond_flag = False
                    quadraticbond_auto_flag = False
                    quarticbond_flag = False
                    morsebond_flag = False
                    morsebond_auto_flag = False
                    
                    # angle types info
                    quadraticangle_flag = False
                    quadraticangle_auto_flag = False
                    quarticangle_flag = False
                    
                    # dihedral types info
                    torsion1_flag = False
                    torsion1_auto_flag = False
                    torsion3_flag = False
                    
                    # improper types info
                    wilson1_flag = False
                    wilson1_auto_flag = False
                    class1_oop_flag = False
                    class1_oop_flag_auto = False
                    
                    # atom types info
                    paircoeff_flag = False
                    pair_coeff_type_9_6 = False
                    pair_coeff_type_12_6 = False
                    
                    # crossterm types info
                    bondbond_flag = False
                    bondbond13_flag = False
                    bondangle_flag = False
                    angleangle_flag = False
                    endbondtorsion_flag = False
                    middlebondtorsion_flag = False
                    angletorsion_flag = False
                    angleangletorsion_flag = False
                
                # atom types info/bondics
                if '#atom_types' in line:
                    atomtypes_flag = True
                elif '#equivalence' in line:
                    equivalence_flag = True
                elif '#auto_equivalence' in line:
                    autoequivalence_flag = True
                elif '#bond_increments' in line:
                    bond_inc_flag = True
                
                # bond types info
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
                
                # angle types info
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
                
                # dihedral types info
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
                # class2 quartic dihedral
                elif '#torsion_3' in line:
                    torsion3_flag = True
                
                # improper types info
                elif '#wilson_out_of_plane' in line and 'auto' not in string:
                    wilson1_flag = True
                elif '#wilson_out_of_plane' in line and 'auto' in string:
                    wilson1_auto_flag = True
                elif '#out_of_plane' in line and 'auto' not in string:
                    class1_oop_flag = True
                elif '#out_of_plane' in line  and 'auto' in string:
                    class1_oop_flag_auto = True
                
                # atom types info/bondics
                elif '#nonbond(9-6)' in line or '#nonbond(12-6)' in line:
                    paircoeff_flag = True
                    if '#nonbond(9-6)' in line:
                        pair_coeff_type_9_6 = True
                    if '#nonbond(12-6)' in line:
                        pair_coeff_type_12_6 = True
                
                # crossterm types info
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

                  
                    
                # Finding atom types information    
                if atomtypes_flag:
                    # Check if len(line) >= 5 and if line[0] is a float
                    if len(line) >= 5 and check_float(line[0]):
                        #print(line)
                        a = Atom_types()
                        ver = float(line[0])
                        ref = float(line[1])
                        type = line[2]
                        a.type = type
                        a.mass = float(line[3])
                        a.element = line[4]
                        # Some .frc files have connections column and others do not
                        try:
                            a.connection = int(line[5])
                            comments = ''
                            for i in range(6, len(line), 1):
                                comments += ' ' + line[i]
                        except:
                            a.connection = 'n.u.' # n.u. stands for not used
                            comments = ''
                            for i in range(5, len(line), 1):
                                comments += ' ' + line[i]
                        a.comment = comments
                        a.ver = ver
                        a.ref = ref
                        self.atom_types[type] = a
                
                
                # Find equivalences information 
                elif equivalence_flag:
                    # Check if len(line) >= 8 and if line[0] is a float
                    if len(line) >= 8 and check_float(line[0]):
                        #print(line)
                        e = Equivalences()
                        ver = float(line[0])
                        ref = float(line[1])
                        type = line[2]
                        e.type = type
                        e.nonb = line[3]
                        e.bond = line[4]
                        e.angle = line[5]
                        e.torsion = line[6]
                        e.oop = line[7]
                        e.ver = ver
                        e.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if type not in self.equivalences:
                            self.equivalences[type] = e
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif type in self.equivalences:
                            logged_version = self.equivalences[type].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['equivalences'].add(type)
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.equivalences[type] = e
                
                
                # Find autoequivalences information    
                elif autoequivalence_flag:
                    # Check if len(line) >= 12 and if line[0] is a float
                    if len(line) >= 12 and check_float(line[0]):
                        #print(line)
                        ae = Auto_equivalences()
                        ver = float(line[0])
                        ref = float(line[1])
                        type = line[2]
                        ae.type = type
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
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if type not in self.auto_equivalences:
                            self.auto_equivalences[type] = ae
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif type in self.auto_equivalences:
                            logged_version = self.auto_equivalences[type].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['auto_equivalences'].add(type)
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.auto_equivalences[type] = ae
                
                
                # Find bond increments information    
                elif bond_inc_flag:
                    # Check if len(line) >= 6 and if line[0] is a float
                    if len(line) >= 6 and check_float(line[0]):
                        #print(line)
                        bi = Bond_incs()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        bi.ij = float(line[4])
                        bi.ji = float(line[5])
                        bi.ver = ver
                        bi.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j) not in self.bond_increments:
                            self.bond_increments[(i,j)] = bi
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j) in self.bond_increments:
                            logged_version = self.bond_increments[(i,j)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['bond_increments'].add(tuple([i,j]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.bond_increments[(i,j)] = bi
                            
                    
                # Find quadratic bond information
                elif quadraticbond_flag:
                    # Check if len(line) >= 6 and if line[0] is a float
                    if len(line) >= 6 and check_float(line[0]):
                        #print(line)
                        qb = Quadratic_bonds()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        qb.r0 = float(line[4])
                        qb.k2 = float(line[5])
                        qb.ver = ver
                        qb.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j) not in self.quadratic_bonds:
                            self.quadratic_bonds[(i,j)] = qb
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j) in self.quadratic_bonds:
                            logged_version = self.quadratic_bonds[(i,j)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['quadratic_bonds'].add(tuple([i,j]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.quadratic_bonds[(i,j)] = qb
                                
                # Find quadratic bond auto information
                elif quadraticbond_auto_flag:
                    # Check if len(line) >= 6 and if line[0] is a float
                    if len(line) >= 6 and check_float(line[0]):
                        #print(line)
                        qb = Quadratic_bonds()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        qb.r0 = float(line[4])
                        qb.k2 = float(line[5])
                        qb.ver = ver
                        qb.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j) not in self.quadratic_bonds_auto:
                            self.quadratic_bonds_auto[(i,j)] = qb
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j) in self.quadratic_bonds_auto:
                            logged_version = self.quadratic_bonds_auto[(i,j)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['quadratic_bonds auto'].add(tuple([i,j]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.quadratic_bonds[(i,j)] = qb
                                
                # Find morse bond information
                elif morsebond_flag:
                    # Check if len(line) >= 6 and if line[0] is a float
                    if len(line) >= 7 and check_float(line[0]):
                        #print(line)
                        mb = Morse_bonds()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        mb.r0 = float(line[4])
                        mb.d = float(line[5])
                        mb.alpha = float(line[6])
                        mb.ver = ver
                        mb.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j) not in self.morse_bonds:
                            self.morse_bonds[(i,j)] = mb
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j) in self.morse_bonds:
                            logged_version = self.morse_bonds[(i,j)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['morse_bonds'].add(tuple([i,j]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.morse_bonds[(i,j)] = mb
                                
                # Find morse bond information
                elif morsebond_auto_flag:
                    # Check if len(line) >= 6 and if line[0] is a float
                    if len(line) >= 7 and check_float(line[0]):
                        #print(line)
                        mb = Morse_bonds()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        mb.r0 = float(line[4])
                        mb.d = float(line[5])
                        mb.alpha = float(line[6])
                        mb.ver = ver
                        mb.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j) not in self.morse_bonds_auto:
                            self.morse_bonds_auto[(i,j)] = mb
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j) in self.morse_bonds_auto:
                            logged_version = self.morse_bonds_auto[(i,j)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['morse_bonds auto'].add(tuple([i,j]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.morse_bonds_auto[(i,j)] = mb
                            
                # Find quartic bond information
                elif quarticbond_flag:
                    # Check if len(line) >= 8 and if line[0] is a float
                    if len(line) >= 8 and check_float(line[0]):
                        #print(line)
                        qb = Quartic_bonds()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        qb.r0 = float(line[4])
                        qb.k2 = float(line[5])
                        qb.k3 = float(line[6])
                        qb.k4 = float(line[7])
                        qb.ver = ver
                        qb.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j) not in self.quartic_bonds:
                            self.quartic_bonds[(i,j)] = qb
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j) in self.quartic_bonds:
                            logged_version = self.quartic_bonds[(i,j)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['quartic_bonds'].add(tuple([i,j]))
    
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.quartic_bonds[(i,j)] = qb                        
                                
                # Find quadratic angle information
                elif quadraticangle_flag:
                    # Check if len(line) >= 7 and if line[0] is a float
                    if len(line) >= 7 and check_float(line[0]):
                        #print(line)
                        qa = Quadratic_angles()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        qa.theta0 = float(line[5])
                        qa.k2 = float(line[6])
                        qa.ver = ver
                        qa.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k) not in self.quadratic_angles:
                            self.quadratic_angles[(i,j,k)] = qa
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k) in self.quadratic_angles:
                            logged_version = self.quadratic_angles[(i,j,k)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['quadratic_angles'].add(tuple([i,j,k]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.quadratic_angles[(i,j,k)] = qa
                                
                # Find quadratic angle information
                elif quadraticangle_auto_flag:
                    # Check if len(line) >= 7 and if line[0] is a float
                    if len(line) >= 7 and check_float(line[0]):
                        #print(line)
                        qa = Quadratic_angles()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        qa.theta0 = float(line[5])
                        qa.k2 = float(line[6])
                        qa.ver = ver
                        qa.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k) not in self.quadratic_angles_auto:
                            self.quadratic_angles_auto[(i,j,k)] = qa
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k) in self.quadratic_angles_auto:
                            logged_version = self.quadratic_angles_auto[(i,j,k)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['quadratic_angles auto'].add(tuple([i,j,k]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.quadratic_angles_auto[(i,j,k)] = qa
                    
                # Find quartic angle information    
                elif quarticangle_flag:
                    # Check if len(line) >= 9 and if line[0] is a float
                    if len(line) >= 9 and check_float(line[0]):
                        #print(line)
                        qa = Quartic_angles()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        qa.theta0 = float(line[5])
                        qa.k2 = float(line[6])
                        qa.k3 = float(line[7])
                        qa.k4 = float(line[8])
                        qa.ver = ver
                        qa.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k) not in self.quartic_angles:
                            self.quartic_angles[(i,j,k)] = qa
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k) in self.quartic_angles:
                            logged_version = self.quartic_angles[(i,j,k)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['quartic_angles'].add(tuple([i,j,k]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.quartic_angles[(i,j,k)] = qa
                    
                # Find torsion 1 information    
                elif torsion1_flag:
                    # Check if len(line) >= 9 and if line[0] is a float
                    if len(line) >= 9 and check_float(line[0]):
                        #print(line)
                        t1 = Torsion1()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
                        t1.kphi = float(line[6])
                        t1.n = int(line[7])
                        t1.phi0 = float(line[8])
                        t1.ver = ver
                        t1.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.torsion_1:
                            self.torsion_1[(i,j,k,l)] = t1
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.torsion_1:
                            logged_version = self.torsion_1[(i,j,k,l)].ver
                            current_version = ver
    
                            # log use of higher versions
                            higher_versions_log['torsion_1'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.torsion_1[(i,j,k,l)] = t1
                                
                # Find torsion 1 auto information
                elif torsion1_auto_flag:
                    # Check if len(line) >= 9 and if line[0] is a float
                    if len(line) >= 9 and check_float(line[0]):
                        #print(line)
                        t1 = Torsion1()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
                        t1.kphi = float(line[6])
                        t1.n = int(line[7])
                        t1.phi0 = float(line[8])
                        t1.ver = ver
                        t1.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.torsion_1_auto:
                            self.torsion_1_auto[(i,j,k,l)] = t1
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.torsion_1_auto:
                            logged_version = self.torsion_1_auto[(i,j,k,l)].ver
                            current_version = ver
    
                            # log use of higher versions
                            higher_versions_log['torsion_1_auto'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.torsion_1_auto[(i,j,k,l)] = t1
                    
                # Find torsion 3 information    
                elif torsion3_flag:
                    # Check if len(line) >= 12 and if line[0] is a float
                    if len(line) >= 12 and check_float(line[0]):
                        #print(line)
                        t3 = Torsion3()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
                        t3.v1 = float(line[6])
                        t3.phi1 = float(line[7])
                        t3.v2 = float(line[8])
                        t3.phi2 = float(line[9])
                        t3.v3 = float(line[10])
                        t3.phi3 = float(line[11])
                        t3.ver = ver
                        t3.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.torsion_3:
                            self.torsion_3[(i,j,k,l)] = t3
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.torsion_3:
                            logged_version = self.torsion_3[(i,j,k,l)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['torsion_3'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.torsion_3[(i,j,k,l)] = t3
                    
                # Find wilson out of plane information    
                elif wilson1_flag:
                    # Check if len(line) >= 8 and if line[0] is a float
                    if len(line) >= 8 and check_float(line[0]):
                        #print(line)
                        w = Wilson_out_of_plane()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
                        w.kchi = float(line[6])
                        w.chi0 = float(line[7])
                        w.ver = ver
                        w.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.wilson_out_of_plane:
                            self.wilson_out_of_plane[(i,j,k,l)] = w
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.wilson_out_of_plane:
                            logged_version = self.wilson_out_of_plane[(i,j,k,l)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['wilson_out_of_plane'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.wilson_out_of_plane[(i,j,k,l)] = w
                
                # Find wilson out of plane auto information    
                elif wilson1_auto_flag:
                    # Check if len(line) >= 8 and if line[0] is a float
                    if len(line) >= 8 and check_float(line[0]):
                        #print(line)
                        w = Wilson_out_of_plane()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
                        w.kchi = float(line[6])
                        w.chi0 = float(line[7])
                        w.ver = ver
                        w.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.wilson_out_of_plane_auto:
                            self.wilson_out_of_plane_auto[(i,j,k,l)] = w
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.wilson_out_of_plane_auto:
                            logged_version = self.wilson_out_of_plane_auto[(i,j,k,l)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['wilson_out_of_plane_auto'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.wilson_out_of_plane_auto[(i,j,k,l)] = w
                                
                # Find class1 oop flag
                elif class1_oop_flag:
                    # Check if len(line) >= 8 and if line[0] is a float
                    if len(line) >= 9 and check_float(line[0]):
                        #print(line)
                        w = Wilson_out_of_plane()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
                        w.kchi = float(line[6])
                        w.n = int(line[7])
                        w.chi0 = float(line[8])
                        w.ver = ver
                        w.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.out_of_plane:
                            self.out_of_plane[(i,j,k,l)] = w
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.out_of_plane:
                            logged_version = self.out_of_plane[(i,j,k,l)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['out_of_plane'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.out_of_plane[(i,j,k,l)] = w
                                
                # Find class1 oop flag
                elif class1_oop_flag_auto:
                    # Check if len(line) >= 8 and if line[0] is a float
                    if len(line) >= 9 and check_float(line[0]):
                        #print(line)
                        w = Wilson_out_of_plane()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
                        w.kchi = float(line[6])
                        w.n = int(line[7])
                        w.chi0 = float(line[8])
                        w.ver = ver
                        w.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.out_of_plane_auto:
                            self.out_of_plane_auto[(i,j,k,l)] = w
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.out_of_plane_auto:
                            logged_version = self.out_of_plane_auto[(i,j,k,l)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['out_of_plane auto'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.out_of_plane_auto[(i,j,k,l)] = w
                    
                # Find pair coeff information    
                elif paircoeff_flag:
                    # Check if len(line) >= 5 and if line[0] is a float and if 9-6 or 12-6
                    if len(line) >= 5 and check_float(line[0]):
                        # Find 9-6
                        if pair_coeff_type_9_6:
                            p = Pair_coeff()
                            ver = float(line[0])
                            ref = float(line[1])
                            type = line[2]
                            p.r = float(line[3])
                            p.eps = float(line[4])
                            p.ver = ver
                            p.ref = ref
                            
                            # If parameter not in dictionary just
                            # add parameter to dictionary
                            if type not in self.pair_coeffs_9_6:
                                self.pair_coeffs_9_6[type] = p
                            
                            # If parameter in dictionary check
                            # version and log the higher version
                            elif type in self.pair_coeffs_9_6:
                                logged_version = self.pair_coeffs_9_6[type].ver
                                current_version = ver
                                
                                # log use of higher versions
                                higher_versions_log['pair_coeffs_9_6'].add(type)
                                
                                # log larger version. If logged_version
                                # is larger then current version let be.
                                # If current version is larger overwrite
                                # logged version.
                                if current_version >= logged_version:
                                    self.pair_coeffs_9_6[type] = p
                                    
                        # Find 12-6
                        if pair_coeff_type_12_6:
                            p = Pair_coeff()
                            ver = float(line[0])
                            ref = float(line[1])
                            type = line[2]
                            p.A = float(line[3])
                            p.B = float(line[4])
                            p.ver = ver
                            p.ref = ref
                            
                            # If parameter not in dictionary just
                            # add parameter to dictionary
                            if type not in self.pair_coeffs_12_6:
                                self.pair_coeffs_12_6[type] = p
                            
                            # If parameter in dictionary check
                            # version and log the higher version
                            elif type in self.pair_coeffs_12_6:
                                logged_version = self.pair_coeffs_12_6[type].ver
                                current_version = ver
                                
                                # log use of higher versions
                                higher_versions_log['pair_coeffs_12_6'].add(type)
                                
                                # log larger version. If logged_version
                                # is larger then current version let be.
                                # If current version is larger overwrite
                                # logged version.
                                if current_version >= logged_version:
                                    self.pair_coeffs_12_6[type] = p
                    
                    
                # Find bondbond information    
                elif bondbond_flag:
                    # Check if len(line) >= 6 and if line[0] is a float
                    if len(line) >= 6 and check_float(line[0]):
                        #print(line)
                        b = Bondbond()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        b.kb_bp = float(line[5])
                        b.ver = ver
                        b.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k) not in self.bondbond:
                            self.bondbond[(i,j,k)] = b
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k) in self.bondbond:
                            logged_version = self.bondbond[(i,j,k)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['bondbond'].add(tuple([i,j,k]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.bondbond[(i,j,k)] = b
                    
                # Find bondbond13 information    
                elif bondbond13_flag:
                    # Check if len(line) >= 7 and if line[0] is a float
                    if len(line) >= 7 and check_float(line[0]):
                        #print(line)
                        b = Bondbond13()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
                        b.kb_bp = float(line[6])
                        b.ver = ver
                        b.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.bondbond13:
                            self.bondbond13[(i,j,k,l)] = b
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.bondbond13:
                            logged_version = self.bondbond13[(i,j,k,l)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['bondbond13'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.bondbond13[(i,j,k,l)] = b
                    
                # Find bondangle information    
                elif bondangle_flag:
                    # Check if len(line) >= 6 and if line[0] is a float
                    if len(line) >= 6 and check_float(line[0]):
                        #print(line)
                        ba = Bondangle()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        ba.kb_theta = float(line[5])
                        ba.ver = ver
                        ba.ref = ref
                        try:
                            ba.kbp_theta = float(line[6])
                        except:
                            ba.kbp_theta = float(line[5])
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k) not in self.bondangle:
                            self.bondangle[(i,j,k)] = ba
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.bondangle:
                            logged_version = self.bondangle[(i,j,k)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['bondangle'].add(tuple([i,j,k]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.bondangle[(i,j,k)] = ba
                    
                # Find angleangle information    
                elif angleangle_flag:
                    # Check if len(line) >= 7 and if line[0] is a float
                    if len(line) >= 7 and check_float(line[0]):
                        #print(line)
                        aa = Angleangle()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
                        aa.k_theta_thetap = float(line[6])
                        aa.ver = ver
                        aa.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.angleangle:
                            self.angleangle[(i,j,k,l)] = aa
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.angleangle:
                            logged_version = self.angleangle[(i,j,k,l)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['angleangle'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.angleangle[(i,j,k,l)] = aa
                    
                # Find end bond torsion information    
                elif endbondtorsion_flag:
                    # Check if len(line) >= 9 and if line[0] is a float
                    if len(line) >= 9 and check_float(line[0]):
                        #print(line)
                        ebt = Endbondtorsion()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
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
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.endbondtorsion:
                            self.endbondtorsion[(i,j,k,l)] = ebt
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.endbondtorsion:
                            logged_version = self.endbondtorsion[(i,j,k,l)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['endbondtorsion'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.endbondtorsion[(i,j,k,l)] = ebt
                    
                    
                # Find middle torsion information    
                elif middlebondtorsion_flag:
                    # Check if len(line) >= 9 and if line[0] is a float
                    if len(line) >= 9 and check_float(line[0]):
                        #print(line)
                        mbt = Middlebondtorsion()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
                        mbt.f1 = float(line[6])
                        mbt.f2 = float(line[7])
                        mbt.f3 = float(line[8])
                        mbt.ver = ver
                        mbt.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.middlebondtorsion:
                            self.middlebondtorsion[(i,j,k,l)] = mbt
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.middlebondtorsion:
                            logged_version = self.middlebondtorsion[(i,j,k,l)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['middlebondtorsion'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.middlebondtorsion[(i,j,k,l)] = mbt
                    
                # Find angle torsion information    
                elif angletorsion_flag:
                    # Check if len(line) >= 9 and if line[0] is a float
                    if len(line) >= 9 and check_float(line[0]):
                        #print(line)
                        at = Angletorsion()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
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
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.angletorsion:
                            self.angletorsion[(i,j,k,l)] = at
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.angletorsion:
                            logged_version = self.angletorsion[(i,j,k,l)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['angletorsion'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.angletorsion[(i,j,k,l)] = at
                
                # Find angle angle torsion information    
                elif angleangletorsion_flag:
                    # Check if len(line) >= 7 and if line[0] is a float
                    if len(line) >= 7 and check_float(line[0]):
                        #print(line)
                        aat = Angleangletorsion()
                        ver = float(line[0])
                        ref = float(line[1])
                        i = line[2]
                        j = line[3]
                        k = line[4]
                        l = line[5]
                        aat.k_ang_ang_tor = float(line[6])
                        aat.ver = ver
                        aat.ref = ref
                        
                        # If parameter not in dictionary just
                        # add parameter to dictionary
                        if (i,j,k,l) not in self.angleangletorsion:
                            self.angleangletorsion[(i,j,k,l)] = aat
                        
                        # If parameter in dictionary check
                        # version and log the higher version
                        elif (i,j,k,l) in self.angleangletorsion:
                            logged_version = self.angleangletorsion[(i,j,k,l)].ver
                            current_version = ver
                            
                            # log use of higher versions
                            higher_versions_log['angleangletorsion'].add(tuple([i,j,k,l]))
                            
                            # log larger version. If logged_version
                            # is larger then current version let be.
                            # If current version is larger overwrite
                            # logged version.
                            if current_version >= logged_version:
                                self.angleangletorsion[(i,j,k,l)] = aat

        
        ###########################################
        # print out use of higher logged versions #
        ###########################################
        def print_higher_versions_log(higher_versions_log):
        
            # print buffer
            print()
            
            # pair_coeffs 9-6
            if higher_versions_log['pair_coeffs_9_6']:
                for type in higher_versions_log['pair_coeffs_9_6']:
                    version = self.pair_coeffs_9_6[type].ver
                    print('Using higher version of parameters for #pair_coeff {} version: {} '.format(type, "{:.2f}".format(version)))
            
            # pair_coeffs 12-6
            if higher_versions_log['pair_coeffs_12_6']:
                for type in higher_versions_log['pair_coeffs_12_6']:
                    version = self.pair_coeffs_12_6[type].ver
                    print('Using higher version of parameters for #pair_coeff {} version: {} '.format(type, "{:.2f}".format(version)))
        
            # equivalences
            if higher_versions_log['equivalences']:
                for type in higher_versions_log['equivalences']:
                    version = self.equivalences[type].ver
                    print('Using higher version of parameters for #equivalence {} version: {} '.format(type, "{:.2f}".format(version)))
                    
            # auto_equivalences
            if higher_versions_log['auto_equivalences']:
                for type in higher_versions_log['auto_equivalences']:
                    version = self.auto_equivalences[type].ver
                    print('Using higher version of parameters for #auto_equivalence {} version: {} '.format(type, "{:.2f}".format(version)))
                    
            # bond increments
            if higher_versions_log['bond_increments']:
                for (i, j) in higher_versions_log['bond_increments']:
                    version = self.bond_increments[(i,j)].ver
                    print('Using higher version of parameters for #bond-inc {} {} version: {} '.format(i, j, "{:.2f}".format(version)))
                    
            # quadratic bonds
            if higher_versions_log['quadratic_bonds']:
                for (i, j) in higher_versions_log['quadratic_bonds']:
                    version = self.quadratic_bonds[(i,j)].ver        
                    print('Using higher version of parameters for #quadratic_bond {} {} version: {} '.format(i, j, "{:.2f}".format(version)))
                    
            # quadratic bonds auto
            if higher_versions_log['quadratic_bonds auto']:
                for (i, j) in higher_versions_log['quadratic_bonds auto']:
                    version = self.quadratic_bonds_auto[(i,j)].ver        
                    print('Using higher version of parameters for #quadratic_bond {} {} version: {} '.format(i, j, "{:.2f}".format(version)))
                    
            # morse bonds
            if higher_versions_log['morse_bonds']:
                for (i, j) in higher_versions_log['morse_bonds']:
                    version = self.morse_bonds[(i,j)].ver        
                    print('Using higher version of parameters for #morse_bonds {} {} version: {} '.format(i, j, "{:.2f}".format(version)))
                    
            # morse bonds auto
            if higher_versions_log['morse_bonds auto']:
                for (i, j) in higher_versions_log['morse_bonds auto']:
                    version = self.morse_bonds_auto[(i,j)].ver        
                    print('Using higher version of parameters for #morse_bonds auto {} {} version: {} '.format(i, j, "{:.2f}".format(version)))
            
            # quartic bonds
            if higher_versions_log['quartic_bonds']:
                for (i, j) in higher_versions_log['quartic_bonds']:
                    version = self.quartic_bonds[(i,j)].ver        
                    print('Using higher version of parameters for #quartic_bond {} {} version: {} '.format(i, j, "{:.2f}".format(version)))
         
            # quadratic angles
            if higher_versions_log['quadratic_angles']:
                for (i, j, k) in higher_versions_log['quadratic_angles']:
                    version = self.quadratic_angles[(i,j,k)].ver                        
                    print('Using higher version of parameters for #quadratic_angle {} {} {} version: {} '.format(i, j, k, "{:.2f}".format(version)))
                    
            # quadratic angles auto
            if higher_versions_log['quadratic_angles auto']:
                for (i, j, k) in higher_versions_log['quadratic_angles auto']:
                    version = self.quadratic_angles_auto[(i,j,k)].ver                        
                    print('Using higher version of parameters for #quadratic_angle {} {} {} version: {} '.format(i, j, k, "{:.2f}".format(version)))
                    
            # quartic angles
            if higher_versions_log['quartic_angles']:
                for (i, j, k) in higher_versions_log['quartic_angles']:
                    version = self.quartic_angles[(i,j,k)].ver                        
                    print('Using higher version of parameters for #quartic_angle {} {} {} version: {} '.format(i, j, k, "{:.2f}".format(version)))
                            
            # torsion_1
            if higher_versions_log['torsion_1']:
                for (i, j, k, l) in higher_versions_log['torsion_1']:
                    version = self.torsion_1[(i,j,k,l)].ver                        
                    print('Using higher version of parameters for #torsion_1 {} {} {} {} version: {} '.format(i, j, k, l, "{:.2f}".format(version)))               
            
            # torsion_1 auto
            if higher_versions_log['torsion_1_auto']:
                for (i, j, k, l) in higher_versions_log['torsion_1_auto']:
                    version = self.torsion_1_auto[(i,j,k,l)].ver                        
                    print('Using higher version of parameters for #torsion_1 {} {} {} {} version: {} '.format(i, j, k, l, "{:.2f}".format(version)))               
                                                                       
            # torsion_3
            if higher_versions_log['torsion_3']:
                for (i, j, k, l) in higher_versions_log['torsion_3']:
                    version = self.torsion_3[(i,j,k,l)].ver                        
                    print('Using higher version of parameters for #torsion_3 {} {} {} {} version: {} '.format(i, j, k, l, "{:.2f}".format(version)))
            
            # wilson_out_of_plane
            if higher_versions_log['wilson_out_of_plane']:
                for (i, j, k, l) in higher_versions_log['wilson_out_of_plane']:
                    version = self.wilson_out_of_plane[(i,j,k,l)].ver                        
                    print('Using higher version of parameters for #wilson_out_of_plane {} {} {} {} version: {} '.format(i, j, k, l, "{:.2f}".format(version)))
                    
            # wilson_out_of_plane_auto
            if higher_versions_log['wilson_out_of_plane_auto']:
                for (i, j, k, l) in higher_versions_log['wilson_out_of_plane_auto']:
                    version = self.wilson_out_of_plane_auto[(i,j,k,l)].ver                        
                    print('Using higher version of parameters for #wilson_out_of_plane_auto {} {} {} {} version: {} '.format(i, j, k, l, "{:.2f}".format(version)))
                    
            # bondbond
            if higher_versions_log['bondbond']:
                for (i, j, k) in higher_versions_log['bondbond']:
                    version = self.bondbond[(i,j,k)].ver                        
                    print('Using higher version of parameters for #bondbond {} {} {} version: {} '.format(i, j, k, "{:.2f}".format(version)))
                    
            # bondbond13
            if higher_versions_log['bondbond13']:
                for (i, j, k, l) in higher_versions_log['bondbond13']:
                    version = self.bondbond13[(i,j,k,l)].ver                        
                    print('Using higher version of parameters for #bondbond13 {} {} {} {} version: {} '.format(i, j, k, l, "{:.2f}".format(version)))
                    
            # bondangle
            if higher_versions_log['bondangle']:
                for (i, j, k) in higher_versions_log['bondangle']:
                    version = self.bondangle[(i,j,k)].ver                        
                    print('Using higher version of parameters for #bondangle {} {} {} version: {} '.format(i, j, k, "{:.2f}".format(version)))
                    
            # angleangle
            if higher_versions_log['angleangle']:
                for (i, j, k, l) in higher_versions_log['angleangle']:
                    version = self.angleangle[(i,j,k,l)].ver                        
                    print('Using higher version of parameters for #angleangle {} {} {} {} version: {} '.format(i, j, k, l, "{:.2f}".format(version)))
                    
            # endbondtorsion
            if higher_versions_log['endbondtorsion']:
                for (i, j, k, l) in higher_versions_log['endbondtorsion']:
                    version = self.endbondtorsion[(i,j,k,l)].ver                        
                    print('Using higher version of parameters for #endbondtorsion {} {} {} {} version: {} '.format(i, j, k, l, "{:.2f}".format(version)))
            
            # middlebondtorsion
            if higher_versions_log['middlebondtorsion']:
                for (i, j, k, l) in higher_versions_log['middlebondtorsion']:
                    version = self.middlebondtorsion[(i,j,k,l)].ver                        
                    print('Using higher version of parameters for #middlebondtorsion {} {} {} {} version: {} '.format(i, j, k, l, "{:.2f}".format(version)))       
                    
            # angletorsion
            if higher_versions_log['angletorsion']:
                for (i, j, k, l) in higher_versions_log['angletorsion']:
                    version = self.angletorsion[(i,j,k,l)].ver                        
                    print('Using higher version of parameters for #angletorsion {} {} {} {} version: {} '.format(i, j, k, l, "{:.2f}".format(version)))
                    
            # angleangletorsion
            if higher_versions_log['angleangletorsion']:
                for (i, j, k, l) in higher_versions_log['angleangletorsion']:
                    version = self.angleangletorsion[(i,j,k,l)].ver                        
                    print('Using higher version of parameters for #angleangletorsion {} {} {} {} version: {} '.format(i, j, k, l, "{:.2f}".format(version)))
                    
            # torsiontorsion
            if higher_versions_log['torsiontorsion']:
                for (i, j, k, l) in higher_versions_log['torsiontorsion']:
                    version = self.torsiontorsion[(i,j,k,l)].ver                        
                    print('Using higher version of parameters for #torsiontorsion {} {} {} {} version: {} '.format(i, j, k, l, "{:.2f}".format(version)))

            # for debuggined
            # print('\n')
            # for i in higher_versions_log:
            #     if higher_versions_log[i]:
            #         for j in higher_versions_log[i]:
            #             print(f'{i}: {j}')
            #     else:
            #         print(f'{i} {higher_versions_log[i]}')
                    
            return
                    
        # print higher versions log info
        #print_higher_versions_log(higher_versions_log)
        