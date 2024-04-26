# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
December 8th, 2021
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
import read_prm
from datetime import date


# frc FF file to read to generate .frc
infile = 'oplsaa.prm' 
outfile = 'all2lmp_oplsaa_tinker.frc'    

# Set version of frc file
version = '1.0'

# Option to set vdw equivs as typeID (default is to use classID)
set_vdw_from_typeID = True

# Find todays date
today = date.today()

# atomic_number to element symbol dict
an2element = {1:'H', 3:'Li', 11:'Na', 19:'K', 37:'Rb', 55:'Cs', 87:'Fr',         # Group 1
              4:'Be', 12:'Mg', 20:'Ca', 38:'Sr', 56:'Ba', 88:'Ra',               # Group 2
              21:'Sc', 39:'Y', 57:'La', 89:'Ac',                                 # Group 3
              22:'Ti', 40: 'Zr', 72:'Hf', 104:'Rf',                              # Group 4
              23:'V', 41: 'Nb', 73:'Ta', 105:'Db',                               # Group 5
              24:'Cr', 42:'Mo', 74:'W', 106:'Sg',                                # Group 6
              25:'Mn', 43:'Tc', 75:'Re', 107:'Bh',                               # Group 7
              26:'Fe', 44:'Ru', 76:'Os', 108:'Hs',                               # Group 8
              27:'Co', 45:'Rh', 77:'Ir', 109:'Mt',                               # Group 9
              28:'Ni', 46:'Pd', 78:'Pt', 110:'Ds',                               # Group 10
              29:'Cu', 47:'Ag', 79:'Au', 111:'Rg',                               # Group 11
              30:'Zn', 48:'Cd', 80:'Hg', 112:'Cn',                               # Group 12
              5:'B', 13:'Al', 31:'Ga', 49:'In', 81:'Ti', 113:'Nh',               # Group 13
              6:'C', 14:'Si', 32:'Ge', 50:'Sn', 82:'Pb', 114:'Fi',               # Group 14
              7:'N', 15:'P', 33:'As', 51:'Sb', 83:'Bi', 115:'Mc',                # Group 15
              8:'O', 16:'S', 34:'Se', 52:'Te', 84:'Po', 116:'Lv',                # Group 16
              9:'F', 17:'Cl', 35:'Br', 53:'I', 85:'At', 117:'Ts',                # Group 17
              2:'He', 10:'Ne', 18:'Ar', 36:'Kr', 54:'Xe', 86:'Rn', 118:'Og',     # Group 18 
              58:'Ce', 59:'Pr', 60:'Nd', 61:'Pm', 62:'Sm', 63:'Eu', 64:'Gd',     # Group 6 Ianthanoid series part1
              65:'Tb', 66:'Dy', 67:'Ho', 68:'Er', 69:'Tm', 70:'Yb', 71:'Lu',     # Group 6 Ianthanoid series part2
              90:'Th', 91:'Pa', 92:'U', 93:'Np', 94:'Pu', 95:'Am', 96:'Cm',      # Group 7 actinoid series part1
              97:'Bk', 98:'Cf', 99:'Es', 100:'Fm', 101:'Md', 102:'No', 103:'Lr', # Group 7 actinoid series part2   
              0:'Unknown'                                                        # If atomic number is unknown                             
              }

# Function to find classID from atom_types
def get_atomtype_from_classID(parms, classID):
    atom_type = '*'
    for i in parms.atom_types:
        a = parms.atom_types[i]
        if a.classID == classID:
            atom_type = '{}_{}'.format(i[0], a.classID)
            break
    return atom_type


# Write outfile .frc file
with open(outfile, 'w') as f:
    
    # Write header
    f.write('! Conversion from Tinker v8.10.5 .prm format to Material Studio BIOSYM .frc format. The "!" character still\n')
    f.write('! acts as a comment character and all2lmp will disregard any information after the "!" character\n\n')
    
    # Write date
    f.write(f'#version {outfile}     {version}     {today.strftime("%b-%d-%Y")}\n')
    f.write('#define all2lmp Tinker .prm conversion to Material Studio .frc file\n\n')
    
    # Read prm file
    parms = read_prm.prm(infile)
    
    ####################
    # Write atom types #
    ####################
    f.write('#atom_types cvff\n\n')
    f.write('> Atom type definitions for any variant of cvff\n')
    f.write('> Masses Tinker .prm file and elements set by Tinkers atomic number information.\n\n')
    
    f.write('!Ver Ref   Type         Mass      Element   Connection   Comment\n')
    f.write('!--- ---   -----      ----------  -------   ----------   ---------------------------\n')
    atom_types = sorted(parms.atom_types.keys())
    for i in atom_types:
        a = parms.atom_types[i]
        try: element = an2element[a.atomic_number]
        except: element = an2element[0]
        newtype = '{}_{}'.format(i[0], i[1])
        comment = a.comment
        try: comment += ' (partial charge: {})'.format(parms.charges[a.classID])
        except: comment += ' (partial charge: failed)' 
        ver = 1.0; ref = 1;
        f.write('{:^3} {:>4}   {:<8}    {:^10.6f}  {:<8}    {:<6}      {}\n'.format(ver, ref, newtype, a.mass, element, a.connection, comment))
        

    ######################
    # Write equivalences #
    ######################
    f.write('\n\n\n')
    f.write('#equivalence cvff\n\n') 

    f.write('> Equivalence table for any variant of cvff\n\n')

    f.write('!		         	                     Equivalences\n')
    f.write('!                 -------------------------------------------------------\n')
    f.write('!Ver  Ref    Type       NonB       Bond       Angle      Torsion      OOP\n')
    f.write('!---- ---    ----       ----       ----       -----      -------     -----\n')
    for i in atom_types:
        a = parms.atom_types[i]
        newtype = '{}_{}'.format(i[0], i[1])
        
        # Find equivs
        vdw_equiv = '{}_{}'.format(i[0], a.classID)
        if set_vdw_from_typeID: vdw_equiv = '{}_{}'.format(i[0], i[1])
        bond_equiv = '{}_{}'.format(i[0], a.classID)
        angle_equiv = '{}_{}'.format(i[0], a.classID)
        torsion_equiv = '{}_{}'.format(i[0], a.classID)
        oop_equiv = '{}_{}'.format(i[0], a.classID)
        ver = 1.0; ref = 1;
        f.write('{:^3} {:>4}     {:<8}   {:<8}   {:<8}   {:<8}   {:<8}    {:<8}\n'.format(ver, ref, newtype, vdw_equiv, bond_equiv, angle_equiv, torsion_equiv, oop_equiv))
        

    ####################
    # Write bond types #
    ####################
    f.write('\n\n\n')
    f.write('#quadratic_bond cvff\n\n') 

    f.write('> Computed with LAMMPS: bond_style harmonic\n')
    f.write('> E = K2 * (R - R0)^2\n\n')

    f.write('!Ver  Ref     I          J          R0         K2\n')  
    f.write('!---- ---   -----      -----     --------   --------\n') 
    for i in parms.quadratic_bonds:
        itype = get_atomtype_from_classID(parms, i[0])
        jtype = get_atomtype_from_classID(parms, i[1])
        r0 = parms.quadratic_bonds[i].r0
        k2 = parms.quadratic_bonds[i].k2
        ver = 1.0; ref = 1;
        f.write('{:^3} {:>4}    {:<8}   {:<8} {:^10.4f} {:^10.4f}\n'.format(ver, ref, itype, jtype, r0, k2))
        
        
    #####################
    # Write angle types #
    #####################
    f.write('\n\n\n')
    f.write('#quadratic_angle cvff\n\n') 

    f.write('> Computed with LAMMPS: angle_style harmonic\n')
    f.write('> E = K2 * (Theta - Theta0)^2\n\n')

    f.write('!Ver  Ref     I          J          K          Theta0         K2\n')       
    f.write('!---- ---   -----      -----      -----       --------     -------\n')
    for i in parms.quadratic_angles:
        itype = get_atomtype_from_classID(parms, i[0])
        jtype = get_atomtype_from_classID(parms, i[1])
        ktype = get_atomtype_from_classID(parms, i[2])
        theta0 = parms.quadratic_angles[i].theta0
        k2 = parms.quadratic_angles[i].k2
        ver = 1.0; ref = 1;
        f.write('{:^3} {:>4}    {:<8}    {:<8}    {:<8} {:^10.4f}   {:^10.4f}\n'.format(ver, ref, itype, jtype, ktype, theta0, k2))
    
    #######################
    # Write torsion types #
    #######################
    f.write('\n\n\n')
    f.write('#torsion_1 opls\n\n') 

    f.write('> Computed with LAMMPS: dihedral_style opls\n')
    f.write('> E = 0.5*K1*[1 + cos(Phi)] + 0.5*K2*[1 + cos(2*Phi)] + 0.5*K3*[1 + cos(3*Phi)] + 0.5*K4*[1 + cos(4*Phi)]\n\n')

    f.write('!Ver  Ref    I        J       K       L             K1         K2         K3         K4  \n')       
    f.write('!---- ---   ----     ----    ----    ----        -------     ------     ------     ------\n')
    for i in parms.torsion_1_opls:
        itype = get_atomtype_from_classID(parms, i[0])
        jtype = get_atomtype_from_classID(parms, i[1])
        ktype = get_atomtype_from_classID(parms, i[2])
        ltype = get_atomtype_from_classID(parms, i[3])
        k1 = parms.torsion_1_opls[i].k1
        k2 = parms.torsion_1_opls[i].k2
        k3 = parms.torsion_1_opls[i].k3
        k4 = parms.torsion_1_opls[i].k4
        ver = 1.0; ref = 1;
        f.write('{:^3} {:>4}    {:<8} {:<8} {:<8} {:<8}  {:<10.4f}  {:<10.4f}  {:<10.4f}  {:<10.4f}\n'.format(ver, ref, itype, jtype, ktype, ltype, k1, k2, k3, k4))
        
    ###################
    # Write OOP types #
    ###################
    f.write('\n\n\n')
    f.write('#out_of_plane cvff\n\n') 

    f.write('> Computed with LAMMPS: improper_style harmonic\n')
    f.write('> E = Kchi * [ 1 + cos(n*Chi - Chi0) ]\n\n')
    
    f.write('!Ver  Ref    I        J       K       L           Kchi       n        chi0\n')       
    f.write('!---- ---   ----     ----    ----    ----        -----     -----     --------\n')
    for i in parms.out_of_plane:
        # Tinker has backwards ordering as compared to LAMMPS ordering of impropers
        itype = get_atomtype_from_classID(parms, i[3])
        jtype = get_atomtype_from_classID(parms, i[2])
        ktype = get_atomtype_from_classID(parms, i[1])
        ltype = get_atomtype_from_classID(parms, i[0])
        kchi = parms.out_of_plane[i].kchi
        n = parms.out_of_plane[i].n
        chi0 = parms.out_of_plane[i].chi0
        ver = 1.0; ref = 1;
        f.write('{:^3} {:>4}    {:<8} {:<8} {:<8} {:<8}  {:<10.4f}  {:<10}  {:<10.4f}\n'.format(ver, ref, itype, jtype, ktype, ltype, kchi, n, chi0))
        

    ######################################
    # Write nonbond parms in cvff format #
    ######################################   
    f.write('\n\n\n')
    f.write('#nonbond(12-6)	cvff \n\n') 

    f.write('>@type A-B\n')
    f.write('>@combination geometric\n\n')
    
    f.write('> Computed with LAMMPS: pair_Style lj/cut/coul/long\n')
    f.write('> E = Aij/r^12 - Bij/r^6\n')
    f.write('> where  Aij = sqrt( Ai * Aj )\n')
    f.write('>        Bij = sqrt( Bi * Bj )\n\n')
    
    
    f.write('!Ver  Ref    I               A                B \n')       
    f.write('!---- ---   ----        -----------      -----------\n')
    for i in parms.pair_coeffs_12_6:
        itype = get_atomtype_from_classID(parms, i)
        if set_vdw_from_typeID:
            atom_type = '*'
            for j in parms.atom_types:
                a = parms.atom_types[j]
                if a.TypeID == i:
                    itype = '{}_{}'.format(j[0], j[1])
                    break
        
        # Compute the sigma value and epsilon
        r = parms.pair_coeffs_12_6[i].r
        epsilon = parms.pair_coeffs_12_6[i].eps
        c=2**(1.0/6.0)
        sigma=r/c
        
        # Find A and B from sigma
        A = 4*epsilon*sigma**12
        B = 4*epsilon*sigma**6
        
        # Write coeffs
        f.write('{:^3} {:>4}    {:<8}    {:<16.6f}    {:<16.6f}\n'.format(ver, ref, itype, A, B))
        
print(f'{outfile} was successfully written. NOTE that * wildcards may have been inserted')
print('in any of the bond, angle, dihedral, and improper coeffs, no checks are performed in')
print('this script for "uniqueness" of the parameters with the * wildcard. It is up to the')
print('user to add numbers to the "backside" of the wildcard operators to ensure a unique')
print('string of types exists (IE * CA_1 CA1 may appear twice in the angle coeffs, where one')
print('can be left as * CA_1 CA_1 and the other should be changed to *1 CA_1 CA_1)')