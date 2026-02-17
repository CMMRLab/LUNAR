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

Will write file as if it was meant for class 1 FF's:
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
import read_dreiding_dff
from datetime import date
import numpy as np


# frc FF file to read to generate .frc
file = 'all2lmp_dreiding.dff'    

# Set version of frc file
version = '1.1'

# Option to use DREIDING/A van der Waals parameters for Implicit Hydrogen atom-types. The written .frc file will contain
# the LJ-parameters for both DREIDNG and DREIDING/A, where the DREIDING will have paramters from Table II for only the 
# atom-types that do not have implicit hydrogens and DREIDING/A from Table VII will be used to generate the LJ-parameters
# for atom-types with implicit Hydrogens. Setting this option as True will set the NonB equivalences as the atom-type itself
# and will keep the atom-types with implicit Hydrogen 4th-character in the written #nonbond(12-6) sections. Setting this as
# False will set the NonB equivalences of the implicit Hydrogen containing atom-types to their nonimplicit Hydrogen atom-type
# counter parts. To have a complete implemention of DREDING AS Layed out in the 1990 paper leave this as False, to have more
# unique LJ-paramters that differ between the atom-type w/implicit H's and the atom-type w/o implicit H's as the older DREIDING/A
# FF, set as True. Josh recommneds the leave as False from his research of the web when trying to figure out how to get a DREIDING
# .frc file for all2lmp. Some .off and .lib files where found of DREIDING and they contained the info as if this setting was set
# to False.
use_DREIDING_A_impHs = False

# Find todays date
today = date.today()


# Function to replace info wisely
def replace_str_index(text,index=0,replacement=''):
    return '%s%s%s'%(text[:index],replacement,text[index+1:])


# Write drieding specifc .frc file
dreiding_frc = 'all2lmp_dreiding.frc'
with open(dreiding_frc, 'w') as f:
    
    # Write header
    f.write('! all2lmp Drieding specific forcefield with cvff class1 format except for inversion terms and energy definitons to be compatable\n')
    f.write('! with all2lmp ff_class=d. The "!" character still acts as a comment character and all2lmp will disregard any information after\n')
    f.write('! the "!" character\n\n')
    
    # Write Dreiding paper info
    f.write('! Implementation of energy coeffs comes from:\n')
    f.write('! Mayo, S. L., Olafson, B. D., & Goddard, W. A. (1990). DREIDING: a generic force field\n')
    f.write('! for molecular simulations. The Journal of Physical Chemistry, 94(26), 8897–8909.\n\n')
    
    # Write date
    f.write(f'#version {dreiding_frc}     {version}     {today.strftime("%b-%d-%Y")}\n')
    f.write('#define all2lmp DRIEDING .frc file\n')
    
    # Read dummy file
    ff = read_dreiding_dff.read(file)
    
    
    ####################
    # Write atom types #
    ####################
    f.write('#atom_types cvff\n\n')
    f.write('> Atom type definitions for any variant of cvff\n')
    f.write('> Masses from CRC 2014/2015 appendix B, using conventional masses.\n\n')

    f.write('!Ver Ref   Type         Mass      Element   Connection   Comment\n')
    f.write('!--- ---   -----      ----------  -------   ----------   ---------------------------\n')
    for i in ff.atom_types:
        a = ff.atom_types[i]
        f.write('{:^3} {:>4}   {:<8}    {:^10.6f}  {:<6}    {:<6}      {}\n'.format(a.ver, int(a.ref), a.type, a.mass, a.element, a.connection, a.comment))
        
    
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
    for i in ff.atom_types:
        a = ff.atom_types[i]
        
        #######################
        # Find nonbond equivs #
        #######################
        nonbond_equiv = '{}'.format(i.replace('Mid', ''))
        
        # If 'RB' in i for bridging between two rinfgs set nonbond equiv
        if 'RB' in i:
            element = '{}{}'.format(i[0], i[1])
            nonbond_equiv = '{}R'.format(element)
            
        # Set nonbond equivs depending on use_DREIDING_A_impHs for the implicit H's
        if not use_DREIDING_A_impHs and len(nonbond_equiv) >= 4:
            if nonbond_equiv[3].isdigit():
                nonbond_equiv = '{}'.format(replace_str_index(nonbond_equiv,index=3,replacement=''))
                
        
        ####################
        # Find bond equivs #
        ####################
        bond_equiv = '{}'.format(i.replace('Mid', ''))
        
        # Reduce bonds down to non-implicit H's
        if len(bond_equiv) >= 4:
            if bond_equiv[3].isdigit():
                bond_equiv = '{}'.format(replace_str_index(bond_equiv,index=3,replacement=''))
                
        # If 'RB' in i for bridging between two rinfgs set bond_equiv
        if 'RB' in i:
            element = '{}{}'.format(i[0], i[1])
            bond_equiv = '{}R'.format(element)
        
        #####################
        # Find angle equivs #
        #####################
        angle_equiv = '{}'.format(i.replace('Mid', ''))
        
        # Reduce angles down to non-implicit H's
        if len(angle_equiv) >= 4:
            if angle_equiv[3].isdigit():
                angle_equiv = '{}'.format(replace_str_index(angle_equiv,index=3,replacement=''))
                
        
        # If 'RB' in i for bridging between two rinfgs set angle_equiv
        if 'RB' in i:
            element = '{}{}'.format(i[0], i[1])
            angle_equiv = '{}R'.format(element)
        
        #######################
        # Find torsion equivs #
        #######################
        if len(i) >= 3:
            # If 3rd character is a digit set spN
            if i[2].isdigit():
                sp = i[2]
            elif i[2] == 'R':
                sp = 'R' # set resonance sp2 state
            else:
                sp = 0
        else:
            sp = 0
        torsion_equiv = 'Sp{}'.format(sp)
        
        # Use Mid keyword if in i
        if 'Mid' in i:
            torsion_equiv = '{}Mid'.format(torsion_equiv)
            
        # If i[0] == 'O' oxygens are treated as more unique and set torsion_equiv accordingly
        if i[0] =='O':
            torsion_equiv = 'O{}'.format(torsion_equiv)
            
        # If 'RB' in i for bridging between two rinfgs set torsion_equiv
        if 'RB' in i:
            element = '{}{}'.format(i[0], i[1])
            #torsion_equiv = '{}R'.format(element)
            torsion_equiv = 'SpRB'
            
            
            
        ###################
        # Find OOP equivs #
        ###################
        oop_equiv = '{}'.format(i.replace('Mid', ''))
        oop_equiv = '{}'.format('NonPlanar')
        
        if len(i) >= 3:
            # If 3rd character is a digit set spN
            if i[2] == 'R' or i[2] == '2':
                oop_equiv = '{}'.format('Planar')
            
        # Write equivs
        f.write('{:^3} {:>4}     {:<8}   {:<8}   {:<8}   {:<8}   {:<8}    {:<8}\n'.format(a.ver, int(a.ref), i, nonbond_equiv, bond_equiv, angle_equiv, torsion_equiv, oop_equiv))
    
    # Write reduction equiv for outer atoms equivs
    f.write('{:^3} {:>4}     {:<8}   {:<8}   {:<8}   {:<8}   {:<8}    {:<8}\n'.format(a.ver, int(a.ref), 'X', 'X', 'X', 'X', 'X', 'X'))
        
        
    ##########################################
    # Generate harmonic bond types and write #
    ##########################################
    f.write('\n\n\n')
    f.write('#quadratic_bond cvff\n\n') 

    f.write('> Computed with LAMMPS: bond_style harmonic\n')
    f.write('> E = K2 * (R - R0)^2\n\n')

    f.write('!Ver  Ref     I          J          R0         K2\n')  
    f.write('!---- ---   -----      -----     --------   --------\n') 
    used_bond_types = []
    for i in ff.bonds:
        for j in ff.bonds:
            # sort to make bond unique
            bond = sorted([i, j])
            ii, jj = bond
            
            # Find bond_i and bond_j
            bond_i = ff.bonds[ii];  bond_j = ff.bonds[jj];
            
            # set ref and ver from bond_i
            ver = bond_i.ver; ref = int(bond_i.ref);
            
            # Set bond radius 
            ri = bond_i.ri; rj = bond_j.ri;
            delta = bond_i.delta_harmonic # use from bond_i, bond_j has same value
            r0 = round(ri + rj - delta, 4)
            
            # Set bond stiffness (apply Drieding rule of thumb to take lowest bond order value).
            # LAMMPS harmonic potential   E=K2*(R-R0)^2 and 
            # Drieding harmonic potential E=0.5*ke*(R-Re)^2 (so K2 must be multiplied be 0.5)
            ke = bond_i.ke_harmonic  # use from bond_i, bond_j has same value
            bo_i = bond_i.n; bo_j = bond_j.n;
            bo = min([bo_i, bo_j]) # Find minimum bond order to find bond stiffness scale factor
            k2 = round(0.5*bo*ke, 4) 
            
            
            # Write bond if bond is not in used_bond_types 
            if bond not in used_bond_types:
                f.write('{:^3} {:>4}    {:<8}   {:<8} {:^10.4f} {:^10.4f}\n'.format(ver, ref, i, j, r0, k2))
                used_bond_types.append([ii, jj])
                
                
    #######################################
    # Generate morse bond types and write #
    #######################################
    f.write('\n\n\n')
    f.write('#morse_bond cvff\n\n') 

    f.write('> Computed with LAMMPS: bond_style morse\n')
    f.write('> E = D * (1 - exp(-ALPHA*(R - R0)))^2\n\n')

    f.write('!Ver  Ref     I          J          R0         D         ALPHA\n')  
    f.write('!---- ---   -----      -----     --------   --------    -------\n') 
    used_bond_types = []
    for i in ff.bonds:
        for j in ff.bonds:
            # sort to make bond unique
            bond = sorted([i, j])
            ii, jj = bond
            
            # Find bond_i and bond_j
            bond_i = ff.bonds[ii];  bond_j = ff.bonds[jj];
            
            # set ref and ver from bond_i
            ver = bond_i.ver; ref = int(bond_i.ref);
            
            # Set bond radius 
            ri = bond_i.ri; rj = bond_j.ri;
            delta = bond_i.delta_morse # use from bond_i, bond_j has same value
            r0 = round(ri + rj - delta, 4)
            
            # Set bond stiffness (apply Drieding rule of thumb to take lowest bond order value).
            # LAMMPS morse potential   E=D*(1 - exp(-ALPHA*(R - R0)))^2 
            # Drieding morse potential E=D*(exp(-ALPHA*n*(R - Re)) - 1)^2 (so ALPHA must be multiplied by n, which is bo)
            ke = bond_i.ke_morse  # use from bond_i, bond_j has same value
            de = bond_i.de_morse # use from bond_i, bond_j has same value
            bo_i = bond_i.n; bo_j = bond_j.n;
            bo = min([bo_i, bo_j]) # Find minimum bond order to find bond stiffness scale factor
            ke = round(bo*ke, 4) 
            d = round(bo*de, 4)
            
            alpha = (ke/(2*de))**(0.5)
            
            # Write bond if bond is not in used_bond_types 
            if bond not in used_bond_types:
                f.write('{:^3} {:>4}    {:<8}   {:<8} {:^10.4f} {:^10.4f} {:^10.4f}\n'.format(ver, ref, i, j, r0, d, alpha))
                used_bond_types.append([ii, jj])
                

    ###########################################
    # Generate harmonic angle types and write #
    ###########################################
    f.write('\n\n\n')
    f.write('#quadratic_angle cvff\n\n') 

    f.write('> Computed with LAMMPS: angle_style harmonic\n')
    f.write('> E = K2 * (Theta - Theta0)^2\n\n')

    f.write('!Ver  Ref     I          J          K          Theta0         K2\n')       
    f.write('!---- ---   -----      -----      -----       --------     -------\n')
    for j in ff.angles:
        angle = ff.angles[j]
        
        # set ref and ver from angle
        ver = angle.ver; ref = int(angle.ref);
        
        # Set anagle potential
        # LAMMPS morse potential   E=K2*(Theta - Theta0)^2
        # Drieding morse potential E=0.5*Kijk(Thetaijk-Theta0)^2 (so K2 must be multiplied be 0.5)
        k2 = round(0.5*angle.ke, 4)
        theta0 = angle.theta
        
        # set i and k as wildcards 
        i = '*'; k = '*'
        
        # Write angle info
        f.write('{:^3} {:>4}    {:<8}    {:<8}    {:<8} {:^10.4f}   {:^10.4f}\n'.format(ver, ref, i, j, k, theta0, k2))
        
    ##############################################
    # Generate harmonic dihedral types and write #
    ##############################################
    f.write('\n\n\n')
    f.write('#torsion_1 cvff\n\n') 

    f.write('> Computed with LAMMPS: dihedral_style harmonic\n')
    f.write('> E = Kphi * [ 1 + cos(n*Phi - Phi0) ]\n\n')

    f.write('!Ver  Ref    I        J       K       L          Kphi        n            Phi0\n')       
    f.write('!---- ---   ----     ----    ----    ----        -------     ------      -------\n')
    for t in ff.torsions:
        tor = ff.torsions[t]
        i = '*'; l = '*'; # set outer atoms as wildcards
        j, k = t # extract j and k inner atoms
        
        # set ref and ver from tor
        ver = tor.ver; ref = int(tor.ref);
        
        # Set torsion potential
        # LAMMPS harmonic potential   E = Kphi *  [ 1 + cos(n*Phi - Phi0) ]
        # Drieding harmonic potential E = 0.5*Vjk*[ 1 - cos(njk{phi - phi0jk}) ] (so Kphi=0.5*vjk)
        n = tor.n
        phi0 = tor.phi0 
        kphi = round((0.5*tor.v), 4)
        
        # Scale kphi by the multiplicity (Iakovos Delasoudas DREIDING multiplicity scaling issue 2/3/2026.)
        nj, nk = 1, 1
        if j[-1].isdigit():
            nj = int(j[-1])
        if j.endswith('R'):
            nj = 2
        if j.endswith('Mid'):
            tmp = j.replace('Mid', '')
            if tmp[-1].isdigit():
                nj = int(tmp[-1])
                
        if k[-1].isdigit():
            nk = int(k[-1])
        if k.endswith('R'):
            nk = 2
        if k.endswith('Mid'):
            tmp = k.replace('Mid', '')
            if tmp[-1].isdigit():
                nk = int(tmp[-1])
        
        kphi = kphi/(nj*nk)

        
        # LAMMPS uses a different equation then drieding so we need to manipulate 0, 90, and 180 accordingly
        if phi0 == 180:
            phi0 = 0
        elif phi0 == 0:
            phi0 = 180
        elif phi0 == 90:
            phi0 = 0
        else:
            phi0 = 0 # Set as safety
        
        # Write torsion info
        f.write('{:^3} {:>4}    {:<8} {:<8} {:<8} {:<8}  {:<10.4f}  {:<10}  {:<10.4f}\n'.format(ver, ref, i, j, k, l, kphi, n,  phi0))
 
        
 
    #############################
    # Write umbrella inversions #
    #############################
    f.write('\n\n\n')
    f.write('#out_of_plane DREIDING\n\n') 

    f.write('> Computed with LAMMPS: improper_style umbrella\n')
    f.write('> E = Kl*(1 - cos(phi)); if phi0==0\n')
    f.write('> E = 0.5*Kl*[(1/sin(phi0)^2]*[(cos(phi) - cos(phi0))^2]; if phi0!=0\n\n')
    
    f.write('!Ver  Ref    I        J          K       L           Kl         phi0\n')       
    f.write('!---- ---   ----     ----        ----    ----        -----      --------\n')
    for j in ff.inversions:
        inv = ff.inversions[j]
        
        # set ref and ver from inv
        ver = inv.ver; ref = int(inv.ref);
        
        # Set i, k, l as wildcards
        i = '*'; k = '*'; l = '*';
        
        # Write inversion info
        f.write('{:^3} {:>4}    {:<8} {:<12} {:<8} {:<8}  {:<10.4f}  {:<10.4f}\n'.format(ver, ref, i, j, k, l, inv.k, inv.phi0))


    ######################################
    # Write nonbond parms in cvff format #
    ######################################   
    f.write('\n\n\n')
    f.write('#nonbond(12-6)	cvff \n\n') 

    f.write('>@type A-B\n')
    f.write('>@combination geometric\n\n')
    
    f.write('> Computed with LAMMPS: pair_Style lj/cut/coul/long\n')
    f.write('> E = Aij/r^12 - Bij/r^6\n')
    f.write('> where  Aij = sqrt( Aii * Ajj )\n')
    f.write('>        Bij = sqrt( Bii * Bjj )\n\n')
    
    
    f.write('!Ver  Ref    I               A                           B \n')       
    f.write('!---- ---   ----        -----------                 -----------\n')
    for i in ff.nonbonds:
        nb = ff.nonbonds[i]
        
        # set ref and ver from nb
        ver = nb.ver; ref = int(nb.ref);
        
        # Find Dreiding r0 and d0
        r0 = nb.r0; d0 = nb.d0;
        
        # Compute the sigma value and epsilon
        epsilon = d0 #d0/4 not sure if the div-by 4 is needed
        c=2**(1.0/6.0)
        sigma=r0/c
        
        # Find A and B from sigma
        A = 4*epsilon*sigma**12
        B = 4*epsilon*sigma**6
        
        # Find it i is an atom-type of implicitH's (intialize as False)
        impH_type = False
        if len(i) >= 4 and i[3].isdigit():
            impH_type = True
        
        # Write nonbond A and B
        if use_DREIDING_A_impHs: # If True write all info to file
            f.write('{:^3} {:>4}    {:<8}    {:<24.8f}    {:<24.8f}\n'.format(ver, ref, i, A, B))
        else:
            if not impH_type:
                f.write('{:^3} {:>4}    {:<8}    {:<24.8f}    {:<24.8f}\n'.format(ver, ref, i, A, B))
                
    ############################################
    # Write nonbond parms in Buckingham format #
    ############################################
    f.write('\n\n\n')
    f.write('#nonbond(Buckingham) DREIDING/LAMMPS \n\n') 

    f.write('>@type A-B-C\n')
    f.write('>@combination rules:\n')
    f.write('> Aij = [Aii*Ajj]^(0.5) \n')
    f.write('> Bij = [Bii*Bjj]^(0.5) \n')
    f.write('> Cij = 0.5*Cii + 0.5*Cjj \n')
    f.write('\n')
    
    f.write('>DREIDING EQN. 32\n')
    f.write('> E = Ae^(-C*r) - B*r^-6\n')
    f.write('> where  Aij = sqrt( Aii * Ajj )\n')
    f.write('>        Bij = sqrt( Bii * Bjj )\n\n')
    
    f.write('>Will be converted to LAMMPS: pair_Style buck/coul/long\n')
    f.write('> E = Ae^(-r/rho) - C/r^6\n')
    f.write('> where  A_LAMMPS   = A_DREDING\n')
    f.write('>        rho_LAMMPS = 1/C_DREIDING\n')
    f.write('>        C_LAMMPS   = B_DREIDING\n\n')
    
    
    f.write('!Ver  Ref    I               A                           B                           C \n')       
    f.write('!---- ---   ----        -----------                 -----------                 -----------\n')
    for i in ff.nonbonds:
        nb = ff.nonbonds[i]
        
        # set ref and ver from nb
        ver = nb.ver; ref = int(nb.ref);
        
        # Find Dreiding r0, d0, and xi
        r0 = nb.r0; d0 = nb.d0; xi = nb.xi
        
        # Convert to A, B, C form for X6 eqn. 32
        A = float( (6*d0/(xi - 6))*np.exp(xi) )
        B = float( xi*d0/(xi - 6)*r0**6 )
        C = float( xi/r0 )
        
        f.write('{:^3} {:>4}    {:<8}    {:<24.8f}    {:<24.8f}    {:<24.8f}\n'.format(ver, ref, i, A, B, C))


    ###########################################
    # Write Hbonding parameters for Gasteiger #
    ###########################################
    f.write('\n\n\n')
    f.write('#nonbond hbonding-Gasteiger \n\n') 
    
    f.write('> E = D_hb*( 5*(R_hb/r)^12 - 6*(R_hb/r)^10 )*cos^4(theta_DHA) \n\n')  
    
    f.write('!Ver  Ref     I          J         D_hb       R_hb      theta_DHA\n')  
    f.write('!---- ---   -----      -----     --------   --------    ---------\n') 
    for i, j in ff.hbond_gasteiger:
        hb = ff.hbond_gasteiger[(i, j)]
        ver = hb.ver; ref = int(hb.ref)       
        f.write('{:^3} {:>4}    {:<8}   {:<8} {:^10.4f} {:^10.4f}  {:^10.4f}\n'.format(ver, ref, i, j, hb.dhb, hb.rhb, hb.theta))
        
    ###########################################
    # Write Hbonding parameters for no-charge #
    ###########################################
    f.write('\n\n\n')
    f.write('#nonbond hbonding-no-charge \n\n') 
    
    f.write('> E = D_hb*( 5*(R_hb/r)^12 - 6*(R_hb/r)^10 )*cos^4(theta_DHA) \n\n')  
    
    f.write('!Ver  Ref     I          J         D_hb       R_hb      theta_DHA\n')  
    f.write('!---- ---   -----      -----     --------   --------    ---------\n') 
    for i, j in ff.hbond_nocharge:
        hb = ff.hbond_nocharge[(i, j)]
        ver = hb.ver; ref = int(hb.ref)       
        f.write('{:^3} {:>4}    {:<8}   {:<8} {:^10.4f} {:^10.4f}  {:^10.4f}\n'.format(ver, ref, i, j, hb.dhb, hb.rhb, hb.theta))

        
        
    ##################################################################
    # Write nonbond parms in sigma/epsilon format (call it DREIDING) #
    ##################################################################   
    # f.write('\n\n\n')
    # f.write('#nonbond(12-6)	DREIDING\n\n') 

    # f.write('>@type Sigma-Epsilon\n')
    # f.write('>@combination geometric\n\n')
    
    # f.write('> Computed with LAMMPS: pair_Style lj/cut/coul/long\n')
    # f.write('> E = 4*Epsilon( (Sigma/R)^12 -  (Sigma/R)^6) \n')
    # f.write('> where  R < Rc\n\n')
    
    
    # f.write('!Ver  Ref    I            Epsilon             Sigma\n')       
    # f.write('!---- ---   ----        -----------      -----------\n')
    # for i in ff.nonbonds:
    #     nb = ff.nonbonds[i]
        
    #     # set ref and ver from nb
    #     ver = nb.ver; ref = int(nb.ref);
        
    #     # Find Dreiding r0 and d0
    #     r0 = nb.r0; d0 = nb.d0;
        
    #     # Compute the sigma value and epsilon
    #     epsilon = d0 #d0/4 not sure if the div-by 4 is needed
    #     c=2**(1.0/6.0)
    #     sigma=r0/c
        
    #     # Find it i is an atom-type of implicitH's (intialize as False)
    #     impH_type = False
    #     if len(i) >= 4 and i[3].isdigit():
    #         impH_type = True
        
    #     # Write nonbond epsilon and sigma
    #     if use_DREIDING_A_impHs: # If True write all info to file
    #         f.write('{:^3} {:>4}    {:<8}    {:<16.6f}    {:<16.6f}\n'.format(ver, ref, i, epsilon, sigma))
    #     else:
    #         if not impH_type:
    #             f.write('{:^3} {:>4}    {:<8}    {:<16.6f}    {:<16.6f}\n'.format(ver, ref, i, epsilon, sigma))

     
        
    ##############
    # Write refs #
    ##############
    f.write('\n\n\n')
    f.write('\n#reference 1\n')
    f.write('@Author Stephen L. Mayo, Barry D. Olafson, and William A. Goddard III*’1\n')
    f.write('@Date Jan-28-2023\n')
    f.write('DREIDING: A Generic Force Field for Molecular Simulations\n')
    
    
    f.write('\n#reference 2\n')
    f.write('@Author tester\n')
    f.write('@Date Jan-28-2023\n')
    f.write('CRC Handbook of Chemistry and Physics, 95th Edition\n')
    
    
print(f'\n\n\nAll done writing {dreiding_frc}\n\n\n')
    
        