#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
July 5th, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

        
class read_assumed_coeffs:
    def __init__(self, file):
    
        self.bond_coeffs = {} # { tuple(element1, element2) : tuple(assumed type1, assumed type2) }
        self.angle_coeffs = {} # { tuple(element1, element2, element3) : tuple(assumed type1, assumed type2, assumed type3) }
        self.dihedral_coeffs = {} # { tuple(element1, element2, element3, element4) : tuple(assumed type1, assumed type2, assumed type3, assumed type4) }
        self.improper_coeffs = {} # { tuple(element1, element2, element3, element4) : tuple(assumed type1, assumed type2, assumed type3, assumed type4) }
        
        
        
        # Opening and reading the assumed coeffs file
        with open(file, 'r') as f:
            
            # Initializing flags    
            bond_flag = False
            angle_flag = False
            dihedral_flag = False
            improper_flag = False
            skip = 0
        
            # Looping through each line of the frc file
            for line in f:
                skip -= 1
                if skip >= 0:
                    continue
                #print(line)
                
                # Strip comment's
                line = line.split('#')[0]
                line = line.rstrip()
                
                # Split line
                line = line.strip()
                line_split = line.split()
                #print(line_split)
                
                
                # Setting flags
                if line == '':
                    bond_flag = False
                    angle_flag = False
                    dihedral_flag = False
                    improper_flag = False
                elif line_split[0] == 'Bond' and line_split[1] == 'Coeffs':
                    bond_flag = True
                    skip = 1
                    continue
                elif line_split[0] == 'Angle' and line_split[1] == 'Coeffs':
                    angle_flag = True
                    skip = 1
                    continue
                elif line_split[0] == 'Dihedral' and line_split[1] == 'Coeffs':
                    dihedral_flag = True
                    skip = 1
                    continue
                elif line_split[0] == 'Improper' and line_split[1] == 'Coeffs':
                    improper_flag = True
                    skip = 1
                    continue            
                
                
                # Finding bond coeffs information    
                if bond_flag:
                    #print(line_split)
                    # Find elements in coeff
                    element_1 = line_split[0]
                    element_2 = line_split[1]
                    
                    # Find assumed types
                    type_1 = line_split[2]
                    type_2 = line_split[3]
                    
                    # Save info into tuples
                    coeff_elements = (element_1, element_2)
                    coeff_types = (type_1, type_2)
                    
                    # Save to dictionary
                    self.bond_coeffs[coeff_elements] = coeff_types
        
                    
                # Finding angle coeffs information
                elif angle_flag:
                    #print(line_split)
                    # Find elements in coeff
                    element_1 = line_split[0]
                    element_2 = line_split[1]
                    element_3 = line_split[2]
                    
                    # Find assumed types
                    type_1 = line_split[3]
                    type_2 = line_split[4]
                    type_3 = line_split[5]
                    
                    # Save info into tuples
                    coeff_elements = (element_1, element_2, element_3)
                    coeff_types = (type_1, type_2, type_3)
                    
                    # Save to dictionary
                    self.angle_coeffs[coeff_elements] = coeff_types
                    
                # Finding dihedral coeffs information    
                elif dihedral_flag:
                    #print(line_split)
                    # Find elements in coeff
                    element_1 = line_split[0]
                    element_2 = line_split[1]
                    element_3 = line_split[2]
                    element_4 = line_split[3]
                    
                    # Find assumed types
                    type_1 = line_split[4]
                    type_2 = line_split[5]
                    type_3 = line_split[6]
                    type_4 = line_split[7]
                    
                    # Save info into tuples
                    coeff_elements = (element_1, element_2, element_3, element_4)
                    coeff_types = (type_1, type_2, type_3, type_4)
                    
                    # Save to dictionary
                    self.dihedral_coeffs[coeff_elements] = coeff_types
                    
                # Finding imroper coeffs information    
                elif improper_flag:
                    # Find elements in coeff
                    element_1 = line_split[0]
                    element_2 = line_split[1]
                    element_3 = line_split[2]
                    element_4 = line_split[3]
                    
                    # Find assumed types
                    type_1 = line_split[4]
                    type_2 = line_split[5]
                    type_3 = line_split[6]
                    type_4 = line_split[7]
                    
                    # Save info into tuples
                    coeff_elements = (element_1, element_2, element_3, element_4)
                    coeff_types = (type_1, type_2, type_3, type_4)
                    
                    # Save to dictionary
                    self.improper_coeffs[coeff_elements] = coeff_types
                    
                    
def check_assumed_auto_fill_types(frc, aafc, frc_file, log):
    
    ########################
    # print start of check #
    ########################
    log.out(f'\n\nChecking that assigned assumed coeff types exists in read in .frc file: {frc_file}')
    
    #############################
    # Check assumed bond coeffs #
    #############################
    flagged_bond_coeff = set([]);
    for i in aafc.bond_coeffs:
        assumed_1, assumed_2 = aafc.bond_coeffs[i]
        
        # Check if the type exists in .frc file
        if (assumed_1, assumed_2) in frc.quartic_bonds: pass
        elif (assumed_2, assumed_1) in frc.quartic_bonds: pass
        elif (assumed_1, assumed_2) in frc.quadratic_bonds: pass
        elif (assumed_2, assumed_1) in frc.quadratic_bonds: pass
        # if it is not in .frc file log fail
        else:
            flagged_bond_coeff.add(i)
            log.warn(f'WARNING assumed bond elements {i[0]} {i[1]} does not have coeff type {assumed_1} {assumed_2} is not in .frc file. Removing from read assumed types to avoid errors.')
            
    ##############################
    # Check assumed angle coeffs #
    ##############################
    flagged_angle_coeff = set([]);
    for i in aafc.angle_coeffs:
        assumed_1, assumed_2, assumed_3 = aafc.angle_coeffs[i]
        
        # Check if the type exists in .frc file
        if (assumed_1, assumed_2, assumed_3) in frc.quartic_angles: pass
        elif (assumed_3, assumed_2, assumed_1) in frc.quartic_angles: pass
        elif (assumed_1, assumed_2, assumed_3) in frc.quadratic_angles: pass
        elif (assumed_3, assumed_2, assumed_1) in frc.quadratic_angles: pass
        # if it is not in .frc file log fail
        else:
            flagged_angle_coeff.add(i)
            log.warn(f'WARNING assumed angle elements {i[0]} {i[1]} {i[2]} does not have coeff type {assumed_1} {assumed_2} {assumed_3} is not in .frc file. Removing from read assumed types to avoid errors.')
            
    #################################
    # Check assumed dihedral coeffs #
    #################################
    flagged_dihedral_coeff = set([]);
    for i in aafc.dihedral_coeffs:
        assumed_1, assumed_2, assumed_3, assumed_4 = aafc.dihedral_coeffs[i]
        
        # Check if the type exists in .frc file
        if (assumed_1, assumed_2, assumed_3, assumed_4) in frc.torsion_3: pass
        elif (assumed_4, assumed_3, assumed_2, assumed_1) in frc.torsion_3: pass
        elif (assumed_1, assumed_2, assumed_3, assumed_4) in frc.torsion_1: pass
        elif (assumed_4, assumed_3, assumed_2, assumed_1) in frc.torsion_1: pass
        # if it is not in .frc file log fail
        else:
            flagged_dihedral_coeff.add(i)
            log.warn(f'WARNING assumed dihedral elements {i[0]} {i[1]} {i[2]} {i[3]} does not have coeff type {assumed_1} {assumed_2} {assumed_3} {assumed_4} is not in .frc file. Removing from read assumed types to avoid errors.')
            
    #################################
    # Check assumed improper coeffs #
    #################################
    flagged_improper_coeff = set([]);
    for i in aafc.improper_coeffs:
        assumed_1, assumed_2, assumed_3, assumed_4 = aafc.improper_coeffs[i]
        
        # Check if the type exists in .frc file
        if (assumed_1, assumed_2, assumed_3, assumed_4) in frc.wilson_out_of_plane: pass
        elif (assumed_4, assumed_3, assumed_2, assumed_1) in frc.wilson_out_of_plane: pass
        elif (assumed_1, assumed_2, assumed_3, assumed_4) in frc.wilson_out_of_plane_auto: pass
        elif (assumed_4, assumed_3, assumed_2, assumed_1) in frc.wilson_out_of_plane_auto: pass
        # if it is not in .frc file log fail
        else:
            flagged_improper_coeff.add(i)
            log.warn(f'WARNING assumed improper elements {i[0]} {i[1]} {i[2]} {i[3]} does not have coeff type {assumed_1} {assumed_2} {assumed_3} {assumed_4} is not in .frc file. Removing from read assumed types to avoid errors.')
    
    ##########################
    # Removing Flagged types #
    # from aafc key from     #
    # dictionary             #
    ##########################

    # Remove bond type coeff that is
    # not in file to stop errors
    for i in flagged_bond_coeff:
        del aafc.bond_coeffs[i]        
    
    # Remove angle type coeff that is
    # not in file to stop errors
    for i in flagged_angle_coeff:
        del aafc.angle_coeffs[i]  
    
    # Remove dihedral type coeff that is
    # not in file to stop errors
    for i in flagged_dihedral_coeff:
        del aafc.dihedral_coeffs[i]  
        
    # Remove improper type coeff that is
    # not in file to stop errors
    for i in flagged_improper_coeff:
        del aafc.improper_coeffs[i]  
        
    return aafc