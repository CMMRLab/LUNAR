# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
Feruary 26th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
from collections import OrderedDict
import os

#################################################
# Function used throughout the rest of the code #
#################################################
# Function to split coeff types into tuple
def split_coeff(types):
    # Find coeff types and split
    types = types.strip()
    types = types.split()
    types = tuple(types)
    return types

# Function for stringing together parameter types
def string_parameter_type(parameter_type):
    string = ''; str_buffer = 6 # Set string buffer size
    for n, i in enumerate(parameter_type):
        if n == 0: string += '{:<{str_buffer}}'.format(i, str_buffer=str_buffer)
        elif n < len(parameter_type)-1: string += '{:^{str_buffer}}'.format(i, str_buffer=str_buffer+2)
        else: string += '{:>{str_buffer}} {:^2}'.format(i, '', str_buffer=str_buffer)
    return string

# Function for stringing together float values for parameters
def string_coeffs(coeff):
    string = ''
    for i in coeff:
        string += '{:^10.4f}'.format(i)
    return string

# Function to find most frequent occurance in list
def most_frequent(List):
    counter = 0; num = List[0];
    for i in List:
        curr_frequency = List.count(i)
        if(curr_frequency> counter):
            counter = curr_frequency
            num = i
    return num

# Function to check for possbile bond_react_merge_prep tuple('N/A's) to set at bottom of lists
def move_NA2bottom(lst, NA_tuple, rm_flag=False):
    if NA_tuple in lst:
        lst.remove(NA_tuple)
        if not rm_flag: lst.append(NA_tuple)
    return lst


###########################################################################
# Function to update m class TypeIDs for atoms, bonds, angles, dihedrals, #
# and impropers and update all energy coeffs based on new merged class    #
###########################################################################
def update_TypeIDs(m, new, log):
    # Update atoms
    for i in m.atoms:
        new_type = 0   
        atom = m.atoms[i]
        str_type = m.pair_coeffs[atom.type].type
        new_type = int(new.atom_types_map[str_type])
        if new_type == 0: log.error(f'ERROR atomID {i} from file {m.filename} failed to update TypeID')
        atom.type = new_type
        
    # Update bonds
    for i in m.bonds:
        new_type = 0 
        bond = m.bonds[i]
        str_type = m.bond_coeffs[bond.type].type
        str_type_tuple = split_coeff(str_type)
        new_type = int(new.bond_types_map[str_type_tuple])
        if new_type == 0: log.error(f'ERROR bondID {i} from file {m.filename} failed to update TypeID')
        bond.type = new_type
        
    # Update angles
    for i in m.angles:
        new_type = 0 
        angle = m.angles[i]
        str_type = m.angle_coeffs[angle.type].type
        str_type_tuple = split_coeff(str_type)
        new_type = int(new.angle_types_map[str_type_tuple])
        if new_type == 0: log.error(f'ERROR angleID {i} from file {m.filename} failed to update TypeID')
        angle.type = new_type
        
    # Update dihedrals
    for i in m.dihedrals:
        new_type = 0 
        dihedral = m.dihedrals[i]
        str_type = m.dihedral_coeffs[dihedral.type].type
        str_type_tuple = split_coeff(str_type)
        new_type = int(new.dihedral_types_map[str_type_tuple])
        if new_type == 0: log.error(f'ERROR dihedralID {i} from file {m.filename} failed to update TypeID')
        dihedral.type = new_type
        
    # Update impropers
    for i in m.impropers:
        new_type = 0 
        improper = m.impropers[i]
        str_type = m.improper_coeffs[improper.type].type
        str_type_tuple = split_coeff(str_type)
        new_type = int(new.improper_types_map[str_type_tuple])
        if new_type == 0: log.error(f'ERROR improperID {i} from file {m.filename} failed to update TypeID')
        improper.type = new_type
        
    # Update energy coeffs
    m.natomtypes = new.natomtypes
    m.nbondtypes = new.nbondtypes
    m.nangletypes = new.nangletypes
    m.ndihedraltypes = new.ndihedraltypes
    m.nimpropertypes = new.nimpropertypes
    
    m.masses = new.masses
    m.pair_coeffs = new.pair_coeffs
    m.bond_coeffs = new.bond_coeffs
    m.angle_coeffs = new.angle_coeffs
    m.dihedral_coeffs = new.dihedral_coeffs
    m.improper_coeffs = new.improper_coeffs
    
    m.bondbond_coeffs = new.bondbond_coeffs
    m.bondangle_coeffs = new.bondangle_coeffs
    m.angleangletorsion_coeffs = new.angleangletorsion_coeffs
    m.endbondtorsion_coeffs = new.endbondtorsion_coeffs
    m.middlebondtorsion_coeffs = new.middlebondtorsion_coeffs
    m.bondbond13_coeffs = new.bondbond13_coeffs
    m.angletorsion_coeffs = new.angletorsion_coeffs
    m.angleangle_coeffs = new.angleangle_coeffs
    return m


################################################
# Class to merge all coeffs read into the code #
################################################
class Coeff_class: pass  # .type .coeffs = [] .consistency = set()
class merged:
    def __init__(self, merge, log):
        #############################################################
        # 1st loop through is to try to find all unique coeff types #
        ############################################################# 
        unique_pair_coeffs = set([]); unique_bond_coeffs = set([]); unique_angle_coeffs = set([])
        unique_dihedral_coeffs = set([]); unique_improper_coeffs = set([])
        
        # Unique style hints logger
        self.unique_style_hints = {'Masses': [], 'Pair_Coeffs': [], 'Bond_Coeffs': [],
                              'Angle_Coeffs': [], 'Dihedral_Coeffs': [],
                              'Improper_Coeffs': [], 'BondBond_Coeffs': [],
                              'BondAngle_Coeffs': [], 'AngleAngleTorsion_Coeffs': [],
                              'EndBondTorsion_Coeffs': [], 'MiddleBondTorsion_Coeffs': [],
                              'BondBond13_Coeffs': [], 'AngleTorsion_Coeffs': [], 'AngleAngle_Coeffs': []}
        
        # nparms logger to set default zeros for each coeff type (force field specific)
        nparms_in_coeff = {'pair':[], 'bond':[], 'angle':[], 'dihedral':[], 'improper':[]}
        
        # log to tally types
        ntypes_log = {'atom':0, 'bond':0, 'angle':0, 'dihedral':0, 'improper':0}
        log.out('-------------------------------------------------------------------------------------------------------------------------------------------')
        log.out('| {:^70} | {:^10} | {:^10} | {:^10} | {:^10} | {:^10} |'.format('filename', 'natom', 'nbond', 'nangle', 'ndihedral', 'improper'))
        log.out('| {:^70} | {:^10} | {:^10} | {:^10} | {:^10} | {:^10} |'.format('(Truncated to 70 characters)', 'types', 'types', 'types', 'types', 'types'))
        log.out('-------------------------------------------------------------------------------------------------------------------------------------------')
        
        # Function to split coeff types into tuple
        def split_coeff(types):
            # Find coeff types and split
            types = types.strip()
            types = types.split()
            types = tuple(types)
            return types
        
        # Function to check comments are in all2lmp style
        def check_comments(types, ntypes, filename, section, log):
            if len(split_coeff(types)) != ntypes or types == 'N/A':
                log.error(f'ERROR {section} comment: {types} in {os.path.basename(filename)} are not in all2lmp.py or bond_react_merge_prep.py style, please address')
            return    
        for i in merge:
            file = merge[i]
        
            # find values/log values and print to table
            filename = file.filename;   natomtypes = file.natomtypes;
            nbondtypes = file.nbondtypes;   nangletypes = file.nangletypes;
            ndihedraltypes = file.ndihedraltypes;   nimpropertypes = file.nimpropertypes;
            ntypes_log['atom'] += natomtypes;   ntypes_log['bond'] += nbondtypes;
            ntypes_log['angle'] += nangletypes;  ntypes_log['dihedral'] += ndihedraltypes;
            ntypes_log['improper'] += nimpropertypes;
            log.out('| {:^70} | {:^10} | {:^10} | {:^10} | {:^10} | {:^10} |'.format(os.path.basename(filename)[0:70], natomtypes, nbondtypes, nangletypes, ndihedraltypes, nimpropertypes))
            
            # Find unique pair coeffs
            pair_coeffs = file.pair_coeffs
            for i in pair_coeffs:
                coeff = pair_coeffs[i]
                types = coeff.type
                unique_pair_coeffs.add(types)
                nparms_in_coeff['pair'].append(len(coeff.coeffs))
                #check_comments(coeff.type, 1, file.filename, 'Pair Coeffs', log)
                #check_comments(file.masses[i].type, 1, file.filename, 'Masses', log)
                
            # Find unique bond coeffs
            bond_coeffs = file.bond_coeffs
            for i in bond_coeffs:
                coeff = bond_coeffs[i]
                types = split_coeff(coeff.type)
                unique_bond_coeffs.add(types)
                nparms_in_coeff['bond'].append(len(coeff.coeffs))
                check_comments(coeff.type, 2, file.filename, 'Bond Coeffs', log)
                
            # Find unique angle coeffs 
            angle_coeffs = file.angle_coeffs
            for i in angle_coeffs:
                coeff = angle_coeffs[i]
                types = split_coeff(coeff.type)
                unique_angle_coeffs.add(types)
                nparms_in_coeff['angle'].append(len(coeff.coeffs))
                check_comments(coeff.type, 3, file.filename, 'Angle Coeffs', log)
                
                try: # Check crossterms
                    check_comments(file.bondbond_coeffs[i].type, 3, file.filename, 'BondBond Coeffs', log)
                    check_comments(file.bondangle_coeffs[i].type, 3, file.filename, 'BondAngle Coeffs', log)
                except: pass
                    
            # Find unique dihedral coeffs
            dihedral_coeffs = file.dihedral_coeffs
            for i in dihedral_coeffs:
                coeff = dihedral_coeffs[i]
                types = split_coeff(coeff.type)
                unique_dihedral_coeffs.add(types)
                nparms_in_coeff['dihedral'].append(len(coeff.coeffs))
                check_comments(coeff.type, 4, file.filename, 'Dihedral Coeffs', log)
                
                try: # Check crossterms
                    check_comments(file.angleangletorsion_coeffs[i].type, 4, file.filename, 'AngleAngleTorsion Coeffs', log)
                    check_comments(file.angleangletorsion_coeffs[i].type, 4, file.filename, 'AngleAngleTorsion Coeffs', log)
                    check_comments(file.endbondtorsion_coeffs[i].type, 4, file.filename, 'EndBondTorsion Coeffs', log)
                    check_comments(file.middlebondtorsion_coeffs[i].type, 4, file.filename, 'MiddleBondTorsion Coeffs', log)
                    check_comments(file.bondbond13_coeffs[i].type, 4, file.filename, 'BondBond13 Coeffs', log)
                    check_comments(file.angletorsion_coeffs[i].type, 4, file.filename, 'AngleTorsion Coeffs', log)
                except: pass
                
            # Find unique improper coeffs 
            improper_coeffs = file.improper_coeffs
            for i in improper_coeffs:
                coeff = improper_coeffs[i]
                types = split_coeff(coeff.type)
                unique_improper_coeffs.add(types)
                nparms_in_coeff['improper'].append(len(coeff.coeffs))
                check_comments(coeff.type, 5, file.filename, 'Improper Coeffs', log)
                
                try: # Check crossterms
                    check_comments(file.angleangle_coeffs[i].type, 5, file.filename, 'AngleAngle Coeffs', log)
                except: pass
    
            # Find unique style hints 
            self.unique_style_hints['Masses'].append(file.mass_coeffs_style_hint)
            self.unique_style_hints['Pair_Coeffs'].append(file.pair_coeffs_style_hint)
            self.unique_style_hints['Bond_Coeffs'].append(file.bond_coeffs_style_hint)
            self.unique_style_hints['Angle_Coeffs'].append(file.angle_coeffs_style_hint)
            self.unique_style_hints['Dihedral_Coeffs'].append(file.dihedral_coeffs_style_hint)
            self.unique_style_hints['Improper_Coeffs'].append(file.improper_coeffs_style_hint)
            self.unique_style_hints['BondBond_Coeffs'].append(file.bondbond_coeffs_style_hint)
            self.unique_style_hints['BondAngle_Coeffs'].append(file.bondangle_coeffs_style_hint)
            self.unique_style_hints['AngleAngleTorsion_Coeffs'].append(file.angleangletorsion_coeffs_style_hint)
            self.unique_style_hints['EndBondTorsion_Coeffs'].append(file.endbondtorsion_coeffs_style_hint)
            self.unique_style_hints['MiddleBondTorsion_Coeffs'].append(file.middlebondtorsion_coeffs_style_hint)
            self.unique_style_hints['BondBond13_Coeffs'].append(file.bondbond13_coeffs_style_hint)
            self.unique_style_hints['AngleTorsion_Coeffs'].append(file.angletorsion_coeffs_style_hint)
            self.unique_style_hints['AngleAngle_Coeffs'].append(file.angleangle_coeffs_style_hint)
            
                
        # print ending of table from intialization
        log.out('-------------------------------------------------------------------------------------------------------------------------------------------')
        log.out('| {:^70} | {:^10} | {:^10} | {:^10} | {:^10} | {:^10} |'.format('Total possible sum: ', ntypes_log['atom'], ntypes_log['bond'], ntypes_log['angle'], ntypes_log['dihedral'], ntypes_log['improper']) )
        log.out('-------------------------------------------------------------------------------------------------------------------------------------------')
        log.out('| {:^70} | {:^10} | {:^10} | {:^10} | {:^10} | {:^10} |'.format('Total unique sum: ', len(unique_pair_coeffs), len(unique_bond_coeffs), len(unique_angle_coeffs), len(unique_dihedral_coeffs), len(unique_improper_coeffs)) )
        log.out('-------------------------------------------------------------------------------------------------------------------------------------------')
        
        
        ##################################################################################
        # Loop through style hints to make sure they are all consistant across all files #
        ##################################################################################
        for coeff in self.unique_style_hints:
            hints = self.unique_style_hints[coeff]
            no_dups = set(hints)
            if 'N/A' not in hints and len(no_dups) > 1:
                log.warn(f'WARNING style hints between some files are different for {coeff} Coeff. Style Hints {str(no_dups)}')
        
        ##########################################################
        # Internally derive an ff_class variable based on number #
        # of parms in bond, angle, dihedral, and improper coeffs #
        ##########################################################
        self.ff_class = 1
        nparms_bond = most_frequent(nparms_in_coeff['bond'])
        nparms_angle = most_frequent(nparms_in_coeff['angle'])
        nparms_dihedral = most_frequent(nparms_in_coeff['dihedral'])
        nparms_improper = most_frequent(nparms_in_coeff['improper'])
        if nparms_bond == 4 and nparms_angle == 4 and nparms_dihedral == 6 and nparms_improper == 2:
            self.ff_class = 2
            log.out('Internally determined force field is class2')
        
        
        #######################################################################
        # Sort each coeff types alphabetically so each coeff is sorted nicely #
        # (all2lmp already sorted each individual type as best as possible)   #
        #######################################################################
        # Atom types
        atom_types_lst = sorted(list(unique_pair_coeffs)); atom_types_dict = {};
        atom_types_lst = move_NA2bottom(atom_types_lst, 'N/A', rm_flag=False)
        
        # Bonds types
        bond_types_lst = list(unique_bond_coeffs); bond_types_dict = {};
        bond_types_lst = sorted(unique_bond_coeffs, key=lambda x: x[1])
        bond_types_lst = sorted(unique_bond_coeffs, key=lambda x: x[0])
        bond_types_lst = move_NA2bottom(bond_types_lst, ('N/A', 'N/A'), rm_flag=False)
        
        # Angles types
        angle_types_lst = list(unique_angle_coeffs); angle_types_dict = {};
        angle_types_lst = sorted(angle_types_lst, key=lambda x: x[2])
        angle_types_lst = sorted(angle_types_lst, key=lambda x: x[1])
        angle_types_lst = sorted(angle_types_lst, key=lambda x: x[0])
        angle_types_lst = move_NA2bottom(angle_types_lst, ('N/A', 'N/A', 'N/A'), rm_flag=False)
        
        # Dihedral types
        dihedral_types_lst = list(unique_dihedral_coeffs); dihedral_types_dict = {};
        dihedral_types_lst = sorted(dihedral_types_lst, key=lambda x: x[3])
        dihedral_types_lst = sorted(dihedral_types_lst, key=lambda x: x[2])
        dihedral_types_lst = sorted(dihedral_types_lst, key=lambda x: x[1])
        dihedral_types_lst = sorted(dihedral_types_lst, key=lambda x: x[0])
        dihedral_types_lst = move_NA2bottom(dihedral_types_lst, ('N/A', 'N/A', 'N/A', 'N/A'), rm_flag=False)
        
        # Improper types Sort by 4th index last since this will be nb==3 or nb!=3 which is 
        # for distinguish between improper coeffs (nb!==3) and angleangle coeffs (nb==3)
        improper_types_lst = list(unique_improper_coeffs); improper_types_dict = {};
        improper_types_lst = sorted(improper_types_lst, key=lambda x: x[3])
        improper_types_lst = sorted(improper_types_lst, key=lambda x: x[2])
        improper_types_lst = sorted(improper_types_lst, key=lambda x: x[1])
        improper_types_lst = sorted(improper_types_lst, key=lambda x: x[0])
        improper_types_lst = sorted(improper_types_lst, key=lambda x: x[4], reverse=True)
        improper_types_lst = move_NA2bottom(improper_types_lst, ('N/A', 'N/A', 'N/A', 'N/A', 'nb!=3'), rm_flag=False)
        
        
        ###########################
        # Set new coeff numbering #
        ###########################      
        for n, i in enumerate(atom_types_lst, 1):
            atom_types_dict[i] = n       
        for n, i in enumerate(bond_types_lst, 1):
            bond_types_dict[i] = n
        for n, i in enumerate(angle_types_lst, 1):
            angle_types_dict[i] = n
        for n, i in enumerate(dihedral_types_lst, 1):
            dihedral_types_dict[i] = n
        for n, i in enumerate(improper_types_lst, 1):
            improper_types_dict[i] = n

        
        #################################################################################################################################
        # Set maps to map convert from string type to newly found numeric type to use on the fly for conversion when writing new files. #
        # Setting them as instances to accesses via merged class to only have to read pass one data struct to file writer functions.    #
        #################################################################################################################################
        self.atom_types_map = atom_types_dict # { string type : equivalent numeric type }
        self.bond_types_map = bond_types_dict # { tuple(string types) : equivalent numeric type }
        self.angle_types_map = angle_types_dict # { tuple(string types) : equivalent numeric type }
        self.dihedral_types_map = dihedral_types_dict # { tuple(string types) : equivalent numeric type }
        self.improper_types_map = improper_types_dict # { tuple(string types) : equivalent numeric type }
        
        
        ############################################################
        # Set ntypes for quick and easy file writing of quantities #
        ############################################################
        self.natomtypes = len(self.atom_types_map)
        self.nbondtypes = len(self.bond_types_map)
        self.nangletypes = len(self.angle_types_map)
        self.ndihedraltypes = len(self.dihedral_types_map)
        self.nimpropertypes = len(self.improper_types_map)
        
        
        ###########################################
        # New coeffs with proper numeric type ids #
        ###########################################
        self.masses = {} # { atom type id : mass }
        self.pair_coeffs = {} # { atom type id : coeffs class }
        self.bond_coeffs = {} # { bond type id : coeffs class }
        self.angle_coeffs = {} # { angle type id : coeffs class }
        self.dihedral_coeffs = {} # { dihedral type id : coeffs class }
        self.improper_coeffs = {} # { improper type id : coeffs class }
        self.bondbond_coeffs = {} # { bondbond type id : coeffs class }
        self.bondangle_coeffs = {} # { bondangle type id : coeffs class }
        self.angleangle_coeffs = {} # { angleangle type id : coeffs class }
        self.angleangletorsion_coeffs = {} # { angleangletorsion type id : coeffs class }
        self.endbondtorsion_coeffs = {} # { endbondtorsion type id : coeffs class }
        self.middlebondtorsion_coeffs = {} # { middlebondtorsion type id : coeffs class }
        self.bondbond13_coeffs = {} # { bondbond13 type id : coeffs class }
        self.angletorsion_coeffs = {} # { angletorsion type id : coeffs class }
        
        
        ########################
        # Style Hints Transfer #
        ########################
        self.mass_coeffs_style_hint = most_frequent(self.unique_style_hints['Masses'])
        self.pair_coeffs_style_hint = most_frequent(self.unique_style_hints['Bond_Coeffs'])
        self.bond_coeffs_style_hint = most_frequent(self.unique_style_hints['Bond_Coeffs'])
        self.angle_coeffs_style_hint = most_frequent(self.unique_style_hints['Angle_Coeffs'])
        self.dihedral_coeffs_style_hint = most_frequent(self.unique_style_hints['Dihedral_Coeffs'])
        self.improper_coeffs_style_hint = most_frequent(self.unique_style_hints['Improper_Coeffs'])
        self.bondbond_coeffs_style_hint = most_frequent(self.unique_style_hints['BondBond_Coeffs'])
        self.bondangle_coeffs_style_hint = most_frequent(self.unique_style_hints['BondAngle_Coeffs'])
        self.angleangle_coeffs_style_hint = most_frequent(self.unique_style_hints['AngleAngle_Coeffs'])
        self.angleangletorsion_coeffs_style_hint = most_frequent(self.unique_style_hints['AngleAngleTorsion_Coeffs'])
        self.endbondtorsion_coeffs_style_hint = most_frequent(self.unique_style_hints['EndBondTorsion_Coeffs'])
        self.middlebondtorsion_coeffs_style_hint = most_frequent(self.unique_style_hints['MiddleBondTorsion_Coeffs'])
        self.bondbond13_coeffs_style_hint = most_frequent(self.unique_style_hints['BondBond13_Coeffs'])
        self.angletorsion_coeffs_style_hint = most_frequent(self.unique_style_hints['AngleTorsion_Coeffs'])
        
        
        ####################################################
        # Intialize all dicts as zeros and update later on #
        ####################################################
        # Intialize all atom types dicts
        for i in atom_types_lst:
            c = Coeff_class()
            c.coeffs = 0; c.type = i; c.consistency = set();
            self.masses[atom_types_dict[i]] = c
            
            c1 = Coeff_class()
            c1.coeffs = [0,0]; c1.type = i; c1.consistency = set();
            self.pair_coeffs[atom_types_dict[i]] = c1
        
        # Intialize all bond types dicts
        for i in bond_types_lst:
            nparms = most_frequent(nparms_in_coeff['bond'])             
            c2 = Coeff_class()
            c2.coeffs = nparms*[0]; c2.type = string_parameter_type(i); c2.consistency = set();
            self.bond_coeffs[bond_types_dict[i]] = c2

        # Intialize all angle types dicts            
        for i in angle_types_lst:
            nparms = most_frequent(nparms_in_coeff['angle'])  
            c3 = Coeff_class()
            c3.coeffs = nparms*[0]; c3.type = string_parameter_type(i); c3.consistency = set();
            self.angle_coeffs[angle_types_dict[i]] = c3
            
            c4 = Coeff_class()
            c4.coeffs = [0,0,0]; c4.type = string_parameter_type(i); c4.consistency = set();
            self.bondbond_coeffs[angle_types_dict[i]] = c4
            
            c5 = Coeff_class()
            c5.coeffs = [0,0,0,0]; c5.type = string_parameter_type(i); c5.consistency = set();
            self.bondangle_coeffs[angle_types_dict[i]] = c5
         
        # Intialize all dihedral types dicts
        for i in dihedral_types_lst:
            nparms = most_frequent(nparms_in_coeff['dihedral'])  
            c6 = Coeff_class()
            c6.coeffs = nparms*[0]; c6.type = string_parameter_type(i); c6.consistency = set();
            self.dihedral_coeffs[dihedral_types_dict[i]] = c6
            
            c7 = Coeff_class()
            c7.coeffs = [0,0,0]; c7.type = string_parameter_type(i); c7.consistency = set();
            self.angleangletorsion_coeffs[dihedral_types_dict[i]] = c7
            
            c8 = Coeff_class()
            c8.coeffs = [0,0,0,0,0,0,0]; c8.type = string_parameter_type(i); c8.consistency = set();
            self.endbondtorsion_coeffs[dihedral_types_dict[i]] = c8
            
            c9 = Coeff_class()
            c9.coeffs = [0,0,0,0,0,0,0]; c9.type = string_parameter_type(i); c9.consistency = set();
            self.middlebondtorsion_coeffs[dihedral_types_dict[i]] = c9
            
            c10 = Coeff_class()
            c10.coeffs = [0,0,0]; c10.type = string_parameter_type(i); c10.consistency = set();
            self.bondbond13_coeffs[dihedral_types_dict[i]] = c10
            
            c11 = Coeff_class()
            c11.coeffs = [0,0,0,0,0,0,0,0]; c11.type = string_parameter_type(i); c11.consistency = set();
            self.angletorsion_coeffs[dihedral_types_dict[i]] = c11

        # Intialize all improper types dicts
        for i in improper_types_lst:
            nparms = most_frequent(nparms_in_coeff['improper'])  
            c12 = Coeff_class()
            c12.coeffs = nparms*[0]; c12.type = string_parameter_type(i); c12.consistency = set();
            self.improper_coeffs[improper_types_dict[i]] = c12
            
            c13 = Coeff_class()
            c13.coeffs = [0,0,0,0,0,0]; c13.type = string_parameter_type(i); c13.consistency = set();
            self.angleangle_coeffs[improper_types_dict[i]] = c13

        
        ################################################
        # Loop through all files and merge coeff types #
        ################################################
        for i in merge:
            file = merge[i]
            filename = os.path.basename(file.filename)
            
            ###############
            # Find masses #
            ###############
            for j in file.masses:
                coeff = file.masses[j]
                mass = coeff.coeffs
                types = coeff.type
                new_type = atom_types_dict[types]
                self.masses[new_type].coeffs = mass
                
                # check for consistency
                self.masses[new_type].consistency.add(tuple(mass))
                if len(self.masses[new_type].consistency) > 1:
                    log.error('\nERROR Masses {:<4} parameter {:<6}\nfrom {} differs from other read in files.'.format(types, mass, filename))

                        
            ####################
            # Find pair coeffs #
            ####################
            for j in file.pair_coeffs:
                coeff = file.pair_coeffs[j]
                coeffs = coeff.coeffs
                types = coeff.type
                new_type = atom_types_dict[types]
                self.pair_coeffs[new_type].coeffs = coeffs
                
                # check for consistency
                self.pair_coeffs[new_type].consistency.add(tuple(coeffs))
                if len(self.pair_coeffs[new_type].consistency) > 1:
                    log.error('\nERROR Pair Coeffs {:<4} parameters {}\nfrom {} differs from other read in files.'.format(types, string_coeffs(coeffs), filename))
                    
                    
            ####################
            # Find bond coeffs #
            ####################
            for j in file.bond_coeffs:
                coeff = file.bond_coeffs[j]
                coeffs = coeff.coeffs
                types = split_coeff(coeff.type)
                new_type = bond_types_dict[types]
                self.bond_coeffs[new_type].coeffs = coeffs
                
                # check for consistency
                self.bond_coeffs[new_type].consistency.add(tuple(coeffs))
                if len(self.bond_coeffs[new_type].consistency) > 1:
                    log.error('\nERROR Bond Coeffs {:<16} parameters {}\nfrom {} differs from other read in files.'.format(coeff.type, string_coeffs(coeffs), filename))


            #####################
            # Find angle coeffs #
            #####################
            for j in file.angle_coeffs:
                coeff = file.angle_coeffs[j]
                coeffs = coeff.coeffs
                types = split_coeff(coeff.type)
                new_type = angle_types_dict[types]
                self.angle_coeffs[new_type].coeffs = coeffs

                # check for consistency
                self.angle_coeffs[new_type].consistency.add(tuple(coeffs))
                if len(self.angle_coeffs[new_type].consistency) > 1:
                    log.error('\nERROR Angle Coeffs {:<16} parameters {}\nfrom {} differs from other read in files.'.format(coeff.type, string_coeffs(coeffs), filename))

            
            ########################
            # Find dihedral coeffs #
            ########################
            for j in file.dihedral_coeffs:
                coeff = file.dihedral_coeffs[j]
                coeffs = coeff.coeffs
                types = split_coeff(coeff.type)
                new_type = dihedral_types_dict[types]
                self.dihedral_coeffs[new_type].coeffs = coeffs

                # check for consistency
                self.dihedral_coeffs[new_type].consistency.add(tuple(coeffs))
                if len(self.dihedral_coeffs[new_type].consistency) > 1:
                    log.error('\nERROR Dihedral Coeffs {:<16} parameters\n{}\nfrom {} differs from other read in files.'.format(coeff.type, string_coeffs(coeffs), filename))

                    
            ########################
            # Find improper coeffs #
            ########################
            for j in file.improper_coeffs:
                coeff = file.improper_coeffs[j]
                coeffs = coeff.coeffs
                types = split_coeff(coeff.type)
                new_type = improper_types_dict[types]
                self.improper_coeffs[new_type].coeffs = coeffs

                # check for consistency
                self.improper_coeffs[new_type].consistency.add(tuple(coeffs))
                if len(self.improper_coeffs[new_type].consistency) > 1:
                    log.error('\nERROR Improper Coeffs {:<16} parameters {}\nfrom {} differs from other read in files.'.format(coeff.type, string_coeffs(coeffs), filename))

            
            ########################
            # Find bondbond coeffs #
            ########################
            if file.bondbond_coeffs:
                for j in file.bondbond_coeffs:
                    coeff = file.bondbond_coeffs[j]
                    coeffs = coeff.coeffs
                    types = split_coeff(coeff.type)
                    new_type = angle_types_dict[types]
                    self.bondbond_coeffs[new_type].coeffs = coeffs
                        
                    # check for consistency
                    self.bondbond_coeffs[new_type].consistency.add(tuple(coeffs))
                    if len(self.bondbond_coeffs[new_type].consistency) > 1:
                        log.error('\nERROR BondBond Coeffs {:<16} parameters {}\nfrom {} differs from other read in files.'.format(coeff.type, string_coeffs(coeffs), filename))

                        
            #########################
            # Find bondangle coeffs #
            #########################
            if file.bondangle_coeffs:
                for j in file.bondangle_coeffs:
                    coeff = file.bondangle_coeffs[j]
                    coeffs = coeff.coeffs
                    types = split_coeff(coeff.type)
                    new_type = angle_types_dict[types]
                    self.bondangle_coeffs[new_type].coeffs = coeffs
                    
                    # check for consistency
                    self.bondangle_coeffs[new_type].consistency.add(tuple(coeffs))
                    if len(self.bondangle_coeffs[new_type].consistency) > 1:
                        log.error('\nERROR BondAngle Coeffs {:<16} parameters {}\nfrom {} differs from other read in files.'.format(coeff.type, string_coeffs(coeffs), filename))
    
                        
            #################################
            # Find angleangletorsion coeffs #
            #################################
            if file.angleangletorsion_coeffs:
                for j in file.angleangletorsion_coeffs:
                    coeff = file.angleangletorsion_coeffs[j]
                    coeffs = coeff.coeffs
                    types = split_coeff(coeff.type)
                    new_type = dihedral_types_dict[types]
                    self.angleangletorsion_coeffs[new_type].coeffs = coeffs
                    
                    # check for consistency
                    self.angleangletorsion_coeffs[new_type].consistency.add(tuple(coeffs))
                    if len(self.angleangletorsion_coeffs[new_type].consistency) > 1:
                        log.error('\nERROR AngleAngleTorsion Coeffs {:<16} parameters {}\nfrom {} differs from other read in files.'.format(coeff.type, string_coeffs(coeffs), filename))
                        

            ##############################
            # Find endbondtorsion coeffs #
            ##############################
            if file.endbondtorsion_coeffs:
                for j in file.endbondtorsion_coeffs:
                    coeff = file.endbondtorsion_coeffs[j]
                    coeffs = coeff.coeffs
                    types = split_coeff(coeff.type)
                    new_type = dihedral_types_dict[types]
                    self.endbondtorsion_coeffs[new_type].coeffs = coeffs
    
                    # check for consistency
                    self.endbondtorsion_coeffs[new_type].consistency.add(tuple(coeffs))
                    if len(self.endbondtorsion_coeffs[new_type].consistency) > 1:
                        log.error('\nERROR EndBondTorsion Coeffs {:<16} parameters\n{}\nfrom {} differs from other read in files.'.format(coeff.type, string_coeffs(coeffs), filename))
                    
                    
            #################################
            # Find middlebondtorsion coeffs #
            #################################
            if file.middlebondtorsion_coeffs:
                for j in file.middlebondtorsion_coeffs:
                    coeff = file.middlebondtorsion_coeffs[j]
                    coeffs = coeff.coeffs
                    types = split_coeff(coeff.type)
                    new_type = dihedral_types_dict[types]
                    self.middlebondtorsion_coeffs[new_type].coeffs = coeffs
                    
                    # check for consistency
                    self.middlebondtorsion_coeffs[new_type].consistency.add(tuple(coeffs))
                    if len(self.middlebondtorsion_coeffs[new_type].consistency) > 1:
                        log.error('\nERROR MiddleBondTorsion Coeffs {:<16} parameters\n{}\nfrom {} differs from other read in files.'.format(coeff.type, string_coeffs(coeffs), filename))            

                    
            ##########################
            # Find bondbond13 coeffs #
            ##########################
            if file.bondbond13_coeffs:
                for j in file.bondbond13_coeffs:
                    coeff = file.bondbond13_coeffs[j]
                    coeffs = coeff.coeffs
                    types = split_coeff(coeff.type)
                    new_type = dihedral_types_dict[types]
                    self.bondbond13_coeffs[new_type].coeffs = coeffs
                                        
                    # check for consistency
                    self.bondbond13_coeffs[new_type].consistency.add(tuple(coeffs))
                    if len(self.bondbond13_coeffs[new_type].consistency) > 1:
                        log.error('\nERROR BondBond13 Coeffs {:<16} parameters\n{}\nfrom {} differs from other read in files.'.format(coeff.type, string_coeffs(coeffs), filename))

                    
            ############################
            # Find angletorsion coeffs #
            ############################
            if file.angletorsion_coeffs:
                for j in file.angletorsion_coeffs:
                    coeff = file.angletorsion_coeffs[j]
                    coeffs = coeff.coeffs
                    types = split_coeff(coeff.type)
                    new_type = dihedral_types_dict[types]
                    self.angletorsion_coeffs[new_type].coeffs = coeffs
                    
                    # check for consistency
                    self.angletorsion_coeffs[new_type].consistency.add(tuple(coeffs))
                    if len(self.angletorsion_coeffs[new_type].consistency) > 1:
                        log.error('\nERROR AngleTorsion Coeffs {:<16} parameters\n{}\nfrom {} differs from other read in files.'.format(coeff.type, string_coeffs(coeffs), filename))
                    
                        
            ##########################
            # Find angleangle coeffs #
            ##########################
            if file.angleangle_coeffs:
                for j in file.angleangle_coeffs:
                    coeff = file.angleangle_coeffs[j]
                    coeffs = coeff.coeffs
                    types = split_coeff(coeff.type)
                    new_type = improper_types_dict[types]
                    self.angleangle_coeffs[new_type].coeffs = coeffs
                    
                    # check for consistency
                    self.angleangle_coeffs[new_type].consistency.add(tuple(coeffs))
                    if len(self.angleangle_coeffs[new_type].consistency) > 1:
                        log.error('\nERROR AngleAngle Coeffs {:<16} parameters\n{}\nfrom {} differs from other read in files.'.format(coeff.type, string_coeffs(coeffs), filename))
                        

        ###########################            
        # Order all created dicts #
        ###########################
        self.masses = dict(OrderedDict(sorted(self.masses.items())))
        self.pair_coeffs = dict(OrderedDict(sorted(self.pair_coeffs.items())))
        self.bond_coeffs = dict(OrderedDict(sorted(self.bond_coeffs.items())))
        self.angle_coeffs = dict(OrderedDict(sorted(self.angle_coeffs.items())))
        self.dihedral_coeffs = dict(OrderedDict(sorted(self.dihedral_coeffs.items())))
        self.improper_coeffs = dict(OrderedDict(sorted(self.improper_coeffs.items())))
        try:
            self.bondbond_coeffs = dict(OrderedDict(sorted(self.bondbond_coeffs.items())))
            self.bondangle_coeffs = dict(OrderedDict(sorted(self.bondangle_coeffs.items())))
            self.angleangletorsion_coeffs = dict(OrderedDict(sorted(self.angleangletorsion_coeffs.items())))
            self.endbondtorsion_coeffs = dict(OrderedDict(sorted(self.endbondtorsion_coeffs.items())))
            self.middlebondtorsion_coeffs = dict(OrderedDict(sorted(self.middlebondtorsion_coeffs.items())))
            self.bondbond13_coeffs = dict(OrderedDict(sorted(self.bondbond13_coeffs.items())))
            self.angletorsion_coeffs = dict(OrderedDict(sorted(self.angletorsion_coeffs.items())))
            self.angleangle_coeffs = dict(OrderedDict(sorted(self.angleangle_coeffs.items())))
        except: pass


# Function to printout merged coeffs
def print_merged(new, log, skip_printing_cross_terms):
    # Print masses
    log.debug('\n\n------------------------------------')
    log.debug('|           Unique Masses          |')
    log.debug('------------------------------------')
    log.debug('| {:^2} | {:^8} | {:^16} |'.format('Id', 'Type', 'Coeffs'))
    log.debug('------------------------------------')
    for i in new.masses:
        coeff = new.masses[i]
        log.debug('| {:^2} | {:^8} | {:^16} |'.format(i, coeff.type, coeff.coeffs[0]))
    log.debug('------------------------------------')
    
    # Print Pair coeffs
    log.debug('\n----------------------------------------')
    log.debug('|           Unique Pair Coeffs         |')
    log.debug('----------------------------------------')
    log.debug('| {:^2} | {:^8} | {:^20} |'.format('Id', 'Type', 'Coeffs'))
    log.debug('----------------------------------------')
    for i in new.pair_coeffs:
        coeff = new.pair_coeffs[i]
        log.debug('| {:^2} | {:^8} | {:^20} |'.format(i, coeff.type, string_coeffs(coeff.coeffs) ) )
    log.debug('----------------------------------------')
    
    # Print Bond coeffs
    log.debug('\n-------------------------------------------------------------------------------')
    log.debug('|                              Unique Bond Coeffs                             |')
    log.debug('-------------------------------------------------------------------------------')
    log.debug('| {:^2} | {:^17} | {:^50} |'.format('Id', 'Type', 'Coeffs'))
    log.debug('-------------------------------------------------------------------------------')
    for i in new.bond_coeffs:
        coeff = new.bond_coeffs[i]
        id1, id2 = split_coeff(coeff.type)
        log.debug('| {:^2} | {:^8} {:^8} | {:^50} |'.format(i, id1, id2, string_coeffs(coeff.coeffs) ) )
    log.debug('-------------------------------------------------------------------------------')
    
    # Print Angle coeffs
    log.debug('\n----------------------------------------------------------------------------------------')
    log.debug('|                                  Unique Angle Coeffs                                 |')
    log.debug('----------------------------------------------------------------------------------------')
    log.debug('| {:^2} | {:^26} | {:^50} |'.format('Id', 'Type', 'Coeffs'))
    log.debug('----------------------------------------------------------------------------------------')
    for i in new.angle_coeffs:
        coeff = new.angle_coeffs[i]
        id1, id2, id3 = split_coeff(coeff.type)
        log.debug('| {:^2} | {:^8} {:^8} {:^8} | {:^50} |'.format(i, id1, id2, id3, string_coeffs(coeff.coeffs) ) )
    log.debug('----------------------------------------------------------------------------------------')
    
    # Print Dihedral coeffs
    log.debug('\n-----------------------------------------------------------------------------------------------------------')
    log.debug('|                                          Unique Dihedral Coeffs                                         |')
    log.debug('-----------------------------------------------------------------------------------------------------------')
    log.debug('| {:^2} | {:^35} | {:^60} |'.format('Id', 'Type', 'Coeffs'))
    log.debug('-----------------------------------------------------------------------------------------------------------')
    for i in new.dihedral_coeffs:
        coeff = new.dihedral_coeffs[i]
        id1, id2, id3, id4 = split_coeff(coeff.type)
        log.debug('| {:^2} | {:^8} {:^8} {:^8} {:^8} | {:^60} |'.format(i, id1, id2, id3, id4, string_coeffs(coeff.coeffs) ) )
    log.debug('-----------------------------------------------------------------------------------------------------------')
    
    # Print Improper coeffs
    log.debug('\n------------------------------------------------------------------------------------------')
    log.debug('|                                Unique Improper Coeffs                                  |')
    log.debug('------------------------------------------------------------------------------------------')
    log.debug('| {:^2} | {:^44} | {:^34} |'.format('Id', 'Type', 'Coeffs'))
    log.debug('------------------------------------------------------------------------------------------')
    for i in new.improper_coeffs:
        coeff = new.improper_coeffs[i]
        id1, id2, id3, id4, id5 = split_coeff(coeff.type)
        log.debug('| {:^2} | {:^8} {:^8} {:^8} {:^8} {:^8} | {:^34} |'.format(i, id1, id2, id3, id4, id5, string_coeffs(coeff.coeffs) ) )
    log.debug('------------------------------------------------------------------------------------------')
    
    # Option to skip crossterm prinouts
    if not skip_printing_cross_terms and new.ff_class == 2:
        # Print Bondbond coeffs
        log.debug('\n-------------------------------------------------------------------------')
        log.debug('|                         Unique Bondbond Coeffs                        |')
        log.debug('-------------------------------------------------------------------------')
        log.debug('| {:^2} | {:^26} | {:^35} |'.format('Id', 'Type', 'Coeffs'))
        log.debug('-------------------------------------------------------------------------')
        for i in new.bondbond_coeffs:
            coeff = new.bondbond_coeffs[i]
            id1, id2, id3 = split_coeff(coeff.type)
            log.debug('| {:^2} | {:^8} {:^8} {:^8} | {:^35} |'.format(i, id1, id2, id3, string_coeffs(coeff.coeffs) ) )
        log.debug('-------------------------------------------------------------------------')
        
        # Print Bondangle coeffs
        log.debug('\n----------------------------------------------------------------------------------------')
        log.debug('|                                Unique Bondangle Coeffs                               |')
        log.debug('----------------------------------------------------------------------------------------')
        log.debug('| {:^2} | {:^26} | {:^50} |'.format('Id', 'Type', 'Coeffs'))
        log.debug('----------------------------------------------------------------------------------------')
        for i in new.bondangle_coeffs:
            coeff = new.bondangle_coeffs[i]
            id1, id2, id3 = split_coeff(coeff.type)
            log.debug('| {:^2} | {:^8} {:^8} {:^8} | {:^50} |'.format(i, id1, id2, id3, string_coeffs(coeff.coeffs) ) )
        log.debug('----------------------------------------------------------------------------------------')
        
        # Print Angleangletorsion coeffs
        log.debug('\n---------------------------------------------------------------------------------------')
        log.debug('|                           Unique Angleangletorsion Coeffs                           |')
        log.debug('---------------------------------------------------------------------------------------')
        log.debug('| {:^2} | {:^35} | {:^40} |'.format('Id', 'Type', 'Coeffs'))
        log.debug('---------------------------------------------------------------------------------------')
        for i in new.angleangletorsion_coeffs:
            coeff = new.angleangletorsion_coeffs[i]
            id1, id2, id3, id4 = split_coeff(coeff.type)
            log.debug('| {:^2} | {:^8} {:^8} {:^8} {:^8} | {:^40} |'.format(i, id1, id2, id3, id4, string_coeffs(coeff.coeffs) ) )
        log.debug('---------------------------------------------------------------------------------------')
        
        # Print Endbondtorsion coeffs
        log.debug('\n----------------------------------------------------------------------------------------------------------------------------------------------')
        log.debug('|                                                        Unique Endbondtorsion Coeffs                                                        |')
        log.debug('----------------------------------------------------------------------------------------------------------------------------------------------')
        log.debug('| {:^2} | {:^35} | {:^95} |'.format('Id', 'Type', 'Coeffs'))
        log.debug('----------------------------------------------------------------------------------------------------------------------------------------------')
        for i in new.endbondtorsion_coeffs:
            coeff = new.endbondtorsion_coeffs[i]
            id1, id2, id3, id4 = split_coeff(coeff.type)
            log.debug('| {:^2} | {:^8} {:^8} {:^8} {:^8} | {:^95} |'.format(i, id1, id2, id3, id4, string_coeffs(coeff.coeffs) ) )
        log.debug('----------------------------------------------------------------------------------------------------------------------------------------------')
        
        # Print Middlebondtorsion coeffs
        log.debug('\n-------------------------------------------------------------------------------------------------')
        log.debug('|                                Unique Middlebondtorsion Coeffs                                |')
        log.debug('-------------------------------------------------------------------------------------------------')
        log.debug('| {:^2} | {:^35} | {:^50} |'.format('Id', 'Type', 'Coeffs'))
        log.debug('-------------------------------------------------------------------------------------------------')
        for i in new.middlebondtorsion_coeffs:
            coeff = new.middlebondtorsion_coeffs[i]
            id1, id2, id3, id4 = split_coeff(coeff.type)
            log.debug('| {:^2} | {:^8} {:^8} {:^8} {:^8} | {:^50} |'.format(i, id1, id2, id3, id4, string_coeffs(coeff.coeffs) ) )
        log.debug('-------------------------------------------------------------------------------------------------')
        
        # Print Bondbond13 coeffs
        log.debug('\n----------------------------------------------------------------------------------')
        log.debug('|                            Unique Bondbond13 Coeffs                            |')
        log.debug('----------------------------------------------------------------------------------')
        log.debug('| {:^2} | {:^35} | {:^35} |'.format('Id', 'Type', 'Coeffs'))
        log.debug('----------------------------------------------------------------------------------')
        for i in new.bondbond13_coeffs:
            coeff = new.bondbond13_coeffs[i]
            id1, id2, id3, id4 = split_coeff(coeff.type)
            log.debug('| {:^2} | {:^8} {:^8} {:^8} {:^8} | {:^35} |'.format(i, id1, id2, id3, id4, string_coeffs(coeff.coeffs) ) )
        log.debug('----------------------------------------------------------------------------------')
        
        # Print Angletorsion coeffs
        log.debug('\n----------------------------------------------------------------------------------------------------------------------------------------------')
        log.debug('|                                                         Unique Angletorsion Coeffs                                                         |')
        log.debug('----------------------------------------------------------------------------------------------------------------------------------------------')
        log.debug('| {:^2} | {:^35} | {:^95} |'.format('Id', 'Type', 'Coeffs'))
        log.debug('----------------------------------------------------------------------------------------------------------------------------------------------')
        for i in new.angletorsion_coeffs:
            coeff = new.angletorsion_coeffs[i]
            id1, id2, id3, id4 = split_coeff(coeff.type)
            log.debug('| {:^2} | {:^8} {:^8} {:^8} {:^8} | {:^95} |'.format(i, id1, id2, id3, id4, string_coeffs(coeff.coeffs) ) )
        log.debug('----------------------------------------------------------------------------------------------------------------------------------------------')
        
        # Print Angleangle coeffs
        log.debug('\n--------------------------------------------------------------------------------------------------------------------')
        log.debug('|                                             Unique Angleangle Coeffs                                             |')
        log.debug('--------------------------------------------------------------------------------------------------------------------')
        log.debug('| {:^2} | {:^44} | {:^60} |'.format('Id', 'Type', 'Coeffs'))
        log.debug('--------------------------------------------------------------------------------------------------------------------')
        for i in new.angleangle_coeffs:
            coeff = new.angleangle_coeffs[i]
            id1, id2, id3, id4, id5 = split_coeff(coeff.type)
            log.debug('| {:^2} | {:^8} {:^8} {:^8} {:^8} {:^8} | {:^60} |'.format(i, id1, id2, id3, id4, id5, string_coeffs(coeff.coeffs) ) )
        log.debug('--------------------------------------------------------------------------------------------------------------------') 
    return