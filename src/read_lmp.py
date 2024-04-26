# -*- coding: utf-8 -*-
"""
@author: Benjamin Jensen, with modifications by Will Pisani, 
         with modification by Josh Kemppainen (image flags,
         type labels, headers, triclinic support, morse/
         harmonic support, moved mass instance to be called
         from Coeff_class, added ability to read style hints)
Revision 3.6
March 7th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49913
"""    
##############################
# Import Necessary Libraries #
##############################
import gzip


# Function to strip comments
def strip_comment(line):
    # remove comments
    end = line.find('#')
    if end >= 0:
        #comment = line[end:]
        line = line[:end]
    return line

# Function to check if variable is an int
def check_integer(variable):
    try:
        int(variable)
        return_boolean = True
    except:
        return_boolean = False
    return return_boolean

# Function get get style hint
def get_style_hint(whole_line):
    if '#' in whole_line:
        try: style_hint = whole_line.split('#')[-1].strip()
        except: style_hint = 'N/A' 
    else: style_hint = 'N/A'
    return style_hint   

# Function for stringing together parameter types
def string_parameter_type(parameter_type):
    string = ''; str_buffer = 6 # Set string buffer size
    for n, i in enumerate(parameter_type):
        if n == 0: string += '{:<{str_buffer}} '.format(i, str_buffer=str_buffer)
        elif n < len(parameter_type)-1: string += ' {:^{str_buffer}} '.format(i, str_buffer=str_buffer+2)
        else: string += ' {:>{str_buffer}} {:^2}'.format(i, '', str_buffer=str_buffer)
    return string

# Function to string type labels (No white spacing)
def string_type_labels(parameter_type):
    string = '';
    for n, i in enumerate(parameter_type):
        string += i
        if n+1 < len(parameter_type):
            string += '-'
    return string

# Function to create multibody comment from type label
def Nbody_comment_from_typelabel(atom_type_labels_forward, typelabel, nbody):
    types = tuple(typelabel.split('-')); string = 'N/A' # Default string
    # Check if types can be split by the '-' character
    if len(types) == nbody: string = string_parameter_type(types)
    else: # else an atom type uses the '-' character and a more advanced approach is needed
        atomtypes = list(atom_type_labels_forward.keys()); indexes = {} # {atomtype:[indexes]}
        atomtypes.extend(['nb==3', 'nb!=3'])
        for i in atomtypes:
            typecount = typelabel.count(i)
            if typecount > 0:
                for j in range(typecount):
                    if i in indexes:
                        indexes[i].append( typelabel.find(i, max(indexes[i])+1) )
                    else:
                        indexes[i] = [typelabel.find(i, 0)]
        typeordering = {} # {index:atomtype}
        for i in indexes:
            for j in indexes[i]:
                typeordering[j] = i
        typeordering = dict(sorted(typeordering.items(), key=lambda item: item[0]))
        types = tuple(typeordering.values())
        string = string_parameter_type(types)
        if string_type_labels(types) != typelabel:
            print(f'WARNING could not resolve {nbody}-body type label "{typelabel}". Best guess is "{string}" and will be used.')
    return string
                    
# Function for label type mapping
def label_type_mapping(Type, forward_map, reverse_map, method, section):
    # forward and reverse map will be used to map type to and from label type
    # method will set which mapping occurs (options: 'forward' or 'reverse')
    # 'forward' mapping will set all types as integers
    # 'reverse' mapping will set all types as sring type labels
    
    # section will determine what to do to each section to make it consistant (options: 'topology' or 'force-feild')
    # 'topology' section will be for atoms, bonds, angles, dihedrals, impropers mapping
    # 'force-feild' section will be for atoms, bonds, angles, dihedrals, impropers, and cross term coeffs
    
    # Topology section mapping (will convert topology paramters to be consistent with force-feild section str(type-label) -> int(type id) )
    if section == 'topology' and forward_map and reverse_map:
        if  method == 'forward' and not check_integer(Type): Type = forward_map[Type]
        if  method == 'reverse': Type = str(Type)
        
    
    # Force-feild section mapping (will convert force-feild paramters to be consistent with topology section int(type id) -> str(type-label) )
    elif section == 'force-feild' and forward_map and reverse_map:
        if  method == 'forward' and check_integer(Type): Type = int(Type) #reverse_map[int(type)]
        if  method == 'reverse' and check_integer(Type): Type = reverse_map[int(Type)]
        
    # else set type as integer if forward and revers map dicts are empty it must be a standard lammps datafile
    else:
        Type = int(Type)
        
    # make sure type is int or string
    try: Type = int(Type)
    except: Type = str(Type)
    return Type           

               
##########################
# Class to read datafile #
##########################
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz .comment
class Bond: pass  # .type .atomids = [atom1id, atom2id]
class Angle: pass  # .type .atomids = [atom1id, atom2id, atom3id]
class Dihedral: pass  # .type .atomids = [atom1id, atom2id, atom3id, atom4id]
class Improper: pass  # .type .atomids = [atom1,atom2,atom3,atom4]
class Coeff_class: pass  # .type .coeffs = []
class Molecule_File:
    def __init__(self, inmolfile, method='forward', sections=['Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers', 'Velocities']):
        # Filename
        self.filename = inmolfile
        
        # Type labels info
        self.atom_type_labels_forward = {}  # {atom type : atom type label}
        self.atom_type_labels_reverse = {}  # {atom type label: atom type}
        self.bond_type_labels_forward = {}  # {bond type: bond type label}
        self.bond_type_labels_reverse = {}  # {bond type label: bond type}
        self.angle_type_labels_forward = {}  # {angle type : angle type label}
        self.angle_type_labels_reverse = {}  # {angle type label : angle type}
        self.dihedral_type_labels_forward = {}  # {dihedral type : dihedral type label}
        self.dihedral_type_labels_reverse = {}  # {dihedral type label : dihedral type}
        self.improper_type_labels_forward = {}  # {improper type : improper type label}
        self.improper_type_labels_reverse = {}  # {improper type label: improper type}
        self.type_labels_flag = False # update if type labels exists
        
        # Structure info
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.angles = {}  # {angle number : angle object}
        self.dihedrals = {}  # {dihedral number : dihedral object}
        self.impropers = {}  # {improper number : improper object}
        self.velocities = {}  # {atom number : tuple of velocities}

        # Parameters
        self.masses = {}  # {atom type : list of coeffs}
        self.pair_coeffs = {}  # {atom type : list of coeffs}
        self.bond_coeffs = {}   # {bond type : list of coeffs}  {1: [340,1.5], 2: [450,1.2], ...}
        self.angle_coeffs = {}  # {angle type : list of coeffs}
        self.dihedral_coeffs = {}  # {dihedral type : list of coeffs}
        self.improper_coeffs = {}  # {improper type : list of coeffs}
        self.bondbond_coeffs = {} # {cross-term number : list of coeffs}, angles
        self.bondangle_coeffs = {} # {cross-term number : list of coeffs}, angles
        self.angleangle_coeffs = {} # {cross-term number : list of coeffs}, impropers
        self.angleangletorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.endbondtorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.middlebondtorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.bondbond13_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.angletorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        
        # Style hints
        self.mass_coeffs_style_hint = 'N/A'
        self.pair_coeffs_style_hint = 'N/A'
        self.bond_coeffs_style_hint = 'N/A'
        self.angle_coeffs_style_hint = 'N/A'
        self.dihedral_coeffs_style_hint = 'N/A'
        self.improper_coeffs_style_hint = 'N/A'
        self.bondbond_coeffs_style_hint = 'N/A'
        self.bondangle_coeffs_style_hint = 'N/A'
        self.angleangle_coeffs_style_hint = 'N/A'
        self.angleangletorsion_coeffs_style_hint = 'N/A'
        self.endbondtorsion_coeffs_style_hint = 'N/A'
        self.middlebondtorsion_coeffs_style_hint = 'N/A'
        self.bondbond13_coeffs_style_hint = 'N/A'
        self.angletorsion_coeffs_style_hint = 'N/A'
        self.atomstyle = 'full' # Default to full and update later if syle hint exists

        # header, box, and misc
        self.total_line = ''
        self.ttype_line = ''
        self.xbox_line = ''
        self.ybox_line = ''
        self.zbox_line = ''
        self.triclinicbox_line = ''
        self.xlo = -0.5
        self.xhi = 0.5
        self.ylo = -0.5
        self.yhi = 0.5
        self.zlo = -0.5
        self.zhi = 0.5
        self.lx = 1
        self.ly = 1
        self.lz = 1
        self.xy = 0
        self.xz = 0
        self.yz = 0
        self.extra_lines = ''
        self.header = ''
        
        # number of quantities
        self.total = 0
        self.natoms = 0
        self.natomtypes = 0
        self.nbonds = 0
        self.nbondtypes = 0
        self.nangles = 0
        self.nangletypes = 0
        self.ndihedrals = 0
        self.ndihedraltypes = 0
        self.nimpropers = 0
        self.nimpropertypes = 0
        self.nbondbond = 0
        self.nbondangle = 0
        self.nangleangle = 0
        self.nangleangletorsion = 0
        self.nendbondtorsion = 0
        self.nmiddlebondtorsion = 0
        self.nbondbond13 = 0
        self.nangletorsion = 0

        
        # # Open and read file
        if inmolfile.endswith('.gz') or '.gz' in inmolfile:
            print('Reading gzipped datafile')
            with gzip.open(inmolfile, 'rt') as file_handle:
                f = file_handle.readlines()
        else:
            with open(inmolfile, 'r') as file_handle:
                f = file_handle.readlines()

        # Initialize flags
        coeff_flag = False
        atomflag = False
        bondflag = False
        angleflag = False
        dihedralflag = False
        improperflag = False
        velocityflag = False 
        atom_type_labels_flag = False
        bond_type_labels_flag = False
        angle_type_labels_flag = False
        dihedral_type_labels_flag = False
        improper_type_labels_flag = False
        skip = 0
        
        # Loop through lines of file
        for n, whole_line in enumerate(f):
            # skip comment lines
            skip -= 1
            if skip >= 0:
                continue

            # remove comments
            line = strip_comment(whole_line)
            line = line.strip()

            # begining of a section, flag the start and skip one line
            if line == '':
                coeff_flag = False
                atomflag = False
                bondflag = False
                angleflag = False
                dihedralflag = False
                improperflag = False
                velocityflag = False
                atom_type_labels_flag = False
                bond_type_labels_flag = False
                angle_type_labels_flag = False
                dihedral_type_labels_flag = False
                improper_type_labels_flag = False
            elif n == 0:
                self.header = line
                continue
            elif 'Atom Type Labels' in line:
                atom_type_labels_flag = True
                self.type_labels_flag = True
                skip = 1; continue;
            elif 'Bond Type Labels' in line:
                bond_type_labels_flag = True
                self.type_labels_flag = True
                skip = 1; continue;
            elif 'Angle Type Labels' in line:
                angle_type_labels_flag = True
                self.type_labels_flag = True
                skip = 1; continue;
            elif 'Dihedral Type Labels' in line:
                dihedral_type_labels_flag = True
                self.type_labels_flag = True
                skip = 1; continue ; 
            elif 'Improper Type Labels' in line:
                improper_type_labels_flag = True
                self.type_labels_flag = True
                skip = 1; continue; 
            elif 'atoms' in line:
                self.total_line = line
                self.total = int(line.split()[0])
                self.natoms = self.total
                continue
            elif 'bonds' in line:
                self.nbonds = int(line.split()[0])
                continue
            elif 'angles' in line:
                self.nangles = int(line.split()[0])
                self.nbondbond = int(line.split()[0])
                self.nbondangle = int(line.split()[0])
                continue
            elif 'dihedrals' in line:
                self.ndihedrals = int(line.split()[0])
                self.nangleangletorsion = int(line.split()[0])
                self.nendbondtorsion = int(line.split()[0])
                self.nmiddlebondtorsion = int(line.split()[0])
                self.nbondbond13 = int(line.split()[0])
                self.nangletorsion = int(line.split()[0])
                continue
            elif 'impropers' in line:
                self.nimpropers = int(line.split()[0])
                self.nangleangle = int(line.split()[0])
                continue
            elif 'atom types' in line:
                self.ttype_line = line
                self.natomtypes = int(line.split()[0])
                continue
            elif 'bond types' in line:
                self.nbondtypes = int(line.split()[0])
                continue
            elif 'angle types' in line:
                self.nangletypes = int(line.split()[0])
                continue
            elif 'dihedral types' in line:
                self.ndihedraltypes = int(line.split()[0])
                continue
            elif 'improper types' in line:
                self.nimpropertypes = int(line.split()[0])
                continue
            elif 'per atom' in line:
                self.extra_lines += line + '\n'
                continue
            elif 'xlo' in line:
                self.xbox_line = line
                line = line.split()
                self.xlo = float(line[0])
                self.xhi = float(line[1])
                self.lx = self.xhi - self.xlo
                continue
            elif 'ylo' in line:
                self.ybox_line = line
                line = line.split()
                self.ylo = float(line[0])
                self.yhi = float(line[1])
                self.ly = self.yhi - self.ylo
                continue
            elif 'zlo' in line:
                self.zbox_line = line
                line = line.split()
                self.zlo = float(line[0])
                self.zhi = float(line[1])
                self.lz = self.zhi - self.zlo
                continue
            elif 'xy' in line and 'xz' in line and 'yz' in line:
                self.triclinicbox_line = line
                line = line.split()
                self.xy = float(line[0])
                self.xz = float(line[1])
                self.yz = float(line[2])
                continue                   
            elif line == 'Masses':
                self.mass_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.atom_type_labels_forward
                reverse_map = self.atom_type_labels_reverse
                coeffs = self.masses
                skip = 1
                continue
            elif line == 'Pair Coeffs':
                self.pair_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.atom_type_labels_forward
                reverse_map = self.atom_type_labels_reverse
                coeffs = self.pair_coeffs
                skip = 1
                continue
            elif line == 'Bond Coeffs':
                self.bond_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.bond_type_labels_forward
                reverse_map = self.bond_type_labels_reverse
                coeffs = self.bond_coeffs
                skip = 1
                continue
            elif line == 'Angle Coeffs':
                self.angle_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.angle_type_labels_forward
                reverse_map = self.angle_type_labels_reverse
                coeffs = self.angle_coeffs
                skip = 1
                continue
            elif line == 'Dihedral Coeffs':
                self.dihedral_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.dihedral_type_labels_forward
                reverse_map = self.dihedral_type_labels_reverse
                coeffs = self.dihedral_coeffs
                skip = 1
                continue
            elif line == 'Improper Coeffs':
                self.improper_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.improper_type_labels_forward
                reverse_map = self.improper_type_labels_reverse
                coeffs = self.improper_coeffs
                skip = 1
                continue
            elif line == 'BondBond Coeffs':
                self.bondbond_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.angle_type_labels_forward
                reverse_map = self.angle_type_labels_reverse
                coeffs = self.bondbond_coeffs
                skip = 1
                continue
            elif line == 'BondAngle Coeffs':
                self.bondangle_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.angle_type_labels_forward
                reverse_map = self.angle_type_labels_reverse
                coeffs = self.bondangle_coeffs
                skip = 1
                continue
            elif line == 'AngleAngle Coeffs':
                self.angleangle_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.improper_type_labels_forward
                reverse_map = self.improper_type_labels_reverse
                coeffs = self.angleangle_coeffs
                skip = 1
                continue
            elif line == 'AngleAngleTorsion Coeffs':
                self.angleangletorsion_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.dihedral_type_labels_forward
                reverse_map = self.dihedral_type_labels_reverse
                coeffs = self.angleangletorsion_coeffs
                skip = 1
                continue
            elif line == 'EndBondTorsion Coeffs':
                self.endbondtorsion_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.dihedral_type_labels_forward
                reverse_map = self.dihedral_type_labels_reverse
                coeffs = self.endbondtorsion_coeffs
                skip = 1
                continue
            elif line == 'MiddleBondTorsion Coeffs':
                self.middlebondtorsion_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.dihedral_type_labels_forward
                reverse_map = self.dihedral_type_labels_reverse
                coeffs = self.middlebondtorsion_coeffs
                skip = 1
                continue
            elif line == 'BondBond13 Coeffs':
                self.bondbond13_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.dihedral_type_labels_forward
                reverse_map = self.dihedral_type_labels_reverse
                coeffs = self.bondbond13_coeffs
                skip = 1
                continue
            elif line == 'AngleTorsion Coeffs':
                self.angletorsion_coeffs_style_hint = get_style_hint(whole_line)
                coeff_flag = True
                forward_map = self.dihedral_type_labels_forward
                reverse_map = self.dihedral_type_labels_reverse
                coeffs = self.angletorsion_coeffs
                skip = 1
                continue
            elif line == 'Atoms' and 'Atoms' in sections:
                atomflag = True
                if '#' in whole_line: self.atomstyle = whole_line.split()[-1]
                skip = 1
                continue
            elif line == 'Bonds' and 'Bonds' in sections:
                bondflag = True
                skip = 1
                continue
            elif line == 'Angles' and 'Angles' in sections:
                angleflag = True
                skip = 1
                continue
            elif line == 'Dihedrals' and 'Dihedrals' in sections:
                dihedralflag = True
                skip = 1
                continue
            elif line == 'Impropers' and 'Impropers' in sections:
                improperflag = True
                skip = 1
                continue
            elif line == 'Velocities' and 'Velocities' in sections:
                velocityflag = True
                skip = 1
                continue

            # Find atom type labels
            if atom_type_labels_flag:
                line = line.split()
                ID = int(line[0])
                Type = line[1]
                self.atom_type_labels_forward[Type] = ID  
                self.atom_type_labels_reverse[ID] = Type     
            # Find bond type labels
            elif bond_type_labels_flag:
                line = line.split()
                ID = int(line[0])
                Type = line[1]
                self.bond_type_labels_forward[Type] = ID
                self.bond_type_labels_reverse[ID] = Type
            # Find angle type labels
            elif angle_type_labels_flag:
                line = line.split()
                ID = int(line[0])
                Type = line[1]
                self.angle_type_labels_forward[Type] = ID
                self.angle_type_labels_reverse[ID] = Type
            # Find dihedral type labels
            elif dihedral_type_labels_flag:
                line = line.split()
                ID = int(line[0])
                Type = line[1]
                self.dihedral_type_labels_forward[Type] = ID
                self.dihedral_type_labels_reverse[ID] = Type
            # Find improper type labels
            elif improper_type_labels_flag:
                line = line.split()
                ID = int(line[0])
                Type = line[1]
                self.improper_type_labels_forward[Type] = ID
                self.improper_type_labels_reverse[ID] = Type


            # Find coeffs    
            if coeff_flag:
                line = line.split()
                ID = label_type_mapping(line[0], forward_map, reverse_map, method, section='force-feild')
                c = Coeff_class(); idcoeffs = [];
                for i in line[1:]:
                    if '.' not in i and 'e' not in i and 'class2' not in i and 'morse' not in i:
                        idcoeffs.append(int(i))
                    elif 'e' not in i and 'class2' not in i and 'morse' not in i:
                        idcoeffs.append(float(i))
                       
                    # for morse and class2 coeffs
                    elif 'class2' in i: idcoeffs.append(i)
                    elif 'morse' in i: idcoeffs.append(i)
                    else: idcoeffs.append(float(i))
                if '#' not in whole_line:
                    c.type = 'N/A'
                else:
                    c.type = whole_line.split('#')[-1].rstrip().lstrip()
                c.coeffs = idcoeffs
                coeffs[ID] = c

            # Find atoms info
            if atomflag:
                line = line.split()
                if 'charge' in self.atomstyle:
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    charge = float(line[2]); x = float(line[3]); y = float(line[4]); z = float(line[5]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[6]); iy = int(line[7]); iz = int(line[8]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.charge = charge
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'molecular' in self.atomstyle:
                    ID = int(line[0]); molid = int(line[1]);
                    Type = label_type_mapping(line[2], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    x = float(line[3]); y = float(line[4]); z = float(line[5]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[6]); iy = int(line[7]); iz = int(line[8]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.molid = molid
                    a.type = Type
                    a.charge = 0
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'full' in self.atomstyle:
                    ID = int(line[0]); molid = int(line[1]);
                    Type = label_type_mapping(line[2], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    charge = float(line[3]); x = float(line[4]); y = float(line[5]); z = float(line[6]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[7]); iy = int(line[8]); iz = int(line[9]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.molid = molid
                    a.type = Type
                    a.charge = charge
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'angle' in self.atomstyle:
                    ID = int(line[0]); molid = int(line[1]);
                    Type = label_type_mapping(line[2], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    x = float(line[3]); y = float(line[4]); z = float(line[5]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[6]); iy = int(line[7]); iz = int(line[8]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.molid = molid
                    a.type = Type
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'atomic' in self.atomstyle:
                    ID = int(line[0]);
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    x = float(line[2]); y = float(line[3]); z = float(line[4]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[5]); iy = int(line[6]); iz = int(line[7]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'body' in self.atomstyle:
                    ID = int(line[0]);
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    bodyflag = int(line[2]); mass = float(line[3]);
                    x = float(line[4]); y = float(line[5]); z = float(line[6]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[7]); iy = int(line[8]); iz = int(line[9]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.bodyflag = bodyflag
                    a.mass = mass
                    a.type = Type
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'bond' in self.atomstyle:
                    ID = int(line[0]); molid = int(line[1]);
                    Type = label_type_mapping(line[2], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    x = float(line[3]); y = float(line[4]); z = float(line[5]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[6]); iy = int(line[7]); iz = int(line[8]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.molid = molid
                    a.type = Type
                    a.charge = 0
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'dielectric' in self.atomstyle:
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    charge = float(line[2]); x = float(line[3]); y = float(line[4]); z = float(line[5]);
                    normx = float(line[6]); normy = float(line[7]); normz = float(line[8]); area = float(line[9]);
                    ed = line[10]; em = line[11]; epsilon = line[12]; curvature = line[13];
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[14]); iy = int(line[15]); iz = int(line[16]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.charge = charge
                    a.x = x
                    a.y = y
                    a.z = z
                    a.normx = normx
                    a.normy = normy
                    a.normz = normz
                    a.area = area
                    a.ed = ed
                    a.em = em
                    a.epsilon = epsilon
                    a.curvature = curvature
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'dipole' in self.atomstyle:
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    charge = float(line[2]); x = float(line[3]); y = float(line[4]); z = float(line[5]);
                    mux = float(line[6]); muy = float(line[7]); muz = float(line[8]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[9]); iy = int(line[10]); iz = int(line[11]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.charge = charge
                    a.x = x
                    a.y = y
                    a.z = z
                    a.mux = mux
                    a.muy = muy
                    a.muz = muz
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'edpd' in self.atomstyle:
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    edpd_temp = float(line[2]); edpd_cv = float(line[3]);
                    x = float(line[4]); y = float(line[5]); z = float(line[6]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[7]); iy = int(line[8]); iz = int(line[9]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.edpd_temp = edpd_temp
                    a.edpd_cv = edpd_cv
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'mdpd' in self.atomstyle:
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    rho = float(line[2]); x = float(line[3]); y = float(line[4]); z = float(line[5]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[6]); iy = int(line[7]); iz = int(line[8]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.rho = rho
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'tdpd' in self.atomstyle:
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    x = float(line[2]); y = float(line[3]); z = float(line[4]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'dpd' in self.atomstyle:
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    theta = float(line[2]); x = float(line[3]); y = float(line[4]); z = float(line[5]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[6]); iy = int(line[7]); iz = int(line[8]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.theta = theta
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'hybrid' in self.atomstyle:
                    ID = int(line[0]);
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    x = float(line[2]); y = float(line[3]); z = float(line[4]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'bpm/sphere' in self.atomstyle:
                    ID = int(line[0]); molid = int(line[1]);
                    Type = label_type_mapping(line[2], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    diameter = float(line[3]); density = float(line[4]);
                    x = float(line[5]); y = float(line[6]); z = float(line[7]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[8]); iy = int(line[9]); iz = int(line[10]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.molid = molid
                    a.type = Type
                    a.charge = charge
                    a.x = x
                    a.y = y
                    a.z = z
                    a.diameter = diameter
                    a.density = density
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'sphere' in self.atomstyle:
                    ID = int(line[0]);
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    diameter = float(line[2]); density = float(line[3]);
                    x = float(line[4]); y = float(line[5]); z = float(line[6]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[7]); iy = int(line[8]); iz = int(line[9]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.charge = charge
                    a.x = x
                    a.y = y
                    a.z = z
                    a.diameter = diameter
                    a.density = density
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'sph' in self.atomstyle:
                    ID = int(line[0]);
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    rho = float(line[2]); esph = float(line[3]); cv = float(line[4]);
                    x = float(line[5]); y = float(line[6]); z = float(line[7]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[8]); iy = int(line[9]); iz = int(line[10]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.charge = charge
                    a.x = x
                    a.y = y
                    a.z = z
                    a.rho = rho
                    a.esph = esph
                    a.cv = cv
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'electron' in self.atomstyle:
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    charge = float(line[2]); espin = float(line[3]); eradius = float(line[4]);
                    x = float(line[5]); y = float(line[6]); z = float(line[7]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[8]); iy = int(line[9]); iz = int(line[10]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.charge = charge
                    a.espin = espin
                    a.eradius = eradius
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'ellipsoid' in self.atomstyle:
                    ID = int(line[0]);
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    ellipsoidflag = int(line[2]); density = float(line[3]);
                    x = float(line[4]); y = float(line[5]); z = float(line[6]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[7]); iy = int(line[8]); iz = int(line[9]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.ellipsoidflag = ellipsoidflag
                    a.density = density
                    a.type = Type
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'line' in self.atomstyle:
                    ID = int(line[0]); molid = int(line[1]);
                    Type = label_type_mapping(line[2], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    lineflag = int(line[3]); density = float(line[4]); 
                    x = float(line[5]); y = float(line[6]); z = float(line[7]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[8]); iy = int(line[9]); iz = int(line[10]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.molid = molid
                    a.type = Type
                    a.lineflag = lineflag
                    a.density = density
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'peri' in self.atomstyle:
                    ID = int(line[0]);
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    volume = float(line[2]); density = float(line[3]);
                    x = float(line[4]); y = float(line[5]); z = float(line[6]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[7]); iy = int(line[8]); iz = int(line[9]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.bodyflag = bodyflag
                    a.volume = volume
                    a.density = density
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'spin' in self.atomstyle:
                    ID = int(line[0]);
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    x = float(line[2]); y = float(line[3]); z = float(line[4]);
                    spx = float(line[5]); spy = float(line[6]); spz = float(line[7]); sp = float(line[8]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[5]); iy = int(line[6]); iz = int(line[7]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.x = x
                    a.y = y
                    a.z = z
                    a.spx = spx
                    a.spy = spy
                    a.spz = spz
                    a.sp = sp
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'smd' in self.atomstyle:
                    ID = int(line[0]);
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    molecule = float(line[2]); volume = float(line[3]); mass = float(line[4]); kradius = float(line[5]);
                    cradius = float(line[6]); x0 = float(line[7]); y0 = float(line[8]); z0 = float(line[9]);
                    x = float(line[10]); y = float(line[11]); z = float(line[12]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[13]); iy = int(line[14]); iz = int(line[15]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.x = x
                    a.y = y
                    a.z = z
                    a.molecule = molecule
                    a.volume = volume
                    a.mass = mass
                    a.kradius = kradius
                    a.cradius = cradius
                    a.x0 = x0
                    a.y0 = y0
                    a.z0 = z0
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'template' in self.atomstyle:
                    ID = int(line[0]);
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    molid = int(line[2]); template_index = int(line[3]); template_atom = int(line[4]);
                    x = float(line[5]); y = float(line[6]); z = float(line[7]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[8]); iy = int(line[9]); iz = int(line[10]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.x = x
                    a.y = y
                    a.z = z
                    a.molid = molid
                    a.template_index = template_index
                    a.template_atom = template_atom
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'tri' in self.atomstyle:
                    ID = int(line[0]); molid = int(line[1]);
                    Type = label_type_mapping(line[2], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    triangleflag = int(line[3]); density = float(line[4]); 
                    x = float(line[5]); y = float(line[6]); z = float(line[7]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[8]); iy = int(line[9]); iz = int(line[10]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.molid = molid
                    a.type = Type
                    a.triangleflag = triangleflag
                    a.density = density
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                elif 'wavepacket' in self.atomstyle:
                    ID = int(line[0])
                    Type = label_type_mapping(line[1], self.atom_type_labels_forward, self.atom_type_labels_reverse, method, section='topology')
                    charge = float(line[2]); espin = float(line[3]); eradius = float(line[4]);
                    etag = float(line[5]); cs_re = float(line[6]); cs_im = float(line[7]);
                    x = float(line[8]); y = float(line[9]); z = float(line[10]);
                    if '#' not in whole_line: comment = 'N/A'
                    else: comment = whole_line.split('#')[-1].rstrip().lstrip()
                    try: ix = int(line[6]); iy = int(line[7]); iz = int(line[8]);
                    except: ix = 0; iy= 0; iz = 0;
                    a = Atom()
                    a.type = Type
                    a.charge = charge
                    a.espin = espin
                    a.eradius = eradius
                    a.etag = etag
                    a.cs_re = cs_re
                    a.cs_im = cs_im
                    a.x = x
                    a.y = y
                    a.z = z
                    a.ix = ix
                    a.iy = iy
                    a.iz = iz
                    a.comment = comment
                    self.atoms[ID] = a
                else: raise Exception('Atom Style not defined as a style hint in the LAMMPS datafile.')

            # Find bonds
            elif bondflag:
                line = line.split()
                ID = int(line[0])
                Type = label_type_mapping(line[1], self.bond_type_labels_forward, self.bond_type_labels_reverse, method, section='topology')
                atom1id = int(line[2]); atom2id = int(line[3]);
                b = Bond()
                b.type = Type
                b.atomids = (atom1id, atom2id)
                self.bonds[ID] = b

            # Find angles
            elif angleflag:
                line = line.split()
                ID = int(line[0])
                Type = label_type_mapping(line[1], self.angle_type_labels_forward, self.angle_type_labels_reverse, method, section='topology')
                atom1id = int(line[2]); atom2id = int(line[3]); atom3id = int(line[4]);
                c = Angle()
                c.type = Type
                c.atomids = (atom1id, atom2id, atom3id)
                self.angles[ID] = c

            # Find dihedrals
            elif dihedralflag:
                line = line.split()
                ID = int(line[0])
                Type = label_type_mapping(line[1], self.dihedral_type_labels_forward, self.dihedral_type_labels_reverse, method, section='topology')
                atom1id = int(line[2]); atom2id = int(line[3]); atom3id = int(line[4]); atom4id = int(line[5]);
                d = Dihedral()
                d.type = Type
                d.atomids = (atom1id, atom2id, atom3id, atom4id)
                self.dihedrals[ID] = d

            # Find impropers
            elif improperflag:
                line = line.split()
                ID = int(line[0])
                Type = label_type_mapping(line[1], self.improper_type_labels_forward, self.improper_type_labels_reverse, method, section='topology')
                atom1id = int(line[2]); atom2id = int(line[3]); atom3id = int(line[4]); atom4id = int(line[5]);
                i = Improper()
                i.type = Type
                i.atomids = (atom1id, atom2id, atom3id, atom4id)
                self.impropers[ID] = i

            # Find velocities
            elif velocityflag:
                line = line.split()
                ID = int(line[0]); vx = float(line[1]); vy = float(line[2]); vz = float(line[3]);
                self.velocities[ID] = (vx, vy, vz)
        
        #-----------------------------------------------------------------------------------#
        # If type labels exist check for 'N/A' comments on atoms, masses, pair coeffs, bond #
        # coeffs, ... etc and try updating comment based on type label (comment will be set #
        # as if all2lmp.py set the comment for compatability with bond_react_merge.py).     #
        #-----------------------------------------------------------------------------------#
        if self.type_labels_flag:
            # Check that forward and reverse dicts are the same length. If they are not the same length the label
            # types are not unique since keys were written to reverse dict more then once and raise exception
            if len(self.atom_type_labels_forward) != len(self.atom_type_labels_reverse):
                raise Exception('Atom labels are not Unique!')
            if len(self.bond_type_labels_forward) != len(self.bond_type_labels_reverse):
                raise Exception('Bond labels are not Unique!')
            if len(self.angle_type_labels_forward) != len(self.angle_type_labels_reverse):
                raise Exception('Angle labels are not Unique!')
            if len(self.dihedral_type_labels_forward) != len(self.dihedral_type_labels_reverse):
                raise Exception('Dihedral labels are not Unique!')   
            if len(self.improper_type_labels_forward) != len(self.improper_type_labels_reverse):
                raise Exception('Improper labels are not Unique!')
            
            # Check if masses have no comments and if so set masses comment from type label
            comments = {self.masses[i].type for i in self.masses}
            if len(comments) == 1 and 'N/A' in comments:
                for i in self.masses:
                    try: self.masses[i].type = self.atom_type_labels_reverse[i]
                    except: print(f'WARNING - file masses have no comments, but has Type Labels, attempted to set masses comment {i}, but failed. Comment left as N/A')
                    
            # Check if pair coeffs have no comments and if so set pair coeffs comment from type label
            comments = {self.pair_coeffs[i].type for i in self.pair_coeffs}
            if len(comments) == 1 and 'N/A' in comments:
                for i in self.pair_coeffs:
                    try: self.pair_coeffs[i].type = self.atom_type_labels_reverse[i]
                    except: print(f'WARNING - file pair coeffs have no comments, but has Type Labels, attempted to set pair coeff comment {i}, but failed. Comment left as N/A')
                    
            # Check if bond coeffs have no comments and if so set bond coeffs comment from type label
            comments = {self.bond_coeffs[i].type for i in self.bond_coeffs}
            if len(comments) == 1 and 'N/A' in comments:
                for i in self.bond_coeffs:
                    try: self.bond_coeffs[i].type = Nbody_comment_from_typelabel(self.atom_type_labels_forward, self.bond_type_labels_reverse[i], 2)
                    except: print(f'WARNING - file bond coeffs have no comments, but has Type Labels, attempted to set bond coeff comment {i}, but failed. Comment left as N/A')
                    
            # Check if angle coeffs have no comments and if so set angle coeffs comment from type label
            comments = {self.angle_coeffs[i].type for i in self.angle_coeffs}
            if len(comments) == 1 and 'N/A' in comments:
                for i in self.angle_coeffs:
                    try:
                        comment = Nbody_comment_from_typelabel(self.atom_type_labels_forward, self.angle_type_labels_reverse[i], 3)
                        self.angle_coeffs[i].type = comment
                        if i in self.bondbond_coeffs: self.bondbond_coeffs[i].type = comment
                        if i in self.bondangle_coeffs: self.bondangle_coeffs[i].type = comment
                    except: print(f'WARNING - file angle coeffs have no comments, but has Type Labels, attempted to set angle coeff comment {i}, but failed. Comment left as N/A')

            # Check if dihedral coeffs have no comments and if so set dihedral coeffs comment from type label
            comments = {self.dihedral_coeffs[i].type for i in self.dihedral_coeffs}
            if len(comments) == 1 and 'N/A' in comments:
                for i in self.dihedral_coeffs:
                    try:
                        comment = Nbody_comment_from_typelabel(self.atom_type_labels_forward, self.dihedral_type_labels_reverse[i], 4)
                        self.dihedral_coeffs[i].type = comment
                        if i in self.angleangletorsion_coeffs: self.angleangletorsion_coeffs[i].type = comment
                        if i in self.endbondtorsion_coeffs: self.endbondtorsion_coeffs[i].type = comment
                        if i in self.middlebondtorsion_coeffs: self.middlebondtorsion_coeffs[i].type = comment
                        if i in self.bondbond13_coeffs: self.bondbond13_coeffs[i].type = comment
                        if i in self.angletorsion_coeffs: self.angletorsion_coeffs[i].type = comment
                    except: print(f'WARNING - file dihedral coeffs have no comments, but has Type Labels, attempted to set dihedral coeff comment {i}, but failed. Comment left as N/A')

            # Check if improper coeffs have no comments and if so set improper coeffs comment from type label
            comments = {self.improper_coeffs[i].type for i in self.improper_coeffs}
            if len(comments) == 1 and 'N/A' in comments:
                for i in self.improper_coeffs:
                    try:
                        comment = Nbody_comment_from_typelabel(self.atom_type_labels_forward, self.improper_type_labels_reverse[i], 5)
                        self.improper_coeffs[i].type = comment
                        if i in self.angleangle_coeffs: self.angleangle_coeffs[i].type = comment
                    except: print(f'WARNING - file improper coeffs have no comments, but has Type Labels, attempted to set improper coeff comment {i}, but failed. Comment left as N/A')

            # Check if atoms have no comments and if so set atoms comment from type label
            comments = {self.atoms[i].comment for i in self.atoms}
            if len(comments) == 1 and 'N/A' in comments:
                for i in self.atoms:
                    try:
                        atomtype = self.atom_type_labels_reverse[self.atoms[i].type]
                        element = atomtype[0].capitalize()
                        comment = '{}/{}'.format(atomtype, element)
                        self.atoms[i].comment = comment
                    except: print(f'WARNING - file atomss have no comments, but has Type Labels, attempted to set atoms comment {i}, but failed. Comment left as N/A')
