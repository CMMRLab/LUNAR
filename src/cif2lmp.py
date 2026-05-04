# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
May 4, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


parse CIF → build lattice → 
apply symmetry ops to asymmetric unit → 
wrap + deduplicate → 
convert to Cartesian → write LAMMPS
"""
##############################
# Import Necessary Libraries #
##############################
from fractions import Fraction
import math
import re
import os


####################################
### Class for reading .cif files ###
####################################
class Atom_mol: pass # .element .x .y .z
class read_cif:    
    def __init__(self, cif_file):
        self.atoms = [] # [atom1, atom2, ...]
        self.symmetry_ops = []
        
        # Unit cell info
        self.a = None
        self.b = None
        self.c = None
        self.alpha = 90
        self.beta  = 90
        self.gamma = 90
        
        # Thermo info
        self.temp = None
    
        #########################################
        ### Opening and reading the .cif file ###
        #########################################
        with open(cif_file, "r") as f:
            lines = [line.strip() for line in f]

        i = 0
        while i < len(lines):
            line = lines[i]
            if not line or line.startswith("#"):
                i += 1
                continue
            
            # Get cell  and thermo data
            data = self.tokenize_cif_line(line)
            if len(data) >= 2:
                key = data[0]
                val = data[1]
                if key == '_cell_length_a':
                    self.a = self.string_to_number(val)
                elif key == '_cell_length_b':
                    self.b = self.string_to_number(val)
                elif key == '_cell_length_c':
                    self.c = self.string_to_number(val)
                elif key == '_cell_angle_alpha':
                    self.alpha = self.string_to_number(val)
                elif key == '_cell_angle_beta':
                    self.beta = self.string_to_number(val)
                elif key == '_cell_angle_gamma':
                    self.gamma = self.string_to_number(val)
                    
                elif key == '_cell_measurement_temperature':
                    self.temp = self.string_to_number(val)
                
            # Loop block
            if line == "loop_":
                headers, rows, i = self.read_loop(lines, i)

                if self.is_atom_loop(headers):
                    self.parse_atom_loop(headers, rows)

                elif self.is_symmetry_loop(headers):
                    self.parse_symmetry_loop(headers, rows)

                continue

            i += 1

        if None in (self.a, self.b, self.c):
            raise ValueError("Missing required cell lengths in CIF file.")
            
        if not self.atoms:
            raise ValueError("No atom fractional-coordinate loop found.")
        
        if not self.symmetry_ops:
            raise ValueError("No explicit symmetry operations found.")

    def is_atom_loop(self, headers):
        return (
            "_atom_site_fract_x" in headers and
            "_atom_site_fract_y" in headers and
            "_atom_site_fract_z" in headers
        )
    
    def parse_atom_loop(self, headers, rows):
        ix = headers.index("_atom_site_fract_x")
        iy = headers.index("_atom_site_fract_y")
        iz = headers.index("_atom_site_fract_z")
    
        if "_atom_site_type_symbol" in headers:
            ielement = headers.index("_atom_site_type_symbol")
        elif "_atom_site_label" in headers:
            ielement = headers.index("_atom_site_label")
        else:
            raise ValueError("Atom loop has positions but no element/label column.")
    
        ioccupancy = headers.index("_atom_site_occupancy") if "_atom_site_occupancy" in headers else None
        iassembly = headers.index("_atom_site_disorder_assembly") if "_atom_site_disorder_assembly" in headers else None
        igroup = headers.index("_atom_site_disorder_group") if "_atom_site_disorder_group" in headers else None
    
        # Find dominant disorder group for each assembly
        assembly_group_occ = {}
    
        for row in rows:
            if len(row) < len(headers):
                continue
    
            if iassembly is None or igroup is None:
                continue
    
            assembly = row[iassembly]
            group = row[igroup]
    
            if assembly in (".", "?") or group in (".", "?"):
                continue
    
            occ = 1.0
            if ioccupancy is not None:
                occ = self.string_to_number(row[ioccupancy])
                if occ is None:
                    occ = 1.0
    
            key = (assembly, group)
            assembly_group_occ[key] = max(assembly_group_occ.get(key, 0.0), occ)
    
        dominant_group = {}
        for (assembly, group), occ in assembly_group_occ.items():
            if assembly not in dominant_group or occ > dominant_group[assembly][1]:
                dominant_group[assembly] = (group, occ)
    
        for row in rows:
            if len(row) < len(headers):
                continue
    
            occupancy = 1.0
            if ioccupancy is not None:
                occupancy = self.string_to_number(row[ioccupancy])
                if occupancy is None:
                    occupancy = 1.0
    
            # Keep nondisordered atoms.
            # For disordered atoms, keep only the dominant group for that assembly.
            if iassembly is not None and igroup is not None:
                assembly = row[iassembly]
                group = row[igroup]
    
                if assembly not in (".", "?") and group not in (".", "?"):
                    keep_group = dominant_group.get(assembly, (group, occupancy))[0]
                    if group != keep_group:
                        continue
    
            atom = Atom_mol()
            atom.element = self.clean_element(row[ielement])
            atom.x = self.string_to_number(row[ix])
            atom.y = self.string_to_number(row[iy])
            atom.z = self.string_to_number(row[iz])
            atom.occupancy = occupancy
    
            self.atoms.append(atom)
        
    def is_symmetry_loop(self, headers):
        return (
            "_symmetry_equiv_pos_as_xyz" in headers or
            "_space_group_symop_operation_xyz" in headers
        )
    
    def parse_symmetry_loop(self, headers, rows):
        if "_symmetry_equiv_pos_as_xyz" in headers:
            iop = headers.index("_symmetry_equiv_pos_as_xyz")
        else:
            iop = headers.index("_space_group_symop_operation_xyz")

        for row in rows:
            if len(row) > iop:
                op = " ".join(row[iop:]).strip()
                op = op.strip("'").strip('"')
                self.symmetry_ops.append(op)
                
    def tokenize_cif_line(self, line):
        return re.findall(r"""(?:'[^']*'|"[^"]*"|\S+)""", line)
                
    def read_loop(self, lines, i):
        """
        Reads a CIF loop beginning at lines[i] == 'loop_'.

        Returns:
            headers : list[str]
            rows    : list[list[str]]
            i       : next unread line index
        """
        i += 1
        headers = []

        # Read loop headers
        while i < len(lines):
            line = lines[i].strip()

            if not line or line.startswith("#"):
                i += 1
                continue

            if line.startswith("_"):
                headers.append(line.split()[0])
                i += 1
            else:
                break

        # Read loop data rows
        rows, tokens = [], []
        while i < len(lines):
            line = lines[i].strip()
        
            if not line or line.startswith("#"):
                i += 1
                continue
        
            if line == "loop_" or line.startswith("data_") or line.startswith("_"):
                break
        
            tokens.extend(self.tokenize_cif_line(line))
        
            while len(tokens) >= len(headers):
                rows.append(tokens[:len(headers)])
                tokens = tokens[len(headers):]
        
            i += 1

        return headers, rows, i
                
    def string_to_number(self, string):
        string = string.strip()
        
        # missing values
        if string in ('?', '.'):
            return None
        
        # strip quotes
        if (string.startswith("'") and string.endswith("'")) or \
           (string.startswith('"') and string.endswith('"')):
            string = string[1:-1]
            
        # remove uncertainty
        string = re.sub(r"\([0-9]+\)$", "", string)
        
        # handle fractions
        if "/" in string:
            return float(Fraction(string))
        
        return float(string)
    
    def clean_element(self, value):
        value = value.strip().strip("'").strip('"')

        # Convert labels like Si1, O2, Fe3 into Si, O, Fe
        match = re.match(r"[A-Za-z]+", value)
        if not match:
            return value

        element = match.group(0)
        if len(element) == 1:
            return element.upper()

        return element[0].upper() + element[1].lower()
    

    
    
######################################################
### Class for converting .cif into a class that is ###
### exactly the same as it comes from read_lmp to  ###
### be able to use .cif files in remaining of the  ###
### codes. *NOTE: after reading you must do the    ###
### following:                                     ###
###     m = cif2lmp.Molecule_File(topofile)        ###
###     m = periodicity.wrap_atoms(m)              ###
###     m = periodicity.reset_image_flags(m)       ###
### To ensure atoms are in the box and to set the  ###
### image flags.                                   ###
######################################################
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz
class Bond: pass  # .type .atomids = [atom1id, atom2id]
class Molecule_File:
    def __init__(self, cif_file):

        # Read cif file
        cif = read_cif(cif_file) 
        self.a = cif.a
        self.b = cif.b
        self.c = cif.c
        self.alpha = cif.alpha
        self.beta  = cif.beta
        self.gamma = cif.gamma
        self.h, self.h_inv = self.cell_to_lammps_h()
        
        self.symmetry_ops = cif.symmetry_ops
        self.asym_natoms = len(cif.atoms)
        self.temp = cif.temp
        self.elements = set()
        

        # Apply symmetry operations            
        coords = [] # [(x, y, z, element)]
        for atom in cif.atoms:
            x0 = atom.x
            y0 = atom.y
            z0 = atom.z
            element = atom.element
            self.elements.add(element)
        
            # Apply sym-ops using nuttered eval
            globals_dict = {'x':x0, 'y':y0, 'z':z0}
            for sym_op in cif.symmetry_ops:
                sym_op = sym_op.split(',')
                if len(sym_op) != 3: continue
                sym_x, sym_y, sym_z = sym_op
                
                # Eval sym-op
                sx = self.nuttered_eval(sym_x, globals_dict=globals_dict)
                sy = self.nuttered_eval(sym_y, globals_dict=globals_dict)
                sz = self.nuttered_eval(sym_z, globals_dict=globals_dict)
                
                # Wrap coords
                wx = self.wrap_frac(sx, tol=1e-8)
                wy = self.wrap_frac(sy, tol=1e-8)
                wz = self.wrap_frac(sz, tol=1e-8)
                coords.append((wx, wy, wz, element))
                
                
        # Remove duplicate symmetry-generated atoms
        self.coords = self.deduplicate_frac_coords(coords, tol=1e-6)
        self.coords = sorted(self.coords, key=lambda x: (x[3], x[0], x[1], x[2])) # sort by element
                             
        # Set system information
        self.filename = cif_file
        self.natoms = len(self.coords) # total atoms in .mol file
        self.nbonds = 0 # total bonds in .mol file
        self.natomtypes = 0
        self.nbondtypes = 0
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.xbox_line = ''
        self.ybox_line = ''
        self.zbox_line = ''
        self.xy = 0;
        self.xz = 0;
        self.yz = 0;
        self.header = '{} {} {}'.format('HEADER, ', os.path.basename(cif_file), ' read w/ cif2lmp')
        self.velocities = {} # { atomid : tuple(velx, vely, velz)}
        self.bond_coeffs = {}
        
        # Fill in atoms
        lx, ly, lz, yz, xz, xy = self.h
        for i, (wx, wy, wz, element) in enumerate(self.coords, start=1):
            x, y, z = self.frac_to_cart(wx, wy, wz)
            
            a = Atom()
            a.type = 1 # Set as 1 this information is not used by all2lmp
            a.molid = 1  # Set as 1 (option in main code to update later on)
            a.charge = 0 # Set as 0 (option in main code to update later on)
            a.element = element
            a.x = x - 0.5*(lx + xy + xz)
            a.y = y - 0.5*(ly + yz)
            a.z = z - 0.5*lz
            a.ix = 0 # box size will be found such that image flags from .mol files will be zero
            a.iy = 0 # box size will be found such that image flags from .mol files will be zero
            a.iz = 0 # box size will be found such that image flags from .mol files will be zero
            self.atoms[i] = a
            self.velocities[i] = (0, 0, 0)
            
        # Set box dimensions string
        xlo = -lx/2
        xhi =  lx/2
        ylo = -ly/2
        yhi =  ly/2
        zlo = -lz/2
        zhi =  lz/2
        self.xbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(xlo, xhi, 'xlo', 'xhi')
        self.ybox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(ylo, yhi, 'ylo', 'yhi')
        self.zbox_line = '{:^15.10f} {:^15.10f} {:^5} {:5}'.format(zlo, zhi, 'zlo', 'zhi')   
        self.xy = xy
        self.xz = xz
        self.yz = yz


    def deduplicate_frac_coords(self, coords, tol=1e-6):
        ndigits = abs(int(math.log10(tol)))
        dedup_coords, seen = [], set()
        for x, y, z, element in coords:
            key = (element,
                   round(x, ndigits),
                   round(y, ndigits),
                   round(z, ndigits))
            if key in seen: continue
    
            seen.add(key)
            dedup_coords.append((x, y, z, element))
    
        return dedup_coords
                
    def wrap_frac(self, x, tol=1e-8):
        x = x % 1.0
        if abs(x) < tol or abs(x - 1.0) < tol:
            x = 0.0
        return x
        
    def nuttered_eval(self, string, globals_dict=None):
        allowed = {'x': None, 'y': None, 'z': None, "__builtins__": {}}
        if globals_dict:
            allowed.update(globals_dict)
        return eval(string, allowed)

    def frac_to_cart(self, fx, fy, fz):
        lx, ly, lz, yz, xz, xy = self.h
    
        x = fx * lx + fy * xy + fz * xz
        y = fy * ly + fz * yz
        z = fz * lz
        return x, y, z
    
    def cart_to_frac(self, x, y, z):
        fx = self.h_inv[0] * x + self.h_inv[5] * y + self.h_inv[4] * z
        fy = self.h_inv[1] * y + self.h_inv[3] * z
        fz = self.h_inv[2] * z
        return fx, fy, fz
    
    def cell_to_lammps_h(self):
        """
        Convert CIF cell parameters:
            a, b, c, alpha, beta, gamma
    
        into LAMMPS triclinic components:
            lx, ly, lz, yz, xz, xy
    
        and sparse h / h_inv arrays.
        """
        a = self.a
        b = self.b
        c = self.c
    
        alpha = math.radians(self.alpha)
        beta  = math.radians(self.beta)
        gamma = math.radians(self.gamma)
    
        # LAMMPS-compatible lower-triangular cell:
        #
        # A = [lx, 0,  0 ]
        # B = [xy, ly, 0 ]
        # C = [xz, yz, lz]
        #
        lx = a
        xy = b * math.cos(gamma)
        ly = b * math.sin(gamma)
        xz = c * math.cos(beta)
        yz = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    
        lz_sq = c*c - xz*xz - yz*yz
        if lz_sq < 0 and abs(lz_sq) < 1.0e-10:
            lz_sq = 0.0
    
        if lz_sq < 0:
            raise ValueError("Invalid cell parameters: computed lz^2 < 0")
    
        lz = math.sqrt(lz_sq)
        h = [lx, ly, lz, yz, xz, xy]
    
        h_inv = 6 * [0.0]
        h_inv[0] = 1.0 / h[0]
        h_inv[1] = 1.0 / h[1]
        h_inv[2] = 1.0 / h[2]
        h_inv[3] = -h[3] / (h[1] * h[2])
        h_inv[4] = (h[3] * h[5] - h[1] * h[4]) / (h[0] * h[1] * h[2])
        h_inv[5] = -h[5] / (h[0] * h[1])
        return h, h_inv
            

################
# TEST READING #
################
if __name__ == "__main__":  
    cif_file = '../EXAMPLES/cif/1511801.cif'
    
    m = Molecule_File(cif_file)
