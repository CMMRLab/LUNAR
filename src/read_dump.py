# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
July 9, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
from itertools import islice


#-------------------------------------------------------------#
# Class to parse and read one portion of a dumpfile at a time #
#-------------------------------------------------------------#
class read_dumpfile:
    def __init__(self, dumpfile):
        self.filename = dumpfile
        self.steps = {} # {step:(step: (start_line, end_line))}
        
        # Generate a bare minimum dump object
        self.parse_timesteps()
    
    # Convert a string to int, float, or leave as str.
    def convert_value(self, string):
        if not isinstance(string, str):
            return string
        
        if any(c in string for c in ('.', 'e', 'E')):
            try:
                return float(string)
            except ValueError:
                return string
    
        try:
            return int(string)
        except ValueError:
            return string
    
    # Sets start and end of line for each section
    def parse_timesteps(self):
        with open(self.filename, 'r') as f:
            # Get the starting section indexes
            steps = [] # [ (step, starting_line_index) ]
            item = ''
            for n, line in enumerate(f):
                if line.startswith('ITEM'):
                    item = line.split('ITEM:')[-1].strip()
                    continue
                    
                if item == 'TIMESTEP':
                    line = line.strip()
                    step = self.convert_value( line )
                    steps.append( (step, n) )
                    
            # Set file max index
            max_index = n + 1
                    
            # Map the starting section indexes to the
            # starting and stopping section indexes
            for i, (step, start) in enumerate(steps):
                if i < len(steps)-1:
                    end = steps[i+1][1] - 1
                else:
                    end = max_index
                self.steps[step] = (start, end)
        return
    
    # Read section into a list
    def read_step(self, step):
        if step not in self.steps:
            raise KeyError(f"Timestep {step} not found in {self.filename}")
        
        start, end = self.steps[step]
        lines = []
        with open(self.filename) as f:
            lines = list(islice(f, start-1, end))
        return lines
    
    # Parse out box dimensions and atom info
    def parse_atoms_and_box(self, step):
        #-----------------------------#
        # Read section and parse info #
        #-----------------------------#
        section = self.read_step(step)
        natoms = 0
        item = ''
        box_lines, atom_lines, atom_header = [], [], ''
        for line in section:
            line = line.strip()
            if line.startswith('ITEM'):
                item = line.split('ITEM:')[-1].strip()
                if 'ATOMS' in item: atom_header = item
                continue

            if 'NUMBER OF ATOMS' in item:
                natoms = int(line)
                
            if 'BOX BOUNDS' in item:
                box_lines.append(line.split())
                
            if 'ATOMS' in item and 'NUMBER OF' not in item:
                atom_lines.append(line)
                
        if len(atom_lines) != natoms:
            raise Exception(f'ERROR parsing step: {step} from {self.filename}. Number of atoms in section does not match "NUMBER OF ATOMS" ITEM.')
        
        #------------------------------------------------------#
        # Convert box lines to cell dimensions like a datafile #
        #------------------------------------------------------#
        # Set defaults
        xlo, xhi = -0.5, 0.5
        ylo, yhi = -0.5, 0.5
        zlo, zhi = -0.5, 0.5
        xy, xz, yz = 0, 0, 0
        
        # Orthogonal box
        if len(box_lines) == 3 and len(box_lines[0]) == 2:
            xlo = float(box_lines[0][0])
            xhi = float(box_lines[0][1])
        
            ylo = float(box_lines[1][0])
            yhi = float(box_lines[1][1])
        
            zlo = float(box_lines[2][0])
            zhi = float(box_lines[2][1])
        
        # Triclinic box from dump bounds
        elif len(box_lines) == 3 and len(box_lines[0]) == 3:
            xlo_bound = float(box_lines[0][0])
            xhi_bound = float(box_lines[0][1])
            xy        = float(box_lines[0][2])
        
            ylo_bound = float(box_lines[1][0])
            yhi_bound = float(box_lines[1][1])
            xz        = float(box_lines[1][2])
        
            zlo_bound = float(box_lines[2][0])
            zhi_bound = float(box_lines[2][1])
            yz        = float(box_lines[2][2])
        
            xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
            xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
        
            ylo = ylo_bound - min(0.0, yz)
            yhi = yhi_bound - max(0.0, yz)
        
            zlo = zlo_bound
            zhi = zhi_bound
        
        # Set final box
        box = {'xlo': xlo, 'xhi': xhi,
               'ylo': ylo, 'yhi': yhi,
               'zlo': zlo, 'zhi': zhi,
                'xy': xy, 'xz': xz, 'yz': yz}
        
        #---------------------------------#
        # Parse atom info from atom_lines #
        #---------------------------------#
        atoms = [] # [ {'column1', value1, 'column2': value2, ....}, ... Natoms ]
        columns = atom_header.split()[1:]
        for line in atom_lines:
            split = line.split()
            if len(columns) != len(split): 
                raise ValueError(f'Expected {len(columns)} atom fields but found {len(split)}. Skipping row in file.')
            atom = {col : self.convert_value(val) for col, val in zip(columns, split)}
            atoms.append( atom )

        return box, atoms
        
    
    
# Test reading dumpfile
if __name__ == "__main__": 
    dumpfile = '../EXAMPLES/free_volume/wildcard_DETDA/detda_rep_1_free_volume_test.dump'
    dump = read_dumpfile(dumpfile)
    
    box, atoms = dump.parse_atoms_and_box(0)
    print(f"{box['xlo']: .16e} {box['xhi']: .16e} xlo xhi")
    print(f"{box['ylo']: .16e} {box['yhi']: .16e} ylo yhi")
    print(f"{box['zlo']: .16e} {box['zhi']: .16e} zlo zhi")
    if box['xy'] != 0.0 or box['xz'] != 0.0 or box['yz'] != 0.0:
        print(f"{box['xy']: .16e} {box['xz']: .16e} {box['yz']: .16e} xy xz yz")
        
    for atom in atoms:
        print(atom)

    
    
    