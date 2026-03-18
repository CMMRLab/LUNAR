# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
November 17, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.read_lmp as read_lmp
import src.masses as masses




class LmpData:
    def __init__(self):
        self.mass_map = masses.mass_map 
        
        
    def _read_lmp(self, topofile):
        return read_lmp.Molecule_File(topofile, method='forward', sections = ['Atoms', 'Bonds', 'Velocities'])
    
    def _get_method_from_string(self, string):
        # Example: 'elements_from_mass(graphite_AB_25x25x25_PCFF.data, 1, Ar)'
        # Output:
        #  method = 'elements_from_mass'
        #  args = ['graphite_AB_25x25x25_PCFF.data', 1, Ar]
        split = string.split('(')
        method = split[0].strip()
        args = [self._string2digit(i.strip()) for i in split[-1].replace(')', '').split(',')]
        
        # Find all methods in current class
        methods = [name for name in dir(self)
                   if callable(getattr(self, name)) and not name.startswith('__')]
        
        # If method is not defined set method as empty string
        if method not in methods:  method = ''
        return method, args
    
    # Function to convert string to float or int and will default to string if needed
    def _string2digit(self, string):
        digit = string
        try: 
            digit = float(string)
            if digit.is_integer():
                digit = int(digit)
        except: pass
        return digit
    
    def elements_from_masses(self, *args):
        # Return if incorrect number of args
        if len(args) < 1: return
        
        # Read the datafile
        topofile = args[0]
        m = self._read_lmp(topofile)
        
        # Attempt getting elements from masses
        elements = {} # {atomTypeID:'element'}
        for i in m.masses:
            mass = m.masses[i].coeffs[0]
            
            # Find closest element based on mass
            diffs = {} # {element: abs(mass - mass_element)}
            for element in self.mass_map:
                masses = self.mass_map[element]
                diffs[element] = min([abs(mass - value) for value in masses])
                
            elements[i] = min(diffs, key=diffs.get)
            
        # Generate element string
        atom_type_ids = sorted(m.masses.keys())
        sorted_elements = [str(elements[i]) for i in atom_type_ids]
        additional_elements = [str(i) for n, i in enumerate(args) if n > 0]
        string = '"{}"'.format( ' '.join(sorted_elements + additional_elements) )
        return string
    
    
if __name__ == "__main__":  
    lmp_data = LmpData()
    elements = lmp_data.elements_from_masses('graphite_AB_25x25x25_PCFF.data')
    print(elements)
    

    # Call method from class based on string
    string = "elements_from_masses(graphite_AB_25x25x25_PCFF.data, Ar)"    
    method_name, args = lmp_data._get_method_from_string(string)
    print(method_name, args)
    if method_name:
        output = getattr(lmp_data, method_name)(*args)
    else:
        output = None
    print(output)
    
        