# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 20, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to convert string to float or int and will default to string if needed
def string2digit(string):
    digit = string
    try: 
        digit = float(string)
        if digit.is_integer():
            digit = int(digit)
    except: pass
    return digit

# Class to read LUNAR/atom_typing log file
class logfile:
    def __init__(self, logfile):
        self.mass = 0 # Will be tallied after file is parsed
        self.filename = logfile
        self.elements = [] # ['element1', 'element2', ...]
        self.rings = {} # { Ring-size : {{'column-name': column-value}}
        self.clusters = {} # { molID : {'column-name': column-value} }
        self.bond_order = {} # { 'element1-elemment2' : {'column-name': column-value} }
        self.bond_length = {} # { 'element1-elemment2' : {'column-name': column-value} }
        self.hybridizations = {} # { 'Atom-type' : {'column-name': column-value} }
        self.ringed_clusters = {} # { 'Ring-type' : {'column-name': column-value} }
        self.atom_bond_order = {} # { 'element' : {'column-name': column-value} }
        self.atom_types_tally = {} # { 'Atom-type' : {'column-name': column-value} }
        self.by_products_tally = {} # { 'Type' : {'column-name': column-value} }
        self.fused_ringed_clusters = {} # { FusedID : {'column-name': column-value} }
        
        # Initialize flags
        columns = []
        table_body = False
        ring_table = False
        elements_list = False
        cluster_table = False
        atom_types_table = False
        by_products_table = False
        bond_order_table = False
        bond_length_table = False
        hybridization_table = False
        ringed_clusters_table = False
        atom_bond_order_table = False
        fused_ringed_clusters_table = False
        relative_line_number = 0
        
        # Parse logfile
        with open(logfile, 'r') as f:
            for linenumber, line in enumerate(f, 1):
                line = line.strip()
                
                # Start setting flags
                if line == '':
                    columns = []
                    table_body = False
                    ring_table = False
                    elements_list = False
                    cluster_table = False
                    atom_types_table = False
                    by_products_table = False
                    bond_order_table = False
                    bond_length_table = False
                    hybridization_table = False
                    ringed_clusters_table = False
                    atom_bond_order_table = False
                    fused_ringed_clusters_table = False
                    relative_line_number = 0
                elif 'Elements found in system:' in line:
                    table_body = False
                    elements_list = True
                    continue
                elif 'Cluster Analysis' in line:
                    table_body = False
                    cluster_table = True
                    continue
                elif 'Atom types Tally' in line:
                    table_body = False
                    atom_types_table = True
                    continue
                elif 'pdb_file' in line:
                    atom_types_table = False
                    continue
                elif 'By Products Tally' in line:
                    table_body = False
                    by_products_table = True
                    continue
                elif 'Bond type bond order statistics' in line:
                    table_body = False
                    bond_order_table = True
                    continue
                elif 'Bond type bond length statistics' in line:
                    table_body = False
                    bond_length_table = True
                    continue
                elif 'Element abo statistics' in line:
                    table_body = False
                    atom_bond_order_table = True
                    continue
                elif 'Hybridization Information' in line:
                    table_body = False
                    hybridization_table = True
                    continue
                elif 'Ringed Clusters' in line and 'Fused' not in line:
                    table_body = False
                    ringed_clusters_table = True
                    continue
                elif 'Fused Ring Clusters' in line:
                    table_body = False
                    fused_ringed_clusters_table = True
                    continue
                elif 'Ring' in line and 'Count' in line and '%Ring' in line:
                    ring_table = True
                    relative_line_number = 0
                
                # Start parsing sections
                if ring_table:
                    columns = self.parse_ring(line, columns, relative_line_number)
                    relative_line_number += 1
                elif elements_list:
                    self.elements.append(line.split()[-1])
                elif cluster_table:
                    columns, table_body = self.parse_table(line, self.clusters, columns, table_body)       
                elif atom_types_table:
                    columns, table_body = self.parse_table(line, self.atom_types_tally, columns, table_body)  
                elif by_products_table:
                    columns, table_body = self.parse_table(line, self.by_products_tally, columns, table_body) 
                elif bond_order_table:
                    columns, table_body = self.parse_table(line, self.bond_order, columns, table_body)
                elif bond_length_table:
                    columns, table_body = self.parse_table(line, self.bond_length, columns, table_body)
                elif atom_bond_order_table:
                    columns, table_body = self.parse_table(line, self.atom_bond_order, columns, table_body)
                elif hybridization_table:
                    columns, table_body = self.parse_table(line, self.hybridizations, columns, table_body)
                elif ringed_clusters_table:
                    columns, table_body = self.parse_table(line, self.ringed_clusters, columns, table_body)
                elif fused_ringed_clusters_table:
                    columns, table_body = self.parse_table(line, self.fused_ringed_clusters, columns, table_body)
                    
        # Go through and sum up all hybridizations to get system mass
        if self.hybridizations:
            for hybrid in self.hybridizations:
                self.mass += self.hybridizations[hybrid]['Mass']
                
    ################################
    # Define methods to this class #
    ################################
    def parse_table(self, line, table_dict, columns, table_body):  
        if '---' in line: table_body = True
        if not table_body:
            if not columns:
                columns = [i.strip() for i in line.split('  ') if i != '']
            else:
                tmp = [i.strip() for i in line.split('  ') if i != '']
                columns = ['{} {}'.format(i, j) for i, j in zip(columns, tmp)]
        elif '----' not in line:
            line_split = line.split()
            key = string2digit(line_split[0])
            table_dict[key] = {key:string2digit(value) for n, (key, value) in enumerate(zip(columns, line_split)) if n > 0}
        return columns, table_body
                
    def parse_ring(self, line, columns, relative_line_number):
        # Get True ring data
        if relative_line_number == 0:
            columns = line.replace('|', '').split()
        elif relative_line_number == 3:
            line_split = line.replace('|', '').split()
            ring = string2digit(line_split[0])
            self.rings[ring] = {key:string2digit(value) for n, (key, value) in enumerate(zip(columns, line_split)) if n > 0}
            
        # Get Partitioned ring data
        elif relative_line_number == 5:
            columns = line.replace('|', '').split()
        elif relative_line_number >= 7 and '----' not in line:
            line_split = line.replace('|', '').split()
            ring = list(self.rings.keys())[-1]
            element = string2digit(line_split[0])
            partitioned_rings = {element:{key:string2digit(value) for n, (key, value) in enumerate(zip(columns, line_split)) if n > 0}}
            self.rings[ring].update(partitioned_rings)
        return columns
                

    

                
# Basic test of the log file reader             
if __name__ == "__main__":
    LUNAR_path = '../../../EXAMPLES/array_processing/atom_typing_logfile_processor/poly_tracking_replicate_1_time_0ps_typed.log.lunar'
    log = logfile(LUNAR_path)
    
    # Show how to access elements list
    print('\n\n')
    print('Elements list: ', log.elements)
    print('System mass  : ', log.mass)
    
    # Show how to access Table: Cluster Data
    if log.clusters:
        print('\nCluster analysis')
        for molID in log.clusters:
            print('{} -> {}'.format(molID, log.clusters[molID]))
        
    # Show how to access Table: Bond Order Data
    if log.bond_order:
        print('\nBond Order')
        for bond_elements in log.bond_order:
            print('{} -> {}'.format(bond_elements, log.bond_order[bond_elements]))
        
    # Show how to access Table: Atom Bond Order Data
    if log.atom_bond_order:
        print('\nAtom Bond Order')
        for element in log.atom_bond_order:
            print('{} -> {}'.format(element, log.atom_bond_order[element]))
            
    # Show how to access Table: Bond Length Data (if bonds_via_distance is used)
    if log.bond_length:
        print('\nBond Length')
        for bond_elements in log.bond_length:
            print('{} -> {}'.format(bond_elements, log.bond_length[bond_elements]))
        
    # Show how to access Table: Hybridization Data
    if log.hybridizations:
        print('\nHybridization analysis')
        for hybrid in log.hybridizations:
            print('{} -> {}'.format(hybrid, log.hybridizations[hybrid]))
        
    # Show how to access Table: Ring Data
    if log.rings:
        print('\nRing analysis')
        for ring in log.rings:
            print('{} -> {}'.format(ring, log.rings[ring]))
        
        print('\nMANUAL Ring analysis')
        print("log.rings[6]['C']['%Mass']", log.rings[6]['C']['%Mass'])
            
    # Show how to access Table: Atom Type Tally
    if log.atom_types_tally:
        print('\nAtom Types Tally')
        for ring in log.atom_types_tally:
            print('{} -> {}'.format(ring, log.atom_types_tally[ring]))
            
    # Show how to access Table: By Products Tally
    if log.by_products_tally:
        print('\nBy Products Tally')
        for by_product in log.by_products_tally:
            print('{} -> {}'.format(by_product, log.by_products_tally[by_product]))
            
    # Show how to access Table: Ringed Cluster Data
    if log.ringed_clusters:
        print('\nRinged Clusters')
        for ring_type in log.ringed_clusters:
            print('{} -> {}'.format(ring_type, log.ringed_clusters[ring_type]))
            
    # Show how to access Table: Fused Ring Cluster Data
    if log.fused_ringed_clusters:
        print('\nFused Ringed Clusters')
        for fusedID in log.fused_ringed_clusters:
            print('{} -> {}'.format(fusedID, log.fused_ringed_clusters[fusedID]))
            
    # Try getting all table values in string format for
    # atom_typing_log_processor.py logger = 'all-table:TABLE-NAME'.
    def table_expand(table, dict_name):
        all_table = []
        for valueID0 in table:
            for valueID1 in table[valueID0]:
                # ring's dictionary may be nested  "3-deep"
                if isinstance(table[valueID0][valueID1], dict):
                    for valueID2 in table[valueID0][valueID1]:
                        if isinstance(valueID0, str):
                            index0 = "'{}'".format(valueID0)
                        else: index0 = "{}".format(valueID0)
                        if isinstance(valueID1, str):
                            index1 = "'{}'".format(valueID1)
                        else: index1 = "{}".format(valueID1)
                        if isinstance(valueID2, str):
                            index2 = "'{}'".format(valueID2)
                        else: index2 = "{}".format(valueID2)
                        all_table.append("{}[{}][{}][{}]".format(dict_name, index0, index1, index2))
                else: # every other dictionary will be nested "2-deep"
                    if isinstance(valueID0, str):
                        index0 = "'{}'".format(valueID0)
                    else: index0 = "{}".format(valueID0)
                    if isinstance(valueID1, str):
                        index1 = "'{}'".format(valueID1)
                    else: index1 = "{}".format(valueID1)
                    all_table.append("{}[{}][{}]".format(dict_name, index0, index1))
        return all_table
    
    print('\n\nTesting all_table: hybridizations')
    all_table = table_expand(log.hybridizations, 'hybridizations')
    print(all_table[:10])
    
    print('\n\nTesting all_table: bond_order')
    all_table = table_expand(log.bond_order, 'bond_order')
    print(all_table[:10])
    
    print('\n\nTesting all_table: ringed_clusters')
    all_table = table_expand(log.ringed_clusters, 'ringed_clusters')
    print(all_table[:10])
    
    print('\n\nTesting all_table: clusters')
    all_table = table_expand(log.clusters, 'clusters')
    print(all_table[:10])

    print('\n\nTesting all_table: rings')
    all_table = table_expand(log.rings, 'rings')
    print(all_table[:10])
    
    # print('\nMANUAL Ring analysis')
    # print("log.rings[3]['Count']", log.rings[3]['Count'])
    # print("log.rings[3]['C']['%Mass']", log.rings[3]['C']['%Mass'])
    
    # print('\n\nTesting all_table: atom_types_tally')
    # all_table = table_expand(log.atom_types_tally, 'atom_types_tally')
    # print(all_table[:10])
    
    print('\n\nTesting all_table: clusters')
    all_table = table_expand(log.clusters, 'clusters')
    print(all_table[:10])
    
    print('\nMANUAL Ring analysis')
    print("clusters[1]['Molecule Formula']['3']", log.clusters[1]['Molecule Formula'])
    
        
        