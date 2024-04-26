# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
April 23rd, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""


#############################################################################################################################
# Python dictionary to set element symbol from the mass of the atom if reading in a LAMMPS .data file or used to set the    #
# mass of an atom when reading in files that already have elemental information such as .mol, .mol2, ... Each key in the    #
# dictionary signifies an element and must have a value of a list with at least one mass in the list. For .mol or .mol2 or  #
# ...files, the 1st mass index will be used for mass calculation in the cluster analysis. No examples are given since this  #
# is self-explanatory, however if a read file does cannot find the information needed from this dictionary and error will   #
# be issued and the code will exit.                                                                                         #
#                                                                                                                           #
# {element symbol: [list of masses to identify element type from datafile. Minimally each list needs 1 mass]}               #
#############################################################################################################################
mass_map = {'Ag': [107.868],
            'Al': [26.98153, 26.981539, 26.98154, 26.982],
            'Ar': [39.944, 39.948],
            'As': [74.922],
            'At': [210],
            'Au': [196.967],
            'B':  [10.81],
            'Ba': [137.327, 137.33],
            'Be': [9.012182, 9.0122],
            'Bi': [208.98],
            'Br': [79.904, 79.909],
            'C':  [12.0000, 12.01115, 12.01100, 12.0112, 10.01115],
            'Ca': [40.078, 40.0798, 40.08],
            'Cl': [35.45, 35.4527, 35.453],
            'Co': [58.9332],
            'Cr': [51.996],
            'Cs': [132.90543, 132.91],
            'Cu': [63.546],
            'D':  [2.01400],
            'F':  [18.998, 18.9984, 18.998403],
            'FI': [289],
            'Fe': [55.847],
            'Fr': [223],
            'Ga': [69.723],
            'Ge': [72.61, 72.63],
            'H':  [1.008, 1.0, 1.00782, 1.00797],
            'He': [4.0026, 4.003],
            'I':  [126.9, 126.9044, 126.90447, 126.9045],
            'In': [114.82],
            'K':  [39.098, 39.0983, 39.1, 39.102],
            'Kr': [83.798, 83.8],
            'L':  [1.0],
            'Li': [6.94, 6.941],
            'Lp': [0.5],
            'Lv': [293],
            'Mc': [290],
            'Mg': [24.305],
            'Mn': [54.93805],
            'Mo': [95.95, 95.94],
            'N':  [14.0, 14.0067, 14.00674, 14.01],
            'Na': [22.98977, 22.9898, 22.99],
            'Ne': [20.18, 20.183],
            'Nh': [286],
            'Ni': [58.69, 58.71],
            'O':  [15.999, 14.9994, 15.99491, 15.9994, 16.0],
            'Og': [294],
            'P':  [30.9738, 30.974, 30.9737],
            'Pb': [207.2],
            'Pd': [106.4, 106.42],
            'Po': [209],
            'Pt': [195.09],
            'Ra': [226],
            'Rb': [85.4678, 85.468],
            'Rn': [222],
            'S':  [32.06, 32.064, 32.066],
            'Sb': [121.76],
            'Se': [78.971],
            'Si': [28.085, 28.0855, 28.086, 28.086],
            'Sn': [118.69, 118.71],
            'Sr': [87.62],
            'TI': [204.38],
            'Te': [127.6],
            'Ti': [47.88],
            'Ts': [294],
            'V':  [50.9415],
            'W':  [183.85],
            'Xe': [131.29, 131.3],
            'mg': [100],
           }

            