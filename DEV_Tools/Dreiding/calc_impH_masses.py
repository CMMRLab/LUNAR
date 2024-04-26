# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 11:42:54 2023

@author: jdkem
"""


element_mass = 15.999000
impH_mass = 1.008
numH = 1



mass = round(element_mass + numH*impH_mass, 6)

print('\n\n\nMass with implicit Hydrogens')
print('{:^10.6f}\n\n\n'.format(mass))


