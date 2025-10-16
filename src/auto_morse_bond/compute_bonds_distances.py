# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 21:09:15 2022

@author: jdkem
"""

####################
# Import Libraries #
####################
import math
 
# Function to find distance
def compute_distance(x1, y1, z1, x2, y2, z2):
    dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
    return math.sqrt(dx*dx + dy*dy + dz*dz)

################################################################
# Statistic's functions to use for analyzing bond stats. *NOTE #
# not using numpy to make this code have zero dependancies*    #
################################################################
def compute_mean(data):
  return sum(data)/len(data)
 
def compute_variance(data):
  mean = compute_mean(data)
  deviations = [(x - mean)**2 for x in data]
  variance = sum(deviations)/len(data)
  return variance
 
def compute_standard_deviation(data):
  variance = compute_variance(data)
  return math.sqrt(variance)

####################################################
# Class to compute bond distances and analyze them #
####################################################
class Types: pass # .symbols .count .avg .min .max .std
class compute:
    def __init__(self, m):
        self.bonds = {} # { bondid : bond distance }
        self.types = {} # { bond type id : type object }
        
        # Find box dimensions to remove periodic boundary conditions
        xline = m.xbox_line.split()
        yline = m.ybox_line.split()
        zline = m.zbox_line.split()
        lx = float(xline[1])-float(xline[0])
        ly = float(yline[1])-float(yline[0])
        lz = float(zline[1])-float(zline[0])
        
        # Initialize self.types with zeros
        bond_type_holder = {i:[] for i in m.bond_coeffs} # { bond type id : [list of all lengths] }
        for i in m.bond_coeffs:
            t = Types()
            t.symbols = m.bond_coeffs[i]
            t.count = 0
            t.avg = 0
            t.min = 0
            t.max = 0
            t.std = 0
            t.lst = []
            self.types[i] = t
       
        
        # Loop through bonds
        for i in m.bonds:
            bond = m.bonds[i]
            id1, id2 = bond.atomids
            
            # Find id1 and id2 x, y, and z
            id1x = m.atoms[id1].x
            id1y = m.atoms[id1].y
            id1z = m.atoms[id1].z
            
            id2x = m.atoms[id2].x
            id2y = m.atoms[id2].y
            id2z = m.atoms[id2].z
            
            # Shift atoms using minimum image convention if bond
            # is periodic (always shifting id2 atom if needed)
            bondlx = abs(id1x - id2x)
            if bondlx > 0.5*lx:
                diff = id1x - id2x
                if diff < 0:
                    id2x = id2x - lx
                elif diff > 0:
                    id2x = id2x + lx
                else:
                    id2x = id2x
                    
            bondly = abs(id1y - id2y)
            if bondly > 0.5*ly:
                diff = id1y - id2y
                if diff < 0:
                    id2y = id2y - ly
                elif diff > 0:
                    id2y = id2y + ly
                else:
                    id2y = id2y

            bondlz = abs(id1z - id2z)
            if bondlz > 0.5*lz:
                diff = id1z - id2z
                if diff < 0:
                    id2z = id2z - lz
                elif diff > 0:
                    id2z = id2z + lz
                else:
                    id2z = id2z

                
            # compute bond distance
            dist = compute_distance(id1x, id1y, id1z, id2x, id2y, id2z)
            
            # Add bond distance to self.bonds and append to holder based on type-id
            self.bonds[i] = dist
            bond_type_holder[bond.type].append(dist)
            
        ###########################
        # Find bonding statistics #
        ###########################
        for i in bond_type_holder:
            lst = bond_type_holder[i] 
            try:
                self.types[i].count = len(lst)
                self.types[i].avg = float('{:.4f}'.format( compute_mean(lst) ))
                self.types[i].min = float('{:.4f}'.format( min(lst) ))
                self.types[i].max = float('{:.4f}'.format( max(lst) ))
                self.types[i].std = float('{:.4f}'.format( compute_standard_deviation(lst) ))
                self.types[i].lst = lst
            except: pass
