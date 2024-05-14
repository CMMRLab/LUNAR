# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
May 17, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    *********************************************************
    * Requirements:                                         *
    *   python numpy module:                                *
    *    - pip3 install numpy (if pip manager is installed) *
    *********************************************************
"""

##############################
# Import Necessary Libraries #
##############################
import numpy as np
import sys
import os

###################################################
# Josh's hand built linear regression model class #
###################################################
class linear_regression:
    def __init__(self, x, y):
        # Length of variables
        n = len(x)  
    
        # Sum of xi and yi 
        sum_xi = sum(x); sum_yi = sum(y);
        
        # xi and yi squared
        xi_2 = [x**2 for x in x]; yi_2 = [y**2 for y in y];
        
        # Sum of xi squared and xi*yi
        sum_xi_2 = sum(xi_2); xi_yi = [x[i]*y[i] for i in range(0, len(x))]; 
                      
        # sum of xi*yi
        sum_xi_yi = sum(xi_yi);
                
        # SSxy
        SSxy = (sum_xi_yi) - ((1/n)*(sum_xi)*(sum_yi));
        
        # SSxx
        SSxx = sum_xi_2 - (1/n)*sum_xi**2;
        
        # y bar and x bar
        y_bar = sum_yi/n; x_bar = sum_xi/n;
    
        # b0 and b1
        self.b1 = SSxy/SSxx; self.b0 = y_bar - self.b1*x_bar;
    
        # SSreg and SS total
        SSreg = (self.b1**2)*(SSxx); SStotal = sum(yi_2) - sum_yi**2/n;
        
        # r2
        self.r2 = SSreg/SStotal;

###################################
# Josh's avg/moving avg functions #
###################################
def moving_average(x, w):
  return np.convolve(x, np.ones(w), 'valid') / w

def avg(lst):
    return sum(lst)/len(lst)

###############################
# Josh's reduce data function #
###############################
def reduce_data(xdata, ydata, xlo, xhi):
    if xlo > xhi:
        print(f'ERROR no data between {xlo}-{xhi}'); sys.exit()
    reducedx = []; reducedy = [];  
    for x, y in zip(xdata, ydata):
        if x >= xlo and x < xhi:
            reducedx.append(x); reducedy.append(y);
    return reducedx, reducedy


#########################################################################
# Function to convert string to float for int, float or fraction. From: #
# https://stackoverflow.com/questions/1806278/convert-fraction-to-float #
#########################################################################
def convert_to_float(frac_str):
    try:
        return float(frac_str)
    except ValueError:
        num, denom = frac_str.split('/')
        try:
            leading, num = num.split(' ')
            whole = float(leading)
        except ValueError:
            whole = 0
        frac = float(num) / float(denom)
        return whole - frac if whole < 0 else whole + frac   


#############################################################
# Function to import file from path useage:                 #
#  module = import_file('path/to/file/default_density.py')  #
#  print(module.mode) # print mode dict                     #
#############################################################
def import_file(path):
    from importlib import util
    root = os.path.basename(path)
    root = root[:root.rfind('.')]
    spec = util.spec_from_file_location(root, path)
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module      