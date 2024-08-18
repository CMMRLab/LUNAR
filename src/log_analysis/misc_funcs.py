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
import matplotlib.pyplot as plt
import numpy as np
import math
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
        
        # lineare regression data
        self.xreg = [min(x), max(x)]
        self.yreg = [i*self.b1 + self.b0 for i in self.xreg]
        

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
        print(f'ERROR (reduce_data) no data between {xlo}-{xhi}')
    reducedx = []; reducedy = [];  
    for x, y in zip(xdata, ydata):
        if x >= xlo and x <= xhi:
            reducedx.append(x); reducedy.append(y);
    return reducedx, reducedy


##########################################################################
# Function to compute vectorized operations on arrays using string math. #
#  Inputs:                                                               #
#     data = {'NAME':[lst], ... }                                        #
#     compute_string = "${A} + ${B}/2 + 100", where A and B are keys in  #
#                       the data dictionary.                             #
#  Outputs:                                                              #
#    computed_list = "${A} + ${B}/2 + 100", using vectorization methods  #
#                    from numpy and Python eval() function.              #
##########################################################################
def compute_thermo_data(data, compute_string):
    #--------------------------------------#
    # d will be used in list(eval(string)) #
    #--------------------------------------#
    d = {i:np.array(data[i]) for i in data}
    nzeros = max([len(d[i]) for i in d])
    message = ''
    
    #------------------------------------------------------------------#
    # The compute_string needs to be transformed, for example:         #
    #  in-string:  "  ${A} + ${B}/2   + 100"                           #
    #  out-string: "d['A'] + d['B']/2 + 100"                           #
    # where A and B are columns in the logfile, stored as keys in the  #
    # dictionary d = {'COLUMN-NAME':np.array()}. Using the Python      #
    # eval() function we can then perform vectorized computes of the   #
    # numpy arrays.                                                    #
    #------------------------------------------------------------------#
    cstring = compute_string
    cstring = cstring.replace("${", "d['")
    cstring = cstring.replace("}", "']")
    
    #---------------------------------------------------------------#
    # Check that each variable in the compute_string exists in the  #
    # dictionary "d". Warn if not in the dictionary.                #
    #---------------------------------------------------------------#
    # Find variables from compute_string
    variables = set(); tmp = ''; dollar_sign = False; bracket0 = False; bracket1 = False;
    for i in compute_string:
        if i == '$': dollar_sign = True
        if i == '{': bracket0 = True
        if i == '}': bracket1 = True
        if dollar_sign and bracket0:
            tmp += str(i)
        if bracket1:
            dollar_sign = False; bracket0 = False; bracket1 = False;
            tmp = tmp.replace('{', '')
            tmp = tmp.replace('}', '')
            variables.add(tmp); tmp = '';
            
    # Check that each variable is in the "d" dictionary
    variable_check = True
    for variable in variables:
        if variable not in d:
            variable_check = False
            message += 'ERROR {} is not in the compute dictionary. Perhaps something was misspelled.\n'.format(variable)     
        
    #---------------------------------------------------------------------#
    # Try making compute, if it fails log a different message and create  #
    # an array of zeros to output.                                        #
    #---------------------------------------------------------------------#
    if variable_check:
        try:
            new_data = list(eval(cstring))
            message += f'{compute_string} was succesfully evaluted.\n'
        except:
            new_data = nzeros*[0]
            message += f'ERROR could not evalute {compute_string} in internal compute.\n'
    else: new_data = nzeros*[0]
    return new_data, message


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
    
    
##############################################
# Function to define of a butterworth filter #
##############################################
def butter_lowpass_filter(xdata, ydata, wn, order):
    from scipy.signal import butter, filtfilt
    
    # Get the filter coefficients and filter ydata
    b, a = butter(order, wn, btype='low', analog=False)
    y = filtfilt(b, a, np.array(ydata))
    return y


#######################################################
# Function to compute integral using Trapezoidal rule #
#######################################################
def compute_integral(x, y):
    ix = []; iy = []; summation = 0
    if len(x) == len(y):
        for i in range(len(x)-1):
            xdiff = x[i+1] - x[i]
            yprod = y[i+1] + y[i]
            summation += (yprod/2)*xdiff
            ix.append( (x[i+1] + x[i])/2 )
            iy.append(  summation )
    else: print('ERROR (compute_integral) inconsistent number of data points between X and Y arrays')
    return ix, iy


#######################################################
# Function to compute derivative using Euler's method #
#######################################################
def compute_derivative(x, y):
    dx = []; dy = [];
    if len(x) == len(y):
        for i in range(len(x)-1):
            xdiff = x[i+1] - x[i]
            ydiff = y[i+1] - y[i]
            dx.append( (x[i+1] + x[i])/2 )
            dy.append( ydiff/xdiff )
    else: print('ERROR (compute_derivative) inconsistent number of data points between X and Y arrays')
    return dx, dy


#########################################################################
# Function to fit y = c0 + c1*x + ... + c_n*x^n to x and y numpy arrays #
#########################################################################
def np_curve_fit(x, y, degree=2, domain=[]):
    # Define polynomial function
    def poly(x, coef):
        yout = np.zeros_like(x)
        for n in range(len(coef)):
            yout += coef[n]*(x**n)
        return yout
    
    # Find domain
    if domain != []:
        indices = np.where( np.logical_and(x >= domain[0], x <= domain[1]) )[0]
        x = x[indices]; y = y[indices];
    else:
        x = x.copy(); y = y.copy();
    
    # Find best fit
    param = np.polynomial.polynomial.polyfit(x, y, degree)
    yfit = poly(x, param)
    return x, y, yfit, param


####################################################
# Function to find peaks and valleys in a data set #
####################################################
def find_peaks_and_valleys(x, y):
    from scipy.signal import find_peaks
    
    # Find peaks
    xdata = np.array(x); ydata = np.array(y);
    peaks, properties = find_peaks(ydata)
    xpeaks = list(xdata[peaks])
    ypeaks = list(ydata[peaks])
    
    # Find valleys
    xvalleys = []; yvalleys = [];
    if len(xpeaks) >= 2 and len(ypeaks) >= 2:
        for i in range(len(peaks)-1):
            lo = peaks[i]; hi = peaks[i+1];
            between_peaksx = xdata[lo:hi]
            between_peaksy = ydata[lo:hi]
            minimum_index = np.min(np.where(between_peaksy == between_peaksy.min())[0])
            xvalleys.append( between_peaksx[minimum_index] )
            yvalleys.append( between_peaksy[minimum_index] )
    return xpeaks, ypeaks, xvalleys, yvalleys


############################################################################
# Function below implements a "piecewise regression" on the data. where    #
# n is the number of breakpoints between the start and end point of the    #
# input data.                                                              #
#                                                                          #
# This method uses both numpy and pwlf. numpy is a global requirement, but #
# pwlf is not, so import within the method to avoid dependancies if this   #
# option is not planned to be used.                                        #
############################################################################
def piecewise_regression(x, y, xlo, xhi, n):
    import pwlf
    npx = np.array(x); npy = np.array(y);
    
    # Perform peicewise regression
    my_pwlf = pwlf.PiecewiseLinFit(npx, npy)
    xbreaks = my_pwlf.fit(n+1)
    ybreaks = my_pwlf.predict(np.array(xbreaks))
    xout = np.linspace(npx.min(), npx.max(), 100)
    yout = my_pwlf.predict(xout)
    
    # Compute the slopes
    def slope(x1, y1, x2, y2):
        m = (y2-y1)/(x2-x1)
        return m
    
    npoints = min([len(xbreaks), len(ybreaks)])-1
    slopes = {(i, i+1):slope(xbreaks[i], ybreaks[i], xbreaks[i+1], ybreaks[i+1]) for i in range(npoints)}
    return xout, yout, xbreaks, ybreaks, slopes


######################################################################################
# The method below implements the Tg/CTE hyperbola fit from the following paper:     #
# Uncertainty quantification in molecular dynamics studies of the glass transition   #
# temperature - Paul N. Patrone, Andrew Dienstfrey, ... - Polymer Volume 87 - 2016   #
#                                                                                    #
# Hyperbola fit for Tg sims. Parameters optimized by scipy.optmize.curve_fit:        #
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html #
#    maxiter sets a cutoff for the maximum number of iterations                      #
#    intial_guess list of parameters 1-N or None, if None scipy intializes the guess #
######################################################################################
def fit_hyperbola(x, y, xlo, xhi, minimum_convergence=None, initial_guess=False, maxiter=10**4):
    from scipy import optimize
      
    # define default outputs (in-case something goes wrong)
    xout = list(x); yout = len(y)*[0];
    params = [0, 0, 0, 0, 0]; center = [0, 0];
    slopes = [0, 0]; transition = [0, 0];
    
    # Convert to float64 and then lists to numpy arrays
    xx = np.array(x); yy = np.array(y);
    
    # Setup intial guess
    p0 = None
    if initial_guess:
        slopeguess = (yy[-1]-yy[0])/(xx[-1]-xx[0])
        p0 = (np.mean(xx),np.mean(yy), slopeguess, slopeguess, np.log((xx[-1]-xx[0])**2/100))

    # Define the hyperbola equation (eqn 1 in paper)
    def hyberbola(t, t0, v0, a, b, c):
        h0 = 0.5*(t-t0) + np.sqrt( ((t-t0)*(t-t0))/4 + np.exp(c) )
        v = v0 + a*(t-t0) + b*h0
        return v
    
    # Find best fit    t0          v0        a        b        c
    parm_bounds = ((np.min(xx), np.min(yy), -np.inf, -np.inf, -np.inf), # lower
                   (np.max(xx), np.max(yy),  np.inf,  np.inf,  np.inf)) # upper
    param, param_cov = optimize.curve_fit(hyberbola, xx, yy, p0=p0, method='trf', bounds=parm_bounds, maxfev=maxiter)
    t0, v0, a, b, c = param # extract out fitting coeffs
    
    # update defaults
    params = list(param)
    yout = list(hyberbola(xout, *param))
    center = [t0, v0]
    slopes = [a, (a+b)]
    
    # Find where function is transitioning (eqn 5 in paper)
    if minimum_convergence is not None:
        dtemp = (np.exp(c/2)*(2*minimum_convergence-1))/(np.sqrt(minimum_convergence*(1-minimum_convergence)))
        transition = [t0-dtemp, t0+dtemp]
    return xout, yout, params, center, slopes, transition


##########################################
# Function to find closest value in list #
##########################################
def closest(lst, value):
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-value))]


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