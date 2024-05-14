# -*- coding: utf-8 -*-
"""
Created by log_analysis.py to store a mode dictionary
called "mode", which will then allow all settings to be
loaded by clicking on this file within log_anaylsis.py.

Date: 05/13/2024
Author: Josh Kemppainen
Purpose: Mode to von misses stress from shear sim
         in XZ-dir


"""
#----------------------------------------------------------------------------------------------------------------#
# This is special mode as it will show how to use strings to perform vectorized manipulation of arrays to compute#
# a new array of data to plot. In general the 'xdata' or 'ydata' keys should be columns in the logfile, however  #
# users may define a compute by setting a string as 'compute: MATH', where MATH is math that the Python eval()   #
# function can handle. A little bit of understanding needs to be understood about how things are setup in the    #
# code to make use of this option detailed below:                                                                #
#   - data is setup in a dictionary as d = {'COLUMN-name': numpy.array(data from logfile)}                       #
#   - numpy can perform vectorized math operations where the + or - or * or / or ** operations can be used       #
#                                                                                                                #
# As an example say a LAMMPS logfile had columns x and y, with x-values = [1,2,3] and y-values = [4,5,6]. The d  #
# dictionary is setup as d = {'x': numpy.array([1,2,3]), 'y': numpy.array([4,5,6])}. Where any math operation    #
# that is compatable with the eval() command and numpy vectorization can be used. You can google around for to   #
# find this, but in general almost any basic math opertation is supported. Here are a few examples of setting    #
# compute strings:                                                                                               #
#    "compute: d['x']    + d['y']"        the code would compute a new array = [5,7,9].                          #
#    "compute: d['x']**2 + d['y']"        the code would compute a new array = [5,9,15].                         #
#    "compute: d['x']**2 / d['y']**2"     the code would compute a new array = [0.0625, 0.16, 0.25].             #
#    "compute: d['x']**2 * d['y']**0.5"   the code would compute a new array = [2.0, 8.94, 22.05].               #
#                                                                                                                #
# You may then set the 'xdata' and/or 'ydata' in the mode dictionaries as these "compute" strings, which will    #
# tell the code to perform those math operations, using the data set by the columns of the logfile.              #
#                                                                                                                #
# This can be useful to compute a new array of data, where the most obvious case is the need to compute the von  #
# mises stress criteria from a shear simulation. As in some research groups shear simulations are performed and  #
# the equivalent tensile yield strength is computed from the von mises stress criteria. Below are some notes     #
# on the von mises stress criteria equation and how the example shear stress log file is setup.                  #
#    svm = sigma_von_mises (computed)                                                                            #
#    s11 = sigma_11 (column in logfile: f_sxx_ave  PLEASE NOTE THIS MAY NOT BE THE SAME IN YOUR LOGFILES)        #
#    s22 = sigma_22 (column in logfile: f_syy_ave  PLEASE NOTE THIS MAY NOT BE THE SAME IN YOUR LOGFILES)        #
#    s33 = sigma_33 (column in logfile: f_szz_ave  PLEASE NOTE THIS MAY NOT BE THE SAME IN YOUR LOGFILES)        #
#    s12 = sigma_12 (column in logfile: f_sxy_ave  PLEASE NOTE THIS MAY NOT BE THE SAME IN YOUR LOGFILES)        #
#    s23 = sigma_23 (column in logfile: f_syz_ave  PLEASE NOTE THIS MAY NOT BE THE SAME IN YOUR LOGFILES)        #
#    s31 = sigma_31 (column in logfile: f_sxz_ave  PLEASE NOTE THIS MAY NOT BE THE SAME IN YOUR LOGFILES)        #
# svm = ( ((s11-s22)**2 + (s22-s33)**2 + (s33-s11)**2 + 6*(s12**2 + s23**2 + s31**2))/2 )**0.5                   #
#----------------------------------------------------------------------------------------------------------------#
# Method 1: Use a "brute force" typing of string (this serves as an example)
svm = "( ((d['f_sxx_ave']-d['f_syy_ave'])**2 + (d['f_syy_ave']-d['f_szz_ave'])**2 + (d['f_szz_ave']-d['f_sxx_ave'])**2 + 6*(d['f_sxy_ave']**2 + d['f_syz_ave']**2 + d['f_sxz_ave']**2))/2 )**0.5"

# Method 2: Using formated strings, this method is prefered as it lets others change column names a bit easier
s11 = "d['{}']".format('f_sxx_ave')
s22 = "d['{}']".format('f_syy_ave')
s33 = "d['{}']".format('f_szz_ave')
s12 = "d['{}']".format('f_sxy_ave')
s23 = "d['{}']".format('f_syz_ave')
s31 = "d['{}']".format('f_sxz_ave')
svm = "( (({s11}-{s22})**2 + ({s22}-{s33})**2 + ({s33}-{s11})**2 + 6*({s23}**2 + {s31}**2 + {s12}**2))/2 )**0.5".format(s11=s11, s22=s22, s33=s33, s12=s12, s23=s23, s31=s31)

# analysis list
analysis = [['moving average', 0, 0.1, 'window=50', 'moving average'],
            ['linear regression', 0, 0.02, 'shift=True', 'Modulus'],
            ['cursor', '', '', 'x=0; y=0;', 'Yield strength']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/property=shear_modulus_xz_strain_rate=2e8.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1',
        'xdata': 'v_etruexz',
        'ydata': f'compute: {svm}', # use variable from above since this is a long equation
        'xlabel': 'True Strain',
        'ylabel': 'Von Mises Stress (MPa)',
        'xscale': '',
        'yscale': '',
        'analysis': analysis
        }

