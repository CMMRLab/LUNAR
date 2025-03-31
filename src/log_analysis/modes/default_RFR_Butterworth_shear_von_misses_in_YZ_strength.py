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
# users may define a compute by setting a string as '${A} + ${B}/2 + 100', where A and B are columns in the log  #
# file. This will perform vectorized math operations on those columns of data from the logfile.                  #
#                                                                                                                #
# As an example say a LAMMPS logfile had columns x and y, with x-values = [1,2,3] and y-values = [4,5,6].        #
#    "compute: ${x}    + ${y}"        the code would compute a new array = [5,7,9].                              #
#    "compute: ${x}**2 + ${y}"        the code would compute a new array = [5,9,15].                             #
#    "compute: ${x}**2 / ${y}**2"     the code would compute a new array = [0.0625, 0.16, 0.25].                 #
#    "compute: ${x}**2 * ${y}**0.5"   the code would compute a new array = [2.0, 8.94, 22.05].                   #
#                                                                                                                #
# You may then set the 'xcompute' and/or 'ycompute' in the mode dictionaries as these "compute" strings, which   #
# will tell the code to perform those math operations, using the data set by the columns of the logfile.         #
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
svm = "( (( ${f_sxx_ave}-${f_syy_ave})**2 + (${f_syy_ave}-${f_szz_ave})**2 + (${f_szz_ave}-${f_sxx_ave})**2 + 6*(${f_sxy_ave}**2 + ${f_syz_ave}**2 + ${f_sxz_ave}**2))/2 )**0.5"

# Method 2: Using formated strings, this method is prefered as it lets others change column names a bit easier
s11 = "{}{}{}".format('${', 'f_sxx_ave', '}')
s22 = "{}{}{}".format('${', 'f_syy_ave', '}')
s33 = "{}{}{}".format('${', 'f_szz_ave', '}')
s12 = "{}{}{}".format('${', 'f_sxy_ave', '}')
s23 = "{}{}{}".format('${', 'f_syz_ave', '}')
s31 = "{}{}{}".format('${', 'f_sxz_ave', '}')
svm = "( (({s11}-{s22})**2 + ({s22}-{s33})**2 + ({s33}-{s11})**2 + 6*({s23}**2 + {s31}**2 + {s12}**2))/2 )**0.5".format(s11=s11, s22=s22, s33=s33, s12=s12, s23=s23, s31=s31)

# analysis list
analysis = [['LAMMPS data (apply Butterworth filter)', '', '', 'qm=msr; psd=False; csv=False; savefig=all; order=2; wn=op', 'LAMMPS Butterworth Filter'],
            ['Regression Fringe Response Modulus', '', '', 'shift=ymin; minxhi=0.01; maxxhi=0.0; xlo=rfs; yp=1; csv=False; savefig=all', 'Modulus']]

# loadable mode
mode = {'logfile': 'EXAMPLES/log_analysis/shear_3_EPON_862_pxld_88.2_replicate_1_FF_PCFF-class2xe.log.lammps',
        'keywords': ['Step', 'Temp'],
        'sections': '1',
        'xdata': 'v_etrueyz',
        'ydata': '',
        'xlabel': 'True Strain',
        'ylabel': 'Von Mises Stress (MPa)',
        'xcompute': '',
        'ycompute': f'{svm}', # use variable from above since this is a long equation,
        'analysis': analysis,
        'nevery': '1',
        'parent_directory': 'logfile'
        }

