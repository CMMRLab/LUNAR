# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
February 14th, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.log_analysis.signal_processing as signal_processing
import src.log_analysis.misc_funcs as misc_funcs
import matplotlib.pyplot as plt
import numpy as np
import math


##################################
# Function to find fringe slopes #
##################################
def compute_fringe_slope(strain, stress, min_strain=None, max_strain=None, direction='forward', stats=False):
    # Set direction
    if direction == 'forward':
        strain = strain.copy()
        stress = stress.copy()
    elif direction == 'reverse':
        strain = strain.copy()
        stress = stress.copy()
        strain.reverse()
        stress.reverse()
    else:
        raise Exception(f'ERROR direction={direction} is not supported. Supported directions are "forward" or "reverse"')
    
    # Set defaults if min_strain or max_strain are None
    if min_strain is None: min_strain = min(strain)
    if max_strain is None: max_strain = max(strain)
    
    # Generate outputs dictionary (only really need 'slopes' and 'fringe', but other statistics
    # can be computed at 25% slower run time, so keep that in mind when using stats=True).
    outputs = {'sum-residuals-squared':[], # cummulative sum of residuals sqaured
               'mean-squared-residual':[], # cummulative sum of residuals sqaured
               'intercepts-variance':[], # cummulative variance of intercepts
               'slopes-variance':[], # cummulative variance of slopes
               'intercepts':[], # cummulative Y-intercepts
               'r-squared-adj': [], # adjusted r squared values
               'r-squared':[], # cummulative r squared values
               'indexes': [], # cummulative indexes of stress-strain lists
               'slopes':[], # cummulative slopes
               'fringe':[]} # strain value
    
    # Start the walked linear regression method
    sum_xi = 0; sum_yi = 0; sum_xi_2 = 0; sum_yi_2 = 0; sum_xi_yi = 0; n = 0;
    for i, (x, y) in enumerate(zip(strain, stress)):
        # Compute cummulative linear regression parameters
        sum_xi += x; sum_yi += y; n += 1
        sum_xi_2 += x*x; sum_yi_2 += y*y; 
        sum_xi_yi += x*y
        
        # Need at least 2 points to perform linear regression
        if n <= 3: continue
        
        # Only compute outputs if x is in the desired range
        if min_strain <= x <= max_strain:
            SSxy = sum_xi_yi - (sum_xi*sum_yi/n)
            SSxx = sum_xi_2 - (sum_xi*sum_xi/n)
            b1 = SSxy/SSxx
            
            outputs['slopes'].append(b1)
            outputs['fringe'].append(x)
            outputs['indexes'].append(i) 
            
            # Additionally statistics that can be computed
            # (creates a 25% slower calculation)
            if stats:
                SStotal = sum_yi_2 - (sum_yi*sum_yi/n)
                MStotal = SStotal/(n - 1)
                x_bar = sum_xi/n
                y_bar = sum_yi/n
                b0 = y_bar - (b1*x_bar)
                SSres = SStotal - (b1*b1*SSxx)
                SSreg = b1*b1*SSxx
                MSres = SSres/(n - 2)
                r2 = SSreg/SStotal
                r2_adj = 1.0 - (MSres/MStotal)
                b0_variance = MSres*( (1/n) + (x_bar*x_bar/SSxx) )
                b1_variance = MSres/SSxx
                
                outputs['mean-squared-residual'].append(MSres)
                outputs['sum-residuals-squared'].append(SSres)
                outputs['intercepts-variance'].append(b0_variance)
                outputs['slopes-variance'].append(b1_variance)
                outputs['intercepts'].append(b0)
                outputs['r-squared'].append(r2)
                outputs['r-squared-adj'].append(r2_adj)
    return outputs


####################################################
# Function to maximize stress within a given bound #
####################################################
def maximize_stress_within_strain_range(strain, stress, xhi, xhi_yield):
    x_yield = None; y_yield = None
    strain_yield, stress_yield = misc_funcs.reduce_data(strain, stress, xhi, xhi_yield)
    if xhi_yield > xhi and strain_yield and stress_yield:
        try:
            y_yield = max(stress_yield)
            x_yield = strain_yield[stress_yield.index(y_yield)]
        except: pass
    return strain_yield, stress_yield, x_yield, y_yield


############################################################################################################
# Function to compute maximum strain for linear region based on maximum slope to machine machine_precision #
############################################################################################################
def maximized_slope(fringe, slopes, machine_precision=1e-8):
    # Find xhi and yhi, dealing with the possibility of having multiple solutions for maximized yhi. If multiply
    # solutions exists, assumed it is a machine precision issue and maximize xhi and yhi to obtain the solution.
    if slopes:
        maximized_slope = max(slopes)
        lo = maximized_slope - machine_precision
        hi = maximized_slope + machine_precision
        indexes = [i for i, y in enumerate(slopes) if lo <= y <= hi]
        maximized = max(indexes)
        xhi = fringe[maximized]
        yhi = slopes[maximized]
    else: xhi = 0; yhi = 0
    return xhi, yhi


#######################################################
# Function to construct perform inverse FFT filtering #
#######################################################
def construct_sine_wave(x, y, Npeaks=2):
    #-----------------------------------#
    # Compute the one-sided FFT and PSD #
    #-----------------------------------#
    # Define sampling rate and number of data points
    N = x.shape[0] # number of data points
    fs = (N-1)/(np.max(x) - np.min(x)) # sampling rate
    d = 1/fs # sampling space

    # Perform one sided FFT
    X = np.fft.rfft(y, axis=0, norm='backward')
    f = np.fft.rfftfreq(N, d=d)
    
    # One sided power spectrum density
    psd = np.real( (X*np.conjugate(X)) )/N
    
    #-----------------------------------------------------------------------------------#
    # Set a scaling factor of 0 or 1 to cancel out (0) or leave (1) certain frequencies #
    #-----------------------------------------------------------------------------------#
    # Find N-largest peaks in PSD to keep and ensure index zero (DC-offset) is in
    # the indices to ensure filtered data oscillates about the correct mean value.
    # The DC-offset maybe very small if data oscillates about a Y-value of zero, 
    # however in such a case, no harm is done in keeping the DC-offset.
    indices = np.argpartition(psd[1:], -Npeaks)[-Npeaks:] + 1
    if 0 not in indices:
        indices = np.append(indices, [0])
    scaling_factors = np.zeros_like(psd)
    scaling_factors[indices] = 1
    psd_peaks = psd[indices]
    psd_freqs = f[indices]
    
    #-------------------------------------------------------------------------------------------#
    # Cancel or leave frequencies in X and then inverse the cleaned fft to get the filterd data #
    #-------------------------------------------------------------------------------------------#
    X_clean = scaling_factors*X
    y_filter = np.fft.irfft(X_clean)
    
    # Since we are performing a one sided FFT, the Nyquist freq may or may not be inlcuded
    # depending on even or odd number of data points, so append a value if Nyquist freq is
    # missing so that y_filter has the same shape as the X-data.
    if y_filter.shape != x.shape:
        y_filter = np.append(y_filter, y_filter[-1])
            
    #-----------------------------------------------------------------------#
    # Compute amplitude and phase of signal based on the dominant frequency #
    #-----------------------------------------------------------------------#
    dominant_frequency_index = np.argmax(np.abs(X[1:])) + 1
    frequency = f[dominant_frequency_index]
    phase = np.angle(X[dominant_frequency_index], deg=False)
    amplitude = (np.max(y_filter) - np.min(y_filter))/2
    return f, psd, psd_freqs, psd_peaks, y_filter, phase, amplitude, frequency



def check_yp(x, y, xlmp, ylmp, yp_derivative):
    if len(x) == len(xlmp):
        residuals = [lmp - clean for lmp, clean in zip(ylmp, y)]
        std = misc_funcs.compute_standard_deviation(residuals)
        mean = misc_funcs.compute_mean(residuals)
        
        f, psd, psd_freqs, psd_peaks, y_filter, phase, amplitude, frequency = construct_sine_wave(np.array(x), np.array(residuals), Npeaks=3)

        
        fig1, ax1 = plt.subplots(1, 1, figsize=(6, 4))
        ax1.plot(x, residuals, 'o', ms=2, color='tab:blue', label='Residuals')
        ax1.plot(x, y_filter, '-', lw=2, color='tab:cyan', label='Standing stress wave')
        
        ax1.axhline(y=mean, lw=3, color='tab:orange', linestyle='--', label='Mean')
        ax1.axhline(y=mean+1*std, lw=3, color='bisque', linestyle='--', label='Mean + 1*sigma')
        ax1.axhline(y=mean+2*std, lw=3, color='tan', linestyle='--', label='Mean + 2*sigma')
        ax1.axhline(y=mean+3*std, lw=3, color='tab:brown', linestyle='--', label='Mean + 3*sigma')

        ax1.axvline(x=yp_derivative[0], lw=3, color='tab:green', linestyle='-', label='Yield Point (X)')        
        ax1.axhline(y=yp_derivative[-1], lw=3, color='tab:green', linestyle='-', label='Yield Point (Y)')
        ax1.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
        
    else:
        print('ERROR could not predict standing wave based on residuals, because cleaned and raw data series are different lengths')
    return


##################################################################################################################
# Function implementing Josh's method for finding linear region in stress vs strain. Inputs explained:           #
#    strain = Python list or tuple of strain data                                                                #
#    stress = Python list or tuple of stress data                                                                #
#    minxi = float or int value to set the minimum xhi value to start computing fringe slopes                    #
#    maxxi = float or int value to set the maximum xhi value to finished computing fringe slopes                 #
#    yp = int value positive or zero (i.e. 1 or 2 or 0 ...). If 0 the calculation is not run. Can also be a      #
#         string if 'min-2d' or 'min-v' or 'max-d'                                                               #
#    offset = float or int value                                                                                 #
#    xlo_method = string or float to control how xlo is determined. The following options are supported:         #
#        'iv'   which uses the vertex of the integrated stress-strain data                                       #
#        'rfs'  which uses the fringe slope method in reverse to maximize the slope                              #
#        float  which uses the float value set by the user, if the float value is less then xhi                  #
#    t1 and t2 = python lists of sublists of length 0 or 1 or 2 (i.e. t1 = [list(raw_data), list(clean_data)] or #
#                t1 = [list(raw_data)]). If the list length is zero for both t1 and t2, no poisson's ratio       #
#                calculations will be performed.                                                                 #
#    write_data = Boolean (True to write data to csv and False to not write data to csv)                         #
#    savefig = String with numbers ranging from 1-3. The numbers can have be seperated or not. Each number       #
#              corresponds to a plot to be saved. Plots are:                                                     #
#                '1' RFR with 4 plots of the RFR calcutions (xhi, xlo, derivaitves, stress-strain)               #
#                '2' Will be generated and saved if t1 or t2 are provided                                        #
#                '3' Will be generated and saved if ffs_stats is True                                            # 
#                '0' Will signify to not save any plots or an empty string                                       #
#                'all' can be passed to plot all three plots                                                     #
#              Examples: savefig = '1' or '1 2' or '1,2' or '1-2' or '1,2,3' or '1 2 3' or 'all' ....            #
#    figname = String to set figure and csv file basename, does not contain extension (i.e. no '.jpeg' or '.csv' #
#              extension, but just the basename name - for example 'Kemppainen-Muzzy-Modulus')                   #
#    dpi = Int to set the dots per inch of the saved figure                                                      #
##################################################################################################################
def compute(strain, stress, minxhi, maxxhi, xlo_method, yp, up, offset, t1, t2, stress_units, write_data, derivative_span, derivative_degree, grid, t_12_avg, savefig, figname, dpi):
    #--------------------------------------------------------------------------------------#
    # Internal plotting options (True to plot, False not to). The following are available: #
    #    ffs_stats, which plots the statistic evolution of the 2nd forward fringe slope.   #
    #    pub_plots, which plots the "regionized" example plots shown in the publication.   #
    #--------------------------------------------------------------------------------------#
    ffs_stats = False
    pub_plots = False
    
    #---------------------------------------------------------------------#
    # Data structs to continually add to through out the analysis process #
    #---------------------------------------------------------------------#
    poisson_xhis = [minxhi] # list to append xhi values for poisson's ratio calculation
    csv_data = {} # { 'name' : [[xdata], [ydata]], ... }
    
    #------------------------------------------------------------------------------------------------------------#
    # This entire function assumes strain increases with time (i.e. tensile or shear tests). Therefore if strain #
    # decreases with time for compression or unloading tests, we need to "reverse" the stress and strain lists.  #
    #------------------------------------------------------------------------------------------------------------#
    if strain[-1] < strain[0]: strain.reverse(); stress.reverse()
    
    #----
    # Update maxxhi is users asked for 'r2' method
    # if maxxhi == 'r2' or 1 == 1:
    #     def maximum_peak(rfr):
    #         # Find peaks
    #         from scipy.signal import find_peaks
    #         xdata = np.array(rfr['fringe']); ydata = np.array(rfr['r-squared']);
    #         peaks, properties = find_peaks(ydata, prominence=None)
    #         xpeaks = xdata[peaks]; ypeaks = ydata[peaks]
            
    #         # Find maximum peak index
    #         maxdex = np.argmax(ypeaks)
    #         max_x = xpeaks[maxdex]
    #         max_y = ypeaks[maxdex]
            
    #         # Related max-peak to rfr-data indexes
    #         rfr_index = np.argmin(np.abs(xdata - max_x))
            
    #         rfr['xpeaks'] = list(xpeaks)
    #         rfr['ypeaks'] = list(ypeaks)
    #         return rfr, rfr_index
        
        
    #     ffs_outputs_r2 = compute_fringe_slope(strain, stress, min_strain=minxhi, max_strain=None, direction='forward', stats=True)
    #     ffs_outputs_r2, maxdex = maximum_peak(ffs_outputs_r2)
        
    #     fig0, ax = plt.subplots(figsize=(6, 4))
        
    #     # Find maximum r2 value
    #     maxxhi = ffs_outputs_r2['fringe'][maxdex]
    #     max_r2 = ffs_outputs_r2['r-squared'][maxdex]
        
        
    #     ax.plot(ffs_outputs_r2['fringe'], ffs_outputs_r2['r-squared'], '-', lw=3, color='tab:blue', label='$r^2$')
    #     ax.plot(ffs_outputs_r2['xpeaks'], ffs_outputs_r2['ypeaks'], '*', ms=8, color='tab:cyan', label='Peaks')
    #     ax.plot(maxxhi, max_r2, 'o', ms=10, color='tab:cyan', label='maxxhi')
    #     ax.set_xlabel('Forward fringe strain')
    #     ax.set_ylabel('Forward fringe $r^2$')
    #     ax.legend()
        
         
    #--------------------------------------------------------------------------------------------#
    # First time: Find xhi by walking from near zero strain to max strain and maximizing modulus #
    #--------------------------------------------------------------------------------------------#   
    # The stress-strain curve may have a "dip" before the linear region (purely an MD thing), where
    # starting from the minimum stress-value elimanates the negative slope initially. As opposed
    # to starting from the minimum strain-value. We just need to be careful the minimum stress is
    # not caused by a fracture event. We will use the "minimum before the maximum" stress to avoid
    # issues in fracture events.
    max_index = stress.index(max(stress))
    try: min_stress = min(stress[:max_index])
    except: min_stress = min(stress)
    starting_strain = strain[stress.index(min_stress)]

    # Set min_strain value based on minxhi
    min_strain_ffs = min(strain); max_strain_ffs = max(strain)
    if minxhi > 0: min_strain_ffs = minxhi 
    if maxxhi > 0: max_strain_ffs = maxxhi
    min_strain_ffs += starting_strain
    
    # Compute forward fringe slopes (ffs) moving 1-data point at a time (starting at the "minimum before the maximum stress")
    ffs_outputs_1st = compute_fringe_slope(strain, stress, min_strain=min_strain_ffs, max_strain=max_strain_ffs, direction='forward', stats=False)
    xhi, yhi = maximized_slope(ffs_outputs_1st['fringe'], ffs_outputs_1st['slopes'], machine_precision=1e-8)
    if maxxhi > 0: # Compute for nicer looking plots
        ffs_outputs_1st = compute_fringe_slope(strain, stress, min_strain=min_strain_ffs, max_strain=max(strain), direction='forward', stats=False) 
    csv_data['forward-fringe-slope-1st-data'] = [ffs_outputs_1st['fringe'], ffs_outputs_1st['slopes']]
    csv_data['forward-fringe-slope-1st-maximum'] = [[xhi], [yhi]]

    #---------------------------------------------------------------------------------------------#
    # Find xlo Method 1: Tristan's method of using the vertex of the strain-energy polynomial fit #
    #---------------------------------------------------------------------------------------------#
    if xlo_method == 'iv':
        # Compute integral to determine strain energy, then compute polynomial fit, and find vertex
        ix, iy = misc_funcs.compute_integral(strain, stress)
        xin, yin, xfit, yfit, param = misc_funcs.np_curve_fit(ix, iy, degree=2, domain=[min(strain), xhi])
        vertex = -param[1]/2/param[2]
        
        # Find xlo and ylo    
        if vertex < 0 or vertex > xhi:
            xlo = min(strain)
        else: xlo = vertex
        lo = misc_funcs.closest(ix, xlo)
        ylo = iy[ix.index(lo)]
        
        # log data for csv
        csv_data['integrated-data'] = [xfit, yfit]
        csv_data['integrated-vertex'] = [[xlo], [ylo]]
        
    #----------------------------------------------------------------------------------#
    # Find xlo Method 2: Josh's method using reversed fringe slope maximization scheme #
    #----------------------------------------------------------------------------------#
    elif xlo_method == 'rfs':
        min_strain_rfs = starting_strain
        max_strain_rfs = xhi - minxhi # use minxhi to set the "span" of the smallest linear region acceptable
        
        # Ensure max_strain_fs is greater then or equal to min_strain_rfs
        if max_strain_rfs <= min_strain_rfs:
            min_index = strain.index(min(strain))
            try: max_strain_rfs = strain[min_index+1]
            except: max_strain_rfs = strain[min_index]
        rstrain, rstress = misc_funcs.reduce_data(strain, stress, min_strain_rfs, xhi)
        
        # Compute reverse fringe slopes (rfs) moving 1-data point at a time for the 1st time
        rfs_outputs_1st = compute_fringe_slope(rstrain, rstress, min_strain=min_strain_rfs, max_strain=max_strain_rfs, direction='reverse', stats=False)
        xlo, ylo = maximized_slope(rfs_outputs_1st['fringe'], rfs_outputs_1st['slopes'], machine_precision=1e-8)
        csv_data['reverse-fringe-slope-1st-data'] = [rfs_outputs_1st['fringe'], rfs_outputs_1st['slopes']]
        csv_data['reverse-fringe-slope-1st-maximum'] = [[xlo], [ylo]]
        
    #----------------------------------------#
    # Find xlo Method 3: Let user decide xlo #
    #----------------------------------------#
    elif isinstance(xlo_method, float) or isinstance(xlo_method, int):
        if float(xlo_method) < xhi:
            xlo = float(xlo_method)
        else: 
            xlo = min(strain)
            print(f'  WARNING user specified xlo "{xlo_method}" is greater then computed xhi "{xhi}". Defaulting to xlo to smallest strain value.')
        
    #---------------------------------------------------------#
    # Raise exception if asking for an unsupported xlo_method #
    #---------------------------------------------------------#
    else: raise Exception(f'ERROR unsupported xlo method "{xlo_method}". Currently supported methods as "iv" or "rfs" or <float-valued>')

    #---------------------------------------------------------------------------------------------#
    # Second time: Find xhi by walking from near zero strain to max strain and maximizing modulus #
    #---------------------------------------------------------------------------------------------#     
    # Re-compute forward fringe slopes (ffs) moving 1-data point at a time for a 2nd time
    min_strain_ffs += xlo
    #max_strain_ffs += xlo
    rstrain, rstress = misc_funcs.reduce_data(strain, stress, xlo, max(strain))
    ffs_outputs_2nd = compute_fringe_slope(rstrain, rstress, min_strain=min_strain_ffs, max_strain=max_strain_ffs, direction='forward', stats=ffs_stats)
    if ffs_outputs_2nd['fringe']:
        xhi, yhi = maximized_slope(ffs_outputs_2nd['fringe'], ffs_outputs_2nd['slopes'], machine_precision=1e-8)
    if maxxhi > 0: # Compute for nicer looking plots
        ffs_outputs_2nd = compute_fringe_slope(rstrain, rstress, min_strain=min_strain_ffs, max_strain=max(strain), direction='forward', stats=ffs_stats)
    poisson_xhis.append(xhi)
    csv_data['forward-fringe-slope-2nd-data'] = [ffs_outputs_2nd['fringe'], ffs_outputs_2nd['slopes']]
    csv_data['forward-fringe-slope-2nd-maximum'] = [[xhi], [yhi]]
            
    #-----------------------------------------------------#
    # Perform final linear regression between xlo and xhi #
    #-----------------------------------------------------#
    rstrain, rstress = misc_funcs.reduce_data(strain, stress, xlo, xhi)
    regression = misc_funcs.linear_regression(rstrain, rstress)
    xreg = regression.xreg; yreg = regression.yreg
    b1 = regression.b1; b0 = regression.b0
    csv_data['stress-strain-data-used-4-regression'] = [rstrain, rstress]
    csv_data['linear-regression-of-stress-strain-data'] = [xreg, yreg]


    #--------------------------------------------------------------------------------------------------------#
    # Find yield point by differentiating fringe slope/strain and finding the valleys of the 2nd derivatives #
    #--------------------------------------------------------------------------------------------------------#  
    # Find where max-stress is to limit some of the strain values for yield point selection (avoids issues
    # in case of fracture events).
    max_stress_index = stress.index(max(stress))
    max_strain = strain[max_stress_index]
    max_stress = stress[max_stress_index]
    if max_strain < xhi: max_strain = max(strain) # Ensure max_strain is larger then xhi
    
    # Find xhi_yield point based on inputs from user. The following options are available:
    #   yp = 0, means do not attempt to find xhi_yield
    #   yp = 1 or  2 or  3 or ..., Means take first or second or third or ... valley as the xhi_yield
    #   yp = 'min-d' means use the minimum of the 2nd derivative bound between the maximum 2nd derivative after the linear region
    #   yp = 'min-v' means use the minimum valley of the 2nd derivatives
    xhi_yield = None; yhi_yield = None; yield_point_derivative = []; strain_yield = []; stress_yield = []; 
    xhi_ultimate = None; yhi_ultimate = None; ultimate_point = []
    if isinstance(yp, int) and yp != 0 or yp in ['min-2d', 'min-v', 'max-d', 'min-r2d2'] or isinstance(up, int) and up != 0:        
        
        # Re-compute forward fringe slopes (ffs) moving 1-data point at a time for a 3rd time (in-case "maxxhi" option was used and then find the fringe and slopes after "xhi")
        rstrain, rstress = misc_funcs.reduce_data(strain, stress, xlo, max(strain))
        ffs_outputs_3rd = compute_fringe_slope(rstrain, rstress, min_strain=min_strain_ffs, max_strain=max(strain), direction='forward', stats=False)
        rfringe, rslopes = misc_funcs.reduce_data(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['slopes'], xhi, max_strain)
        
        # Update derivative_span if user supplies the 'hlr' keyword for half of the linear region span
        if derivative_span == 'hlr':
            tmp_strain, tmp_stress = misc_funcs.reduce_data(strain, stress, xlo, xhi)
            derivative_span = math.floor(len(tmp_strain)/2)
            print(f'  RFR: "ds" = hlr, updating span to half the linear region data points = {derivative_span}')
        
        # Compute some derivatives using user settings
        if derivative_span == 0 or derivative_degree < 2:
            dfringe, dslopes1, dslopes2 = misc_funcs.compute_derivative(rfringe, rslopes)
            dfringe_full, dslopes1_full, dslopes2_full = misc_funcs.compute_derivative(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['slopes'])
            print('  RFR: Using central difference derivatives')
        elif derivative_span > 0 and derivative_degree >= 2:
            dfringe, dslopes1, dslopes2 = misc_funcs.poly_derivative(rfringe, rslopes, span=derivative_span, deg=derivative_degree)
            dfringe_full, dslopes1_full, dslopes2_full = misc_funcs.poly_derivative(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['slopes'], span=derivative_span, deg=derivative_degree)
            print(f'  RFR: Using polynomial fitting derivatives w/ "ds" = {derivative_span} and "dd" = {derivative_degree}')
        else:
            print(f'  WARNING "ds" = {derivative_span} and "dd" = {derivative_degree}. To compute 2nd derivatives "dd" must be greater then 2, defaulting to central difference')
            
        # log the data
        csv_data['reduced-fringe-slope-data-4-differentiating'] = [rfringe, rslopes]
        csv_data['reduced-fringe-slope-1st-derivative'] = [dfringe, dslopes1]
        csv_data['reduced-fringe-slope-2nd-derivative'] = [dfringe, dslopes2]
        csv_data['reduced-fringe-slope-1st-derivative-FULL'] = [dfringe_full, dslopes1_full]
        csv_data['reduced-fringe-slope-2nd-derivative-FULL'] = [dfringe_full, dslopes2_full]
       
        # Find peaks and valleys (slowly decreasing prominence - std_dev/3 works for MD models, but for expeirmental
        # stress-strain curves it does not. The slow decreasing of promence is meant for expeirmental stress-strain
        # cases as a "Catch-all" to allow for the rest of the "yp" method to function properly)
        std_dev =  misc_funcs.compute_standard_deviation(dslopes2)
        incremented_search = True
        if incremented_search:
            for i in range(1, 11):
                std_div = i*3
                prominence = std_dev/std_div
                xpeaks, ypeaks, xvalleys, yvalleys, valley_depths = signal_processing.find_peaks_and_valleys(dfringe, dslopes2, prominence)
                if len(xvalleys) >= 4: break
        
        # Find peaks and valleys based on a fixed std_dev/3
        else:
            std_div = 3
            prominence = std_dev/std_div
            xpeaks, ypeaks, xvalleys, yvalleys, valley_depths = signal_processing.find_peaks_and_valleys(dfringe, dslopes2, prominence)
            
        # Perform one more check for valleys and if not set prominence as None to maximize the possibilty of finding valleys
        if len(xvalleys) < 1 and len(yvalleys) < 1:
            std_div = None
            xpeaks, ypeaks, xvalleys, yvalleys, valley_depths = signal_processing.find_peaks_and_valleys(dfringe, dslopes2, None)
        
        # Find valleys
        valleys = [(x, y, d[0], d[1], d[2]) for x, y, d in zip(xvalleys, yvalleys, valley_depths)] # [ (x-value, y-value, avg-depth, small-depth, large-depth) ]
        
        # Set yield point based on yp index, if there are any valleys, else default to minimum 2nd derivative
        if valleys and isinstance(yp, int) or isinstance(up, int):
            csv_data['reduced-fringe-slope-2nd-derivative-peaks'] = [xpeaks, ypeaks]
            csv_data['reduced-fringe-slope-2nd-derivative-valleys'] = [xvalleys, yvalleys]
            # valleys = sorted(valleys, key=lambda x: x[2]) # Sort by valley depth
            # valleys = sorted(valleys, key=lambda x: x[1]) # Sort by yvalley value
            if isinstance(yp, int):
                try:
                    yield_index = yp - 1
                    indexes = [i for i in range(len(valleys))]
                    if yield_index in indexes: 
                        xhi_yield, yhi_yield, avg_depth, small_depth, large_depth = valleys[yield_index]
                        csv_data['fringe-slope-2nd-derivative-valley-index={}'.format(yp)] = [[xhi_yield], [yhi_yield], 1]
                    else:
                        print(f'  WARNING could not find xhi yield point. Input yp = {yp}, attempting to find nearest valley.')
                        closes_index = misc_funcs.closest(indexes, yield_index)
                        xhi_yield, yhi_yield = valleys[closes_index]
                        csv_data['fringe-slope-2nd-derivative-valley-index={}'.format(closes_index+1)] = [[xhi_yield], [yhi_yield]]
                        print(f'    yp of {closes_index+1} found yield point')
                        
                    # Find closests strain value to min_2d_fringe
                    x_yield = misc_funcs.closest(strain, xhi_yield)
                    y_yield = stress[strain.index(x_yield)]
                    poisson_xhis.append(x_yield)
                    yield_point_derivative = [x_yield, y_yield]
                except: print('  WARNING linear region is near end of stress-strain. Cant compute yield strength.')
            
            # Testing ultimate predictions
            if isinstance(up, int):                
                # Find index in values
                indexes = [i for i in range(len(valleys))]
                if up > 0: ultimate_index = up - 1
                else: ultimate_index = len(indexes) - abs(up)
                if ultimate_index in indexes: 
                    xhi_ultimate, yhi_ultimate, avg_depth, small_depth, large_depth = valleys[ultimate_index]
                    csv_data['fringe-slope-2nd-derivative-valley-index={}'.format(up)] = [[xhi_ultimate], [yhi_ultimate], 1]
                else:
                    print(f'  WARNING could not find xhi yield point. Input yp = {yp}, attempting to find nearest valley.')
                    closes_index = misc_funcs.closest(indexes, ultimate_index)
                    xhi_ultimate, yhi_ultimate = valleys[closes_index]
                    csv_data['fringe-slope-2nd-derivative-valley-index={}'.format(closes_index+1)] = [[xhi_ultimate], [yhi_ultimate]]
                    print(f'    up of {closes_index+1} found ultimate point')
                    
                # Find closests strain value to min_2d_fringe
                x_ultimate = misc_funcs.closest(strain, xhi_ultimate)
                y_ultimate = stress[strain.index(x_ultimate)]
                ultimate_point = [x_ultimate, y_ultimate]
                    
        # Set yield point based on maximum valley depth
        elif valleys and yp == 'max-d':
            try:
                # Find deepest valley depth
                xhi_yield, yhi_yield, avg_depth, small_depth, large_depth = max(valleys, key=lambda x: abs(x[2]))
                csv_data['fringe-slope-2nd-derivative-max-valley-depth'] = [[xhi_yield], [yhi_yield]]
                
                # Find closests strain value to min_2d_fringe
                x_yield = misc_funcs.closest(strain, xhi_yield)
                y_yield = stress[strain.index(x_yield)]
                poisson_xhis.append(x_yield)
                yield_point_derivative = [x_yield, y_yield]
            except: print('  WARNING linear region is near end of stress-strain. Cant compute yield strength.')
            
        # Set yield point based on minimum valley
        elif valleys and yp == 'min-v':
            try:
                yhi_yield = min(yvalleys)
                xhi_yield = xvalleys[yvalleys.index(yhi_yield)]
                csv_data['fringe-slope-2nd-derivative-min-valley'] = [[xhi_yield], [yhi_yield]]
                
                # Find closests strain value to min_2d_fringe
                x_yield = misc_funcs.closest(strain, xhi_yield)
                y_yield = stress[strain.index(x_yield)]
                poisson_xhis.append(x_yield)
                yield_point_derivative = [x_yield, y_yield]
            except: print('  WARNING linear region is near end of stress-strain. Cant compute yield strength.')
        
        
        # Set yield point based on minimum of 2nd derivative in a given bound
        elif not valleys or yp == 'min-2d':
            # If yp is an int and no valleys are found this code will default to using the minimum of 2nd derivative in a given bound
            if isinstance(yp, int):
                print('  WARNING no valleys where found, defaulting to using minimum of 2nd derivative of fringe slope.')
            
            # Find the maximum 2nd derivative to bound the minimum afterwards
            try:
                max_2d_index = dslopes2.index(max(dslopes2))
                max_2d_fringe = dfringe[max_2d_index]
                
                # Find the minimum 2nd derivative bound between the max_2d_index and the end of the data set
                min_2d_slope = min(dslopes2)#[max_2d_index:])
                min_2d_fringe = dfringe[dslopes2.index(min_2d_slope)] 
                csv_data['fringe-slope-2nd-derivative-minimum'] = [[min_2d_fringe], [min_2d_slope]]
                
                # Find closests strain value to min_2d_fringe
                x_yield = misc_funcs.closest(strain, min_2d_fringe)
                y_yield = stress[strain.index(x_yield)]
                poisson_xhis.append(x_yield)
                yield_point_derivative = [x_yield, y_yield]
            except: 
                min_2d_fringe = None; min_2d_slope = None
                print('  WARNING linear region is near end of stress-strain. Cant compute yield strength.')
        
    # Set yield point based on minimum of 2nd derivative in a given bound
    if yp == 'min-r2d2':
        rstrain, rstress = misc_funcs.reduce_data(strain, stress, xlo, max(strain))
        ffs_outputs_3rd = compute_fringe_slope(rstrain, rstress, min_strain=xlo, max_strain=max(strain), direction='forward', stats=True)
        
        # Update derivative_span if user supplies the 'hlr' keyword for half of the linear region span
        if derivative_span == 'hlr':
            tmp_strain, tmp_stress = misc_funcs.reduce_data(strain, stress, xlo, xhi)
            derivative_span = math.floor(len(tmp_strain)/2)
            print(f'  RFR: "ds" = hlr, updating span to half the linear region data points = {derivative_span}')
        
        # Compute some derivatives using user settings
        if derivative_span == 0 or derivative_degree < 2:
            drf, dr2_d1, dr2_d2 = misc_funcs.compute_derivative(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['r-squared'])
            print('  RFR: Using central difference derivatives')
        elif derivative_span > 0 and derivative_degree >= 2:
            drf, dr2_d1, dr2_d2 = misc_funcs.poly_derivative(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['r-squared'], span=derivative_span, deg=derivative_degree)
            print(f'  RFR: Using polynomial fitting derivatives w/ "ds" = {derivative_span} and "dd" = {derivative_degree}')
        else:
            print(f'  WARNING "ds" = {derivative_span} and "dd" = {derivative_degree}. To compute 2nd derivatives "dd" must be greater then 2, defaulting to central difference')
        
        drf, dr2_d1, dr2_d2 = misc_funcs.compute_derivative(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['r-squared'])
        
        drf_bounded, dr2_d2_bounded = misc_funcs.reduce_data(drf, dr2_d2, xhi, max_strain)
        upper_bound = drf_bounded[dr2_d2_bounded.index(max(dr2_d2_bounded))]
        drf_bounded, dr2_d2_bounded = misc_funcs.reduce_data(drf, dr2_d2, xhi, upper_bound)
        
        min_dr2_d2 = min(dr2_d2_bounded)
        x_yield = drf_bounded[dr2_d2_bounded.index(min_dr2_d2)]
        y_yield = stress[strain.index(x_yield)]
        poisson_xhis.append(x_yield)
        yield_point_derivative = [x_yield, y_yield]
        
    # Set yield point based on r^2 crossing of mean +- 3*sigma or r^2 in linear region
    elif yp in ['mean(r2)-3s', 'max(r2)-3s']:     
        rstrain, rstress = misc_funcs.reduce_data(strain, stress, xlo, max(strain))
        ffs_outputs_3rd = compute_fringe_slope(rstrain, rstress, min_strain=xlo, max_strain=max(strain), direction='forward', stats=True)
        r2_linear_x, r2_linear_y = misc_funcs.reduce_data(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['r-squared'], xlo, xhi)
        r2_sigma = misc_funcs.compute_standard_deviation(r2_linear_y)
        if yp == 'mean(r2)-3s':
            r2_reference = misc_funcs.compute_mean(r2_linear_y)
        if yp == 'max(r2)-3s':
            r2_reference = max(r2_linear_y)
        r2_lo = r2_reference - 3*r2_sigma
        r2_hi = r2_reference + 3*r2_sigma
        
        # Find when value drops below the r2_lo threshold
        r2_nonlinear_x, r2_nonlinear_y = misc_funcs.reduce_data(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['r-squared'], xhi, max(strain))
        x_crossing, y_crossing = misc_funcs.value_crossing(r2_nonlinear_x, r2_nonlinear_y, yvalue=r2_lo, style='high')
        x_cross = r2_linear_x[-1]; y_cross = r2_linear_y[-1]
        for x_cross, y_cross in zip(x_crossing, y_crossing):
            if y_cross <= r2_lo: break
        x_yield = x_cross
        y_yield = stress[strain.index(x_yield)]
        poisson_xhis.append(x_yield)
        yield_point_derivative = [x_yield, y_yield]
        
    # Set yield point based on the maximum of the 3rd forward fringe pass
    elif yp == 'max-3ffs':
        rstrain, rstress = misc_funcs.reduce_data(strain, stress, xhi, max(strain))
        start = xhi #+ minxhi
        ffs_outputs_3rd = compute_fringe_slope(rstrain, rstress, min_strain=start, max_strain=max(strain), direction='forward', stats=False)
        max_3rd_y = max(ffs_outputs_3rd['slopes'])
        max_3rd_x =  ffs_outputs_3rd['fringe'][ffs_outputs_3rd['slopes'].index(max_3rd_y)]
        
        yp_index = strain.index(max_3rd_x)
        x_yield = strain[yp_index]
        y_yield = stress[yp_index]
        poisson_xhis.append(x_yield)
        yield_point_derivative = [x_yield, y_yield]

    # Set yield point based on the maximum of the 3rd and 4th forward fringe pass
    elif yp == 'min-4ffs':    
        rstrain, rstress = misc_funcs.reduce_data(strain, stress, xhi, max(strain))    
        start = xhi + minxhi
        ffs_outputs_3rd = compute_fringe_slope(rstrain, rstress, min_strain=start, max_strain=max(strain), direction='forward', stats=False)   
        ffs_outputs_4th = compute_fringe_slope(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['slopes'], min_strain=start, max_strain=max(strain), direction='forward', stats=False)
        
        ffs_outputs_5th = compute_fringe_slope(ffs_outputs_4th['fringe'], ffs_outputs_4th['slopes'], min_strain=start, max_strain=max(strain), direction='forward', stats=False)
        
        ffs_x_bounded, ffs_y_bounded = misc_funcs.reduce_data(ffs_outputs_4th['fringe'], ffs_outputs_4th['slopes'], xhi, max_strain)
        yhi_yield = min(ffs_y_bounded)
        xhi_yield = ffs_x_bounded[ffs_y_bounded.index(yhi_yield)]
        
        yp_index = strain.index(xhi_yield)
        x_yield = strain[yp_index]
        y_yield = stress[yp_index]
        poisson_xhis.append(x_yield)
        yield_point_derivative = [x_yield, y_yield]


    #-------------------------------------#
    # Find yield point by using N% offset #
    #-------------------------------------#
    offset_x = []; offset_y = []; yield_point_offset = [];
    if offset > 0:
        # Find offset line
        dx = abs(misc_funcs.compute_mean([strain[i+1]-strain[i] for i in range(len(strain)-1)]))
        dy = abs(misc_funcs.compute_mean([stress[i+1]-stress[i] for i in range(len(stress)-1)]))
        if dx < 0.0001: dx = 0.0001
        maxy = max(stress) + 2*dy; x = min(strain); y = min(stress);
        while y <= maxy:
            x += dx; y = b1*x + b0
            offset_x.append(x + offset); offset_y.append(y)
        csv_data['offset-line={}'.format(offset)] = [offset_x, offset_y]
            
        # Find two closests points to use as the intersection
        distances = {} # {(stress-strain-index, offset-line-index):distance}
        xskip = xhi + offset
        for n1, (x1, y1) in enumerate(zip(strain, stress)):
            if x1 < xskip: continue
            xclosest = misc_funcs.closest(offset_x, x1)   
            n2 = offset_x.index(xclosest)
            x2 = offset_x[n2]; y2 = offset_y[n2];
            deltax = x2 - x1; deltay = y2 - y1;
            distance = deltax*deltax + deltay*deltay # not using square root for speed-up (its not required anyways, since this is relative)
            distances[(n1, n2)] = distance
        try:
            intersection = min(distances, key=distances.get)
            xyield0 = strain[intersection[0]]
            xyield1 = offset_x[intersection[1]]
            xyield2 = (xyield0 + xyield1)/2
            xyield = misc_funcs.closest(strain, xyield2)
            yyield = stress[strain.index(xyield)]
            yield_point_offset = [xyield, yyield]
            poisson_xhis.append(xyield)
            csv_data['offset-line_yield_point ({})'.format(offset)] = [[xyield], [yyield]]
        except: print(f'  ERROR RFR-modulus could not find yield point using offset method. Likely due to lack of data at {offset} offset')
    
    #-------------------------------------------------#
    # Plot results for fringe-slope and stress-strain #
    #-------------------------------------------------#
    # Set xlimits
    delta = 0.01
    xlimits = (min(strain)-delta, max(strain)+delta)
    
    # Create new plot
    if yp != 0 or yield_point_offset:
        fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 9))
    else:
        fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 7))
    
    # Plot fringe modulus method
    # if yp != 0 or yp in ['min-2d', 'min-v', 'max-d', 'min-c']:
    #     if xhi_yield is not None:
    #         ax1.axvline(xhi_yield, color='tab:cyan', ls='--', lw=2, label='Yield point ({})'.format(xhi_yield))
    ax1.plot(ffs_outputs_1st['fringe'], ffs_outputs_1st['slopes'], '.', ms=4.0, color='tab:blue', label='Fringe slope 1st')
    ax1.plot(ffs_outputs_2nd['fringe'], ffs_outputs_2nd['slopes'], '.', ms=4.0, color='tab:green', label='Fringe slope 2nd')
    if yp == 'min-4ffs': 
        ax1.plot(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['slopes'], '.', ms=4.0, color='tab:cyan', label='Fringe slope 3rd')
    ax1.plot(xhi, yhi, 'o', ms=10.0, markeredgecolor='black', color='tab:green', label='Upper break point (x,y):\n ({:.4f}, {:.4f})'.format(xhi, yhi))
    ax1.set_xlabel('Forward fringe strain')
    ax1.set_ylabel('Forward fringe slope {}'.format(stress_units))
    ax1.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
    ax1.set_title('Forward fringe slope to find xhi', fontsize=10)
    ax1.set_xlim(xlimits)
    if grid == 'on':
        ax1.grid()
    
    # Plot integral method for finding xlo
    if xlo_method == 'iv':
        ax2.plot(ix, iy, '.', ms=4.0, color='tab:blue', label='Strain energy')
        ax2.plot(xin, yfit, '.', ms=4.0, color='tab:orange', label='Polynomial fit')
        ax2.plot(xlo, ylo, 'o', ms=10.0, markeredgecolor='black', color='tab:green', label='Lower break point (x,y):\n ({:.4f}, {:.4f})'.format(xlo, ylo))
        ax2.set_xlabel('Strain')
        ax2.set_ylabel('Strain Energy')
        ax2.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
        ax2.set_title('Strain energy to find xlo for linear regression', fontsize=10)
        ax2.set_xlim(xlimits)
        if grid == 'on':
            ax2.grid()
    if xlo_method == 'rfs':
        ax2.plot(rfs_outputs_1st['fringe'], rfs_outputs_1st['slopes'], '.', ms=4.0, color='tab:blue', label='Reversed fringe slope')
        ax2.plot(xlo, ylo, 'o', ms=10.0, markeredgecolor='black', color='tab:green', label='Lower break point (x,y):\n ({:.4f}, {:.4f})'.format(xlo, ylo))
        ax2.set_xlabel('Reversed fringe strain')
        ax2.set_ylabel('Reversed fringe slope {}'.format(stress_units))
        ax2.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
        ax2.set_title('Reversed fringe slope from 1st xhi to find xlo', fontsize=10)
        #ax2.set_xlim(xlimits)
        if grid == 'on':
            ax2.grid()
        
    # Plot yield strength (if using derivative option)
    if yp != 0 or yp in ['min-2d', 'min-v', 'max-d', 'min-r2d2', 'mean(r2)-3s', 'max(r2)-3s', 'max-3ffs', 'min-4ffs']:
        # Plot yield strength method
        if yp in ['min-2d', 'min-v', 'max-d'] or isinstance(yp, int):
            ax3.plot(dfringe_full, dslopes2_full, '-', color='lime', lw=3.0, label='$(d^2slope)/(dfringe^2)$ unbounded')
            ax3.plot(dfringe, dslopes2, '-', color='tab:blue', lw=3.0, label='$(d^2slope)/(dfringe^2)$ bounded')
            ax3.axvline(xlo, color='tab:red', ls='--', lw=2, label='Linear-region xlo')
            ax3.axvline(xhi, color='tab:red', ls='--', lw=2, label='Linear-region xhi')
            ax3.axvline(max_strain, color='tab:green', ls='--', lw=2, label='Max stress')
            if xvalleys and yvalleys and yp != 'min-2d':
                ax3.plot(xhi_yield, yhi_yield, 'o', ms=10.0, markeredgecolor='black', color='tab:green', label='Yield point (x,y):\n ({:.4f}, {:.4f})'.format(xhi_yield, yhi_yield))
                if std_div is not None:
                    ax3.plot(xvalleys, yvalleys, 'o', ms=6.0, markeredgecolor='black', color='tab:orange', label='Valleys (w/ prominence=STD-DEV/{})'.format(std_div))
                else:
                    ax3.plot(xvalleys, yvalleys, 'o', ms=6.0, markeredgecolor='black', color='tab:orange', label='Valleys (w/ prominence={})'.format(std_div))
                ax3.axvline(xhi_yield, color='tab:cyan', ls='--', lw=2, label='Yield point ({})'.format(xhi_yield))
            if yp == 'min-2d' or not xvalleys and not yvalleys and  min_2d_fringe is not None and min_2d_slope is not None:
                ax3.plot(min_2d_fringe, min_2d_slope, 'o', ms=10.0, markeredgecolor='black', color='tab:green', label='Yield point (x,y):\n ({:.4f}, {:.4f})'.format(min_2d_fringe, min_2d_slope))
                ax3.axvline(min_2d_fringe, color='tab:cyan', ls='--', lw=2, label='Yield point ({})'.format(min_2d_fringe))
            ax3.set_xlabel('Strain')
            ax3.set_ylabel('$d^2strain$ / $dslope^2$ {}'.format(stress_units))
            ax3.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
            ax3.set_xlim(xlimits)
            if grid == 'on':
                ax3.grid()
            
        elif yp == 'max-3ffs':
            ax3.plot(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['slopes'], '.', ms=4.0, color='tab:blue', label='Fringe slope 3rd')
            # ax3.axvline(xlo, color='tab:red', ls='--', lw=2, label='Linear-region xlo')
            # ax3.axvline(xhi, color='tab:red', ls='--', lw=2, label='Linear-region xhi')
            ax3.plot(max_3rd_x, max_3rd_y, 'o', ms=10.0, markeredgecolor='black', color='tab:green', label='Yeild point (x,y):\n ({:.4f}, {:.4f})'.format(max_3rd_x, max_3rd_y))
            
            ax3.set_xlabel('Forward fringe strain')
            ax3.set_ylabel('Forward fringe slope {}'.format(stress_units))
            ax3.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
            ax3.set_xlim(xlimits)
            if grid == 'on':
                ax3.grid()
            
        elif yp == 'min-4ffs':  
            ax3_color = 'tab:blue' 
            ax5_color = 'tab:green'
            
            ax3.plot(ffs_outputs_4th['fringe'], ffs_outputs_4th['slopes'], '.', ms=4.0, color=ax3_color, label='Fringe slope 1 of 3rd')
            ax3.plot(xhi_yield, yhi_yield, 'o', ms=10.0, markeredgecolor='black', color=ax3_color, label='Yield point (x,y):\n ({:.4f}, {:.4f})'.format(xhi_yield, yhi_yield))
            ax3.axvline(xlo, color='tab:red', ls='--', lw=2, label='Linear-region xlo')
            ax3.axvline(xhi, color='tab:red', ls='--', lw=2, label='Linear-region xhi')
            ax3.tick_params(axis='y', colors=ax3_color)
            ax3.spines['right'].set_color(ax3_color)
            ax3.set_xlabel('Forward fringe strain')
            ax3.set_ylabel('Forward fringe $r^2$', color=ax3_color)
            ax3.set_xlim(xlimits)
            
            ax5 = ax3.twinx()
            ax5.plot(ffs_outputs_5th['fringe'], ffs_outputs_5th['slopes'], '.', ms=4.0, color=ax5_color, label='Fringe slope 2 of 3rd')
            #ax5.plot(xhi_yield, yhi_yield, 'o', ms=10.0, markeredgecolor='black', color=ax5_color, label='Yield point (x,y):\n ({:.4f}, {:.4f})'.format(xhi_yield, yhi_yield))
            ax5.axvline(max_strain, color='lime', ls='--', lw=2, label='Max stress')
            ax5.tick_params(axis='y', colors=ax5_color)
            ax5.spines['right'].set_color(ax5_color)
            ax5.set_ylabel('$r2^2$ / $dslope^2$', color=ax5_color)
            
            ax3.legend(loc='lower right', bbox_to_anchor=(1.0, 0.2), fancybox=True, ncol=1, fontsize=8)
            ax5.legend(loc='lower right', bbox_to_anchor=(1.0, 0.0), fancybox=True, ncol=1, fontsize=8)
            if grid == 'on':
                ax3.grid()

        elif yp == 'min-r2d2':
            ax3_color = 'tab:blue' 
            ax5_color = 'tab:green'
            
            ax3.plot(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['r-squared'], '-', color=ax3_color, lw=3, label='$r^2$')
            ax3.axvline(xlo, color='tab:red', ls='--', lw=2, label='Linear-region xlo')
            ax3.axvline(xhi, color='tab:red', ls='--', lw=2, label='Linear-region xhi')

            ax3.tick_params(axis='y', colors=ax3_color)
            ax3.spines['right'].set_color(ax3_color)
            ax3.set_xlabel('Forward fringe strain')
            ax3.set_ylabel('Forward fringe $r^2$', color=ax3_color)
            ax3.set_xlim(xlimits)

            ax5 = ax3.twinx()
            #ax5.plot(drf, dr2_d2, '-', color='lime', lw=3, label='$(dr2^2)/(dfringe^2)$ unbounded')
            ax5.plot(drf_bounded, dr2_d2_bounded, '-', color=ax5_color, lw=3, label='$(dr2^2)/(dfringe^2)$ bounded')
            ax5.plot(x_yield, min_dr2_d2, 'o', ms=10.0, markeredgecolor='black', color=ax5_color, label='Yield point (x,y):\n ({:.4f}, {:.4f})'.format(x_yield, min_dr2_d2))
            ax5.axvline(upper_bound, color=ax5_color, ls='--', lw=2, label='Max $(dr2^2)/(dfringe^2)$')
            ax5.axvline(max_strain, color='lime', ls='--', lw=2, label='Max stress')
            #ax5.axvline(x_yield, color=ax5_color, ls='--', lw=2, label='Yield point')
            ax5.tick_params(axis='y', colors=ax5_color)
            ax5.spines['right'].set_color(ax5_color)
            ax5.set_ylabel('$r2^2$ / $dslope^2$', color=ax5_color)
            
            ax3.legend(loc='lower right', bbox_to_anchor=(1.00, 0.80), fancybox=True, ncol=1, fontsize=8)
            ax5.legend(loc='lower right', bbox_to_anchor=(0.42, 0.70), fancybox=True, ncol=1, fontsize=8)
            if grid == 'on':
                ax3.grid()
            
        elif yp in ['mean(r2)-3s', 'max(r2)-3s']:
            ax3.plot(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['r-squared'], color='tab:blue', lw=3.0, label='2nd forward fringe $r^2$')
            ax3.plot(r2_linear_x, r2_linear_y, color='lime', lw=4.0, label='2nd forward fringe linear region $r^2$')
            
            ax3.axvline(xlo, color='tab:red', ls='--', lw=1.5, label='Linear-region xlo')
            ax3.axvline(xhi, color='tab:red', ls='--', lw=1.5, label='Linear-region xhi')
            
            if yp == 'mean(r2)-3s':
                ax3.axhline(r2_lo, color='tab:orange', ls='--', lw=1.5, label='Linear $r^2$ ($\\mu$ - 3*$\\sigma$) w/ $\\mu$={:.4f}'.format(r2_reference))
                ax3.axhline(r2_hi, color='tab:orange', ls='--', lw=1.5, label='Linear $r^2$ ($\\mu$ + 3*$\\sigma$) w/ $\\sigma$={:.4f}'.format(r2_sigma))
            if yp == 'max(r2)-3s':
                ax3.axhline(r2_lo, color='tab:orange', ls='--', lw=1.5, label='Linear $r^2$ (max - 3*$\\sigma$) w/ max={:.4f}'.format(r2_reference))
                ax3.axhline(r2_hi, color='tab:orange', ls='--', lw=1.5, label='Linear $r^2$ (max + 3*$\\sigma$) w/ $\\sigma$={:.4f}'.format(r2_sigma))
            
            ax3.plot(x_cross, y_cross, 'o', ms=10.0, markeredgecolor='black', color='tab:green', label='Yield point (x,y):\n ({:.4f}, {:.4f})'.format(x_cross, y_cross))
            
            ax3.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=10)
            if grid == 'on':
                ax3.grid()

        # Plot stress-vs-strain, yield region, and yield strength
        ax4.plot(strain, stress, '-', lw=4, color='tab:blue', label='Stress')
        if strain_yield and stress_yield:
            ax4.plot(strain_yield, stress_yield, '.', mfc='white', ms=5, markeredgecolor='tab:orange', lw=0.01, label='Elastic region')
        ax4.plot(xreg, yreg, '-', lw=4.0, color='lime', label='Linear regression (slope={:.4f})'.format(b1))
        if yield_point_derivative:
            ax4.plot(yield_point_derivative[0], yield_point_derivative[1], 'o', ms=10.0, markeredgecolor='black', color='tab:green', label='Yield point from "yp" (x,y):\n ({:.4f}, {:.4f}) NOT SHIFTED YET'.format(yield_point_derivative[0], yield_point_derivative[1]))
        ax4.set_xlabel('Strain')
        ax4.set_ylabel('Stress {}'.format(stress_units))
        ax4.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
        ax4.set_title('Yield region and yield strength', fontsize=10)
        ax4.set_xlim(xlimits)
        if grid == 'on':
            ax4.grid()
        
    # Plot yield strength (if using offset option)
    if yield_point_offset:
        if xhi_yield is None and yhi_yield is None and yp == 0:
            ax4.plot(strain, stress, '.', mfc='white', ms=5, markeredgecolor='tab:blue', lw=0.01, label='Stress')
            ax4.plot(xreg, yreg, '-', lw=4.0, color='lime', label='Linear regression')
        ax4.plot(offset_x, offset_y, '.', ms=3.0, color='tab:cyan', label=f'Linear regression offset ({offset})')
        if yield_point_offset:
            ax4.plot(yield_point_offset[0], yield_point_offset[1], 'o', ms=10.0, markeredgecolor='black', color='teal', label='Yield point from "offset {}" (x,y):\n ({:.4f}, {:.4f}) NOT SHIFTED YET'.format(offset, yield_point_offset[0], yield_point_offset[1]))
        ax4.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
        ax4.set_xlabel('Strain')
        ax4.set_ylabel('Stress')
        ax4.set_xlim(xlimits)
        if grid == 'on':
            ax4.grid()
            
    # Plot yield strength (if using derivative option)
    if up != 0 and xvalleys and yvalleys:
        ax3.plot(xhi_ultimate, yhi_ultimate, 'o', ms=10.0, markeredgecolor='black', color='tab:purple', label='Ultimate point (x,y):\n ({:.4f}, {:.4f})'.format(xhi_ultimate, yhi_ultimate))
        ax3.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
        if ultimate_point:
            ax4.plot(ultimate_point[0], ultimate_point[1], 'o', ms=10.0, markeredgecolor='black', color='tab:purple', label='Ultimate point from "up" (x,y):\n ({:.4f}, {:.4f}) NOT SHIFTED YET'.format(*ultimate_point))
            ax4.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
            
    # Apply tight layout to avoid overlapping plots and save plot if user wants
    fig1.tight_layout()
    if '0' not in str(savefig) and '1' in str(savefig) or 'all' in str(savefig):
        fig1.savefig(figname+'_1.jpeg', dpi=dpi)
        
    
    #----------------------------------------------------------#
    # Compute poisson's ratio if t1 and t2 are not empty lists #
    #----------------------------------------------------------#
    nu1 = None; nu2 = None; nu12 = None
    if len(t1) >= 1 and len(t2) >= 1:
        # Find data sets to analyze
        if len(t1) == 1:
            t1_raw = []; t1_reg = t1[0]
        elif len(t1) >= 2:
            t1_raw = t1[0]
            t1_reg = t1[-1] # use the last logged entry
        else:
            raise Exception(f'  ERROR length of t1 can only be 1 or 2. Current length of t1 is {len(t1)}')
            
        if len(t2) == 1:
            t2_raw = []; t2_reg = t2[0]
        elif len(t2) >= 2:
            t2_raw = t2[0]
            t2_reg = t2[-1]  # use the last logged entry
        else:
            raise Exception(f'  ERROR length of t2 can only be 1 or 2. Current length of t2 is {len(t2)}')
            
        if len(t1) == 1 and len(t2) == 1:
            t12_raw = []
            t12_reg = [(i + j)/2 for i, j in zip(t1_reg, t2_reg)]
        elif len(t1) >= 2 and len(t2) >= 2:
            t12_raw = [(i + j)/2 for i, j in zip(t1_raw, t2_raw)]
            t12_reg = [(i + j)/2 for i, j in zip(t1_reg, t2_reg)]
        else:
            raise Exception(f'  ERROR length of t1 and t2 can only be 1 or 2. Current length of t1 is {len(t1)} and t2 is {len(t2)}')
        
        # Find linear regions
        poisson_xhi = max(poisson_xhis)
        axial1, trans1 = misc_funcs.reduce_data(strain, t1_reg, xlo, poisson_xhi)
        axial2, trans2 = misc_funcs.reduce_data(strain, t2_reg, xlo, poisson_xhi)
        axial12, trans12 = misc_funcs.reduce_data(strain, t12_reg, xlo, poisson_xhi)
    
        # Perform linear regression on linear regions
        nu1_lr = misc_funcs.linear_regression(axial1, trans1)
        nu2_lr = misc_funcs.linear_regression(axial2, trans2)
        nu12_lr = misc_funcs.linear_regression(axial12, trans12)
        nu1 = abs(nu1_lr.b1)
        nu2 = abs(nu2_lr.b1)
        nu12 = abs(nu12_lr.b1)
        
        # Plot results
        if t_12_avg:
            fig2, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 5))
        else:
            fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
        if len(t1) == 1:
            ax1.plot(strain, t1_reg, '.', mfc='white', ms=5, markeredgecolor='tab:blue', lw=0.01, label='LAMMPS "raw" data')
            label1 = '{}\ny = {:.6f}x + {:.6f}'.format('Linear regression on "raw" data', nu1_lr.b1, nu1_lr.b0)
        if len(t1) >= 2:
            ax1.plot(strain, t1_raw, '.', mfc='white', ms=5, markeredgecolor='tab:cyan', lw=0.01, label='LAMMPS "raw" data')
            ax1.plot(strain, t1_reg, '.', mfc='white', ms=5, markeredgecolor='tab:blue', lw=0.01, label='LAMMPS "cleaned" data')
            label1 = '{}\ny = {:.6f}x + {:.6f}'.format('Linear regression on "clean" data', nu1_lr.b1, nu1_lr.b0)
        ax1.plot(axial1, trans1, '.', mfc='white', ms=5, markeredgecolor='tab:orange', lw=0.01, label='Linear region')
        ax1.plot(nu1_lr.xreg, nu1_lr.yreg, '--', lw=2, color='tab:green', label=label1)
        ax1.set_xlabel('Axial strain')
        ax1.set_ylabel('Transverse strain (t1)')
        ax1.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
        ax1.set_title('Axial vs Transverse-1 Strain', fontsize=10)
        ax1.set_xlim(xlimits)
        if grid == 'on':
            ax1.grid()
        
        if len(t2) == 1:
            ax2.plot(strain, t2_reg, '.', mfc='white', ms=5, markeredgecolor='tab:blue', lw=0.01, label='LAMMPS "raw" data')
            label2 = '{}\ny = {:.6f}x + {:.6f}'.format('Linear regression on "raw" data', nu2_lr.b1, nu2_lr.b0)
        if len(t2) >= 2:
            ax2.plot(strain, t2_raw, '.', mfc='white', ms=5, markeredgecolor='tab:cyan', lw=0.01, label='LAMMPS "raw" data')
            ax2.plot(strain, t2_reg, '.', mfc='white', ms=5, markeredgecolor='tab:blue', lw=0.01, label='LAMMPS "clean" data')
            label2 = '{}\ny = {:.6f}x + {:.6f}'.format('Linear regression on "clean" data', nu2_lr.b1, nu2_lr.b0)
        ax2.plot(axial2, trans2, '.', mfc='white', ms=5, markeredgecolor='tab:orange', lw=0.01, label='Linear region')
        ax2.plot(nu2_lr.xreg, nu2_lr.yreg, '--', lw=2, color='tab:green', label=label2)
        ax2.set_xlabel('Axial strain')
        ax2.set_ylabel('Transverse strain (t2)')
        ax2.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
        ax2.set_title('Axial vs Transverse-2 Strain', fontsize=10)
        ax2.set_xlim(xlimits)
        
        if t_12_avg:
            if len(t1) == 1 and len(t2) == 1:
                ax3.plot(strain, t12_reg, '.', mfc='white', ms=5, markeredgecolor='tab:blue', lw=0.01, label='LAMMPS "raw" data')
                label3 = '{}\ny = {:.6f}x + {:.6f}'.format('Linear regression on "raw" data', nu12_lr.b1, nu12_lr.b0)
            if len(t1) >= 2 and len(t2) >= 2:
                ax3.plot(strain, t12_raw, '.', mfc='white', ms=5, markeredgecolor='tab:cyan', lw=0.01, label='LAMMPS "raw" data')
                ax3.plot(strain, t12_reg, '.', mfc='white', ms=5, markeredgecolor='tab:blue', lw=0.01, label='LAMMPS "clean" data')
                label3 = '{}\ny = {:.6f}x + {:.6f}'.format('Linear regression on "clean" data', nu12_lr.b1, nu12_lr.b0)
            ax3.plot(axial12, trans12, '.', mfc='white', ms=5, markeredgecolor='tab:orange', lw=0.01, label='Linear region')
            ax3.plot(nu12_lr.xreg, nu12_lr.yreg, '--', lw=2, color='tab:green', label=label3)
            ax3.set_xlabel('Axial strain')
            ax3.set_ylabel('Transverse strain (t1 + t2)/2')
            ax3.legend(loc='lower right', bbox_to_anchor=(1, 0), fancybox=True, ncol=1, fontsize=8)
            ax3.set_title('Axial vs Transverse-average Strain', fontsize=10)
            ax3.set_xlim(xlimits)
        
        
        
        if grid == 'on':
            ax2.grid()
        
        # Apply tight layout to avoid overlapping plots and save plot if user wants
        fig2.tight_layout()
        if '0' not in str(savefig) and '2' in str(savefig) or 'all' in str(savefig):
            fig2.savefig(figname+'_2.jpeg', dpi=dpi)
            
    #-------------------------------#
    # Write logged data to csv file #
    #-------------------------------#
    if write_data:
        misc_funcs.savedata_to_csv(csv_data, figname+'.csv')

        
    ######################################
    # Additional plots (if user desires) #
    ######################################
    #----------------------------------------------------#
    # Plot forward fringe slope statistics if user wants #
    #----------------------------------------------------#
    if ffs_stats:
        fig3, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(12, 9))
        
        ax1.plot(ffs_outputs_2nd['fringe'], ffs_outputs_2nd['r-squared'], '-', lw=3, label='$r^2$')
        ax1.axvline(xlo, color='tab:red', ls='--', lw=2, label='Linear-region xlo')
        ax1.axvline(xhi, color='tab:red', ls='--', lw=2, label='Linear-region xhi')
        ax1.set_xlabel('Forward fringe strain')
        ax1.set_ylabel('Forward fringe $r^2$')
        ax1.legend()
        
        ax2.plot(ffs_outputs_2nd['fringe'], ffs_outputs_2nd['intercepts'], '-', lw=3, label='Y-intercepts')
        ax2.axvline(xlo, color='tab:red', ls='--', lw=2, label='Linear-region xlo')
        ax2.axvline(xhi, color='tab:red', ls='--', lw=2, label='Linear-region xhi')
        ax2.set_xlabel('Forward fringe strain')
        ax2.set_ylabel('Forward fringe intercepts')
        ax2.legend()
        
        ax3.plot(ffs_outputs_2nd['fringe'], ffs_outputs_2nd['sum-residuals-squared'], '-', lw=3, label='sum($residuals^2$)')
        ax3.axvline(xlo, color='tab:red', ls='--', lw=2, label='Linear-region xlo')
        ax3.axvline(xhi, color='tab:red', ls='--', lw=2, label='Linear-region xhi')
        ax3.set_xlabel('Forward fringe strain')
        ax3.set_ylabel('Forward fringe sum($residuals^2$)')
        ax3.legend()
        
        ax4.plot(ffs_outputs_2nd['fringe'], [abs(i)**0.5 for i in ffs_outputs_2nd['slopes-variance']], '-', lw=3, label='slope std-dev')
        ax4.axvline(xlo, color='tab:red', ls='--', lw=2, label='Linear-region xlo')
        ax4.axvline(xhi, color='tab:red', ls='--', lw=2, label='Linear-region xhi')
        ax4.set_xlabel('Forward fringe strain')
        ax4.set_ylabel('Forward fringe slope std-dev')
        ax4.legend()
        
        ax5.plot(ffs_outputs_2nd['fringe'], [abs(i)**0.5 for i in ffs_outputs_2nd['intercepts-variance']], '-', lw=3, label='intercept std-dev')
        ax5.axvline(xlo, color='tab:red', ls='--', lw=2, label='Linear-region xlo')
        ax5.axvline(xhi, color='tab:red', ls='--', lw=2, label='Linear-region xhi')
        ax5.set_xlabel('Forward fringe strain')
        ax5.set_ylabel('Forward fringe Y-intercept std-dev')
        ax5.legend()
        
        ax6.plot(ffs_outputs_2nd['fringe'], ffs_outputs_2nd['mean-squared-residual'], '-', lw=3, label='mean($residuals^2$)')
        ax6.axvline(xlo, color='tab:red', ls='--', lw=2, label='Linear-region xlo')
        ax6.axvline(xhi, color='tab:red', ls='--', lw=2, label='Linear-region xhi')
        ax6.set_xlabel('Forward fringe strain')
        ax6.set_ylabel('Forward fringe mean($residuals^2$)')
        ax6.legend()
        
        # Apply tight layout to avoid overlapping plots and save plot if user wants
        fig3.tight_layout()
        if '0' not in str(savefig) and '3' in str(savefig) or 'all' in str(savefig):
            fig3.savefig(figname+'_3.jpeg', dpi=dpi)
            
    ###########################################################
    # Plot figures found in the RFR publication if user wants #
    ###########################################################
    # https://stackoverflow.com/questions/45597092/expanded-legend-over-2-subplots
    # https://matplotlib.org/stable/users/explain/axes/legend_guide.html
    if pub_plots and yield_point_derivative:
        # Adjustable plotting parameters
        font = {'family': 'sans-serif',
                'color':  'black',
                'weight': 'normal',
                'size': 12}
        axis_thickness = 3
        yield_color = 'tab:purple'
        regression_color = '#ff9d3aff'
        fixed_bound_color = 'black'
        moving_bound_color = '#bbbbbbff'
        line_width = 6
        legends = False # Show legends
        dim = 4.0 # set subplot dimensiions as 4in x 4in
        
        # Define regions to "regionize" plots
        #          'name':       [lo-X, hi-X, color]
        # regions = {'Toe region': [min(strain), xlo, 'tab:blue'],
        #             'Linear elastic region': [xlo, xhi, 'darkblue'],
        #             'Non-Linear elastic region': [xhi, yield_point_derivative[0], 'lime'],
        #             'Plastic region': [yield_point_derivative[0], max(strain), 'mediumaquamarine']}
        
        regions = {'Toe region': [min(strain), xlo, '#253494ff'],
                    'Linear elastic region': [xlo, xhi, '#2c7fb8ff'],
                    'Non-Linear elastic region': [xhi, yield_point_derivative[0], '#41b6c4ff'],
                    'Plastic region': [yield_point_derivative[0], max(strain), '#a1dab4ff']}
        
            
        # Define function to "regionize" X-, Y-data
        def regionize(x, y, regions):
            regionized = {} # {'name': [[regionize X-data], [regionized Y-data], color]}
            for region in regions:
                lo, hi, color = regions[region]
                rx, ry = misc_funcs.reduce_data(x, y, lo, hi)
                regionized[region] = [rx[:-1], ry[:-1], color]
            return regionized
        
        # Define function to set axis line thickness
        def set_axis_thickness(axis, thickness):
            axis.spines['left'].set_linewidth(thickness)
            axis.spines['right'].set_linewidth(thickness)
            axis.spines['top'].set_linewidth(thickness)
            axis.spines['bottom'].set_linewidth(thickness)
            return
        
        #---------------------------------#
        # "Regionized" stress-strain plot #
        #---------------------------------#
        fig4, ax1 = plt.subplots(figsize=(dim, dim))
        
        # Generate noisy data
        noise = (35/100)*max(stress)*np.random.normal(loc=0.0, scale=1/3, size=len(strain))
        noisy_stress = np.array(stress) + noise
        noisy_strain = np.array(strain)
        
        # Fill in data with some interpolations
        xvals = np.linspace(min(noisy_strain), max(noisy_strain), 3*len(noisy_strain))
        yinterp = np.interp(xvals, noisy_strain, noisy_stress)
        
        
        ax1.plot(xvals, yinterp, '.', ms=8, color='#bbbbbbff', label='MD simulation')
        ax1.plot(strain, stress, '-', lw=line_width, color='#2c7fb8ff', label='Filtered Curve')
        set_axis_thickness(ax1, axis_thickness)
        ax1.set_xlim(xlimits)
        ax1.set_xticks([])
        ax1.set_yticks([])
        #ax1.axis('off')
        if legends:
            ax1.legend()
        fig4.tight_layout()
        if '0' not in str(savefig):
            fig4.savefig(figname+'_published_0.jpeg', dpi=dpi)
            fig4.savefig(figname+'_published_0.eps', dpi=dpi, format='eps')
        
        
        fig4, ax1 = plt.subplots(figsize=(dim, dim))
        regionized_strain_stress = regionize(strain, stress, regions)
        for region in regionized_strain_stress:
            x, y, color = regionized_strain_stress[region]
            ax1.plot(x, y, '-', lw=line_width, color=color, label=region)
        set_axis_thickness(ax1, axis_thickness)
        ax1.set_xlim(xlimits)
        ax1.set_xticks([])
        ax1.set_yticks([])
        #ax1.axis('off')
        if legends:
            ax1.legend()
        fig4.tight_layout()
        if '0' not in str(savefig):
            fig4.savefig(figname+'_published_1.jpeg', dpi=dpi)
            fig4.savefig(figname+'_published_1.eps', dpi=dpi, format='eps')
            
        #---------------------------------------#
        # "Regionized" fringe slope explanation #
        #---------------------------------------#
        fig5, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, figsize=(4*dim, 2*dim))
        regionized_ffs_1st = regionize(ffs_outputs_1st['fringe'], ffs_outputs_1st['slopes'], regions)
        plus_minus_ylimts = (6/100)*(max(ffs_outputs_1st['slopes']) - min(ffs_outputs_1st['slopes'])) 
        fringe_slope_ylimts = (min(ffs_outputs_1st['slopes'])-plus_minus_ylimts, max(ffs_outputs_1st['slopes'])+plus_minus_ylimts)
        
        # Segment 1
        segment_1 = [min(strain), 0.25*max(strain)]
        segment_1_strain, segment_1_stress = misc_funcs.reduce_data(strain, stress, min(segment_1), max(segment_1))
        segment_1_x = misc_funcs.closest(ffs_outputs_1st['fringe'], max(segment_1))
        segment_1_y = ffs_outputs_1st['slopes'][ffs_outputs_1st['fringe'].index(segment_1_x)]
        segment_1_lr = misc_funcs.linear_regression(segment_1_strain, segment_1_stress)
        for region in regionized_strain_stress:
            x, y, color = regionized_strain_stress[region]
            ax1.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax1.plot(segment_1_lr.xreg, segment_1_lr.yreg, '--', lw=line_width/2, color=regression_color, label='Linear regression')
        ax1.axvline(min(segment_1), color=fixed_bound_color, ls='--', lw=line_width/2, label='Fixed bound') 
        ax1.axvline(max(segment_1), color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
            
        for region in regionized_ffs_1st:
            x, y, color = regionized_ffs_1st[region]
            x, y = misc_funcs.reduce_data(x, y, min(segment_1), max(segment_1))
            ax5.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax5.axvline(min(segment_1), color=fixed_bound_color, ls='--', lw=line_width/2, label='Fixed bound') 
        ax5.axvline(max(segment_1), color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
        ax5.plot(segment_1_x, segment_1_y, 'o', ms=12, color=regression_color, label='Linear regression')
        
        set_axis_thickness(ax1, axis_thickness)
        set_axis_thickness(ax5, axis_thickness)
        ax1.set_xlim(xlimits)
        ax5.set_xlim(xlimits)
        ax5.set_ylim(fringe_slope_ylimts)
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax5.set_xticks([])
        ax5.set_yticks([])
            
        # Segment 2
        segment_2 = [min(strain), 0.5*max(strain)]
        segment_2_strain, segment_2_stress = misc_funcs.reduce_data(strain, stress, min(segment_2), max(segment_2))
        segment_2_x = misc_funcs.closest(ffs_outputs_1st['fringe'], max(segment_2))
        segment_2_y = ffs_outputs_1st['slopes'][ffs_outputs_1st['fringe'].index(segment_2_x)]
        segment_2_lr = misc_funcs.linear_regression(segment_2_strain, segment_2_stress)
        for region in regionized_strain_stress:
            x, y, color = regionized_strain_stress[region]
            ax2.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax2.plot(segment_2_lr.xreg, segment_2_lr.yreg, '--', lw=line_width/2, color=regression_color, label='Linear regression')
        ax2.axvline(min(segment_2), color=fixed_bound_color, ls='--', lw=line_width/2, label='Fixed bound') 
        ax2.axvline(max(segment_2), color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
            
        for region in regionized_ffs_1st:
            x, y, color = regionized_ffs_1st[region]
            x, y = misc_funcs.reduce_data(x, y, min(segment_2), max(segment_2))
            ax6.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax6.axvline(min(segment_2), color=fixed_bound_color, ls='--', lw=line_width/2, label='Fixed bound') 
        ax6.axvline(max(segment_2), color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
        ax6.plot(segment_2_x, segment_2_y, 'o', ms=12, color=regression_color, label='Linear regression')
        
        set_axis_thickness(ax2, axis_thickness)
        set_axis_thickness(ax6, axis_thickness)
        ax2.set_xlim(xlimits)
        ax6.set_xlim(xlimits)
        ax6.set_ylim(fringe_slope_ylimts)
        ax2.set_xticks([])
        ax2.set_yticks([])
        ax6.set_xticks([])
        ax6.set_yticks([])
            
        # Segment 3
        segment_3 = [min(strain), 0.75*max(strain)]
        segment_3_strain, segment_3_stress = misc_funcs.reduce_data(strain, stress, min(segment_3), max(segment_3))
        segment_3_x = misc_funcs.closest(ffs_outputs_1st['fringe'], max(segment_3))
        segment_3_y = ffs_outputs_1st['slopes'][ffs_outputs_1st['fringe'].index(segment_3_x)]
        segment_3_lr = misc_funcs.linear_regression(segment_3_strain, segment_3_stress)
        for region in regionized_strain_stress:
            x, y, color = regionized_strain_stress[region]
            ax3.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax3.plot(segment_3_lr.xreg, segment_3_lr.yreg, '--', lw=line_width/2, color=regression_color, label='Linear regression')
        ax3.axvline(min(segment_3), color=fixed_bound_color, ls='--', lw=line_width/2, label='Fixed bound') 
        ax3.axvline(max(segment_3), color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
            
        for region in regionized_ffs_1st:
            x, y, color = regionized_ffs_1st[region]
            x, y = misc_funcs.reduce_data(x, y, min(segment_3), max(segment_3))
            ax7.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax7.axvline(min(segment_3), color=fixed_bound_color, ls='--', lw=line_width/2, label='Fixed bound') 
        ax7.axvline(max(segment_3), color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
        ax7.plot(segment_3_x, segment_3_y, 'o', ms=12, color=regression_color, label='Linear regression')
        
        set_axis_thickness(ax3, axis_thickness)
        set_axis_thickness(ax7, axis_thickness)
        ax3.set_xlim(xlimits)
        ax7.set_xlim(xlimits)
        ax7.set_ylim(fringe_slope_ylimts)
        ax3.set_xticks([])
        ax3.set_yticks([])
        ax7.set_xticks([])
        ax7.set_yticks([])
        
        # Segment 4
        segment_4 = [min(strain), 1.0*max(strain)]
        segment_4_strain, segment_4_stress = misc_funcs.reduce_data(strain, stress, min(segment_4), max(segment_4))
        segment_4_x = misc_funcs.closest(ffs_outputs_1st['fringe'], max(segment_4))
        segment_4_y = ffs_outputs_1st['slopes'][ffs_outputs_1st['fringe'].index(segment_4_x)]
        segment_4_lr = misc_funcs.linear_regression(segment_4_strain, segment_4_stress)
        for region in regionized_strain_stress:
            x, y, color = regionized_strain_stress[region]
            x, y = misc_funcs.reduce_data(x, y, min(segment_4), max(segment_4))
            ax4.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax4.plot(segment_4_lr.xreg, segment_4_lr.yreg, '--', lw=line_width/2, color=regression_color, label='Linear regression')
        ax4.axvline(min(segment_4), color=fixed_bound_color, ls='--', lw=line_width/2, label='Fixed bound') 
        ax4.axvline(max(segment_4), color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
            
        for region in regionized_ffs_1st:
            x, y, color = regionized_ffs_1st[region]
            ax8.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax8.axvline(min(segment_4), color=fixed_bound_color, ls='--', lw=line_width/2, label='Fixed bound') 
        ax8.axvline(max(segment_4), color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
        ax8.plot(segment_4_x, segment_4_y, 'o', ms=12, color=regression_color, label='Linear regression')
            
        set_axis_thickness(ax4, axis_thickness)
        set_axis_thickness(ax8, axis_thickness)
        ax4.set_xlim(xlimits)
        ax8.set_xlim(xlimits)
        ax8.set_ylim(fringe_slope_ylimts)
        ax4.set_xticks([])
        ax4.set_yticks([])
        ax8.set_xticks([])
        ax8.set_yticks([])
            
        fig5.tight_layout()
        if '0' not in str(savefig):
            fig5.savefig(figname+'_published_2.jpeg', dpi=dpi)
            fig5.savefig(figname+'_published_2.eps', dpi=dpi, format='eps')
            
            
        #-----------------------#
        # "Regionized" FRF plot #
        #-----------------------#
        fig6, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(2*dim, 3*dim))
        fringe_slopes = ffs_outputs_1st['slopes'] + rfs_outputs_1st['slopes'] + ffs_outputs_2nd['slopes']
        plus_minus_ylimts = (6/100)*(max(fringe_slopes) - min(fringe_slopes)) 
        fringe_slope_ylimts = (min(fringe_slopes)-plus_minus_ylimts, max(fringe_slopes)+plus_minus_ylimts)
        
        # 1st forward fringe
        ffs_1st_xlo = min(strain) 
        ffs_1st_xhi = csv_data['forward-fringe-slope-1st-maximum'][0][0]
        regionized_ffs_1st = regionize(ffs_outputs_1st['fringe'], ffs_outputs_1st['slopes'], regions)
        ffs_1st_strain, ffs_1st_stress = misc_funcs.reduce_data(strain, stress, ffs_1st_xlo, ffs_1st_xhi)
        ffs_1st_lr = misc_funcs.linear_regression(ffs_1st_strain, ffs_1st_stress)
        for region in regionized_strain_stress:
            x, y, color = regionized_strain_stress[region]
            ax1.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax1.axvline(ffs_1st_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Fixed bound') 
        ax1.axvline(ffs_1st_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
        ax1.plot(ffs_1st_lr.xreg, ffs_1st_lr.yreg, '--', lw=line_width/2, color=regression_color, label='Linear regression')
        
        ffs_1st_xhi_tmp, ffs_1st_yhi = maximized_slope(ffs_outputs_1st['fringe'], ffs_outputs_1st['slopes'], machine_precision=1e-8)
        for region in regionized_ffs_1st:
            x, y, color = regionized_ffs_1st[region]
            ax2.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax2.axvline(ffs_1st_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Fixed bound') 
        ax2.axvline(ffs_1st_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
        #ax2.axhline(ffs_1st_yhi, color=regression_color, ls='--', lw=line_width/2, label='Maximized slope')
        ax2.plot(ffs_1st_xhi_tmp, ffs_1st_yhi, 'o', ms=12, color=regression_color, label='Maximized slope')

        #ax1.set_ylabel('Stress', fontdict=font)
        #ax2.set_ylabel('1st Forward Fringe Slope', fontdict=font)
        set_axis_thickness(ax1, axis_thickness)
        set_axis_thickness(ax2, axis_thickness)
        ax1.set_xlim(xlimits)
        ax2.set_xlim(xlimits)
        ax2.set_ylim(fringe_slope_ylimts)
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax2.set_xticks([])
        ax2.set_yticks([])
        if legends:
            ax1.legend(bbox_to_anchor=(0.0, 1.2, 2, 0.12), loc=2, ncol=3, mode='expand', fontsize=font['size'])

        # 1st reverse fringe
        rfs_1st_xlo = xlo
        rfs_1st_xhi = csv_data['forward-fringe-slope-1st-maximum'][0][0]
        regionized_rfs_1st = regionize(rfs_outputs_1st['fringe'], rfs_outputs_1st['slopes'], regions) 
        rfs_1st_strain, rfs_1st_stress = misc_funcs.reduce_data(strain, stress, rfs_1st_xlo, rfs_1st_xhi)
        rfs_1st_lr = misc_funcs.linear_regression(rfs_1st_strain, rfs_1st_stress)
        for region in regionized_strain_stress:
            x, y, color = regionized_strain_stress[region]
            ax3.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax3.axvline(rfs_1st_xlo, color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
        ax3.axvline(rfs_1st_xhi, color=fixed_bound_color, ls='--', lw=line_width/2, label='Fix bound') 
        ax3.plot(rfs_1st_lr.xreg, rfs_1st_lr.yreg, '--', lw=line_width/2, color=regression_color, label='Linear regression')
        
        rfs_1st_xhi_tmp, rfs_1st_yhi = maximized_slope(rfs_outputs_1st['fringe'], rfs_outputs_1st['slopes'], machine_precision=1e-8)
        for region in regionized_rfs_1st:
            x, y, color = regionized_rfs_1st[region]
            ax4.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax4.axvline(rfs_1st_xlo, color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
        ax4.axvline(rfs_1st_xhi, color=fixed_bound_color, ls='--', lw=line_width/2, label='Fix bound') 
        #ax4.axhline(rfs_1st_yhi, color=regression_color, ls='--', lw=line_width/2, label='Maximized slope')
        ax4.plot(rfs_1st_xhi_tmp, rfs_1st_yhi, 'o', ms=12, color=regression_color, label='Maximized slope')

        #ax3.set_ylabel('Stress', fontdict=font)
        #ax4.set_ylabel('1st Reverse Fringe Slope', fontdict=font)  
        set_axis_thickness(ax3, axis_thickness)
        set_axis_thickness(ax4, axis_thickness)
        ax3.set_xlim(xlimits)
        ax4.set_xlim(xlimits)
        ax4.set_ylim(fringe_slope_ylimts)
        ax3.set_xticks([])
        ax3.set_yticks([])
        ax4.set_xticks([])
        ax4.set_yticks([])
        
        # 2nd forward fringe
        ffs_2nd_xlo = xlo 
        ffs_2nd_xhi = csv_data['forward-fringe-slope-2nd-maximum'][0][0]
        regionized_ffs_2nd = regionize(ffs_outputs_2nd['fringe'], ffs_outputs_2nd['slopes'], regions)  
        ffs_2nd_strain, ffs_2nd_stress = misc_funcs.reduce_data(strain, stress, ffs_2nd_xlo, ffs_2nd_xhi)
        ffs_2nd_lr = misc_funcs.linear_regression(ffs_2nd_strain, ffs_2nd_stress)
        for region in regionized_strain_stress:
            x, y, color = regionized_strain_stress[region]
            ax5.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax5.axvline(ffs_2nd_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Fixed bound') 
        ax5.axvline(ffs_2nd_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
        ax5.plot(ffs_2nd_lr.xreg, ffs_2nd_lr.yreg, '--', lw=line_width/2, color=regression_color, label='Linear regression')
        
        ffs_2nd_xhi_tmp, ffs_2nd_yhi = maximized_slope(ffs_outputs_2nd['fringe'], ffs_outputs_2nd['slopes'], machine_precision=1e-8)
        for region in regionized_ffs_2nd:
            x, y, color = regionized_ffs_2nd[region]
            
            # # For the linear region - interpolate down to xlo (incase minxhi is used)
            # if region == 'Linear elastic region' and minxhi > 0:
            #     dxs = [abs(x[i+1]-x[i]) for i in range(len(x)-1)]
            #     dx = sum(dxs)/len(dxs)
            #     spanx = min(x) - xlo
            #     nx = math.floor(spanx/dx)
            #     xinsert = []; yinsert = []
            #     for n in range(nx):
            #         xinsert.append(min(x) - n*dx)
            #         yinsert.append(y[0])
            #     xinsert.reverse()
            #     xtmp = x.copy(); ytmp = y.copy()
            #     x = xinsert + xtmp
            #     y = yinsert + ytmp
            
            ax6.plot(x, y, '-', lw=line_width, color=color, label=region)
        ax6.axvline(ffs_2nd_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Fixed bound') 
        ax6.axvline(ffs_2nd_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Moving bound') 
        #ax6.axhline(ffs_2nd_yhi, color=regression_color, ls='--', lw=line_width/2, label='Maximized slope')
        ax6.plot(ffs_2nd_xhi_tmp, ffs_2nd_yhi, 'o', ms=12, color=regression_color, label='Maximized slope')

        #ax5.set_xlabel('Strain', fontdict=font)
        #ax6.set_xlabel('Strain', fontdict=font)
        #ax5.set_ylabel('Stress', fontdict=font)
        #ax6.set_ylabel('2nd Forward Fringe Slope', fontdict=font)
        set_axis_thickness(ax5, axis_thickness)
        set_axis_thickness(ax6, axis_thickness)
        ax5.set_xlim(xlimits)
        ax6.set_xlim(xlimits)
        ax6.set_ylim(fringe_slope_ylimts)
        ax5.set_xticks([])
        ax5.set_yticks([])
        ax6.set_xticks([])
        ax6.set_yticks([])
        
        fig6.tight_layout()
        if '0' not in str(savefig):
            fig6.savefig(figname+'_published_3.jpeg', dpi=dpi)
            fig6.savefig(figname+'_published_3.eps', dpi=dpi, format='eps')
            
        #-----------------------------------------------------#
        # "Regionized" yield strength: 2nd derivarive methods #
        #-----------------------------------------------------#
        if isinstance(yp, int) and yp != 0 or yp in ['min-2d', 'min-v', 'max-d', 'min-r2d2']:
            fig7, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(2*dim, 2*dim))
            regionized_ffs_2nd_dervative2 = regionize(dfringe_full, dslopes2_full, regions)  
            
            # Stress-strain
            for region in regionized_strain_stress:
                x, y, color = regionized_strain_stress[region]
                ax1.plot(x, y, '-', lw=line_width, color=color, label=region)
            ax1.axvline(ffs_2nd_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Lower linear elastic region') 
            ax1.axvline(ffs_2nd_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Upper linear elastic region') 
            ax1.axvline(yield_point_derivative[0], color=yield_color, ls='--', lw=line_width/2, label='Yield strain') 
            ax1.plot(yield_point_derivative[0], yield_point_derivative[1], 'o', ms=12, color=yield_color, label='Yield Point')
            
            # 2nd forward fringe
            for region in regionized_ffs_2nd:
                x, y, color = regionized_ffs_2nd[region]
                ax2.plot(x, y, '-', lw=line_width, color=color, label=region)
            ax2.axvline(ffs_2nd_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Lower linear elastic region') 
            ax2.axvline(ffs_2nd_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Upper linear elastic region') 
            ax2.axvline(yield_point_derivative[0], color=yield_color, ls='--', lw=line_width/2, label='Yield strain') 
                
            # 2nd dervative of 2nd forward fringe
            ymins = []; ymaxs = []
            for region in regionized_ffs_2nd_dervative2:
                x, y, color = regionized_ffs_2nd_dervative2[region]
                ax3.plot(x, y, '-', lw=line_width, color=color, label=region)
                if y: ymins.append(min(y)); ymaxs.append(max(y))
            ax3.axvline(ffs_2nd_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Lower linear elastic region') 
            ax3.axvline(ffs_2nd_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Upper linear elastic region') 
            ax3.axvline(yield_point_derivative[0], color=yield_color, ls='--', lw=line_width/2, label='Yield strain') 
            if xvalleys and yvalleys and yp != 'min-2d':
                ax3.plot(xhi_yield, yhi_yield, 'o', ms=12, color=yield_color, label='Yield Point')
            if yp == 'min-2d' or not xvalleys and not yvalleys and  min_2d_fringe is not None and min_2d_slope is not None:
                ax3.plot(min_2d_fringe, min_2d_slope, 'o', ms=12, color=yield_color, label='Yield Point')
            #ax2.set_ylim(fringe_slope_ylimts)
                
            # Define y-limits (may not need for every data set, but the published example would benefit from this)
            ylimits_override = True # over ride default or using all y-value's for maximum y-limit
            if ylimits_override:
                ymins = regionized_ffs_2nd_dervative2['Non-Linear elastic region'][1] + regionized_ffs_2nd_dervative2['Plastic region'][1]
                ymaxs = regionized_ffs_2nd_dervative2['Non-Linear elastic region'][1] + regionized_ffs_2nd_dervative2['Plastic region'][1]
            ymin = min(ymins) # Default
            ymax = max(ymaxs) # Default
            plus_minus = 35/100 # 30% buffer
            span = plus_minus*(ymax-ymin)
            ylimits = (ymin-span, ymax+span)
            
            # Poisson's ratio
            trans1 = None
            if len(t1) >= 1:
                # Find data sets to analyze
                if len(t1) == 1:
                    trans1 = t1[0]
                elif len(t1) >= 2:
                    trans1 = t1[-1] # use the last logged entry
            if trans1 is not None:
                trans = regionize(strain, trans1, regions)
                nu_strain, nu_trans = misc_funcs.reduce_data(strain, trans1, ffs_2nd_xlo, yield_point_derivative[0])
                nu_lr = misc_funcs.linear_regression(nu_strain, nu_trans)
                for region in trans:
                    x, y, color = trans[region]
                    ax4.plot(x, y, '-', lw=line_width, color=color, label=region)
                
                ax4.plot(nu_lr.xreg, nu_lr.yreg, '--', lw=line_width/2, color=regression_color, label='Linear regression')
                ax4.axvline(ffs_2nd_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Lower linear elastic region') 
                ax4.axvline(ffs_2nd_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Upper linear elastic region') 
                ax4.axvline(yield_point_derivative[0], color=yield_color, ls='--', lw=line_width/2, label='Yield strain') 
                set_axis_thickness(ax4, axis_thickness)
                ax4.set_xlim(xlimits)
                ax4.set_xticks([])
                ax4.set_yticks([])
                #ax1.axis('off')
                if legends:
                    ax2.legend()
            
            set_axis_thickness(ax1, axis_thickness)
            set_axis_thickness(ax2, axis_thickness)
            set_axis_thickness(ax3, axis_thickness)
            set_axis_thickness(ax4, axis_thickness)
            ax1.set_xlim(xlimits)
            ax2.set_xlim(xlimits)
            ax2.set_ylim(fringe_slope_ylimts)
            ax3.set_xlim(xlimits) 
            ax3.set_ylim(ylimits)  
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax3.set_xticks([])
            ax3.set_yticks([])
            ax4.set_xticks([])
            ax4.set_yticks([])
            if legends:
                ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.75), fancybox=True, ncol=1)
                
            fig7.tight_layout()
            if savefig:
                fig7.savefig(figname+'_published_4.jpeg', dpi=dpi)
                fig7.savefig(figname+'_published_4.eps', dpi=dpi, format='eps')
                
        #---------------------------------------------#
        # "Regionized" yield strength: ffs 3rd method #
        #---------------------------------------------#
        if yp == 'max-3ffs':
            fig7, (ax1, ax2) = plt.subplots(2, 1, figsize=(1*dim, 2*dim))
            
            # Stress-strain
            for region in regionized_strain_stress:
                x, y, color = regionized_strain_stress[region]
                ax1.plot(x, y, '-', lw=line_width, color=color, label=region)
            ax1.axvline(ffs_2nd_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Lower linear elastic region') 
            ax1.axvline(ffs_2nd_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Upper linear elastic region') 
            ax1.axvline(yield_point_derivative[0], color=regression_color, ls='--', lw=line_width/2, label='Yield strain') 
            ax1.plot(yield_point_derivative[0], yield_point_derivative[1], 'o', ms=12, color=regression_color, label='Yield Point')
            
            # 3rd forward fringe
            regionized_ffs_3rd = regionize(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['slopes'], regions)  
            for region in regionized_ffs_3rd:
                x, y, color = regionized_ffs_3rd[region]
                ax2.plot(x, y, '-', lw=line_width, color=color, label=region)
            ax2.axvline(ffs_2nd_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Lower linear elastic region') 
            ax2.axvline(ffs_2nd_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Upper linear elastic region') 
            ax2.axvline(yield_point_derivative[0], color=regression_color, ls='--', lw=line_width/2, label='Yield strain') 
            
            set_axis_thickness(ax1, axis_thickness)
            set_axis_thickness(ax2, axis_thickness)
            ax1.set_xlim(xlimits)
            ax2.set_xlim(xlimits)
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax2.set_xticks([])
            ax2.set_yticks([])
            if legends:
                ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.75), fancybox=True, ncol=1)
                
            fig7.tight_layout()
            if savefig:
                fig7.savefig(figname+'_published_4.jpeg', dpi=dpi)
                fig7.savefig(figname+'_published_4.eps', dpi=dpi, format='eps')
            
        #------------------------------------------#
        # "Regionized" yield strength: r^2 methods #
        #------------------------------------------#    
        if yp in ['mean(r2)-3s', 'max(r2)-3s']:  
            fig7, (ax1, ax2) = plt.subplots(2, 1, figsize=(1*dim, 2*dim))
            
            # Stress-strain
            for region in regionized_strain_stress:
                x, y, color = regionized_strain_stress[region]
                ax1.plot(x, y, '-', lw=line_width, color=color, label=region)
            ax1.axvline(ffs_2nd_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Lower linear elastic region') 
            ax1.axvline(ffs_2nd_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Upper linear elastic region') 
            ax1.axvline(yield_point_derivative[0], color=regression_color, ls='--', lw=line_width/2, label='Yield strain') 
            ax1.plot(yield_point_derivative[0], yield_point_derivative[1], 'o', ms=12, color=regression_color, label='Yield Point')
            
            
            # 3rd forward fringe
            regionized_ffs_3rd = regionize(ffs_outputs_3rd['fringe'], ffs_outputs_3rd['r-squared'], regions)  
            for region in regionized_ffs_3rd:
                x, y, color = regionized_ffs_3rd[region]
                ax2.plot(x, y, '-', lw=line_width, color=color, label=region)
            ax2.axvline(ffs_2nd_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Lower linear elastic region') 
            ax2.axvline(ffs_2nd_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Upper linear elastic region') 
            ax2.axvline(yield_point_derivative[0], color=regression_color, ls='--', lw=line_width/2, label='Yield strain') 
            
            if yp == 'mean(r2)-3s':
                ax2.axhline(r2_lo, color='tab:purple', ls='--', lw=1.5, label='Linear $r^2$ ($\\mu$ - 3*$\\sigma$)')
                ax2.axhline(r2_hi, color='tab:purple', ls='--', lw=1.5, label='Linear $r^2$ ($\\mu$ + 3*$\\sigma$)')
            if yp == 'max(r2)-3s':
                ax2.axhline(r2_lo, color='tab:purple', ls='--', lw=1.5, label='Linear $r^2$ (max - 3*$\\sigma$)')
                ax2.axhline(r2_hi, color='tab:purple', ls='--', lw=1.5, label='Linear $r^2$ (max + 3*$\\sigma$)')
                        
            set_axis_thickness(ax1, axis_thickness)
            set_axis_thickness(ax2, axis_thickness)
            ax1.set_xlim(xlimits)
            ax2.set_xlim(xlimits)
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax2.set_xticks([])
            ax2.set_yticks([])
            if legends:
                ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.75), fancybox=True, ncol=1)
                #ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.75), fancybox=True, ncol=1)
                
            fig7.tight_layout()
            if savefig:
                fig7.savefig(figname+'_published_4.jpeg', dpi=dpi)
                fig7.savefig(figname+'_published_4.eps', dpi=dpi, format='eps')
            
        #------------------------------------------#
        # Generate axial vs transverse strain plot #
        #------------------------------------------#
        trans1 = None
        if len(t1) >= 1:
            # Find data sets to analyze
            if len(t1) == 1:
                trans1 = t1[0]
            elif len(t1) >= 2:
                trans1 = t1[-1] # use the last logged entry
        if trans1 is not None:
            fig8, (ax1, ax2) = plt.subplots(2, 1, figsize=(1*dim, 2*dim))
            for region in regionized_strain_stress:
                x, y, color = regionized_strain_stress[region]
                ax1.plot(x, y, '-', lw=line_width, color=color, label=region)
            ax1.axvline(ffs_2nd_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Lower linear elastic region') 
            ax1.axvline(ffs_2nd_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Upper linear elastic region') 
            ax1.axvline(yield_point_derivative[0], color=yield_color, ls='--', lw=line_width/2, label='Yield strain') 
            set_axis_thickness(ax1, axis_thickness)
            ax1.set_xlim(xlimits)
            ax1.set_xticks([])
            ax1.set_yticks([])
            #ax1.axis('off')
            
            
            trans = regionize(strain, trans1, regions)
            nu_strain, nu_trans = misc_funcs.reduce_data(strain, trans1, ffs_2nd_xlo, yield_point_derivative[0])
            nu_lr = misc_funcs.linear_regression(nu_strain, nu_trans)
            for region in trans:
                x, y, color = trans[region]
                ax2.plot(x, y, '-', lw=line_width, color=color, label=region)
            
            ax2.plot(nu_lr.xreg, nu_lr.yreg, '--', lw=line_width/2, color=regression_color, label='Linear regression')
            ax2.axvline(ffs_2nd_xlo, color=fixed_bound_color, ls='--', lw=line_width/2, label='Lower linear elastic region') 
            ax2.axvline(ffs_2nd_xhi, color=moving_bound_color, ls='--', lw=line_width/2, label='Upper linear elastic region') 
            ax2.axvline(yield_point_derivative[0], color=yield_color, ls='--', lw=line_width/2, label='Yield strain') 
            set_axis_thickness(ax2, axis_thickness)
            ax2.set_xlim(xlimits)
            ax2.set_xticks([])
            ax2.set_yticks([])
            #ax1.axis('off')
            if legends:
                ax2.legend()
                
            fig8.tight_layout()
            if '0' not in str(savefig):
                fig8.savefig(figname+'_published_5.jpeg', dpi=dpi)
                fig8.savefig(figname+'_published_5.eps', dpi=dpi, format='eps')
            
            
        #-----------------------------------------------------------------#
        # Option to "linearize" the "Linear elastic region and write a    #
        # new logfile, as this region may not be perfectly linear. This   #
        # is meant to generate a nice quality example for the manuscript. #
        #-----------------------------------------------------------------#
        linearize = False
        if linearize:
            # Apply a linear regression to the linear region and "rebuild Y-data"
            b0 = 0; b1 = 0
            for region in regionized_strain_stress:
                x, y, color = regionized_strain_stress[region]
                if region == 'Linear elastic region':
                    elastic_regression = misc_funcs.linear_regression(x, y)
                    b0 = elastic_regression.b0
                    b1 = elastic_regression.b1
                    y = [i*b1 + b0 for i in x]
                    regionized_strain_stress[region] = [x, y, color]

            # Write data to a log file
            with open(figname+'_linearized.txt', 'w') as f:
                f.write('{} {}\n'.format('strain', 'stress'))
                for region in regionized_strain_stress:
                    x, y, color = regionized_strain_stress[region]
                    for i, j in zip(x, y):
                        f.write('{} {}\n'.format(i, j))
        
    return xlo, xhi, yield_point_derivative, yield_point_offset, nu1, nu2, nu12, ultimate_point