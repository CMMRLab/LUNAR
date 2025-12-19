# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
November 18th, 2024
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
import src.log_analysis.misc_funcs as misc_funcs
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import numpy as np
import math




################################################################
# Function to find peaks and valleys in a data set using scipy #
################################################################
def find_peaks_and_valleys(x, y, prominence=None):
    from scipy.signal import find_peaks
    
    # Find peaks
    xdata = np.array(x); ydata = np.array(y);
    peaks, properties = find_peaks(ydata, prominence=prominence)
    xpeaks = list(xdata[peaks]); ypeaks = list(ydata[peaks])
    
    # Find valleys
    xvalleys = []; yvalleys = []; valley_depths = [] # [ (avg, small, large) ]
    if len(xpeaks) >= 2 and len(ypeaks) >= 2:
        for i in range(len(peaks)-1):
            lo = peaks[i]; hi = peaks[i+1];
            between_peaksx = xdata[lo:hi]
            between_peaksy = ydata[lo:hi]
            minimum_index = np.min(np.where(between_peaksy == between_peaksy.min())[0])
            xvalleys.append( between_peaksx[minimum_index] )
            yvalleys.append( between_peaksy[minimum_index] ) 
            
            # Compute percent difference in y-direction
            y_peak_avg = (ydata[lo] + ydata[hi])/2
            y_valley = between_peaksy[minimum_index]
            avg = y_valley - y_peak_avg
            small = y_valley - min(lo, hi)
            large = y_valley - max(lo, hi)
            valley_depths.append( (avg, small, large) )
    return xpeaks, ypeaks, xvalleys, yvalleys, valley_depths



##################################################################################
# Function to automatically determine the best mirroring options for iFFT filter #
##################################################################################
def iFFT_quadrant_mirroring(xdata, ydata, threshold, quadrant_mirror):
    #------------------------------------------------------------------------#
    # Set default settings for calling butter_lowpass_filter (since multiple #
    # tests will be run, we DO NOT WANT to save data or figures.)            #
    #------------------------------------------------------------------------#
    savefig = False; figname = ''; dpi = 300
    
    #-------------------------------------------------#
    # Determine half_data to only check for residuals #
    # either from lo-half_data or half_data-hi        #
    #-------------------------------------------------#
    half_data = int(xdata.shape[0]/2)
    
    #------------------------------#
    # First: Optimize the "lo" end #
    #------------------------------#
    lo_quads2test = [1, 2, 3, 4]; lo_summed_residuals2 = {} # {quadrant_mirror:sum-of-residuals-squared}
    for quad in lo_quads2test: 
        quadrants = '{},{}'.format(quad, 1) # hold hi constant at 1
        y_filter, qm = iFFT_filter(xdata, ydata, threshold, quadrants, savefig, figname, dpi, plot_PSD=False)
        residuals = ydata - y_filter
        residuals = residuals[:half_data] # we only care about the first half fit
        lo_summed_residuals2[quad] = np.sum(residuals**2)
        
    # Find minimized sum of residuals squared
    lo = min(lo_summed_residuals2, key=lo_summed_residuals2.get)
    
    #-------------------------------#
    # Second: Optimize the "hi" end #
    #-------------------------------#
    hi_quads2test = [1, 2, 3, 4]; hi_summed_residuals2 = {} # {quadrant_mirror:sum-of-residuals-squared}
    for quad in hi_quads2test: 
        quadrants = '{},{}'.format(1, quad) # hold lo constant at 1
        y_filter, qm = iFFT_filter(xdata, ydata, threshold, quadrants, savefig, figname, dpi, plot_PSD=False)
        residuals = ydata - y_filter
        residuals = residuals[half_data:] # we only care about the last half fit
        hi_summed_residuals2[quad] = np.sum(residuals**2)
        
    # Find minimized sum of residuals squared
    hi = min(hi_summed_residuals2, key=hi_summed_residuals2.get)
    
    #------------------------------------#
    # Set optimal quadrant_mirror string #
    #------------------------------------#
    optimal_quadrant_mirror = '{},{}'.format(lo, hi)
    if '-p' in str(quadrant_mirror):
        optimal_quadrant_mirror = '{},{}-p'.format(lo, hi)
    return optimal_quadrant_mirror
    


###############################################
# Function to implement an inverse FFT filter #
###############################################
def iFFT_filter(x, y, threshold, quadrant_mirror, savefig, figname, dpi, plot_PSD=True):
    xdata = np.array(x); ydata = np.array(y)

    #-----------------------------------------------------------------------#
    # Automagically detect which are the best quandrant mirroring locations #
    #-----------------------------------------------------------------------#
    if 'msr' in str(quadrant_mirror):
        quadrant_mirror = iFFT_quadrant_mirroring(xdata, ydata, threshold, quadrant_mirror)
        print(f'    Optimized quadrant_mirror for iFFT filter was found to be {quadrant_mirror}.')   
        
    #----------------------------------------------------------------------------#
    # Perfrom quadrant mirroring operations if two digits are in quadrant_mirror #
    #----------------------------------------------------------------------------#
    digits = [int(i) for i in str(quadrant_mirror) if i.isdigit()]
    if len(digits) == 2:
        lo, hi = digits
        xdata, ydata, lo_trim, hi_trim, lo_xdata, lo_ydata, hi_xdata, hi_ydata = data_extension(xdata, ydata, lo=lo, hi=hi)
        
        # Plot the data if the '-p' flag is at the end of quadrant_mirror
        if str(quadrant_mirror).endswith('-p'):
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(lo_xdata, lo_ydata, 'o', ms=4, color='tab:green', label='lo data in quadrant {}'.format(lo))
            ax.plot(xdata[lo_trim:hi_trim], ydata[lo_trim:hi_trim], 'o', ms=4, color='tab:blue', label='Orignal data w/o mirror')
            ax.plot(hi_xdata, hi_ydata, 'o', ms=4, color='tab:cyan', label='hi data in quadrant {}'.format(hi))
            ax.axhline(0, ls='--', color='black', lw=1)
            ax.axvline(0, ls='--', color='black', lw=1)
            ax.legend()
    
    #--------------------------#
    # Compute and plot the PSD #
    #--------------------------#
    x_fft, y_fft, x_psd, y_psd, fs, N = compute_FFT(xdata, ydata)
    if 'mean' in str(threshold):
        # Attempt getting scaling factor as threshold could be:
        #   'mean'
        #   'mean:SF', where 'SF' is the optional scaling factor (e.g.
        #              threshold='mean:0.25', means scale mean by 0.25)
        try: scale_factor = float(threshold.split(':')[-1])
        except: scale_factor = 1.0
        
        # Compute threshold
        threshold = scale_factor*np.mean(y_psd)
    elif isinstance(threshold, (float, int)):
        threshold = threshold
    else:
        threshold = np.mean(y_psd)
        print(f'WARNING threshold={threshold} is not supported. Defaulting to mean of PSD ({threshold})')

    #-----------------------------------------------------------------------------------#
    # Set a scaling factor of 0 or 1 to cancel out (0) or leave (1) certain frequencies #
    #-----------------------------------------------------------------------------------#
    scaling_factors = np.zeros_like(y_psd)
    for i, mag in enumerate(y_psd):
        # Kill small Fourier Coeffs
        if mag < threshold:
            scaling_factor = 0
        else:
            scaling_factor = 1
        scaling_factors[i] = scaling_factor
    
    #-----------------------------------------------------------------------------------------------#
    # Cancel or leave frequencies in y_ftt and then inverse the cleaned fft to get the filterd data #
    #-----------------------------------------------------------------------------------------------#
    y_clean = scaling_factors*y_fft
    y_filter = np.fft.irfft(y_clean)
    if y_filter.shape != xdata.shape:
        y_filter = np.append(y_filter, y_filter[-1])
    
    # If quadrant mirroring and plotting, append the filter response to the plot
    if str(quadrant_mirror).endswith('-p') and len(digits) == 2:
        ax.plot(xdata, y_filter, '-', lw=4, color='tab:orange', label='Filtered data')
        ax.legend()
        fig.tight_layout()
        if '0' not in str(savefig) and '2' in str(savefig) or 'all' in str(savefig):
            fig.savefig(figname+'_qm.jpeg', dpi=dpi)
    
    # If quadrant mirroring was used, get orginal length of data and
    if len(digits) == 2: y_filter = y_filter[lo_trim:hi_trim] 
    
    #------------------------------#
    # Plot and save PSD if desired #
    #------------------------------#
    if plot_PSD:
        fig, ax = plt.subplots(1, figsize=(6, 4))
        ax.stem(x_psd, y_psd, linefmt='tab:blue', markerfmt='.', label='$|X(f)|^2/N$')
        ax.axhline(threshold, ls='--', color='tab:olive', lw=1, label='threshold={:.8f}'.format(threshold))
        ax.set_xlabel('Frequency (1/X-units)', fontsize=12)
        ax.set_ylabel('Power', fontsize=12)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), fancybox=True, ncol=1, fontsize=12)
        fig.tight_layout()
        if '0' not in str(savefig) and '1' in str(savefig) or 'all' in str(savefig):
            fig.savefig(figname+'_PSD.jpeg', dpi=dpi)
    return y_filter, quadrant_mirror



##############################################################################
# Function to extend data at the "lo" and "hi" end. This function assumes    #
# xdata to be increasing and ydata to be the dependant variable. Parameters  #
# and their meanings:                                                        #
#    xdata = numpy array of X-data                                           #
#    ydata = numpy array of Y-data                                           #
#                                                                            #
#    lo = integer of 1, 2, or 3 to set "lo" xdata end (or start of xdata)    #
#    hi = integer of 1, 2, 3 or 4 to set "hi" xdata end (or end of xdata)    #
#      lo/hi integer meanings can be thought of as quadrant numbers for "lo" #
#      and have the same meaning for "hi", but at the tail end of the data   #
#      1 = means shut off data mirroring.                                    #
#      2 = means mirror data (lo=mirrored into quadrant 2)                   #
#      3 = means apply 180 degree rotation (lo=mirrored into quadrant 3)     #
#      4 = means flip the Y-data and append at the end of the Y-data         #
#                                                                            #
#    index = integer to set "symmetry" index. For example zero means         #
#            lo-symmetry index is zero (i.e. xdata[0] and ydata[0]) and      #
#            hi-symmetry index is 0-1 (i.e xdata[-1] and ydata[-1])          #
##############################################################################
def data_extension(xdata, ydata, lo=1, hi=1, index=0):
    # Setup number of data points based on index mirroring
    ndata = xdata.shape[0] - (index+1)
    
    # Perform lo padding operations
    if lo == 2:
        lo_xdata = min(xdata) + xdata[index] - xdata[::-1][index+1:]
        lo_ydata = ydata[::-1][index+1:]
    elif lo == 3:
        lo_xdata = min(xdata) + xdata[index] - xdata[::-1][index+1:]
        lo_ydata = ydata[index] - ydata[::-1][index+1:]
    elif lo == 4:
        lo_xdata = min(xdata) + xdata[index] - xdata[::-1][index+1:]
        lo_ydata = ydata[index+1:] - (ydata[-(index+1)] - ydata[index]) 
    else:
        lo_xdata = np.array([])
        lo_ydata = np.array([])
    
    # Perform hi padding operations
    if hi == 2:
        hi_xdata = -min(xdata) + max(xdata) + xdata[index+1:]  
        hi_ydata = ydata[::-1][index+1:]
    elif hi == 3:
        hi_xdata = -min(xdata) + max(xdata) + xdata[index+1:]  
        hi_ydata = ydata[index] - ydata[::-1][index+1:] + 2*ydata[-(index+1)] + ydata[index]
    elif hi == 4:
        hi_xdata = -min(xdata) + max(xdata) + xdata[index+1:]  
        hi_ydata = ydata[-(index+1)] + ydata[index+1:]
    else:
        hi_xdata = np.array([])
        hi_ydata = np.array([])
        
    # Assemble data
    if lo in [2, 3, 4] and hi in [2, 3, 4]:
        xdata = np.concatenate((lo_xdata, xdata, hi_xdata), axis=0)
        ydata = np.concatenate((lo_ydata, ydata, hi_ydata), axis=0)  
        lo_trim = ndata
        hi_trim = -ndata
    elif lo in [2, 3, 4] and hi == 1:
        xdata = np.concatenate((lo_xdata, xdata), axis=0)
        ydata = np.concatenate((lo_ydata, ydata), axis=0)
        lo_trim = ndata
        hi_trim = xdata.shape[0]
    elif lo == 1 and hi in [2, 3, 4]:
        xdata = np.concatenate((xdata, hi_xdata), axis=0)
        ydata = np.concatenate((ydata, hi_ydata), axis=0)  
        lo_trim = 0
        hi_trim = -ndata
    else:
        lo_trim = 0
        hi_trim = xdata.shape[0]
    return xdata, ydata, lo_trim, hi_trim, lo_xdata, lo_ydata, hi_xdata, hi_ydata


##########################################################################
# Function to automatically determine the optimal data extension integer #
##########################################################################
def determine_mirroring_locations(xdata, ydata, wn, order, quadrant_mirror):
    #------------------------------------------------------------------------#
    # Set default settings for calling butter_lowpass_filter (since multiple #
    # tests will be run, we DO NOT WANT to save data or figures.)            #
    #------------------------------------------------------------------------#
    write_data = False; savefig = False; figname = ''; dpi = 300
    
    #-------------------------------------------------#
    # Determine half_data to only check for residuals #
    # either from lo-half_data or half_data-hi        #
    #-------------------------------------------------#
    half_data = int(xdata.shape[0]/2)
    
    #------------------------------#
    # First: Optimize the "lo" end #
    #------------------------------#
    lo_quads2test = [1, 2, 3, 4]; lo_summed_residuals2 = {} # {quadrant_mirror:sum-of-residuals-squared}
    for quad in lo_quads2test: 
        quadrants = '{},{}'.format(quad, 1) # hold hi constant at 1
        ybutter, qm = butter_lowpass_filter(xdata, ydata, wn, order, quadrants, write_data, savefig, figname, dpi, pflag=False)
        residuals = ydata - ybutter
        residuals = residuals[:half_data] # we only care about the first half fit
        lo_summed_residuals2[quad] = np.sum(residuals**2)
        
    # Find minimized sum of residuals squared
    lo = min(lo_summed_residuals2, key=lo_summed_residuals2.get)
    
    #-------------------------------#
    # Second: Optimize the "hi" end #
    #-------------------------------#
    hi_quads2test = [1, 2, 3, 4]; hi_summed_residuals2 = {} # {quadrant_mirror:sum-of-residuals-squared}
    for quad in hi_quads2test: 
        quadrants = '{},{}'.format(1, quad) # hold lo constant at 1
        ybutter, qm = butter_lowpass_filter(xdata, ydata, wn, order, quadrants, write_data, savefig, figname, dpi, pflag=False)
        residuals = ydata - ybutter
        residuals = residuals[half_data:] # we only care about the last half fit
        hi_summed_residuals2[quad] = np.sum(residuals**2)
        
    # Find minimized sum of residuals squared
    hi = min(hi_summed_residuals2, key=hi_summed_residuals2.get)
    
    #------------------------------------#
    # Set optimal quadrant_mirror string #
    #------------------------------------#
    optimal_quadrant_mirror = '{},{}'.format(lo, hi)
    if '-p' in str(quadrant_mirror):
        optimal_quadrant_mirror = '{},{}-p'.format(lo, hi)
    return optimal_quadrant_mirror


##############################################
# Function to define of a butterworth filter #
##############################################
def butter_lowpass_filter(xdata, ydata, wn, order, quadrant_mirror, write_data, savefig, figname, dpi, pflag=True):
    from scipy.signal import butter, sosfiltfilt
    
    #--------------------------------------------------------------#
    # Setup csv data structure to append to for quadrant mirroring #
    #--------------------------------------------------------------#
    csv_data = {} # { 'name' : [[xdata], [ydata]], ... }
    
    #-------------------------------------------------------------------------------------#
    # Check to see if xdata is increassing or decreasing. If xdata is decreasing, we will #
    # flip before filtering and then flip the filtered data back to the orginal direcion. #
    #-------------------------------------------------------------------------------------#
    decreasing_x = False
    if xdata[0] >= xdata[-1]:
        decreasing_x = True
        xdata = xdata[::-1]
        ydata = ydata[::-1]
    
    #-----------------------------------------------------------------------#
    # Automagically detect which are the best quandrant mirroring locations #
    #-----------------------------------------------------------------------#
    if 'msr' in str(quadrant_mirror):
        quadrant_mirror = determine_mirroring_locations(xdata, ydata, wn, order, quadrant_mirror)
        if pflag: print(f'    Optimized quadrant_mirror for butterworth filter was found to be {quadrant_mirror}.')
        
    #----------------------------------------------------------------------------#
    # Perfrom quadrant mirroring operations if two digits are in quadrant_mirror #
    #----------------------------------------------------------------------------#
    digits = [int(i) for i in str(quadrant_mirror) if i.isdigit()]
    if len(digits) == 2:
        lo, hi = digits
        xdata, ydata, lo_trim, hi_trim, lo_xdata, lo_ydata, hi_xdata, hi_ydata = data_extension(xdata, ydata, lo=lo, hi=hi)
            
        # Plot the data if the '-p' flag is at the end of quadrant_mirror
        if str(quadrant_mirror).endswith('-p'):
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(lo_xdata, lo_ydata, 'o', ms=4, color='tab:green', label='lo data in quadrant {}'.format(lo))
            ax.plot(xdata[lo_trim:hi_trim], ydata[lo_trim:hi_trim], 'o', ms=4, color='tab:blue', label='Orignal data w/o mirror')
            ax.plot(hi_xdata, hi_ydata, 'o', ms=4, color='tab:cyan', label='hi data in quadrant {}'.format(hi))
            ax.axhline(0, ls='--', color='black', lw=1)
            ax.axvline(0, ls='--', color='black', lw=1)
            ax.legend()
            
            csv_data['lo-QM={}'.format(lo)] = [list(lo_xdata), list(lo_ydata)]
            csv_data['hi-QM={}'.format(hi)] = [list(hi_xdata), list(hi_ydata)]
            csv_data['original-data'] = [list(xdata[lo_trim:hi_trim]), list(ydata[lo_trim:hi_trim])]
    
    #-----------------------------------------------------------#
    # Get the filter coefficients and filter the prepared ydata #
    #-----------------------------------------------------------#
    sos = butter(order, wn, btype='low', analog=False, output='sos', fs=None)
    y = sosfiltfilt(sos, ydata, axis=-1, padtype=None)
    csv_data['filtered-data'] = [list(xdata), list(y)]
    
    # If quadrant mirroring and plotting, append the filter response to the plot
    if str(quadrant_mirror).endswith('-p') and len(digits) == 2:
        ax.plot(xdata, y, '-', lw=4, color='tab:orange', label='Filtered data')
        ax.legend()
        
        # Comment/uncomment (Josh's tmp solution to labeling plot)
        # fs=12
        # ax.set_xlabel('True Strain', fontsize=fs)
        # ax.set_ylabel('True Stress (MPa)', fontsize=fs)
        # ax.tick_params(axis='both', which='major', labelsize=fs)
        if '0' not in str(savefig) and '2' in str(savefig) or 'all' in str(savefig):
            fig.tight_layout()
            fig.savefig(figname+'.jpeg', dpi=dpi)

    # If users want data to be written, write the data 
    if write_data:
        misc_funcs.savedata_to_csv(csv_data, figname+'.csv')
    
    # If quadrant mirroring was used, get orginal length of data and
    if len(digits) == 2: y = y[lo_trim:hi_trim] 
    if decreasing_x: y = y[::-1]
    return y, quadrant_mirror


##############################################
# Function to convert power to decibels (db) #
##############################################
def power_to_db(power, ref_power=1):
  return 10*np.log10(power/ref_power)


##################################################
# Function to compute X and Y componets of a FFT #
##################################################
def compute_FFT(x, y):
    # Convert to numpy arrays
    xx = np.array(x)
    yy = np.array(y)
    
    # Define sampling rate and number of data points
    N = xx.shape[0]                      # number of data points
    fs = (N-1)/(np.max(xx) - np.min(xx)) # sampling rate
    d = 1/fs                             # sampling space

    # Perform one sided FFT
    y_fft = np.fft.rfft(yy, axis=0, norm='backward')
    freq = np.fft.rfftfreq(N, d=d)
    
    # Compute one sided power spectral density from FFT
    y_psd = np.real( (y_fft*np.conjugate(y_fft))/N )
    
    # Compute onde sided amplitudes from FFT
    y_amp = np.abs(y_fft)/N 
    y_amp[1:-1] *= 2
    if N % 2 == 0:
        y_amp[-1] /= 2
    return freq, y_fft, y_amp, y_psd, fs, N


###############################################################################
# Function to plot FFT, Power spectrum, and Wn; if Wn is supplied by the user #
###############################################################################
def generate_and_plot_fft_wn(xdata, ydata, wn, write_data, savefig, figname, dpi):
    # Compute power spectrum and find wns
    freq, y_fft, y_amp, y_psd, fs, N = compute_FFT(xdata, ydata)
    wns = freq/(0.5*fs)
    
    # Plot the FFT, Power Spectrum, and nearest Wn    
    plot_fft_and_wn(freq, y_fft, y_amp, y_psd, wns, fs, wn=wn, threshold=None, mode='psd', figname=figname, write_data=False, savefig=savefig, dpi=dpi)
    return


#############################################################
# Function to plot psd or amplitude from fft and a wn value #
#############################################################
def plot_fft_and_wn(freq, y_fft, y_amp, y_psd, wns, fs, wn=None, threshold=None, mode='psd', figname='FFT_PSD_amp', write_data=False, savefig=False, dpi=300):
    #-----------------------------------------#
    # Setup settings and values based on mode #
    #-----------------------------------------#
    # Check for a provided wn
    if wn is not None:
        try: index = np.absolute(wns-wn).argmin()
        except: index = None
    else: index = None
    
    # Set up values
    if mode == 'psd':
        y = y_psd
        ylabel = r'Power (Y-units$^2$)'
        spectral_label = r'$|X(f)|^2/N$'
    elif mode == 'amp':
        y = y_amp
        ylabel = r'Amplitude (Y-units)'
        spectral_label = r'$||X(f)||/N$'
    else:
        raise Exception(f'ERROR mode={mode} is not supported. Supported modes are "psd" or "amp"')
    
    
    #---------------------------#
    # Plot data on linear-scale #
    #---------------------------#
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    ax1.stem(freq, y, linefmt='tab:blue', markerfmt='.', label=spectral_label)
    if threshold is not None:
        ax1.axhline(threshold, ls='--', color='tab:olive', lw=1, label='Threshold={:.8f}'.format(threshold))
    if wn is not None and index is not None:
        ax1.plot(freq[index], y[index],  'o', ms=5, color='tab:orange', label='$W_n$ (x={:.6f}; y={:.6f})'.format(freq[index], y[index]))
        ax1.axvline(freq[index], ls='--', color='tab:red', lw=1, label='Cutoff={:.8f}'.format(freq[index]))
    ax1.set_xlabel('Frequency (1/X-units)', fontsize=12)
    ax1.set_ylabel(ylabel, fontsize=12)
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), fancybox=True, ncol=1, fontsize=12)
    
    # Generate second x-axis for axis1
    def forward_conversion1(freq):
        return freq/(0.5*fs)
    def reverse_conversion1(freq):
        return freq
    second_x1 = ax1.secondary_xaxis('top', functions=(forward_conversion1, reverse_conversion1))
    second_x1.set_xlabel('$W_n$', fontsize=12)
    
    
    #------------------------#
    # Plot data on log-scale #
    #------------------------#
    ax2.stem(wns, y, linefmt='tab:blue', markerfmt='.', label=spectral_label)
    ax2.set_yscale('log')
    if threshold:
        ax2.axhline(threshold, ls='--', color='tab:olive', lw=1, label='Threshold={:.8f}'.format(threshold))
    if wn is not None and index is not None:
        ax2.plot(wns[index], y[index],  'o', ms=5, color='tab:orange', label='$W_n$ (x={:.6f}; y={:.6f})'.format(wns[index], y[index]))
        ax2.axvline(wns[index], ls='--', color='tab:red', lw=1, label='Cutoff={:.8f}'.format(wns[index]))
    ax2.set_xlabel('$W_n$', fontsize=12)
    ax2.set_ylabel(ylabel, fontsize=12)
    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), fancybox=True, ncol=1, fontsize=12)
    ax2.xaxis.set_major_formatter(ScalarFormatter())
    
    # Generate second x-axis for axis2
    def forward_conversion2(wns):
        return (0.5*fs)*wns
    def reverse_conversion2(wns):
        return wns
    second_x2 = ax2.secondary_xaxis('top', functions=(forward_conversion2, reverse_conversion2))
    second_x2.set_xlabel('Frequency (1/X-units)', fontsize=12)
    
    # Finalize plot
    fig.tight_layout()
    if savefig:
        fig.savefig(figname+'.jpeg', dpi=dpi)
        
    # Write data to csv
    if write_data:
        csv_data = {} # { 'name' : [[xdata], [ydata]], ... }
        csv_data['FFT'] = [freq, y_fft]
        csv_data['PSD'] = [freq, y_psd]
        csv_data['amp'] = [freq, y_amp]
        misc_funcs.savedata_to_csv(csv_data, figname+'.csv')
    return


#####################################################################
# Function to optimize wn for butter_lowpass_filter using residuals #
#####################################################################
def butter_optimize_wn_with_FFT(xdata, ydata, wn_method, write_data, savefig, figname, dpi):
    # Compute power spectrum and find wns
    freq, y_fft, y_amp, y_psd, fs, N = compute_FFT(xdata, ydata)
    wns = freq/(0.5*fs)
    
    # Find threshold value if user supplies an approiapte string: wn=op<thres>-p or wn=oa<thres>-p
    threshold = 'mean'; wn_error = False
    if '<' in wn_method and '>' in wn_method:
        value = wn_method.split('<')[-1].split('>')[0].split(',')[0]
        print(value)
        
        try:
            values = [float(i) for i in wn_method.split('<')[-1].split('>')[0].split(',')]
            if len(values) == 1:
                threshold = values[0]
            else: wn_error = True
        except: wn_error = True
        if wn_error:
            print(f'ERROR could not parse wn={wn_method} string, will use mean. Correct syntax: wn=op<thres>-p or wn=oa<thres>-p')

    
    # Find first time psd or amp crosses the threshold
    wn = None; mode = 'psd'
    if wn_method.startswith('op'):
        if threshold == 'mean': threshold = np.mean(y_psd)
        index = np.min(np.where(y_psd <= threshold)[0])
        if index == 0: index = 1 # wn cant be the DC-offset
        wn = wns[index]
        mode = 'psd'
    if wn_method.startswith('oa'):
        if threshold == 'mean': threshold = np.mean(y_amp)
        index = np.min(np.where(y_amp <= threshold)[0])
        if index == 0: index = 1 # wn cant be the DC-offset
        wn = wns[index]
        mode = 'amp'
    if wn is None:
        wn = 0.01
        print(f'ERROR something went wrong during wn optimztion, using default wn={wn}')

    
    # Create new plot
    if wn_method.endswith('-p'):
        plot_fft_and_wn(freq, y_fft, y_amp, y_psd, wns, fs, wn=wn, threshold=threshold, mode=mode, 
                        figname=figname, write_data=write_data, savefig=savefig, dpi=dpi)
    
    return wn


#####################################################################
# Function to optimize wn for butter_lowpass_filter using residuals #
#####################################################################
def butter_optimize_wn_with_residuals(xdata, ydata, order, qm, wn_method, savefig, figname, dpi, delta=0.001):
    # Remove '-p' flag from qm or lrs otherwise a lot of plots will be generated
    if '-p' in str(qm): qm = str(qm).replace('-p', '')
    
    # Set default settings for calling butter_lowpass_filter (since multiple
    # tests will be run, we DO NOT WANT to save data or figures.)
    write_data = False; savefig = False; figname = ''; dpi = 300
    
    # Find "spectra" of wns
    wns = np.arange(2*delta, 1.0, delta)
    summed_residuals = np.zeros_like(wns)
    for n, wn in enumerate(wns):
        ytemp, qm_temp = butter_lowpass_filter(xdata, ydata, wn, order, qm, write_data, savefig, figname, dpi, pflag=False)
        residuals = ydata - ytemp
        summed_residuals[n] = np.sum(residuals)
        
    # Define zero point (could be zero or mean of summed_residuals, etc ...)
    zero_point = np.mean(summed_residuals)
    
    # Compute first zero crossing point  
    zero_crossing_x = None; zero_crossing_y = None
    x_cross, y_cross = misc_funcs.value_crossing(wns, summed_residuals, yvalue=zero_point, style='low')
    if len(x_cross) >= 1:
        zero_crossing_index = np.min(np.where(wns == x_cross[0])[0])
    else: # Set default if there is no zero crossing
        default_wn = 0.01
        zero_crossing_index = np.absolute(wn-default_wn).argmin()
    
    # Compute optimized wn
    decimal_places = len(str(delta).split('.')[-1]) + 1
    zero_crossing_x = wns[zero_crossing_index] + delta/2
    zero_crossing_x = round(zero_crossing_x, decimal_places)
    zero_crossing_y = zero_point #summed_residuals[zero_crossing_index]
    optimized_wn = zero_crossing_x
    
    decimal_places = len(str(delta).split('.')[-1]) + 1
    optimized_wn = round(optimized_wn, decimal_places)
    
    # Create new plot
    if wn_method.endswith('-p'):
        fig_res, ax = plt.subplots()
        plt.plot(wns, summed_residuals,  '.', ms=5, color='tab:orange', label='sum(residuals)')
        plt.axhline(zero_point, ls='--', color='red', lw=2, label='zero point')
        plt.plot(zero_crossing_x, zero_crossing_y,  'o', ms=5, color='tab:blue', label='x={:.4f}; y={:.4f}'.format(zero_crossing_x, zero_crossing_y))
        ax.set_xlabel('W$_n$', fontsize=12)
        ax.set_ylabel('sum(residuals)', fontsize=12)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, ncol=3, fontsize=12)
        if savefig:
            fig_res.savefig(figname+'.jpeg', dpi=dpi)
    return optimized_wn


#####################################################################
# Function to optimize wn for butter_lowpass_filter using residuals #
#####################################################################
def butter_optimize_wn_with_CVE(xdata, ydata, n, wn_method, order, delta, stop, qm='1,1', savefig=False, figname='BW_CVE', dpi=300):
    from tqdm import tqdm
    BWF_write_data, BWF_savefig = False, False
    fold = math.floor(xdata.size/n)
    folds = np.array([int(i*n) for i in range(fold-1)])
    if stop > 1: stop = 1
    if delta < 0: delta = 1/1000
    wns = np.arange(delta, stop, delta)
    cves = np.zeros_like(wns)    
    print('\n\nStarting wn optimization using cross-validation error ...')
    for n, wn in tqdm(enumerate(wns), total=wns.size):
        cve = 0
        for i in folds:
            px = np.delete(xdata, i)
            py = np.delete(ydata, i) 
            ytemp, qm_temp = butter_lowpass_filter(px, py, wn, order, qm, BWF_write_data, BWF_savefig, figname, dpi, pflag=False)
            dy = ydata[i] - ytemp[i]
            cve += dy*dy
        cves[n] = math.sqrt(cve/px.size)
        
    min_index = np.argmin(cves)
    wn = wns[min_index]
    cve = cves[min_index]

    if wn_method.endswith('-p'):
        fig, ax = plt.subplots()
        ax.plot(wns, cves,  '-o', lw=2, ms=5, color='tab:blue', label='')
        ax.plot(wn, cve,  'o', ms=7, color='tab:orange', label='Minimum CVE ({:.4f}, {:.4f})'.format(wn, cve))
        ax.set_xlabel('W$_n$', fontsize=12)
        ax.set_ylabel('CVE', fontsize=12)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, ncol=3, fontsize=12)
        if savefig:
            fig.savefig(figname+'.jpeg', dpi=dpi)
    return wn

############################
### Perform some testing ###
############################
if __name__ == "__main__":  
    # %matplotlib qt 
    # %matplotlib inline
    
    # Generate stress and strain 
    max_strain = 0.2
    max_stress = 150
    strain = np.linspace(0, max_strain, 2000)
    stress = max_stress*(1 - np.exp(-strain/max_strain**2))
    
    # Generate noisy stress-strain data
    np.random.seed(73149)
    noise = np.random.normal(loc=0.0, scale=1.0, size=strain.size)
    mag = [0, 25, 50, 100]
    
    noisy0 = stress + (mag[0]/300)*max_stress*noise
    noisy1 = stress + (mag[1]/300)*max_stress*noise
    noisy2 = stress + (mag[2]/300)*max_stress*noise
    noisy3 = stress + (mag[3]/300)*max_stress*noise
    
    
    # Setup wn string
    wn_method = 'cve<0.001, 0.05, 1>-p'
    wn_error = False
    delta = 0.001
    stop = 1.0
    n = 1
    if '<' in wn_method and '>' in wn_method:
        try:
            values = [float(i) for i in wn_method.split('<')[-1].split('>')[0].split(',')]
            if len(values) == 3:
                delta, stop, n = values
            else: wn_error = True
        except: wn_error = True
        if wn_error:
            print(f'ERROR could not parse wn={wn_method} string. Correct syntax: wn=cve<delta, stop, n>-p')
            print('Using defaults of delta={}, stop={}, n={} for cross-validation calculations'.format(delta, stop, n))
    
    
    # Filter data with settings from above
    order = 2
    wn_noisy2 = butter_optimize_wn_with_CVE(strain, noisy2, n, wn_method, order, delta, stop, qm='1,1', savefig=False, figname='BW_CVE', dpi=300)
    filt2, qm_noisy2 = butter_lowpass_filter(strain, noisy2, wn_noisy2, order, '3,2', False, False, 'figname', 300, pflag=False)
    
    
    # Plot the results
    fig, ax = plt.subplots()
    ax.plot(strain, stress, '-', lw=2, color='tab:blue', label='True signal')
    ax.plot(strain, noisy2, 'o', ms=2, color='tab:cyan', label='{}% Noisy'.format(mag[2]))
    ax.plot(strain, filt2, 'o', ms=2, color='tab:orange', label='{}% Filtered Noisy'.format(mag[2]))
    ax.set_xlabel('Strain', fontsize=12)
    ax.set_ylabel('Stress', fontsize=12)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, ncol=3, fontsize=12)
    
    
    # Testing leave-one-out
    xdata = np.array([i for i in range(5)])
    ydata = xdata**2
    n = 1
    
    fold = math.floor(xdata.size/n)
    folds = np.array([int(i*n) for i in range(fold-1)])
    for i in folds:
        px = np.delete(xdata, i)
        py = np.delete(ydata, i) 
        
        print()
        print(i)
        #print(xdata[i], px[i])
        print(xdata, xdata[i])
        print(px, px[i])
    

