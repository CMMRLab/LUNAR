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
import matplotlib.pyplot as plt
import numpy as np




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
            y_peak_avg = (lo + hi)/2
            y_valley = between_peaksy[minimum_index]
            avg = y_valley - y_peak_avg
            small = y_valley - min(lo, hi)
            large = y_valley - max(lo, hi)
            valley_depths.append( (avg, small, large) )
    return xpeaks, ypeaks, xvalleys, yvalleys, valley_depths


##################################################
# Function to compute X and Y componets of a FFT #
##################################################
def compute_FFT(x, y):
    # Convert to numpy arrays
    xx = np.array(x)
    yy = np.array(y)
    
    # Define sampling rate and number of data points
    dx = np.mean(np.abs(np.diff(xx)))
    if dx != 0: 
        fs = 1/dx # sampling rate
    else: fs = xx.shape[0]/(np.max(xx) - np.min(xx))
    N = xx.shape[0] # number of data points
    d = 1/fs # sampling space

    # Perform one sided FFT
    fhat = np.fft.rfft(yy, axis=0, norm='backward')
    x_fft = np.fft.rfftfreq(N, d=d)
    y_fft = fhat 
    return x_fft, y_fft, fs, N


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
    x_fft, y_fft, x_psd, y_psd, fs, N = compute_power_spectral_density(xdata, ydata)
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


####################################################
# Function to compute Power Spectrum Density (PSD) #
####################################################
def compute_power_spectral_density(x, y):
    # Compute one sided FFT
    x_fft, y_fft, fs, N = compute_FFT(x, y)
    
    # One sided power spectrum density
    x_psd = x_fft.copy()
    y_psd = np.real( (y_fft*np.conjugate(y_fft))/N )
    return x_fft, y_fft, x_psd, y_psd, fs, N


##############################################
# Function to convert power to decibels (db) #
##############################################
def power_to_db(power, ref_power=1):
  return 10*np.log10(power/ref_power)


##############################################################
# Function to plot fft, power spectrum and selected wn value #
##############################################################
def plot_fft_power_and_wn(x_fft, y_fft, fs, wns, x_psd, y_psd, y_psd_dB, figname, write_data, savefig, dpi, wn_x=None, wn_y=None, cutoff_frequency=None, mean_psd=True):
    fig_psd, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    
    # Power spectrum plot in natural units
    ax1.stem(x_psd, y_psd, linefmt='tab:blue', markerfmt='.', label='$|X(f)|^2/N$')
    if mean_psd:
        y_psd_mean = np.mean(np.array(y_psd))
        ax1.axhline(y_psd_mean, ls='--', color='tab:olive', lw=1, label='mean(Power)={:.8f}'.format(y_psd_mean))
    if wn_x is not None and wn_y is not None and cutoff_frequency is not None:
        ax1.axvline(cutoff_frequency, ls='--', color='tab:red', lw=1, label='Cutoff={:.8f}'.format(cutoff_frequency))
    ax1.set_xlabel('Frequency (1/X-units)', fontsize=12)
    ax1.set_ylabel('Power', fontsize=12)
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), fancybox=True, ncol=1, fontsize=12)
    
    
    # Generate second x-axis for axis1
    def forward_conversion1(x_psd):
        return x_psd/(0.5*fs)
    def reverse_conversion1(x_psd):
        return x_psd
    second_x1 = ax1.secondary_xaxis('top', functions=(forward_conversion1, reverse_conversion1))
    second_x1.set_xlabel('$W_n$', fontsize=12)
    
    # Power spectrum plot in dB
    ax2.stem(wns, y_psd_dB, linefmt='tab:blue', markerfmt='.', label='$|X(f)|^2/N$')
    if wn_x is not None and wn_y is not None and cutoff_frequency is not None:
        ax2.plot(wn_x, wn_y,  'o', ms=5, color='tab:orange', label='$W_n$ (x={:.6f}; y={:.6f})'.format(wn_x, wn_y))
        ax2.axvline(wn_x, ls='--', color='tab:red', lw=1, label='Cutoff={:.8f}'.format(cutoff_frequency))
    ax2.set_xlabel('$W_n$', fontsize=12)
    ax2.set_ylabel('Power (dB - $ref_{power}$ at mean)', fontsize=12)
    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0), fancybox=True, ncol=1, fontsize=12)
    
    # Generate second x-axis for axis2
    def forward_conversion2(wns):
        return (0.5*fs)*wns
    def reverse_conversion2(wns):
        return wns
    second_x2 = ax2.secondary_xaxis('top', functions=(forward_conversion2, reverse_conversion2))
    second_x2.set_xlabel('Frequency (1/X-units)', fontsize=12)
    
    # finalize plot
    fig_psd.tight_layout()
    if savefig:
        fig_psd.savefig(figname+'.jpeg', dpi=dpi)
        
    # Write data to csv
    if write_data:
        csv_data = {} # { 'name' : [[xdata], [ydata]], ... }
        csv_data['FFT'] = [x_fft, y_fft]
        csv_data['PSD-natural'] = [x_psd, y_psd]
        csv_data['PSD-dB-at-mean-freq'] = [x_psd, y_psd_dB]
        csv_data['PSD-dB-at-mean-wns'] = [wns, y_psd_dB]
        misc_funcs.savedata_to_csv(csv_data, figname+'.csv')
    return


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


###############################################################################
# Function to plot FFT, Power spectrum, and Wn; if Wn is supplied by the user #
###############################################################################
def generate_and_plot_fft_power_and_wn(xdata, ydata, wn, write_data, savefig, figname, dpi):
    # Compute power spectrum
    x_fft, y_fft, x_psd, y_psd, fs, N = compute_power_spectral_density(xdata, ydata)
    
    # Compute corrected y_fft
    y_fft = np.abs(np.real(y_fft/(N/2))) 
    
    # Convert power spectrum to dB with the mean power as the reference
    psd_mean = np.mean(y_psd)
    y_psd_dB = power_to_db(y_psd, psd_mean)

    # Find corresponding wns
    wns = x_psd/(0.5*fs)
    
    # Find nearest index for supplied wn, and compute nearest wn_x and wn_y
    nearest_index = np.absolute(wns-wn).argmin()
    wn_x = wns[nearest_index]
    wn_y = y_psd_dB[nearest_index]
    
    # Find cutoff frequency
    cutoff_index = np.min(np.where(wns == wn_x)[0])
    cutoff_frequency = x_fft[cutoff_index]
    
    # Plot the FFT, Power Spectrum, and nearest Wn
    plot_fft_power_and_wn(x_fft, y_fft, fs, wns, x_psd, y_psd, y_psd_dB, figname, write_data, savefig, dpi, wn_x, wn_y, cutoff_frequency, False)
    return


#####################################################################
# Function to optimize wn for butter_lowpass_filter using residuals #
#####################################################################
def butter_optimize_wn_with_power_spectrum(xdata, ydata, wn_method, write_data, savefig, figname, dpi):
    # Compute power spectrum
    x_fft, y_fft, x_psd, y_psd, fs, N = compute_power_spectral_density(xdata, ydata)
    
    # Compute corrected y_fft
    y_fft = np.abs(np.real(y_fft/(N/2))) 
    
    # Convert power spectrum to dB with the mean power as the reference
    psd_mean = np.mean(y_psd)
    y_psd_dB = power_to_db(y_psd, psd_mean)

    # Find corresponding wns
    wns = x_psd/(0.5*fs)
    
    # Define zero point and find first zero crossing point
    wn_x = None; wn_y = None; default_wn = 0.01
    x_cross, y_cross = misc_funcs.value_crossing(x_psd, y_psd_dB, yvalue=0, style='high')
    if len(x_cross) >= 1:
        zero_crossing_index = np.absolute(x_psd-x_cross[0]).argmin()
    else: # Set default if there is no zero crossing
        zero_crossing_index = np.absolute(wns-default_wn).argmin()
    
    # Find optimzed Wn value for the Butterworth filter
    try:
        wn_x = wns[zero_crossing_index]
        wn_y = y_psd_dB[zero_crossing_index]
        optimized_wn = round(wn_x, 6)
        
        # Find cutoff frequency
        cutoff_index = np.min(np.where(wns == wn_x)[0])
        cutoff_frequency = x_fft[cutoff_index]
    except:
        optimized_wn = default_wn 
        cutoff_frequency = 1
    
    # Create new plot
    if wn_method.endswith('-p'):
        plot_fft_power_and_wn(x_fft, y_fft, fs, wns, x_psd, y_psd, y_psd_dB, figname, write_data, savefig, dpi, wn_x, wn_y, cutoff_frequency, True)
        
    # # Try looking at derivatives
    # dxn_psd, dy1_psd, dy2_psd = misc_funcs.compute_derivative(wns, y_psd)
    # ix1_psd, iy1_psd = misc_funcs.compute_integral(wns, y_psd)
    
    # critical_x, critical_y = misc_funcs.value_crossing(dxn_psd, dy1_psd, yvalue=0, style='low')
    # inflection_x, inflection_y = misc_funcs.value_crossing(dxn_psd, dy2_psd, yvalue=0, style='low')
    # ymin = min(dy2_psd)    # Lower bound for veritical lines
    # ymax = max(dy2_psd)    # Lower bound for veritical lines
        
    # fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 9))
    # ax1.plot(x_psd, y_psd)
    # ax1.set_xlabel('X - PSD')
    # ax1.set_ylabel('Y - PSD')        
            
    # ax2.plot(dxn_psd, dy1_psd)
    # ax2.set_xlabel('d1 X - PSD')
    # ax2.set_ylabel('d1 Y - PSD') 
    
    # ax3.plot(dxn_psd, dy2_psd)
    # ax3.set_xlabel('d2 X - PSD')
    # ax3.set_ylabel('d2 Y - PSD') 
    
    # # if critical_x:
    # #     ax2.vlines(critical_x[0], ymin, ymax, color='tab:red', label='Critical Points')
    # # if inflection_x:
    # #     ax3.vlines(inflection_x[0], ymin, ymax, color='tab:red', label='Inflection Points')
    
    
    # ax4.plot(ix1_psd, iy1_psd)
    # ax4.set_xlabel('i1 X - PSD')
    # ax4.set_ylabel('i1 Y - PSD') 
    
    # fig1.tight_layout()
    
    return optimized_wn


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