# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
August 25, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

PEKK Tg ~= 196-206 C -> 469-479
"""
##############################
# Import Necessary Libraries #
##############################
import src.log_analysis.misc_funcs as misc_funcs
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import numpy as np


##################################
# Function to find fringe slopes #
##################################
def compute_fringe_slope(xdata, ydata, min_span=0, direction='forward', stats=False):
    # Set direction
    if direction == 'forward':
        xdata = xdata.copy()
        ydata = ydata.copy()
    elif direction == 'reverse':
        xdata = xdata.copy()
        ydata = ydata.copy()
        xdata.reverse()
        ydata.reverse()
    else:
        raise Exception(f'ERROR direction={direction} is not supported. Supported directions are "forward" or "reverse"')
    
    # Generate outputs dictionary (only really need 'slopes' and 'fringe', but other statistics
    # can be computed at 25% slower run time, so keep that in mind when using stats=True).
    outputs = {'sum-residuals-squared':[], # cummulative sum of residuals sqaured
               'mean-squared-residual':[], # cummulative sum of residuals sqaured
               'intercepts-variance':[], # cummulative variance of intercepts
               'slopes-variance':[], # cummulative variance of slopes
               'intercepts':[], # cummulative Y-intercepts
               'r-squared-adj': [], # adjusted r squared values
               'r-squared':[], # cummulative r squared values
               'r-squared-std':[], # rolling standard deviation
               'indexes': [], # cummulative indexes of ydata-strain lists
               'slopes':[], # cummulative slopes
               'fringe':[]} # strain value
    
    # Start the walked linear regression method
    sum_xi = 0; sum_yi = 0; sum_xi_2 = 0; sum_yi_2 = 0; sum_xi_yi = 0; n = 0;
    start = xdata[0]
    for i, (x, y) in enumerate(zip(xdata, ydata)):
        # Compute cummulative linear regression parameters
        sum_xi += x; sum_yi += y; n += 1
        sum_xi_2 += x*x; sum_yi_2 += y*y; 
        sum_xi_yi += x*y
        span = abs(start - x)
        
        # Need at least 2 points to perform linear regression
        if n <= 3: continue
        
        # Only compute outputs if span is in the desired range
        if span >= min_span:
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
                
                # try: 
                #     #r2_std = statistics.stdev(outputs['r-squared'])
                #     r2_std = np.std(np.array(outputs['r-squared']))
                # except: r2_std = 0
                # outputs['r-squared-std'].append(r2_std)
    return outputs

def maximum_peak(rfr):
    # Find peaks
    xdata = np.array(rfr['fringe']); ydata = np.array(rfr['r-squared']);
    peaks, properties = find_peaks(ydata, prominence=None)
    xpeaks = xdata[peaks]; ypeaks = ydata[peaks]
    
    # Find maximum peak index
    maxdex = np.argmax(ypeaks)
    max_x = xpeaks[maxdex]
    max_y = ypeaks[maxdex]
    
    # Related max-peak to rfr-data indexes
    rfr_index = np.argmin(np.abs(xdata - max_x))
    
    rfr['xpeaks'] = list(xpeaks)
    rfr['ypeaks'] = list(ypeaks)
    return rfr, rfr_index


#####################################################
# Function for computing the RFR-thermal properties #
#####################################################
def compute(xdata, ydata, min_span, rtemp, x_units, y_units, figname, dpi):
    #------------------------------------------------------------------------------------------------------------#
    # This entire function assumes strain increases with time (i.e. tensile or shear tests). Therefore if strain #
    # decreases with time for compression or unloading tests, we need to "reverse" the stress and strain lists.  #
    #------------------------------------------------------------------------------------------------------------#
    if xdata[-1] < xdata[0]: xdata.reverse(); ydata.reverse()
    
    
    #-------------------------------------------------#
    # Start computing regression frige response stats #
    #-------------------------------------------------#
    # First forward pass
    ffs = compute_fringe_slope(xdata, ydata, min_span=min_span, direction='forward', stats=True)
    ffs, max_r2_ffs = maximum_peak(ffs)
    
    #max_r2_ffs = ffs['r-squared'].index(max(ffs['r-squared']))
    lo_bounds = sorted([xdata[0], ffs['fringe'][max_r2_ffs]])
    xdata_lo, ydata_lo = misc_funcs.reduce_data(xdata, ydata, min(lo_bounds), max(lo_bounds))
    lo_lr = misc_funcs.linear_regression(xdata_lo, ydata_lo)
    
    # Find the minimum after the maximum
    ymin = min(ffs['r-squared'][max_r2_ffs:])
    mindex = ffs['r-squared'].index(ymin)
    xmin = ffs['fringe'][mindex]

    # First reverse pass (truncated to correct location)
    xlo = max(lo_bounds)  # this pass should not be be allowed to go below the lower linear-region
    xhi = xmin + min_span # need to allow for min_span properly
    rxdata, rydata = misc_funcs.reduce_data(xdata, ydata, xlo, xhi)
    rfs = compute_fringe_slope(rxdata, rydata, min_span=min_span, direction='reverse', stats=True)
    rfs, max_r2_rfs = maximum_peak(rfs)
    
    #max_r2_rfs = rfs['r-squared'].index(max(rfs['r-squared']))
    hi_bounds = sorted([rfs['fringe'][max_r2_rfs], rxdata[-1]])
    xdata_hi, ydata_hi = misc_funcs.reduce_data(xdata, ydata, min(hi_bounds), max(hi_bounds))
    hi_lr = misc_funcs.linear_regression(xdata_hi, ydata_hi)
    
    # Compute intersection
    m1 = lo_lr.b1; b1 = lo_lr.b0
    m2 = hi_lr.b1; b2 = hi_lr.b0
    if m1 != m2:
        x_intersection = (b2 - b1) / (m1 - m2)
        y_intersection = m1*x_intersection + b1
        intersection = [x_intersection, y_intersection]
    else:
        center_index = int(len(xdata)/2)
        intersection = [xdata[center_index], ydata[center_index]]
        print('  ERROR RFR-thermal could not find the intersection between the lo-CTE and the hi-CTE linear regression b/c they are parallel. Setting intersction at the center')
    
    # Set transition region
    transition = [max(lo_bounds), min(hi_bounds)]
    transition_y = [min(ydata) for _ in transition]
    
    # Compute room temp volume if user provides rtemp range
    rx, ry, vol_avg = None, None, None
    if len(rtemp) == 2:
        try: 
            rtemp = tuple([float(i) for i in rtemp])
            digits = True
        except: digits = False
        if digits:
            xlo = min(rtemp)
            xhi = max(rtemp)
            rx, ry = misc_funcs.reduce_data(xdata, ydata, xlo, xhi)
            vol_avg = sum(ry)/len(ry)
    
    # Compute CTE if user provides rtemp range
    lo_cte, hi_cte = None, None
    if vol_avg is not None:
        lo_cte = (1/(3*vol_avg))*lo_lr.b1
        hi_cte = (1/(3*vol_avg))*hi_lr.b1
    
    #-------------------#
    # setup output data #
    #-------------------#
    outputs = {'lo_bounds': lo_bounds,
               'hi_bounds': hi_bounds,
               'lo_line': [lo_lr.b0, lo_lr.b1, lo_lr.r2],
               'hi_line': [hi_lr.b0, hi_lr.b1, hi_lr.r2],
               'tg': intersection,
               'lo_cte': lo_cte,
               'hi_cte': hi_cte,
               'vol_avg': vol_avg,
               'transition': transition,
               'lo_reg_x': list(lo_lr.xreg), 
               'lo_reg_y': list(lo_lr.yreg), 
               'hi_reg_x': list(hi_lr.xreg), 
               'hi_reg_y': list(hi_lr.yreg)}
    
    
    #-----------------------------------------------------#
    # Generate the visualizations for this anaysis method #
    #-----------------------------------------------------#
    # Set a map from known Y-units to Y-values
    y_units_map = {'(A^3)': 'Volume', '(A$^3$)':'Volume', '(g/cc)':'Density', '($g/cm^3$)':'Density'}
    
    # plot the results
    fs = 12
    legend_fs_scale = 0.95
    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))

    
    ax1.plot(xdata, ydata, '.', ms=6, color='tab:grey', label='LAMMPS data')
    ax1.plot(lo_lr.xreg, lo_lr.yreg, marker='|', ls='--', ms=14, lw=4, color='tab:blue', label='Lo-slope ({:.6f})'.format(lo_lr.b1))
    ax1.plot(hi_lr.xreg, hi_lr.yreg, marker='|', ls='--', ms=14, lw=4, color='tab:red', label='Hi-slope ({:.6f})'.format(hi_lr.b1))
    if rx is not None and ry is not None:
        ax1.plot(rx, ry, '.', ms=6, color='cyan', label='{} ({:.2f})'.format('$Vol_{avg}$', vol_avg))
    
    ax1.plot(*intersection, 'p', ms=12, color='tab:purple', label='$T_g$ ({:.6f})'.format(intersection[0]))
    
    label = 'Transition Region\n  - xlo={:.6f}\n  - xhi={:.6f}\n  - span={:.6f}'.format(transition[0], transition[1], max(transition)-min(transition))
    label = 'Transition Region\n  - xlo={:.6f}\n  - xhi={:.6f}'.format(transition[0], transition[1])
    ax1.plot(transition, transition_y, marker='|', ls='-', ms=14, lw=4, color='black', label=label)
    ax1.legend(loc='upper left', bbox_to_anchor=(-0.35, 1.4), fancybox=True, ncol=2, fontsize=legend_fs_scale*fs)
    y_values = y_units_map.get(y_units, 'Y-value')
    ax1.set_xlabel('Temperature {}'.format(x_units), fontsize=fs)
    ax1.set_ylabel('{} {}'.format(y_values, y_units), fontsize=fs)
    ax1.tick_params(axis='both', labelsize=fs)
    
    if lo_cte is not None and hi_cte is not None:
        title = 'Lo-CTE={:.6e}; Hi-CTE={:.6e}'.format(lo_cte, hi_cte)
    else: title = 'Lo-CTE={}; Hi-CTE={}'.format(lo_cte, hi_cte)
    ax1.set_title(title, loc='right', y=1.35, pad=15)
    
    
    ax2.plot(ffs['fringe'], ffs['r-squared'], color='tab:blue', label='Lo-Foward')
    ax2.plot(ffs['xpeaks'], ffs['ypeaks'], '*', ms=8, color='darkblue', label='Lo-peaks')
    ax2.plot(ffs['fringe'][max_r2_ffs], ffs['r-squared'][max_r2_ffs], 'o', ms=10, color='tab:blue', label='max(Foward)')
    ax2.axvline(xmin, color='tab:cyan', ls='--', lw=2, label='min-after-max\n(Forward)')
    ax2.axhline(ymin, color='tab:cyan', ls='--', lw=2)
    ax2.plot(xmin, ymin, 'o', ms=10, color='tab:cyan')
    ax2.legend(loc='upper left', bbox_to_anchor=(-0.025, 1.4), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs)
    ax2.set_xlabel('Temperature {}'.format(x_units), fontsize=fs, color='black')
    ax2.set_ylabel('Forward-fringe $r^2$', fontsize=fs, color='tab:blue')
    ax2.tick_params(axis='x', labelsize=fs, colors='black')
    ax2.tick_params(axis='y', labelsize=fs, colors='tab:blue')
    #ax2.set_yscale('log')
    
    ax3 = ax2.twinx()
    ax3.plot(rfs['fringe'], rfs['r-squared'], color='tab:red', label='Hi-Reverse')
    ax3.plot(rfs['xpeaks'], rfs['ypeaks'], '*', ms=8, color='darkred', label='Hi-peaks')
    
    ax3.plot(rfs['fringe'][max_r2_rfs], rfs['r-squared'][max_r2_rfs], 'o', ms=10, color='tab:red', label='max(Reverse)')
    
    ax3.axvline(max(lo_bounds), color='salmon', ls='--', lw=2, label='max(Forward)\nlower bound')
    ax3.legend(loc='upper left', bbox_to_anchor=(0.525, 1.4), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs)
    ax3.set_ylabel('Reverse-fringe $r^2$', fontsize=fs, color='tab:red')
    ax3.tick_params(axis='x', labelsize=fs, colors='black')
    ax3.tick_params(axis='y', labelsize=fs, colors='tab:red')
    
    fig1.tight_layout()
    fig1.savefig(figname+'_1.jpeg', dpi=dpi)

    return outputs