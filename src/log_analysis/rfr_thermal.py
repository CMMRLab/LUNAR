# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
July 7, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html#scipy.stats.gaussian_kde
https://docs.scipy.org/doc/scipy/tutorial/stats/kernel_density_estimation.html
"""
##############################
# Import Necessary Libraries #
##############################
import src.log_analysis.misc_funcs as misc_funcs
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks, peak_prominences
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
    return outputs


#################################
# Function to order legend grid #
#################################
def order_legend_grid(legend_grid, labels, ncols=2):
    nrows = len(legend_grid)
    desired_order = []
    for col in range(ncols):
        for row in range(nrows):
            item = legend_grid[row][col]
            if item not in labels:
                continue
            if item is not None:
                desired_order.append(item)
    return desired_order

    
##########################################################################
# Run multi pass FBF or BFB search to find optimal line segments.        #
# use IC as initial condition, will default to edge of data              #
# Will work for F, FB, FBF, FBFB... non-alternating raises an exception. #
##########################################################################
def rfr_multi_pass(xdata, ydata, min_span, sequence='FBF', xlo=None, xhi=None, ic=None):
    # Trim data
    if xlo is not None or xhi is not None:
        if xlo is None: xlo = min(xdata)
        if xhi is None: xhi = max(xdata)
        xdata, ydata = misc_funcs.reduce_data(xdata, ydata, xlo, xhi)
        
    # Set starting point (initial condition) for multi-pass
    min_xdata = min(xdata)
    max_xdata = max(xdata)
    if sequence.upper().startswith('F'):
        if ic is None: ic = min_xdata
    elif sequence.upper().startswith('B'):
        if ic is None: ic = max_xdata
    else:
        raise Exception(f'sequence={sequence} not supported, must be FBF or BFB')
        
    # Start multi-pass sequence
    passes = [] # [RFR_dict1, RFR_dict2, ... Npasses]
    prior = None
    for direction in sequence:
        # Ensure direction is alternating
        if prior is not None and prior == direction.upper():
            raise Exception(f'sequence={sequence} does not alternate F and B')
        
        # Forward pass
        if direction.upper() == 'F':
            # Compute forward pass 
            xreduce, yreduce = misc_funcs.reduce_data(xdata, ydata, ic, max_xdata)
            fwd = compute_fringe_slope(xreduce, yreduce, min_span=min_span, direction='forward', stats=True)
            xseries = fwd['fringe']
            yseries = fwd['r-squared']
            if not xseries: continue
            
            # Set minimum and maximum values
            max_index  = yseries.index( max(yseries) )
            fwd['max_r2']     = yseries[max_index]
            fwd['max_fringe'] = xseries[max_index]
            
            # min_index  = yseries.index( min(yseries) )
            # fwd['min_r2']     = yseries[min_index]
            # fwd['min_fringe'] = xseries[min_index]
        
            # Log pass and reset ic
            passes.append(fwd)
            prior = direction
            ic = fwd['max_fringe']
            
        # Reverse pass
        elif direction.upper() == 'B':
            # Compute forward pass 
            xreduce, yreduce = misc_funcs.reduce_data(xdata, ydata, min_xdata, ic)
            rev = compute_fringe_slope(xreduce, yreduce, min_span=min_span, direction='reverse', stats=True)
            xseries = rev['fringe']
            yseries = rev['r-squared']
            if not xseries: continue
            
            # Set minimum and maximum values
            max_index  = yseries.index( max(yseries) )
            rev['max_r2']     = yseries[max_index]
            rev['max_fringe'] = xseries[max_index]
            
            # min_index  = yseries.index( min(yseries) )
            # rev['min_r2']     = yseries[min_index]
            # rev['min_fringe'] = xseries[min_index]
            
            # Log pass and reset ic
            passes.append(rev)
            prior = direction
            ic = rev['max_fringe']

        else:
            raise Exception(f'direction={direction} not supported, must be F or B')
            
    # Set bounds based on number of passes
    if len(passes) > 1:  # Use last two FB or BF
        bounds = sorted([passes[-1]['max_fringe'], passes[-2]['max_fringe']])
    elif len(passes) == 1:
        bounds = sorted([passes[0]['max_fringe'], ic])
    else: # Default to min and max of series
        bounds = sorted([min_xdata, max_xdata])
        
    #  Compute linear regression within bounds
    xdata_linear, ydata_linear = misc_funcs.reduce_data(xdata, ydata, min(bounds), max(bounds))
    regression = misc_funcs.linear_regression(xdata_linear, ydata_linear)
            
    return passes, bounds, regression


##############################################################
# Function to compute intersection of two linear regressions #
##############################################################
def regression_intersection(line1, line2, xdata, ydata):
    # Compute intersection
    m1 = line1.b1; b1 = line1.b0
    m2 = line2.b1; b2 = line2.b0
    if m1 != m2:
        x_intersection = (b2 - b1) / (m1 - m2)
        y_intersection = m1*x_intersection + b1
        intersection = [x_intersection, y_intersection]
        intersects = True
    else:
        center_index = int(len(xdata)/2)
        intersection = [xdata[center_index], ydata[center_index]]
        #print('  ERROR RFR-thermal could not find the intersection between the lo-CTE and the hi-CTE linear regression b/c they are parallel. Setting intersction at the center')
        #print('  X-region: ', line1.xreg)
        intersects = False
    return intersection, intersects


#######################################################################################
# Function to compute the estimate PDF and all of its properties using a Gaussian KDE #
#######################################################################################
def compute_gaussian_kde(xdata, ydata, tgs, slopes, lowers, uppers, widths, bound_weights, distance_from_glassy, ic_increment, glassy_slope, settings):
    xpeaks, ypeaks, peak_info = None, None, None
    tg_kde, temp, kernel_effective_width = None, None, None
    selected_tg, tg_selection, selected_label = None, None, None
    if len(tgs) >= 2 and np.std(tgs, ddof=1) != 0:
        # Set the effective width of the Gaussian kernel used for kernel density
        # estimation (KDE) of the candidate Tg values. The kernel width is defined as
        # the total span containing approximately 99.7% of the Gaussian probability
        # mass (±3σ). This allows the smoothing scale to be specified directly in
        # temperature units (K), rather than by the Gaussian standard deviation.
        #
        # Since candidate Tg values are generated from a discrete IC temperature grid,
        # choosing a kernel width comparable to the IC increment helps suppress
        # artificial multimodal structure while preserving physically meaningful Tg
        # populations.
        #
        # SciPy's gaussian_kde expects bw_method as a dimensionless scaling factor
        # multiplied by the sample standard deviation. Therefore, convert the desired
        # kernel width into the corresponding Gaussian standard deviation (σ =
        # kernel_width / 6), then compute the required bw_method factor.
        tg_std = np.std(tgs, ddof=1)
        
        kernel_effective_width = settings.get('kernel_width', ic_increment)
        kernel_std = kernel_effective_width / 6.0
        bw_method = kernel_std / tg_std
        
        # Generate evalution points for the PDF (with proper bounds
        # w/ max padding of 5% of tgs or +-4*std of Gaussian kernel bandwidth)
        xmin = min(tgs)
        xmax = max(tgs)
        padding = max(0.05*(xmax - xmin), 4*kernel_std) # +- 4*std should cover most area in the KDE
        if padding == 0: padding = ic_increment
        
        # Clamp data to physically meaningful bounds (xdata span) if we overshoot
        xmin_pdf = max(min(xdata), xmin - padding)
        xmax_pdf = min(max(xdata), xmax + padding)
        dx_pdf   = 1/100 # Temperature increment (K)
        temp     = np.arange(xmin_pdf, xmax_pdf + dx_pdf, dx_pdf)

        
        # Compute optional weights for the Gaussian KDE. The user may supply an
        # arbitrary NumPy-compatible mathematical expression through:
        #
        #     weights = "<expression>"
        #
        # where the expression is evaluated element-wise and must return one weight
        # per candidate Tg. The following variables are available:
        #
        #   rw : Normalized rubbery-region width
        #        rw = widths / max(widths)
        #
        #   rs : Rubbery-region slope normalized by the glassy slope
        #        rs = ctes / glassy_slope
        #
        # Supported NumPy functions are available through the np namespace.
        #
        # Example expressions:
        #
        #     weights = "rw"
        #     weights = "rs"
        #     weights = "0.5*rw + 0.5*rs"
        #     weights = "rw**2"
        #     weights = "np.sqrt(rw)"
        #     weights = "np.exp(-rs)"
        weighting = settings.get('weights', None)
        weights = None
        if weighting is not None:
            # Compute possible weighting options, where the following are supported
            #  rw = rubbery width normalized to one (divided by maximum width)
            #  rs = rubbery slope divided by glassy slope
            rw = widths/max(widths)   # Set width-weights via a normalized width of regions
            rs = slopes/glassy_slope  # Set slope-weights via a normalized cte relative to glassy slope
            bw = bound_weights        # Set bound penalty where no-penetration=1 and larger-penetraions approach zero
            
            kw = kernel_effective_width # Set scalar value for kernel width
            gd = distance_from_glassy/max(distance_from_glassy) # Set distance-from-glassy weight normalized by max distances
            
            # Setup the globals namespace to limit scope of what eval() can do
            allowed_builtins = ['min','max','sum','abs','len','map','range','reversed']
            copied_builtins = globals()['__builtins__'].copy()
            globals_dict = {'rw':rw, 'rs':rs, 'kw': kw, 'bw':bw, 'gd':gd}
            globals_dict['__builtins__'] = {key:copied_builtins[key] for key in allowed_builtins}
            globals_dict['np'] = globals()['np']

            # Try computing weights with vectorized numpy and eval
            try:
                weights = np.asarray(eval(weighting, globals_dict), dtype=float)
                weights = np.maximum(weights, 0)   # clamp negatives to zero
                weights = np.minimum(weights, 1e4) # clamp large values to 1e4
                weights = np.nan_to_num(weights, nan=1.0) # set Nan's to 1.0
                
                # If all weights evaluate to zero reset weights to ones (this is
                # useful to shut-off weighting with zeros like "0.0*rw + 0.0*rs",
                # but still keep the syntax for setting weights for testing)
                weight_set = set(weights)
                if len(weight_set) == 1 and sum(weight_set) == 0:
                    weights = np.ones_like(weights)
                    print(weights)
                
                if weights.ndim != 1:
                    raise ValueError("weights expression must return a 1-D array")
                
                if len(weights) != len(tgs):
                    raise ValueError("weights expression returned the wrong length")
                
                if np.any(weights < 0):
                    raise ValueError("weights must be non-negative")
                    
                print('Successfully evaluated weighting compute = "{}"'.format(weighting))
            except:
                print('Failed to evaulate weighting compute = "{}". Using weights=None'.format(weighting))
        
        # Compute KDE and PDF (estimated PDF, but area under curve is ~1.0) for Tg data
        tg_kde_func = gaussian_kde(tgs, bw_method=bw_method, weights=weights)
        tg_kde = tg_kde_func(temp)           
        
        # Require detected peaks to be separated by at least one kernel width.
        # Convert the desired temperature separation into the corresponding number
        # of points on the PDF evaluation grid. A minimum of one point is required
        # by scipy.signal.find_peaks().
        dtemp = temp[1] - temp[0]
        min_peak_separation = kernel_effective_width/2
        distance = max(1, int(np.ceil(min_peak_separation / dtemp)))
        
        # Require peaks to stand-out more then the y-jump between two adjacent points
        prominence = max(np.diff(tg_kde))

        # Find the peaks in the PDF
        peaks, props = find_peaks(tg_kde, prominence=prominence)#, distance=distance)
        xpeaks = temp[peaks]
        ypeaks = tg_kde[peaks]
        if len(peaks) == 0:
            selected_tg = temp[np.argmax(tg_kde)]
            selected_label = 'PDF(max)'
            tg_selection = 'max'
            return temp, tg_kde, tg_selection, xpeaks, ypeaks, peak_info, selected_tg, selected_label, kernel_effective_width
        
        # Find peak densities, prominences, and bases
        peak_density = tg_kde[peaks]
        peak_prom, prom_left_bases, prom_right_bases = peak_prominences(tg_kde, peaks)
        
        # Find peak probability bounds using minima between neighboring peaks.
        # This partitions the PDF into local peak regions instead of using
        # prominence bases, which can sometimes span too far.
        peak_left_bases  = np.zeros(len(peaks), dtype=int)
        peak_right_bases = np.zeros(len(peaks), dtype=int)
        for n, pidx in enumerate(peaks):
            pidx = int(pidx)
        
            # Left bound
            if n == 0:
                lb = 0
            else:
                prev_peak = int(peaks[n - 1])
                local_min = np.argmin(tg_kde[prev_peak:pidx + 1])
                lb = prev_peak + local_min
        
            # Right bound
            if n == len(peaks) - 1:
                rb = len(tg_kde) - 1
            else:
                next_peak = int(peaks[n + 1])
                local_min = np.argmin(tg_kde[pidx:next_peak + 1])
                rb = pidx + local_min
        
            peak_left_bases[n] = lb
            peak_right_bases[n] = rb
        
        # Find probability mass associated with each peak using left/right bases
        peak_prob = np.zeros(len(peaks), dtype=float)
        for n, (lb, rb) in enumerate(zip(peak_left_bases, peak_right_bases)):
            lb = int(lb)
            rb = int(rb)
        
            if rb > lb:
                peak_prob[n] = np.trapz(tg_kde[lb:rb + 1], temp[lb:rb + 1])
            else:
                peak_prob[n] = 0.0
        
        # Peak properties for labeling/selection
        peak_info = []
        for n, (pidx, xt, yd, prom, prob) in enumerate(zip(peaks, xpeaks, peak_density, peak_prom, peak_prob)):
            tmp = {'peak_index': n,
                   'pdf_index': pidx,
                   'temp': xt,
                   'density': yd,
                   'prom': prom,
                   'prob': prob,
                   'left_base': int(peak_left_bases[n]),
                   'right_base': int(peak_right_bases[n]),
                   'selected': False,
                   'label': ''}
            peak_info.append(tmp)
        
        # Tg selection rule syntax where the 'tg' key in settings can have a valued string like:
        #
        # "PDF(max)"
        #     Select the Tg corresponding to the global maximum of the estimated
        #     probability density function (PDF).
        #
        # "PDF(temp, index)"
        #     Select the PDF peak by index after sorting peaks by increasing temperature.
        #     index=0 is the lowest-temperature peak, index=1 is the next lowest, etc.
        #
        # "PDF(density, index)"
        #     Select the PDF peak by index after sorting peaks by decreasing probability
        #     density. index=0 is the highest-density peak, index=1 is the second
        #     highest-density peak, etc.
        #
        # "PDF(prom, index)"
        #     Select the PDF peak by index after sorting peaks by decreasing peak
        #     prominence. index=0 is the most prominent peak, index=1 is the second
        #     most prominent peak, etc.
        #
        # "PDF(prob, index)"
        #     Partition the PDF into Tg populations using the local minima between
        #     neighboring peaks, compute the probability associated with each
        #     population by integrating the PDF over its region, and select the peak
        #     by index after sorting by decreasing probability. index=0 is the most
        #     probable Tg population, index=1 is the second most probable, etc.
        tg_selection = settings.get('tg', 'PDF(max)')
        
        # Strip out un-needed characters
        tg_selection = tg_selection.replace('PDF', '')
        tg_selection = tg_selection.replace('(', '')
        tg_selection = tg_selection.replace(')', '')
        tg_selection = tg_selection.strip()
        
        # Get style and index (if applicable)
        tg_index = None
        if ',' in tg_selection:
            options = tg_selection.split(',')
            if len(options) != 2:
                raise Exception('ERROR Tg selection rule doesnt follow PDF(style, index) syntax')
            tg_selection, tg_index = options
            tg_selection = tg_selection.strip()
            tg_index = tg_index.strip()
            try: 
                tg_index = int(tg_index)
            except: 
                raise Exception('ERROR Tg selection index is not an int.')
                
        # Check that the selection method and index are supported
        if tg_selection not in ['max', 'temp', 'density', 'prom', 'prob']:
            raise Exception(f'ERROR Tg selection style="{tg_selection}" is not supported. Currently supported styles: "max" or "temp" or "density"')
        if tg_selection in ['temp', 'density', 'prom', 'prob'] and tg_index is None:
            raise Exception('ERROR Tg selection style requires an index of a group to select using the syntax "PDF(style, index)"')
                
        # Ensure index is in bounds, if not round to nearest index
        if isinstance(tg_index, int):
            if tg_index >= 0: # Support positve index rounding
                possible_indexes = [i for i in range(len(peaks))]
                diff_indexes     = [abs(tg_index - i) for i in possible_indexes]
                min_diff_index   = diff_indexes.index(min(diff_indexes))
                bounded_index    = possible_indexes[min_diff_index]
                if bounded_index != tg_index:
                    print(f'WARNING - user supplied index={tg_index}, but there are not enough peaks to support this index. Rounding to nearest index={bounded_index}')
                    tg_index = bounded_index
            else: # Support negative index rounding
                possible_indexes = [-(i + 1) for i in range(len(peaks))]
                diff_indexes     = [abs(tg_index - i) for i in possible_indexes]
                min_diff_index   = diff_indexes.index(min(diff_indexes))
                bounded_index    = possible_indexes[min_diff_index]
                if bounded_index != tg_index:
                    print(f'WARNING - user supplied index={tg_index}, but there are not enough peaks to support this index. Rounding to nearest index={bounded_index}')
                    tg_index = bounded_index
        
        # Find Tg based on maximum of PDF
        if tg_selection == 'max':
            selected_tg = temp[np.argmax(tg_kde)]
            
            # Label peaks by density ranking for display
            ranked = sorted(peak_info, key=lambda d: d['density'], reverse=True)
            for rank, info in enumerate(ranked):
                info['label_rank'] = rank
                info['label'] = 'density[{}]\n{:.3g}'.format(rank, info['density'])
            
            selected_peak = min(peak_info, key=lambda d: abs(d['temp'] - selected_tg))
            selected_peak['selected'] = True
            selected_label = 'PDF(max)'
            
        # Sort peaks in order of increasing temp (alreay in this order) and select user defined index
        elif tg_selection == 'temp':
            selected_tg = xpeaks[tg_index]
            
            # Label peaks by density ranking for display
            ranked = sorted(peak_info, key=lambda d: d['temp'])
            for rank, info in enumerate(ranked):
                info['label_rank'] = rank
                info['label'] = 'temp[{}]\n{:.2f}'.format(rank, info['temp'])
            
            selected_peak = ranked[tg_index]
            selected_peak['selected'] = True
            selected_tg = selected_peak['temp']
            selected_label = 'PDF(temp,{})'.format(tg_index)
        
        # Sort peaks in order of decreasing density and select user defined index
        elif tg_selection == 'density':
            ranked = sorted(peak_info, key=lambda d: d['density'], reverse=True)
            
            # Label peaks by density ranking for display
            for rank, info in enumerate(ranked):
                info['label_rank'] = rank
                info['label'] = 'density[{}]\n{:.3g}'.format(rank, info['density'])
            
            # Select the Tg
            selected_peak = ranked[tg_index]
            selected_peak['selected'] = True
            selected_tg = selected_peak['temp']
            selected_label = 'PDF(density,{})'.format(tg_index)
            
        # Sort peaks in order of decreasing prominences and select user defined index
        elif tg_selection == 'prom':
            ranked = sorted(peak_info, key=lambda d: d['prom'], reverse=True)
            
            # Label peaks by density ranking for display
            for rank, info in enumerate(ranked):
                info['label_rank'] = rank
                info['label'] = 'prom[{}]\n{:.3g}'.format(rank, info['prom'])
            
            # Select the Tg
            selected_peak = ranked[tg_index]
            selected_peak['selected'] = True
            selected_tg = selected_peak['temp']
            selected_label = 'PDF(prom,{})'.format(tg_index)
            
        # Sort peaks in order of decreasing probability mass and select user defined index
        elif tg_selection == 'prob':
            ranked = sorted(peak_info, key=lambda d: d['prob'], reverse=True)
        
            # Label peaks by probability mass ranking for display
            for rank, info in enumerate(ranked):
                info['label_rank'] = rank
                info['label'] = 'prob[{}]\n{:.3g}'.format(rank, info['prob'])
        
            # Select the Tg
            selected_peak = ranked[tg_index]
            selected_peak['selected'] = True
            selected_tg = selected_peak['temp']
            selected_label = 'PDF(prob,{})'.format(tg_index)
            
    else:
        if tgs.any():
            selected_tg = tgs[0]
            selected_label = ''
            print('WARNING - there are not enough possible Tg - values to make a Gaussian kernel')
            print('OR all Tg values are the same making the standard deviation zero and not possible to compute.')
            print('Using lowest temp Tg - IC.')
        else:
            selected_tg = (max(xdata) - min(xdata))/2
            selected_label = ''
            print('WARNING - the generation of the PDF failed, setting the selected Tg value at the center of the X-data')
        
    return temp, tg_kde, tg_selection, xpeaks, ypeaks, peak_info, selected_tg, selected_label, kernel_effective_width
    

    
#####################################################
# Function for computing the RFR-thermal properties #
#####################################################
def compute(xdata, ydata, min_span, rtemp, settings, x_label, y_label, x_units, y_units, figname, dpi):
    #------------------------------------------------------------------------------------#
    # This entire function assumes temp increases with time. Therefore if temp decreases # 
    # with time for cooling, we need to "reverse" the temp and volume/density lists.     #
    #------------------------------------------------------------------------------------#
    if xdata[-1] < xdata[0]: xdata.reverse(); ydata.reverse()
    fig1 = None
    
    # Log results to csv struct
    csv_data = {} # { label: [[xdata], [ydata]], ... }
    csv_data['LAMMPS_data'] = [xdata, ydata]


    #---------------------------------#
    # Glassy Forward-backward-forward #
    #---------------------------------#
    # Glassy-forward-first pass (gf1). Starting at lowest temp and walking to the right.
    # Glassy-reverse-first pass (gr1). Starting at max(gf1) and walking to the left.
    # Glassy-forward-first pass (gf2). Starting at max(gr1) and walking to the right.
    passes, glassy_bounds, glassy_lr = rfr_multi_pass(xdata, ydata, min_span, sequence='FBF', xlo=None, xhi=None)
    gf1, gr1, gf2 = passes
    gf1['color'] = 'tab:cyan'
    gr1['color'] = 'tab:blue'
    gf2['color'] = 'darkblue'
    #print( list(gf1.keys()) )
    csv_data['Glassy_Forward_1'] = [gf1['fringe'], gf1['r-squared']]
    csv_data['Glassy_Reverse_1'] = [gr1['fringe'], gr1['r-squared']]
    csv_data['Glassy_Forward_2'] = [gf2['fringe'], gf2['r-squared']]
    

    #------------------------------------------------#
    # Rubbery IC Backward-forward-backward IC checks #
    #------------------------------------------------#
    # Setup IC grid to search over
    ic_increment = min_span/2 # couple ic_increment to min_span parameter
    ic_start     = (max(glassy_bounds) // ic_increment + 1) * ic_increment  # round up on glassy upper, 1st will fail and be ignored
    ic_end       = max(xdata)
    ic_steps     = int( (ic_end - ic_start) // ic_increment + 1 )
    ic_list      = [float(ic_start + n*ic_increment) for n in range(ic_steps)]
    
    # Get any user settings
    tau = settings.get('tau', min_span/4)
    min_tg = settings.get('min_tg', max(glassy_bounds)) # Default to max of glassy bounds if not provided
    
    # Start guess linear sections from different ICs
    xlo = max(glassy_bounds)
    guesses = {} # { float(ic): [ [rfr_pass1, ... rfr_passN], [lo_bound, hi_bound], regression, guessed_tg], ... Nic's }
    lowers, uppers, tgs, ics, slopes, bound_weights, distance_from_glassy, bound_penetrations = [], [], [], [], [], [], [], []
    for ic in ic_list:
        # Compute FBF or BFB from IC
        passes, rubbery_bounds, guessed_lr = rfr_multi_pass(xdata, ydata, min_span, sequence='BFB', ic=ic, xlo=xlo, xhi=None)
        if len(passes) != 3: continue # skip incomplete FBF or BFB searches
        
        # Get the intersection between this IC guess and the glassy region
        tg_guess, intersects = regression_intersection(glassy_lr, guessed_lr, xdata, ydata)
        if not intersects: continue # skip over non-intersecting guesses
        if tg_guess[0] < min(xdata): continue # Tg should be within xdata-span (lower-side)
        if tg_guess[0] > max(xdata): continue # Tg should be within xdata-span (upper-side)
        if tg_guess[0] < min_tg: continue # Skip over too small Tg's based on user inputs 

        # Compute bound penertaion penality
        bound_penetration = 0.0
        #bound_penetration += max(0.0, min(xdata) - tg_guess[0])         # Tg should be within xdata-span (lower-side). Currently handle with hard continues above.
        #bound_penetration += max(0.0, tg_guess[0] - max(xdata))         # Tg should be within xdata-span (upper-side). Currently handle with hard continues above.
        bound_penetration += max(0.0, max(glassy_bounds) - tg_guess[0])  # Tg should not be within the glassy region (upper bound)
        bound_penetration += max(0.0, tg_guess[0] - min(rubbery_bounds)) # Tg should not be within the rubbery region (lower bound)
        
        # Compute bound weight
        bound_weight = np.exp(-bound_penetration / tau)
        
        # Compute distance from glassy
        distance_from_glassy.append( abs(tg_guess[0] - max(glassy_bounds)) )
        
        guesses[ic] = [passes, rubbery_bounds, guessed_lr, tg_guess]
        bound_penetrations.append( bound_penetration )
        bound_weights.append( bound_weight )
        lowers.append( min(rubbery_bounds) )
        uppers.append( max(rubbery_bounds) )
        slopes.append( guessed_lr.b1 )
        tgs.append( tg_guess[0] )
        ics.append( ic )
    
    # Default fail save
    if not guesses:
        print('WARNING - No valid ICs found - setting IC at maximum of glassy bound')
        
        # Compute FBF or BFB from IC
        ic = ( max(glassy_bounds) + min(glassy_bounds) ) / 2
        ic_list.append(ic)
        passes, rubbery_bounds, guessed_lr = rfr_multi_pass(xdata, ydata, min_span, sequence='BFB', ic=ic, xlo=xlo, xhi=None)
        tg_guess, intersects = regression_intersection(glassy_lr, guessed_lr, xdata, ydata)

        # Compute bound penertaion penality
        bound_penetration = 0.0
        #bound_penetration += max(0.0, min(xdata) - tg_guess[0])         # Tg should be within xdata-span (lower-side). Currently handle with hard continues above.
        #bound_penetration += max(0.0, tg_guess[0] - max(xdata))         # Tg should be within xdata-span (upper-side). Currently handle with hard continues above.
        bound_penetration += max(0.0, max(glassy_bounds) - tg_guess[0])  # Tg should not be within the glassy region (upper bound)
        bound_penetration += max(0.0, tg_guess[0] - min(rubbery_bounds)) # Tg should not be within the rubbery region (lower bound)
        
        # Compute bound weight
        bound_weight = np.exp(-bound_penetration / tau)
        
        # Compute distance from glassy
        distance_from_glassy.append( abs(tg_guess[0] - max(glassy_bounds)) )
        
        guesses[ic] = [passes, rubbery_bounds, guessed_lr, tg_guess]
        bound_penetrations.append( bound_penetration )
        bound_weights.append( bound_weight )
        lowers.append( min(rubbery_bounds) )
        uppers.append( max(rubbery_bounds) )
        slopes.append( guessed_lr.b1 )
        tgs.append( tg_guess[0] )
        ics.append( ic )


    #---------------------------------------------------------#
    # Rubbery Backward-forward-backward (Josh's PDF matching) #
    #---------------------------------------------------------#    
    # Convert list to numpy arrays
    distance_from_glassy= np.array( distance_from_glassy )
    bound_weights = np.array( bound_weights )
    uppers   = np.array( uppers )
    lowers   = np.array( lowers )
    widths   = uppers - lowers
    tg_vals  = np.array( tgs )
    slope_vals = np.array( slopes )
    glassy_slope = glassy_lr.b1
    try:
        # Compute the PDF via the KDE
        temp, tg_kde, tg_selection, xpeaks, ypeaks, peak_info, selected_tg, selected_label, kernel_width = compute_gaussian_kde(xdata, ydata, tg_vals, slope_vals, lowers, uppers, widths, bound_weights, distance_from_glassy, ic_increment, glassy_slope, settings)               
        csv_data['Tg_KDE'] = [temp, tg_kde]
        csv_data['Tg_KDE_peaks'] = [xpeaks, ypeaks]
    
    except: # it is unlikely above will fail, but if it does set Tg at the center and warn the user
        temp, tg_kde, tg_selection, xpeaks, ypeaks, peak_info = None, None, None, None, None, None
        selected_tg = (max(xdata) - min(xdata))/2
        selected_label = ''
        print('WARNING - the generation of the PDF failed, setting the selected Tg value at the center of the X-data')
        
    # Go and find closests Tg in guesses (we will update based on widest in group if possible)
    closests_ic = {ic: float( abs(selected_tg - tg_guess[0]) ) for ic, (passes, bounds, guessed_lr, tg_guess) in guesses.items()}
    rubbery_ic = min(closests_ic, key=closests_ic.get) # Find key of lowest value (minimized match of selected Tg)

    # Use closests IC from Tristans group elimation for rubbery domain
    passes, rubbery_bounds, rubbery_lr, intersection = guesses[rubbery_ic]
    rr1, rf1, rr2 = None, None, None
    if len(passes) == 1:
        rr1 = passes
    elif len(passes) == 2:
        rr1, rf1 = passes
    elif len(passes) == 3:
        rr1, rf1, rr2 = passes 

    if rr1 is not None: 
        rr1['color'] = 'tab:orange'
        csv_data['Rubbery_Reverse_1'] = [rr1['fringe'], rr1['r-squared']]
    if rf1 is not None: 
        rf1['color'] = 'tab:red'
        csv_data['Rubbery_Forward_1'] = [rf1['fringe'], rf1['r-squared']]
    if rr2 is not None: 
        rr2['color'] = 'darkred'
        csv_data['Rubbery_Reverse_2'] = [rr2['fringe'], rr2['r-squared']]
    
    csv_data['Rubbery_Guessed_ICs']     = ics
    csv_data['Rubbery_Possible_Lowers'] = lowers
    csv_data['Rubbery_Possible_Uppers'] = uppers
    csv_data['Rubbery_Possible_Slopes'] = slopes
    csv_data['Rubbery_Possible_Tgs']    = tgs
    csv_data['Rubbery_Tg_penetrations'] = bound_penetrations
            

    
    #-------------------------------#
    # Thermal property calculations #
    #-------------------------------#    
    # Set transition region
    transition = [max(glassy_bounds), min(rubbery_bounds)]
    transition_y = [min(ydata) for _ in transition]
    
    # Compute reference temp volume if user provides rtemp range. Possible rtemp inputs
    #  rtemp=xlo,xhi;   # V_room = avg(V(xlo, xhi)); # Volume evaluated by averaging the MD-data within xlo and xhi. Same volume used for glassy and rubbery.
    #  rtemp=tg;        # V_tg   = V(tg);            # Volume evaluated at Tg. Same volume used for glassy and rubbery.
    #  rtemp=tg,tol;    # V_tg   = V(tg +- tol);     # Volume evaluated at Tg +- a tolerance. Different volume used for glassy and rubbery.
    #  rtemp=mid;       # V_mid  = (V2 - V1))/2      # Volume evaluated at mid-point of each region's linear regression. Different volume used for glassy and rubbery.
    #  rtemp=avg;       # V_avg  = avg(V(T1, T2))    # Volume evaluated by averaging the MD-data within each region. Different volume used for glassy and rubbery.
    #  rtemp=log(mean); # V_log  = (V2 - V1) / ( ln(V2/V1) ); # Most accurate average over entire region. Different volume used for glassy and rubbery.
    rx, ry, glassy_volume, rubbery_volume = None, None, None, None # Defaults for rtemp=xlo,xhi
    glassy_temp, rubbery_temp = None, None  # Defaults for rtemp=tg, tol
    
    # rtemp=xlo,xhi; because both indexes converted to digits. Average volume in range of xlo, xhi
    if len(rtemp) == 2:
        try: 
            rtemp = tuple([float(i) for i in rtemp])
            digits = True
        except: digits = False
        if digits: # Only operate on if both indexes are digits
            xlo = min(rtemp)
            xhi = max(rtemp)
            rx, ry = misc_funcs.reduce_data(xdata, ydata, xlo, xhi)
            if rx and ry:
                vol_avg = sum(ry)/len(ry)
                glassy_volume  = vol_avg
                rubbery_volume = vol_avg
            else: rx, ry = None, None
            
    # All other rtemp, will be strings
    if isinstance(rtemp, str):
        rtemp = rtemp.split(',')
        
        # rtemp=tg,tol; where tol is +- about Tg for glassy and rubbery CTE
        if 'tg' in rtemp and len(rtemp) >= 2:
            try:
                tol = float(rtemp[1])
                digit = True
            except: digit = False
            if digit: # Only operate on if index 1 is a digit
                # Each linear regression is: y = b0 + b1*x; b0=intercept and b1=slope
                glassy_temp  = intersection[0] - tol
                rubbery_temp = intersection[0] + tol
                
                glassy_volume  = glassy_lr.b0 + glassy_lr.b1*glassy_temp
                rubbery_volume = rubbery_lr.b0 + rubbery_lr.b1*rubbery_temp
                #print(glassy_temp, intersection[0], rubbery_temp)
                
    
        import sys
        sys.exit()
    
    # Compute CTE if user provides rtemp range
    glassy_cte, rubbery_cte = None, None
    if glassy_volume is not None:
        glassy_cte  = (1/(3*glassy_volume))*glassy_lr.b1
    if rubbery_volume is not None:
        rubbery_cte = (1/(3*rubbery_volume))*rubbery_lr.b1
    
    
    #-------------------#
    # setup output data #
    #-------------------#
    outputs = {}
    outputs = {'glassy_bounds':  glassy_bounds,
               'rubbery_bounds': rubbery_bounds,
               'glassy_line':  [glassy_lr.b0,  glassy_lr.b1,  glassy_lr.r2],
               'rubbery_line': [rubbery_lr.b0, rubbery_lr.b1, rubbery_lr.r2],
               'tg': intersection,
               'glassy_cte':  glassy_cte,
               'rubbery_cte': rubbery_cte,
               
               # 'glassy_slope':   glassy_lr.b1,
               # 'rubbery_slope':  rubbery_lr.b1,
               'glassy_volume':  glassy_volume,
               'rubbery_volume': rubbery_volume,
               
               
               'transition':  transition,
               'glassy_reg_x':  list(glassy_lr.xreg), 
               'glassy_reg_y':  list(glassy_lr.yreg), 
               'rubbery_reg_x': list(rubbery_lr.xreg), 
               'rubbery_reg_y': list(rubbery_lr.yreg)}
    
    
    #-----------------------------------------------------#
    # Generate the visualizations for this anaysis method #
    #-----------------------------------------------------#
    # Set a map from known Y-units to Y-values
    y_units_map = {'(A^3)': 'Volume', '(A$^3$)':'Volume', '(g/cc)':'Density', '($g/cm^3$)':'Density',
                   '(A)': 'Length'}
    
    # If these are empty set defaults
    if not x_label:
        x_label = 'Temperature'
    if not y_label:
        y_label = y_units_map.get(y_units, 'Y-value')
    
    # Ensure this is the largest dx-span in data to reduce number of points in plots for speed of rendering
    plot_dx = 0.5
    
    # plot the results
    lw = 2
    ms = 6
    fs = 10
    legend_fs_scale = 0.8
    fig1 = plt.figure(figsize=(3*5, 9), constrained_layout=True)
    gs = fig1.add_gridspec(nrows=3, ncols=3, height_ratios=[0.4, 1.0, 1.0], hspace=0.15)
    
    # Top legend axes
    leg_ax1 = fig1.add_subplot(gs[0, 0])
    leg_ax2 = fig1.add_subplot(gs[0, 1])
    leg_ax3 = fig1.add_subplot(gs[0, 2])
    
    # Middle plotting axes
    ax1 = fig1.add_subplot(gs[1, 0])
    ax2 = fig1.add_subplot(gs[1, 1])
    ax3 = fig1.add_subplot(gs[1, 2])
    
    # Bottom plotting axes
    ax4 = fig1.add_subplot(gs[2, 0])
    ax5 = fig1.add_subplot(gs[2, 1])
    ax6 = fig1.add_subplot(gs[2, 2])
    
    for lax in [leg_ax1, leg_ax2, leg_ax3]:
        lax.axis('off')
        

    # LAMMPS data / things we care about
    #-----------------------------------
    xplot, yplot = reduce_xy_for_plot(xdata, ydata, plot_dx)
    ax1.plot(xplot, yplot, '.', ms=ms, color='tab:grey', label='LAMMPS data')
    ax1.plot(glassy_lr.xreg,  glassy_lr.yreg,  marker='|', ls='--', ms=14, lw=4, color='tab:blue', label='Glassy-slope ({:.6f})'.format(glassy_lr.b1))
    ax1.plot(rubbery_lr.xreg, rubbery_lr.yreg, marker='|', ls='--', ms=14, lw=4, color='tab:red',  label='Rubbery-slope ({:.6f})'.format(rubbery_lr.b1))
    csv_data['Glassy_line']  = [glassy_lr.xreg,  glassy_lr.yreg]
    csv_data['Rubbery_line'] = [rubbery_lr.xreg, rubbery_lr.yreg]
    csv_data['Tg'] = [ [intersection[0]], [intersection[1]] ]
    if rx is not None and ry is not None:
        ax1.plot(rx, ry, '.', ms=6, color='cyan', label='{} ({:.2f})'.format('$Vol_{avg}$', vol_avg))
    
    ax1.plot(*intersection, 'p', ms=12, color='tab:purple', label='$T_g$ ({:.6f})'.format(intersection[0]))
    
    label = 'Transition Region\n  - xlo={:.6f}\n  - xhi={:.6f}\n  - span={:.6f}'.format(transition[0], transition[1], max(transition)-min(transition))
    label = 'Transition Region\n  - xlo={:.6f}\n  - xhi={:.6f}'.format(transition[0], transition[1])
    ax1.plot(transition, transition_y, marker='|', ls='-', ms=14, lw=4, color='black', label=label)

    ax1.set_xlabel('{} {}'.format(x_label, x_units), fontsize=fs)
    ax1.set_ylabel('{} {}'.format(y_label, y_units), fontsize=fs)
    ax1.tick_params(axis='both', labelsize=fs)
    
    
    handles, labels = ax1.get_legend_handles_labels()
    leg_ax1.legend(handles, labels, loc='lower center', fancybox=True, ncol=2, fontsize=legend_fs_scale*fs)

    if glassy_cte is not None and rubbery_cte is not None:
        title = 'Glassy-CTE={:.6e}\nRubbery-CTE={:.6e}'.format(glassy_cte, rubbery_cte)
    else: title = 'Glassy-CTE={}\nRubbery-CTE={}'.format(glassy_cte, rubbery_cte)
    ax1.set_title(title, loc='left', x=0.05, y=0.8, pad=15)

    
    # Glassy RFR-thermal
    #-------------------
    ax2.plot(gf1['fringe'], gf1['r-squared'], '-', lw=lw, color=gf1['color'], label='1st pass (forward)')
    ax2.plot(gr1['fringe'], gr1['r-squared'], '-', lw=lw, color=gr1['color'], label='2nd pass (reverse)')
    ax2.plot(gf2['fringe'], gf2['r-squared'], '-', lw=lw, color=gf2['color'], label='3rd pass (forward)')
    
    ax2.plot(gf1['max_fringe'], gf1['max_r2'], '<', ms=ms, mec=gf1['color'], mfc='white', mew=1.25*lw, label='max(1st pass)')
    ax2.plot(gr1['max_fringe'], gr1['max_r2'], '>', ms=ms, mec=gr1['color'], mfc='white', mew=1.25*lw, label='max(2nd pass)')
    ax2.plot(gf2['max_fringe'], gf2['max_r2'], '<', ms=ms, mec=gf2['color'], mfc='white', mew=1.25*lw, label='max(3rd pass)')
    
    ax2.set_xlabel('{} {}'.format(x_label, x_units), fontsize=fs, color='black')
    ax2.set_ylabel('Glassy $R^2$', fontsize=fs, color='black')
    ax2.tick_params(axis='x', labelsize=fs, colors='black')
    ax2.tick_params(axis='y', labelsize=fs, colors='black')
    
    handles, labels = ax2.get_legend_handles_labels()
    legend_grid = (['1st pass (forward)', 'max(1st pass)'],
                   ['2nd pass (reverse)', 'max(2nd pass)'],
                   ['3rd pass (forward)', 'max(3rd pass)'])
    desired_order = order_legend_grid(legend_grid, labels, ncols=2)
    lookup = dict(zip(labels, handles))
    leg_ax2.legend([lookup[l] for l in desired_order], desired_order, loc='lower center',
                    fancybox=True, ncol=2, fontsize=legend_fs_scale*fs)
    
    # Rubbery RFR-thermal
    #--------------------
    if rr1 is not None: ax3.plot(rr1['fringe'], rr1['r-squared'], '-', lw=lw, color=rr1['color'], label='1st pass (reverse)')
    if rf1 is not None: ax3.plot(rf1['fringe'], rf1['r-squared'], '-', lw=lw, color=rf1['color'], label='2nd pass (forward)')
    if rr2 is not None: ax3.plot(rr2['fringe'], rr2['r-squared'], '-', lw=lw, color=rr2['color'], label='3rd pass (reverse)')
    
    if rr1 is not None: ax3.plot(rr1['max_fringe'], rr1['max_r2'], '>', ms=ms, mec=rr1['color'], mfc='white', mew=1.25*lw, label='max(1st pass)')
    if rf1 is not None: ax3.plot(rf1['max_fringe'], rf1['max_r2'], '<', ms=ms, mec=rf1['color'], mfc='white', mew=1.25*lw, label='max(2nd pass)')
    if rr2 is not None: ax3.plot(rr2['max_fringe'], rr2['max_r2'], '>', ms=ms, mec=rr2['color'], mfc='white', mew=1.25*lw, label='max(3rd pass)')
    
    ax3.set_xlabel('{} {}'.format(x_label, x_units), fontsize=fs, color='black')
    ax3.set_ylabel('Rubbery $R^2$', fontsize=fs, color='black')
    ax3.tick_params(axis='x', labelsize=fs, colors='black')
    ax3.tick_params(axis='y', labelsize=fs, colors='black')
    
    handles, labels = ax3.get_legend_handles_labels()
    legend_grid = (['1st pass (reverse)', 'max(1st pass)'],
                   ['2nd pass (forward)', 'max(2nd pass)'],
                   ['3rd pass (reverse)', 'max(3rd pass)'])
    desired_order = order_legend_grid(legend_grid, labels, ncols=2)
    lookup = dict(zip(labels, handles))
    leg_ax3.legend([lookup[l] for l in desired_order], desired_order, loc='lower center',
                    fancybox=True, ncol=2, fontsize=legend_fs_scale*fs)
    
    # IC guessing and predictions
    #----------------------------
    ax4.plot(ics, lowers, '-v', lw=lw, ms=ms, color='tab:blue',   label='Lower bound')
    ax4.plot(ics, uppers, '-^', lw=lw, ms=ms, color='tab:red',    label='Upper bound')
    ax4.plot(ics, tgs,    '-p', lw=lw, ms=ms, color='tab:purple', label='Tg')
    
    ax4.legend(loc='upper left', bbox_to_anchor=(0, 1), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs)
    ax4.set_xlabel('Guessed IC Temperature {}'.format(x_units), fontsize=fs, color='black')
    ax4.set_ylabel('Possible Rubbery states {}'.format(x_units), fontsize=fs, color='black')
    ax4.tick_params(axis='x', labelsize=fs, colors='black')
    ax4.tick_params(axis='y', labelsize=fs, colors='black')
    
    twin_color = 'tab:cyan'
    ax4_twin = ax4.twinx()
    ax4_twin.plot(ics, slopes, '-d', lw=lw, ms=ms, color=twin_color)
    ax4_twin.set_ylabel('Possible rubbery slopes ({}/{})'.format(y_units.replace('(', '').replace(')', ''), x_units.replace('(', '').replace(')', '')), fontsize=fs, color=twin_color)
    ax4_twin.tick_params(axis='y', labelsize=fs, colors=twin_color)
    
    # Probability distribution of Tg guesses
    #---------------------------------------   
    if tg_kde is not None:
        # Resample the KDE onto a uniform plotting grid that is smaller in dx to reduce size of plot
        temp_plot = np.arange(temp[0], temp[-1] + plot_dx, plot_dx)
        tg_kde_plot = np.interp(temp_plot, temp, tg_kde)
        
        ax5.plot(temp_plot, tg_kde_plot, lw=2, color='tab:purple', label='Convergence Kernel Density')
        ax5.fill_between(temp_plot, tg_kde_plot, alpha=0.25, color='tab:purple')
    if xpeaks is not None:
        ax5.plot(xpeaks, ypeaks, '*', ms=2*ms, mfc='tab:purple', mec='black', mew=1.2, label='PDF - Peaks')
        
        # Annotate peaks with selection criterion labels
        ymax = max(tg_kde)
        xrange = max(temp) - min(temp)
        ytext_hi = 1.40*ymax
        ytext_lo = 1.20*ymax
        max_labeled_peak = 3
        for n, info in enumerate(peak_info):
            label_rank = info.get('label_rank', n)
            if label_rank > max_labeled_peak and not info['selected']:
                continue
            
            xt = info['temp']
            yd = info['density']

            # Alternate small left/right x-offsets
            dx = (-1)**n * 0.02*xrange
            xtxt = xt + dx

            # Alternate between two label heights
            ytext = ytext_hi if n % 2 == 0 else ytext_lo

            if info['selected']:
                bbox = dict(boxstyle='round,pad=0.25', fc='white', ec='black', lw=1.2)
                arrow_color = 'black'
                text = '{}\n{}'.format(selected_label, info['label'])
                weight = 'bold'
                csv_data['selected_peak'] = [[xt], [yd]]
            else:
                bbox = dict(boxstyle='round,pad=0.25', fc='white', ec='black', lw=0.8, alpha=0.9)
                arrow_color = '0.35'
                text = info['label']
                weight = 'normal'

            arrowprops = dict(arrowstyle='->', color=arrow_color, lw=1.0, shrinkA=0, shrinkB=4)
            ax5.annotate(text, xy=(xt, yd), xytext=(xtxt, ytext), textcoords='data', ha='center', va='bottom',
                         fontsize=0.65*fs, fontweight=weight, bbox=bbox, arrowprops=arrowprops, zorder=10)
            
            
        # Plot left/right bases used for peak probability integration
        if tg_selection == 'prob':
            for n, info in enumerate(peak_info):
                lb = info['left_base']
                rb = info['right_base']
    
                xlb = temp[lb]
                xrb = temp[rb]
                ylb = tg_kde[lb]
                yrb = tg_kde[rb]
                key = 'Tg_KDE_base_{}'.format(n)
                csv_data[key] = [[xlb, xrb], [ylb, yrb]]
    
                if info['selected']:
                    base_color = 'black'
                    base_lw = 1.2
                    base_alpha = 0.8
                    shade_alpha = 0.55
                else:
                    base_color = '0.35'
                    base_lw = 0.8
                    base_alpha = 0.45
                    shade_alpha = 0.03
    
                ax5.axvline(xlb, color=base_color, ls=':', lw=base_lw, alpha=base_alpha, zorder=3)
                ax5.axvline(xrb, color=base_color, ls=':', lw=base_lw, alpha=base_alpha, zorder=3)
    
                ax5.fill_between(temp[lb:rb+1], tg_kde[lb:rb+1], alpha=shade_alpha, color='tab:purple', zorder=2)
                ax5.plot(xlb, ylb, '|', ms=1.8*ms, color=base_color, alpha=base_alpha, zorder=6)
                ax5.plot(xrb, yrb, '|', ms=1.8*ms, color=base_color, alpha=base_alpha, zorder=6)

        ax5.set_ylim(top=2.25*ymax)
    
    ax5.axvline(selected_tg, color='black', ls='--', lw=lw, label='Selected Tg = {:.2f}'.format(selected_tg))
    ax5.set_xlabel('{} {}'.format(x_label, x_units), fontsize=fs)
    ax5.set_ylabel('Estimated Probability Density (1/{})'.format(x_units.replace('(', '').replace(')', '')))
    ax5.legend(loc='upper right', bbox_to_anchor=(1, 1), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs)
        
    # Guessed linear regressions at each IC
    #--------------------------------------
    ax6.plot(xplot, yplot, '.', ms=ms, color='tab:grey', alpha=0.35, label='LAMMPS data')
    ax6.plot(glassy_lr.xreg,  glassy_lr.yreg,  marker='|', ls='--', ms=14, lw=4, color='tab:blue', label='Glassy region')
    
    norm = mpl.colors.Normalize(vmin=min(ic_list), vmax=max(ic_list))
    cmap = plt.cm.coolwarm  # blue -> red
    for ic, (_, _, guess_lr, tg_guess) in guesses.items():
        color = cmap(norm(ic))
        ax6.plot(guess_lr.xreg, guess_lr.yreg,  marker='|', ls='-', ms=1.5*ms, lw=lw, color=color, alpha=1.0)
        ax6.plot(*tg_guess, 'p', ms=1.5*ms, color=color)

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = fig1.colorbar(sm, ax=ax6)
    cbar.set_label('Guessed IC Temperature {}'.format(x_units), fontsize=fs)
    cbar.ax.tick_params(labelsize=fs)
    ax6.set_xlabel('{} {}'.format(x_label, x_units), fontsize=fs, color='black')
    ax6.set_ylabel('{} {}'.format(y_label, y_units), fontsize=fs, color='black')
    ax6.tick_params(axis='both', labelsize=fs, colors='black')
    #ax6.set_title('Guess LR by IC', fontsize=fs)
    ax6.legend(loc='upper left', bbox_to_anchor=(0, 1), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs)

    
    #fig1.tight_layout()
    #fig1.tight_layout(rect=[0, 0, 1, 0.95])
    if figname:
        fig1.set_constrained_layout_pads(hspace=0.05, wspace=0.05)
        fig1.savefig(figname+'_1.jpeg', dpi=dpi)
        
    # Write data to csv file
    write_data = True
    if write_data:
        csvname = '{}.csv'.format(figname)
        print(f'\n\nWriting {csvname} ...')
        misc_funcs.savedata_to_csv(csv_data, csvname)
        
    # Generate publication figures
    pub_figures = False
    #pub_figures = True
    if pub_figures:
        publication_figures(outputs, csv_data, figname, ic_list, guesses, min_span, peak_info, settings, ic_increment)
        

    return outputs, fig1


###############################################
# Function for generating publication figures #
###############################################
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

def interpolate_constant_dx(x, y, dx):
    """
    Linearly interpolate paired x-y data onto a uniformly spaced x-grid.

    Parameters
    ----------
    x, y : array-like
        One-dimensional input data.
    dx : float
        Desired positive x-spacing.

    Returns
    -------
    x_uniform, y_uniform : numpy.ndarray
        Interpolated data with constant x-spacing.

    Notes
    -----
    - Nonfinite data are removed.
    - Data are sorted by x.
    - Repeated x-values are combined by averaging their y-values.
    - Extrapolation is not performed.
    - The final original x-value is included only when it lies on the
      constant-dx grid.
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()

    if x.size != y.size:
        raise ValueError("x and y must have the same number of values.")

    if not np.isfinite(dx) or dx <= 0:
        raise ValueError("dx must be a positive finite number.")

    # Remove NaN and infinite values.
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if x.size < 2:
        raise ValueError("At least two finite data points are required.")

    # Sort by increasing x.
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    # Average y-values associated with duplicate x-values.
    x_unique, inverse = np.unique(x, return_inverse=True)

    y_sums = np.bincount(inverse, weights=y)
    y_counts = np.bincount(inverse)
    y_unique = y_sums/y_counts

    if x_unique.size < 2:
        raise ValueError("At least two unique x-values are required.")

    # Construct an exactly uniform grid without extrapolation.
    span = x_unique[-1] - x_unique[0]
    ratio = span/dx
    tolerance = 10*np.finfo(float).eps*max(1.0, abs(ratio))
    number_of_steps = int(np.floor(ratio + tolerance))

    x_uniform = (
        x_unique[0]
        + dx*np.arange(number_of_steps + 1)
    )

    # Linear interpolation.
    y_uniform = np.interp(
        x_uniform,
        x_unique,
        y_unique,
    )

    return x_uniform, y_uniform

# Function to annotate a curve
def annotate_curve(ax, xdata, ydata, text, point_fraction, text_offset,
                   color='black', fontsize=10, ha='left'):
    """
    Add a label with a leader line pointing to a curve.

    point_fraction
        Position along the curve's x-range where the line should point:
        0 is the minimum x-value and 1 is the maximum x-value.
    text_offset
        Label offset from the selected curve point in display points: (x, y).
    """
    xdata = np.asarray(xdata, dtype=float).ravel()
    ydata = np.asarray(ydata, dtype=float).ravel()
    
    # Interpolate to at least dx=1
    if xdata.size >= 2:
        dx = 1
        xdata, ydata = interpolate_constant_dx(xdata, ydata, dx)

    # Keep only paired, finite coordinates.
    npoints = min(xdata.size, ydata.size)
    xdata = xdata[:npoints]
    ydata = ydata[:npoints]
    valid = np.isfinite(xdata) & np.isfinite(ydata)
    xdata = xdata[valid]
    ydata = ydata[valid]
    if xdata.size == 0:
        return

    point_fraction = np.clip(point_fraction, 0.0, 1.0)
    target_x = xdata.min() + point_fraction*(xdata.max() - xdata.min())
    index = np.argmin(np.abs(xdata - target_x))

    ax.annotate(
        text,
        xy=(xdata[index], ydata[index]),
        xycoords='data',
        xytext=text_offset,
        textcoords='offset points',
        ha=ha,
        va='center',
        color=color,
        fontsize=fontsize,
        zorder=10,
        bbox=dict(boxstyle='round,pad=0.2',
                  facecolor='white', edgecolor=color, alpha=0.8),
        arrowprops=dict(arrowstyle='-',
                        color=color, lw=1.25,
                        shrinkA=3, shrinkB=3),
    )
    
    return None

def find_intersection(xpoints1, ypoints1, xpoints2, ypoints2, tolerance=1e-12):
    """
    Find the intersection of two lines and construct a three-point
    V-shaped line for plotting.

    Returns
    -------
    intersection : tuple
        (x, y) coordinates of the intersection.
    x_v, y_v : lists
        Three-point lists with the intersection as the middle point.
    """
    x1, x2 = xpoints1
    y1, y2 = ypoints1
    x3, x4 = xpoints2
    y3, y4 = ypoints2

    rx, ry = x2 - x1, y2 - y1
    sx, sy = x4 - x3, y4 - y3

    denominator = rx * sy - ry * sx

    if np.isclose(denominator, 0.0, atol=tolerance):
        print("The lines are parallel or coincident.")
        return None, None, None

    t = ((x3 - x1) * sy - (y3 - y1) * sx) / denominator

    x_intersection = x1 + t * rx
    y_intersection = y1 + t * ry

    # Use the leftmost endpoint from either line for the left branch
    endpoints = [
        (x1, y1),
        (x2, y2),
        (x3, y3),
        (x4, y4),
    ]

    left_point = min(endpoints, key=lambda point: point[0])
    right_point = max(endpoints, key=lambda point: point[0])

    vx = [
        left_point[0],
        x_intersection,
        right_point[0],
    ]

    vy = [
        left_point[1],
        y_intersection,
        right_point[1],
    ]

    intersection = (x_intersection, y_intersection)
    vx = (x2, x_intersection, x3)
    vy = (y2, y_intersection, y3)

    return intersection, vx, vy

# Function to Generate publication figures
def publication_figures(outputs, csv_data, figname, ic_list, guesses, min_span, peak_info, settings, ic_increment):
    from matplotlib.lines import Line2D
    from matplotlib.legend_handler import HandlerTuple
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter
    from matplotlib.ticker import ScalarFormatter
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.patches import ConnectionPatch
    from matplotlib.transforms import blended_transform_factory
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    plt.rcParams.update({
        # Main text
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7,

        # Plot text
        'axes.labelsize': 7,
        'axes.titlesize': 7,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'legend.fontsize': 6,

        # Mathematical text
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold',

        # Keep fonts editable when exported
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
    })

    def add_moving_bound_arrows(ax, fixed_x, moving_x, y_positions=(0.20, 0.50, 0.80), length_fraction=0.15,
                                color='#00A9C6', lw=3, mutation_scale=18, zorder=0.5):
    
        # Tail lies between the fixed and moving bounds;
        # arrowhead lands on the moving boundary.
        # x_tail = moving_x + length_fraction*(fixed_x - moving_x)
        # x_head = moving_x
        
        
        # Constant arrow width as a fraction of the x-axis range.
        x_min, x_max = ax.get_xlim()
        arrow_width = length_fraction * abs(x_max - x_min)
        
        # Point from the fixed boundary toward the moving boundary.
        direction = np.sign(moving_x - fixed_x)
        if direction == 0:
            return
        
        x_head = moving_x
        x_tail = moving_x - direction * arrow_width

    
        # x uses data coordinates; y uses axes-fraction coordinates.
        transform = blended_transform_factory(ax.transData, ax.transAxes)
        for y_position in y_positions:
            ax.annotate(
                '',
                xy=(x_head, y_position),
                xytext=(x_tail, y_position),
                xycoords=transform,
                textcoords=transform,
                arrowprops={
                    'arrowstyle': '-|>',
                    'color': color,
                    'lw': lw,
                    'mutation_scale': mutation_scale,
                    'shrinkA': 0,
                    'shrinkB': 0,
                },
                annotation_clip=False,
                zorder=zorder,  # Behind plotted lines and markers
            )
        return
    
    #-------------------------#
    # General figure settings #
    #-------------------------#
    lw = 6
    ms = 12
    fs = 12
    dpi = 1000
    dim = 4
    kde_color         = '#89C653' #'tab:purple'
    kde_color         = '#AA3377'  # muted magenta
    #kde_color         = '#6A3D9A'  # darker purple
    
    regression_color  = '#ff9d3aff'
    fixed_bound_color  = 'black'
    moving_bound_color = '#bbbbbbff'
    legend_fs_scale  = 0.75
    axis_thickness   = 3
    temp_label   = 'Temperature'
    volume_label = 'Volume'
    
    print('\n\n\nCSV-data:')
    for key, value in csv_data.items():
        print(key, len(value))
    print('\n\n\n')
        
    #---------------#
    # Setup regions #
    #---------------#
    # Normal regions
    glassy_line_x,  glassy_line_y  = csv_data['Glassy_line']
    rubbery_line_x, rubbery_line_y = csv_data['Rubbery_line']  
    
    # https://www.color-hex.com/color-palettes/
    #marroonish   = '#640017'
    #dark_red     = '#B2182B'
    #cyanish_blue = '#41b6c4ff'
    #light_blue   = '#2c7fb8ff'
    #medium_blue  = '#253494ff'
    regions = {'Glassy region':           [min(glassy_line_x),  max(glassy_line_x),  '#0097A7'],
               'Glass transition region': [max(glassy_line_x),  min(rubbery_line_x), '#773786'],
               'Rubbery region':          [min(rubbery_line_x), max(rubbery_line_x), '#76001c'],#'#640017'],
              }
    
    # Beyound normal regions
    xpeaks, ypeaks  = csv_data['Tg_KDE_peaks']
    selected_peak_x = csv_data['selected_peak'][0][0]
    selected_peak_y = csv_data['selected_peak'][1][0] 
    closests_index  = np.nanargmin(np.abs(xpeaks - selected_peak_x))
    try:
        # Go and find closests Tg in guesses (we will update based on widest in group if possible)
        selected_transition = xpeaks[closests_index+1]
        selected_transition = xpeaks[-1]
        closests_ic  = {ic: float( abs(selected_transition - tg_guess[0]) ) for ic, (passes, bounds, guessed_lr, tg_guess) in guesses.items()}
        softening_ic = min(closests_ic, key=closests_ic.get) # Find key of lowest value (minimized match of selected transition)
        
        passes, softening_bounds, rubbery_lr, intersection = guesses[softening_ic]
        regions['High-temperature transition region'] = [max(rubbery_line_x),   min(softening_bounds),  '#777777']#'slategray']
        regions['High-temperature region']            = [min(softening_bounds), max(softening_bounds),  '#333333']#'dimgray']
    except: pass
    
    print('\n\n\nRegions:')
    for key, value in regions.items():
        print(key, value)
    print('\n\n\n')
    

    
    #-------------------#
    # Figure3: Namining #
    #-------------------#
    fs = 12
    temperature, volume = csv_data['LAMMPS_data']
    regionized_temp_volume = regionize(temperature, volume, regions)
    
    xlo = min(temperature)
    xhi = max(temperature)
    ylo = min(volume)
    yhi = max(volume)
    print('TEMPERATURE RANGE: ', xlo, xhi)
    
    dx  = xhi - xlo
    dy  = yhi - ylo
    sf  = 0.2
    
    xlo -= sf*dx
    xhi += sf*dx
    ylo -= sf*dy
    yhi += sf*dy
    
    fig, ax = plt.subplots(figsize=(dim, dim))
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        ax.plot(x, y, '-', lw=lw, color=color, label=region)
        
        location = (-40, 40)
        if 'High' in region:
            location = (-90, 35)
            
        
        region = ' \n'.join(region.split())
        annotate_curve(ax, x, y, region, 0.30, location,
                       color=color, fontsize=legend_fs_scale*fs)
        
    rfr_tg, rfr_vx, rfr_vy = find_intersection(glassy_line_x,  glassy_line_y, rubbery_line_x, rubbery_line_y, tolerance=1e-12)
    
    ax.plot(glassy_line_x,  glassy_line_y,  '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression')
    ax.plot(rubbery_line_x, rubbery_line_y, '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression')
    ax.plot(rfr_vx, rfr_vy, '--', lw=lw/6, color=regression_color)
    ax.plot(*rfr_tg, 'p', ms=ms, mfc=regression_color, mec='black', mew=1.2, label="$T_g$")
    
    annotate_curve(ax, (rfr_tg[0]), (rfr_tg[1]), "$T_g$", 0.30, (20, -20),
                   color=regression_color, fontsize=legend_fs_scale*fs)
    
    annotate_curve(ax, glassy_line_x,  glassy_line_y, "Glassy\nlinear\nregression", 0.50, (30, -20),
                   color=regression_color, fontsize=legend_fs_scale*fs)
    
    annotate_curve(ax, rubbery_line_x, rubbery_line_y, "Rubbery\nlinear\nregression", 0.50, (20, -20),
                   color=regression_color, fontsize=legend_fs_scale*fs)
        
    ax.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    ax.set_xlim( (xlo, xhi) )
    ax.set_ylim( (ylo, yhi) )
    
    ax.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax, axis_thickness)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig(figname+'_published_1.jpeg', dpi=dpi)
    
    
    #---------------------#
    # Figure4: RFR-method #
    #---------------------#
    # Take an RFR-pass
    passes, bounds, regression = rfr_multi_pass(temperature, volume, min_span, sequence='F', xlo=None, xhi=None)
    forward_pass = passes[0]
    fringe   = forward_pass['fringe']
    rsquared = forward_pass['r-squared']
    regionized_fringe_r2 = regionize(fringe, rsquared, regions)

    legend_fs_scale  = 0.71
    fs = 23
    zlo = min(rsquared)
    zhi = max(rsquared)
    dz  = zhi - zlo
    sf  = 0.1
    zlo -= sf*dz
    zhi += sf*dz
    fringe_point_label = r'Fringe $R^2$ at moving bound'
    y_positions = (0.10, 0.30, 0.50, 0.70, 0.90)
    y_positions = (0.20, 0.40, 0.60, 0.80)
    
    # Generate figure layout's
    nrows, ncols = 2, 4
    
    left   = 0.065
    right  = 0.995
    bottom = 0.075
    top    = 0.82
    
    fig_width = 4 * dim
    
    # Make each subplot box physically square.
    fig_height = (
        fig_width
        * (right - left)
        * nrows
        / ((top - bottom) * ncols)
    )
    
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(fig_width, fig_height),
        gridspec_kw={'wspace': 0, 'hspace': 0},
    )
    
    ((ax1, ax2, ax3, ax4),
     (ax5, ax6, ax7, ax8)) = axes
    
    fig.subplots_adjust(
        left=left,
        right=right,
        bottom=bottom,
        top=top,
        wspace=0,
        hspace=0,
    )
    
    # Explicitly preserve square subplot boxes.
    for ax in axes.flat:
        ax.set_box_aspect(1)
    
    # Segment 1
    dtemp = max(temperature) - min(temperature)
    segment_1 = [min(temperature), min(temperature) + 0.25*dtemp]
    segment_1_x, segment_1_y = misc_funcs.reduce_data(temperature, volume, min(segment_1), max(segment_1))
    segment_1_lr = misc_funcs.linear_regression(segment_1_x, segment_1_y)
    segment_1_fsx = misc_funcs.closest(fringe, max(segment_1))
    segment_1_fsy = rsquared[fringe.index(segment_1_fsx)]
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        ax1.plot(x, y, '-', lw=lw, color=color, label=region)

    ax1.axvline(min(segment_1), color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax1.axvline(max(segment_1), color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax1.plot(segment_1_lr.xreg, segment_1_lr.yreg, '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression')
    annotate_curve(ax1, segment_1_lr.xreg, segment_1_lr.yreg, 'Linear regression', 1.0, (40, -10), 
                   color=regression_color, fontsize=legend_fs_scale*fs)

    # Fringe r^2
    for region in regionized_fringe_r2:
        x, y, color = regionized_fringe_r2[region]
        x, y = misc_funcs.reduce_data(x, y, min(segment_1), max(segment_1))
        ax5.plot(x, y, '-', lw=lw, color=color, label=region)

    ax5.axvline(min(segment_1), color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax5.axvline(max(segment_1), color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax5.plot(segment_1_fsx, segment_1_fsy, 'o', ms=ms, mfc=regression_color, mec='black', label=fringe_point_label)
    annotate_curve(ax5, [segment_1_fsx], [segment_1_fsy], 'Fringe $R^2$', 1.0, (40, -10), 
                   color=regression_color, fontsize=legend_fs_scale*fs)
    
    #ax1.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax1.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    
    ax5.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax5.set_ylabel('Fringe $R^2$', fontsize=fs, color='black')
    
    ax1.set_xlim( (xlo, xhi) )
    ax1.set_ylim( (ylo, yhi) )
    ax5.set_xlim( (xlo, xhi) )
    ax5.set_ylim( (zlo, zhi) )
    
    set_axis_thickness(ax1, axis_thickness)
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    ax5.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax5, axis_thickness)
    ax5.set_xticks([])
    ax5.set_yticks([0.9, 1.0])
    add_moving_bound_arrows(ax1, fixed_x=min(segment_1), moving_x=max(segment_1), y_positions=y_positions, 
                            color=moving_bound_color, lw=6, mutation_scale=22, length_fraction=0.1, zorder=5)
    add_moving_bound_arrows(ax5, fixed_x=min(segment_1), moving_x=max(segment_1), y_positions=y_positions, 
                            color=moving_bound_color, lw=6, mutation_scale=22, length_fraction=0.1, zorder=5)
    
    # Segment 2
    segment_2 = [min(temperature), min(temperature) + 0.50*dtemp]
    segment_2_x, segment_2_y = misc_funcs.reduce_data(temperature, volume, min(segment_2), max(segment_2))
    segment_2_lr = misc_funcs.linear_regression(segment_2_x, segment_2_y)
    segment_2_fsx = misc_funcs.closest(fringe, max(segment_2))
    segment_2_fsy = rsquared[fringe.index(segment_2_fsx)]
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        ax2.plot(x, y, '-', lw=lw, color=color, label=region)
    
    ax2.axvline(min(segment_2), color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax2.axvline(max(segment_2), color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax2.plot(segment_2_lr.xreg, segment_2_lr.yreg, '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression')
    
    # Fringe r^2
    for region in regionized_fringe_r2:
        x, y, color = regionized_fringe_r2[region]
        x, y = misc_funcs.reduce_data(x, y, min(segment_2), max(segment_2))
        ax6.plot(x, y, '-', lw=lw, color=color, label=region)
    
    ax6.axvline(min(segment_2), color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax6.axvline(max(segment_2), color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax6.plot(segment_2_fsx, segment_2_fsy, 'o', ms=ms, mfc=regression_color, mec='black', label=fringe_point_label)
    
    #ax2.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    #ax2.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    
    ax6.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    #ax6.set_ylabel('Fringe $R^2$', fontsize=fs, color='black')
    
    ax2.set_xlim( (xlo, xhi) )
    ax2.set_ylim( (ylo, yhi) )
    ax6.set_xlim( (xlo, xhi) )
    ax6.set_ylim( (zlo, zhi) )
    
    set_axis_thickness(ax2, axis_thickness)
    ax2.set_xticks([])
    ax2.set_yticks([])

    
    set_axis_thickness(ax6, axis_thickness)
    ax6.set_xticks([])
    ax6.set_yticks([])
    add_moving_bound_arrows(ax2, fixed_x=min(segment_2), moving_x=max(segment_2), y_positions=y_positions, 
                            color=moving_bound_color, lw=6, mutation_scale=22, length_fraction=0.1, zorder=5)
    add_moving_bound_arrows(ax6, fixed_x=min(segment_2), moving_x=max(segment_2), y_positions=y_positions, 
                            color=moving_bound_color, lw=6, mutation_scale=22, length_fraction=0.1, zorder=5)
    
    # Segment 3
    segment_3 = [min(temperature), min(temperature) + 0.75*dtemp]
    segment_3_x, segment_3_y = misc_funcs.reduce_data(temperature, volume, min(segment_3), max(segment_3))
    segment_3_lr = misc_funcs.linear_regression(segment_3_x, segment_3_y)
    segment_3_fsx = misc_funcs.closest(fringe, max(segment_3))
    segment_3_fsy = rsquared[fringe.index(segment_3_fsx)]
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        ax3.plot(x, y, '-', lw=lw, color=color, label=region)
    
    ax3.axvline(min(segment_3), color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax3.axvline(max(segment_3), color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax3.plot(segment_3_lr.xreg, segment_3_lr.yreg, '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression')
    
    # Fringe r^2
    for region in regionized_fringe_r2:
        x, y, color = regionized_fringe_r2[region]
        x, y = misc_funcs.reduce_data(x, y, min(segment_3), max(segment_3))
        ax7.plot(x, y, '-', lw=lw, color=color, label=region)
    
    ax7.axvline(min(segment_3), color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax7.axvline(max(segment_3), color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax7.plot(segment_3_fsx, segment_3_fsy, 'o', ms=ms, mfc=regression_color, mec='black', label=fringe_point_label)
    
    #ax3.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    #ax3.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    
    ax7.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    #ax7.set_ylabel('Fringe $R^2$', fontsize=fs, color='black')
    
    ax3.set_xlim( (xlo, xhi) )
    ax3.set_ylim( (ylo, yhi) )
    ax7.set_xlim( (xlo, xhi) )
    ax7.set_ylim( (zlo, zhi) )
    
    set_axis_thickness(ax3, axis_thickness)
    ax3.set_xticks([])
    ax3.set_yticks([])
    
    set_axis_thickness(ax7, axis_thickness)
    ax7.set_xticks([])
    ax7.set_yticks([])
    add_moving_bound_arrows(ax3, fixed_x=min(segment_3), moving_x=max(segment_3), y_positions=y_positions, 
                            color=moving_bound_color, lw=6, mutation_scale=22, length_fraction=0.1, zorder=5)
    add_moving_bound_arrows(ax7, fixed_x=min(segment_3), moving_x=max(segment_3), y_positions=y_positions, 
                            color=moving_bound_color, lw=6, mutation_scale=22, length_fraction=0.1, zorder=5)
    
    # Segment 4
    segment_4 = [min(temperature), min(temperature) + 1.00*dtemp]
    segment_4_x, segment_4_y = misc_funcs.reduce_data(temperature, volume, min(segment_4), max(segment_4))
    segment_4_lr = misc_funcs.linear_regression(segment_4_x, segment_4_y)
    segment_4_fsx = misc_funcs.closest(fringe, max(segment_4))
    segment_4_fsy = rsquared[fringe.index(segment_4_fsx)]
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        ax4.plot(x, y, '-', lw=lw, color=color, label=region)
    
    ax4.axvline(min(segment_4), color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax4.axvline(max(segment_4), color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax4.plot(segment_4_lr.xreg, segment_4_lr.yreg, '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression')
    
    # Fringe r^2
    for region in regionized_fringe_r2:
        x, y, color = regionized_fringe_r2[region]
        x, y = misc_funcs.reduce_data(x, y, min(segment_4), max(segment_4))
        ax8.plot(x, y, '-', lw=lw, color=color, label=region)
    
    ax8.axvline(min(segment_4), color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax8.axvline(max(segment_4), color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax8.plot(segment_4_fsx, segment_4_fsy, 'o', ms=ms, mfc=regression_color, mec='black', label=fringe_point_label)
    
    #ax4.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    #ax4.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    
    ax8.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    #ax8.set_ylabel('Fringe $R^2$', fontsize=fs, color='black')
    
    ax4.set_xlim( (xlo, xhi) )
    ax4.set_ylim( (ylo, yhi) )
    ax8.set_xlim( (xlo, xhi) )
    ax8.set_ylim( (zlo, zhi) )
    
    set_axis_thickness(ax4, axis_thickness)
    ax4.set_xticks([])
    ax4.set_yticks([])
    ax4.set_xticks([])
    ax4.set_yticks([])
    
    set_axis_thickness(ax8, axis_thickness)
    ax8.set_xticks([])
    ax8.set_yticks([])
    ax8.set_xticks([])
    ax8.set_yticks([])
    add_moving_bound_arrows(ax4, fixed_x=min(segment_4), moving_x=max(segment_4), y_positions=y_positions, 
                            color=moving_bound_color, lw=6, mutation_scale=22, length_fraction=0.1, zorder=5)
    add_moving_bound_arrows(ax8, fixed_x=min(segment_4), moving_x=max(segment_4), y_positions=y_positions, 
                            color=moving_bound_color, lw=6, mutation_scale=22, length_fraction=0.1, zorder=5)
    
    # Add legend and annotations
    all_axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

    # Collect each legend entry only once.
    handles_by_label = {}
    for axis in all_axes:
        axis_handles, axis_labels = axis.get_legend_handles_labels()
        for handle, label in zip(axis_handles, axis_labels):
            if label and not label.startswith('_'):
                handles_by_label.setdefault(label, handle)
    
    # Control the legend order.
    region_keys = list(regions.keys())
    other_keys = ['Fixed bound', 'Moving bound', 'Linear regression', fringe_point_label]
    legend_order = (region_keys + other_keys)
    
    legend_handles = [
        handles_by_label[label]
        for label in legend_order
        if label in handles_by_label
    ]
    
    legend_labels = [
        label
        for label in legend_order
        if label in handles_by_label
    ]
    
    # Rename entries in the displayed legend, if desired.
    display_names = {
        'Fixed bound': 'Fixed bound for linear regression',
        'Moving bound': 'Moving bound for linear regression',
        'Linear regression': 'Bounded linear regression',
    }
    
    displayed_labels = [
        display_names.get(label, label)
        for label in legend_labels
    ]
    
    # Panel boundaries in figure coordinates.
    grid_left = ax1.get_position().x0
    grid_right = ax4.get_position().x1
    grid_top = ax1.get_position().y1
    grid_width = grid_right - grid_left
    
    # Choose enough columns to produce two legend rows.
    legend_ncols = (len(legend_handles) + 1)//2
    legend_ncols = 4
    
    figure_legend = fig.legend(
        legend_handles,
        displayed_labels,
        loc='upper left',
        bbox_to_anchor=(grid_left, 0.89, grid_width, 0.09),
        bbox_transform=fig.transFigure,
        ncol=legend_ncols,
        mode='expand',
        frameon=True,
        fancybox=True,
        framealpha=1.0,
        edgecolor='black',
        borderaxespad=0,
        borderpad=0.6,
        columnspacing=1.2,
        labelspacing=0.35,
        handlelength=1.8,
        handletextpad=0.4,
        fontsize=legend_fs_scale*fs,
    )
    figure_legend.get_frame().set_linewidth(1.5)
    
    # Heading immediately above the panels.
    heading_y = grid_top + 0.032
    
    heading = fig.text(
        grid_left,
        heading_y,
        'Fringe Response Calculation',
        ha='left',
        va='center',
        fontsize=1.15*fs,
        fontweight='bold',
        color='black',
    )
    
    # Draw once so the heading's actual width can be measured.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    
    heading_bbox = (
        heading.get_window_extent(renderer=renderer)
        .transformed(fig.transFigure.inverted())
    )
    
    # Start the arrow just after the heading.
    arrow_start = heading_bbox.x1 + 0.007
    
    heading_arrow = FancyArrowPatch(
        (arrow_start, heading_y),
        (grid_right, heading_y),
        transform=fig.transFigure,
        arrowstyle='-|>',
        mutation_scale=25,
        linewidth=4,
        color='black',
        capstyle='round',
        joinstyle='round',
        clip_on=False,
        zorder=20,
    )
    fig.add_artist(heading_arrow)
    
    # Add Lin example points (https://doi.org/10.1038/s41524-021-00647-w)
    # Splitting via minimum and selecting maximum on either side
    # temperature, volume
    # fringe, rsquared
    use_lin = False
    lin_scale = 0.7
    lin_color = 'crimson'
    fringe   = np.array( fringe )
    rsquared = np.array( rsquared ) 
    mindex = np.argmin(rsquared) # Location of minimum R^2
    zerdex = 0                   # Minimum possible index (zerdex like zero)
    maxdex = int(len(rsquared))  # Maximum possible index

    if use_lin:
        # To the left of the minimum R^2
        lo_fringe = fringe[zerdex:mindex]
        lo_rsquared = rsquared[zerdex:mindex]
        lo_maxdex   = np.argmax(lo_rsquared)
        
        # To the right of the minimum R^2
        hi_fringe   = fringe[mindex:maxdex]
        hi_rsquared = rsquared[mindex:maxdex]
        hi_maxdex   = np.argmax(hi_rsquared)
        
        ax8.plot(lo_fringe[lo_maxdex], lo_rsquared[lo_maxdex], 'o', ms=ms, mfc=lin_color, mec='black', zorder=0.1)
        annotate_curve(ax8, lo_fringe[lo_maxdex], lo_rsquared[lo_maxdex], 'Lin et al. glassy start\ngoing to the right', 1.0, (-10, -30), 
                       color=lin_color, fontsize=legend_fs_scale*fs*lin_scale)
    
        ax8.plot(hi_fringe[hi_maxdex], hi_rsquared[hi_maxdex], 'o', ms=ms, mfc=lin_color, mec='black', zorder=0.1)
        annotate_curve(ax8, hi_fringe[hi_maxdex], hi_rsquared[hi_maxdex], 'Lin et al. rubbery start\ngoing to the right', 1.0, (-100, 30), 
                       color=lin_color, fontsize=legend_fs_scale*fs*lin_scale)
    
        # Lin et al. example on glassy end 
        window  = 200
        lin_xlo = lo_fringe[lo_maxdex]
        lin_xhi = lin_xlo + window
        lin_x, lin_y = misc_funcs.reduce_data( list(temperature), list(volume), lin_xlo, lin_xhi)
        lin_lr = misc_funcs.linear_regression(lin_x, lin_y)
        ax4.plot(lin_lr.xreg, lin_lr.yreg, '--', lw=lw/2, color=lin_color, zorder=8)
        annotate_curve(ax4, lin_lr.xreg, lin_lr.yreg, 'Lin et al.\nglassy fit (a={})'.format(window), 0.5, (-50, -50), 
                       color=lin_color, fontsize=legend_fs_scale*fs*lin_scale)
        
        # Lin et al. example on rubbery end
        lin_xlo = hi_fringe[hi_maxdex]
        lin_xhi = lin_xlo + window
        lin_x, lin_y = misc_funcs.reduce_data( list(temperature), list(volume), lin_xlo, lin_xhi)
        lin_lr = misc_funcs.linear_regression(lin_x, lin_y)
        ax4.plot(lin_lr.xreg, lin_lr.yreg, '--', lw=lw/2, color=lin_color, zorder=8)
        annotate_curve(ax4, lin_lr.xreg, lin_lr.yreg, 'Lin et al.\nrubbery fit (a={})'.format(window), 0.5, (-150, 0), 
                       color=lin_color, fontsize=legend_fs_scale*fs*lin_scale)
    
    fig.savefig(figname+'_published_2.jpeg', dpi=dpi)
    
    #-------------------------------------#
    # "Regionized" FRF plot -> Glassy FBF #
    #-------------------------------------#
    # Generate a BFB Pass
    ic = min(temperature)
    passes, rubbery_bounds, guessed_lr = rfr_multi_pass(temperature, volume, min_span, sequence='FBF', ic=ic, xlo=None, xhi=None)
    gf1, gr1, gf2 = passes
    gf1, gr1, gf2 = passes
    
    gf1_fringe, gf1_r2 = [gf1['fringe'], gf1['r-squared']]
    gr1_fringe, gr1_r2 = [gr1['fringe'], gr1['r-squared']]
    gf2_fringe, gf2_r2 = [gf2['fringe'], gf2['r-squared']]
    
    zlo = min( [min(gf1_r2), min(gr1_r2), min(gf2_r2)] )
    zhi = max( [max(gf1_r2), max(gr1_r2), max(gf2_r2)] )
    dz  = zhi - zlo
    sf  = 0.1
    zlo -= sf*dz
    zhi += sf*dz
    
    
    fs = 14
    y_positions = (0.10, 0.30, 0.50, 0.70, 0.90)
    y_positions = (0.20, 0.40, 0.60, 0.80)
    fig, axes = plt.subplots(
        3, 2,
        figsize=(2*dim, 3.3*dim),  # Extra height for row headings
        gridspec_kw={'wspace': 0.40, 'hspace': 0.30})
    
    ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = axes
    
    # Forward: 1
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        ax1.plot(x, y, '-', lw=lw, color=color, label=region)
        
    ax1.axvline(ic, color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax1.axvline(gf1['max_fringe'], color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    
    segment_1 = [ic, gf1['max_fringe']]
    segment_1_x, segment_1_y = misc_funcs.reduce_data(temperature, volume, min(segment_1), max(segment_1))
    segment_1_lr = misc_funcs.linear_regression(segment_1_x, segment_1_y)
    ax1.plot(segment_1_lr.xreg, segment_1_lr.yreg, '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression within bounds')
        
    regionized_fringe_r2 = regionize(gf1_fringe, gf1_r2, regions)
    for region in regionized_fringe_r2:
        x, y, color = regionized_fringe_r2[region]
        ax2.plot(x, y, '-', lw=lw, color=color, label=region)
        
    ax2.axvline(ic, color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax2.axvline(gf1['max_fringe'], color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax2.plot(gf1['max_fringe'], gf1['max_r2'], 'o', ms=ms, mec='black', mfc=regression_color, mew=1, label='Point of maximum $R^2$')
    
    max_r2_x = gf1['max_fringe']
    if max_r2_x == 387:# TODO: JDK override for cleaner publication plots
        max_r2_x = 375
    annotate_curve(ax2, (gf1['max_fringe']), (gf1['max_r2']), '{:.1f}'.format(max_r2_x), 1.0, (50, 0), 
                   color=regression_color, fontsize=legend_fs_scale*fs)
    
    ax1.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax1.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    ax1.set_xlim( (xlo, xhi) )
    ax1.set_ylim( (ylo, yhi) )
    set_axis_thickness(ax1, axis_thickness)
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    add_moving_bound_arrows(ax1, fixed_x=ic, moving_x=gf1['max_fringe'], y_positions=y_positions, 
                            color=moving_bound_color, lw=3, mutation_scale=18, length_fraction=0.075, zorder=5)
    
    ax2.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax2.set_ylabel('Fringe $R^2$', fontsize=fs, color='black')
    ax2.set_xlim( (xlo, xhi) )
    ax2.set_ylim( (zlo, zhi) )
    
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax2, axis_thickness)
    ax2.set_xticks([])
    ax2.set_yticks([0.9, 1.0])
    
    add_moving_bound_arrows(ax2, fixed_x=ic, moving_x=gf1['max_fringe'], y_positions=y_positions,
                            color=moving_bound_color, lw=3, mutation_scale=18, length_fraction=0.075, zorder=5)
    
    # Update ic
    ic = max(segment_1)
    
    
    # Backward: 1
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        ax3.plot(x, y, '-', lw=lw, color=color, label=region)
        
    ax3.axvline(ic, color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax3.axvline(gr1['max_fringe'], color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    
    segment_2 = [ic, gr1['max_fringe']]
    segment_2_x, segment_2_y = misc_funcs.reduce_data(temperature, volume, min(segment_2), max(segment_2))
    segment_2_lr = misc_funcs.linear_regression(segment_2_x, segment_2_y)
    ax3.plot(segment_2_lr.xreg, segment_2_lr.yreg, '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression within bounds')
    
    regionized_fringe_r2 = regionize(gr1_fringe, gr1_r2, regions)
    for region in regionized_fringe_r2:
        x, y, color = regionized_fringe_r2[region]
        ax4.plot(x, y, '-', lw=lw, color=color, label=region)
        
    ax4.axvline(ic, color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax4.axvline(gr1['max_fringe'], color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax4.plot(gr1['max_fringe'], gr1['max_r2'], 'o', ms=ms, mec='black', mfc=regression_color, mew=1, label='Point of maximum $R^2$')
    annotate_curve(ax4, (gr1['max_fringe']), (gr1['max_r2']), '{:.1f}'.format(gr1['max_fringe']), 1.0, (70, 0), 
                   color=regression_color, fontsize=legend_fs_scale*fs)
    
    ax3.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax3.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    ax3.set_xlim( (xlo, xhi) )
    ax3.set_ylim( (ylo, yhi) )
    set_axis_thickness(ax3, axis_thickness)
    ax3.set_xticks([])
    ax3.set_yticks([])
    
    add_moving_bound_arrows(ax3, fixed_x=ic, moving_x=gr1['max_fringe'], y_positions=y_positions, 
                            color=moving_bound_color, lw=3, mutation_scale=18, length_fraction=0.075, zorder=5)
    
    ax4.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax4.set_ylabel('Fringe $R^2$', fontsize=fs, color='black')
    ax4.set_xlim( (xlo, xhi) )
    ax4.set_ylim( (zlo, zhi) )
    
    ax4.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax4, axis_thickness)
    ax4.set_xticks([])
    ax4.set_yticks([0.9, 1.0])
    
    add_moving_bound_arrows(ax4, fixed_x=ic, moving_x=gr1['max_fringe'], y_positions=y_positions, 
                            color=moving_bound_color, lw=3, mutation_scale=18, length_fraction=0.075, zorder=5)
    
    # Update ic
    ic = min(segment_2)
    
    # Forward: 2
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        ax5.plot(x, y, '-', lw=lw, color=color, label=region)
        
    ax5.axvline(ic, color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax5.axvline(gf2['max_fringe'], color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    
    segment_3 = [ic, gf2['max_fringe']]
    segment_3_x, segment_3_y = misc_funcs.reduce_data(temperature, volume, min(segment_3), max(segment_3))
    segment_3_lr = misc_funcs.linear_regression(segment_3_x, segment_3_y)
    ax5.plot(segment_3_lr.xreg, segment_3_lr.yreg, '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression within bounds')
    
    regionized_fringe_r2 = regionize(gf2_fringe, gf2_r2, regions)
    for region in regionized_fringe_r2:
        x, y, color = regionized_fringe_r2[region]
        ax6.plot(x, y, '-', lw=lw, color=color, label=region)
        
    ax6.axvline(ic, color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax6.axvline(gf2['max_fringe'], color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax6.plot(gf2['max_fringe'], gf2['max_r2'], 'o', ms=ms, mec='black', mfc=regression_color, mew=1, label='Point of maximum $R^2$')
    
    max_r2_x = gf2['max_fringe']
    if max_r2_x == 387:# TODO: JDK override for cleaner publication plots
        max_r2_x = 375
    annotate_curve(ax6, (gf2['max_fringe']), (gf2['max_r2']), '{:.1f}'.format(max_r2_x), 1.0, (50, 0), 
                   color=regression_color, fontsize=legend_fs_scale*fs)
    
    ax5.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax5.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    ax5.set_xlim( (xlo, xhi) )
    ax5.set_ylim( (ylo, yhi) )
    set_axis_thickness(ax5, axis_thickness)
    ax5.set_xticks([])
    ax5.set_yticks([])
    
    add_moving_bound_arrows(ax5, fixed_x=ic, moving_x=gf2['max_fringe'], y_positions=y_positions, 
                            color=moving_bound_color, lw=3, mutation_scale=18, length_fraction=0.075, zorder=5)
    
    ax6.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax6.set_ylabel('Fringe $R^2$', fontsize=fs, color='black')
    ax6.set_xlim( (xlo, xhi) )
    ax6.set_ylim( (zlo, zhi) )
    
    ax6.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax6, axis_thickness)
    ax6.set_xticks([])
    ax6.set_yticks([0.9, 1.0])
    
    add_moving_bound_arrows(ax6, fixed_x=ic, moving_x=gf2['max_fringe'], y_positions=y_positions,
                            color=moving_bound_color, lw=3, mutation_scale=18, length_fraction=0.075, zorder=5)
    
    
    # Collect legend entries from ax1 and ax2.
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    
    # Remove duplicate labels while preserving their order.
    unique_entries = {}
    for handle, label in zip(handles1 + handles2, labels1 + labels2):
        if label and not label.startswith('_'):
            unique_entries.setdefault(label, handle)

    
    row_headings = [
        (ax1, '(a) Forward Fringe R^2: First guess on upper bound'),
        (ax3, '(b) Backward Fringe R^2: Find final lower bound'),
        (ax5, '(c) Forward Fringe R^2: Find final upper bound'),
    ]
    
        
    # Reserve space above the top row.
    fig.tight_layout(rect=[0, 0, 1, 0.86])
    
    # Must come after tight_layout()/subplots_adjust().
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    
    for left_ax, heading in row_headings:
        position = left_ax.get_position()
        heading_y = position.y1 + 0.008
    
        # Split at the first colon and keep the colon in the bold portion.
        bold_text, colon, normal_text = heading.partition(':')
        bold_text += colon
    
        bold_artist = fig.text(
            position.x0,
            heading_y,
            bold_text,
            ha='left',
            va='bottom',
            fontsize=fs,
            fontweight='bold',
            fontstyle='italic',
            color='black',
        )
    
        # Render and measure the bold portion.
        fig.canvas.draw()
        bold_bbox = (
            bold_artist.get_window_extent(renderer=renderer)
            .transformed(fig.transFigure.inverted())
        )
    
        # Add the remaining text immediately after the bold portion.
        fig.text(
            bold_bbox.x1,
            heading_y,
            normal_text,
            ha='left',
            va='bottom',
            fontsize=fs,
            fontweight='normal',
            #fontstyle='italic',
            color='black',
        )
        
    # Legend centered above ax1 and ax2.
    fig.legend(
        unique_entries.values(),
        unique_entries.keys(),
        loc='upper center',
        bbox_to_anchor=(0.5, 1.05),
        ncol=2,
        fontsize=fs*0.90,
        frameon=True,
        edgecolor='black',
        fancybox=True,
        columnspacing=1.2,
        handletextpad=0.5,
    )
    
    fig.savefig(
        figname + '_published_3.jpeg',
        dpi=dpi,
        bbox_inches='tight',
    )
    
    #--------------------------------------#
    # "Regionized" FRF plot -> Rubbery BFB #
    #--------------------------------------#
    # Generate a BFB Pass
    rubbery_span =  max(rubbery_line_x) - min(rubbery_line_x) 
    ic = min(rubbery_line_x) + 0.6*rubbery_span # Start at 60%
    ic = max(rubbery_line_x) + 0.6*rubbery_span # Start at 60%
    # try:
    #     high_temp_transition = regions['High-temperature transition region']
    #     #high_temp_region     = regions['High-temperature region']
    #     ic = high_temp_transition[1] #+ min_span/100 #0.1*( high_temp_region[1] - high_temp_region[0] )
    #     print(ic)
    # except: pass
    passes, rubbery_bounds, guessed_lr = rfr_multi_pass(temperature, volume, min_span, sequence='BFB', ic=ic, xlo=None, xhi=None)
    rr1, rf1, rr2 = passes
    rr1_fringe, rr1_r2 = [rr1['fringe'], rr1['r-squared']]
    rf1_fringe, rf1_r2 = [rf1['fringe'], rf1['r-squared']]
    rr2_fringe, rr2_r2 = [rr2['fringe'], rr2['r-squared']]
    
    zlo = min( [min(rr1_r2), min(rf1_r2), min(rr2_r2)] )
    zhi = max( [max(rr1_r2), max(rf1_r2), max(rr2_r2)] )
    dz  = zhi - zlo
    sf  = 0.1
    zlo -= sf*dz
    zhi += sf*dz
    
    
    fs = 14
    y_positions = (0.10, 0.30, 0.50, 0.70, 0.90)
    y_positions = (0.20, 0.40, 0.60, 0.80)
    fig, axes = plt.subplots(
        3, 2,
        figsize=(2*dim, 3.3*dim),  # Extra height for row headings
        gridspec_kw={'wspace': 0.40, 'hspace': 0.30})
    
    ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = axes
    # for pair in axes:
    #     for ax in pair:
    #         ax.set_box_aspect(1)
    
    # Backward: 1
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        ax1.plot(x, y, '-', lw=lw, color=color, label=region)
        
    ax1.axvline(ic, color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax1.axvline(rr1['max_fringe'], color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    
    segment_1 = [ic, rr1['max_fringe']]
    segment_1_x, segment_1_y = misc_funcs.reduce_data(temperature, volume, min(segment_1), max(segment_1))
    segment_1_lr = misc_funcs.linear_regression(segment_1_x, segment_1_y)
    ax1.plot(segment_1_lr.xreg, segment_1_lr.yreg, '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression within bounds')
        
    regionized_fringe_r2 = regionize(rr1_fringe, rr1_r2, regions)
    for region in regionized_fringe_r2:
        x, y, color = regionized_fringe_r2[region]
        ax2.plot(x, y, '-', lw=lw, color=color, label=region)
        
    ax2.axvline(ic, color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax2.axvline(rr1['max_fringe'], color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax2.plot(rr1['max_fringe'], rr1['max_r2'], 'o', ms=ms, mec='black', mfc=regression_color, mew=1, label='Point of maximum $R^2$')
    annotate_curve(ax2, (rr1['max_fringe']), (rr1['max_r2']), '{:.1f}'.format(rr1['max_fringe']), 1.0, (-40, 0), 
                   color=regression_color, fontsize=legend_fs_scale*fs)
    
    ax1.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax1.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    ax1.set_xlim( (xlo, xhi) )
    ax1.set_ylim( (ylo, yhi) )
    set_axis_thickness(ax1, axis_thickness)
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    add_moving_bound_arrows(ax1, fixed_x=ic, moving_x=rr1['max_fringe'], y_positions=y_positions, 
                            color=moving_bound_color, lw=3, mutation_scale=18, length_fraction=0.075, zorder=5)
    
    ax2.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax2.set_ylabel('Fringe $R^2$', fontsize=fs, color='black')
    ax2.set_xlim( (xlo, xhi) )
    ax2.set_ylim( (zlo, zhi) )
    
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax2, axis_thickness)
    ax2.set_xticks([])
    ax2.set_yticks([0.9, 1.0])
    
    add_moving_bound_arrows(ax2, fixed_x=ic, moving_x=rr1['max_fringe'], y_positions=y_positions,
                            color=moving_bound_color, lw=3, mutation_scale=18, length_fraction=0.075, zorder=5)

    
    # Update ic
    ic = min(segment_1)
    
    
    # Forward: 1
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        ax3.plot(x, y, '-', lw=lw, color=color, label=region)
        
    ax3.axvline(ic, color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax3.axvline(rf1['max_fringe'], color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    
    segment_2 = [ic, rf1['max_fringe']]
    segment_2_x, segment_2_y = misc_funcs.reduce_data(temperature, volume, min(segment_2), max(segment_2))
    segment_2_lr = misc_funcs.linear_regression(segment_2_x, segment_2_y)
    ax3.plot(segment_2_lr.xreg, segment_2_lr.yreg, '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression within bounds')
    
    regionized_fringe_r2 = regionize(rf1_fringe, rf1_r2, regions)
    for region in regionized_fringe_r2:
        x, y, color = regionized_fringe_r2[region]
        ax4.plot(x, y, '-', lw=lw, color=color, label=region)
        
    ax4.axvline(ic, color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax4.axvline(rf1['max_fringe'], color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax4.plot(rf1['max_fringe'], rf1['max_r2'], 'o', ms=ms, mec='black', mfc=regression_color, mew=1, label='Point of maximum $R^2$')
    annotate_curve(ax4, (rf1['max_fringe']), (rf1['max_r2']), '{:.1f}'.format(rf1['max_fringe']), 1.0, (-70, 0), 
                   color=regression_color, fontsize=legend_fs_scale*fs)
    
    ax3.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax3.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    ax3.set_xlim( (xlo, xhi) )
    ax3.set_ylim( (ylo, yhi) )
    set_axis_thickness(ax3, axis_thickness)
    ax3.set_xticks([])
    ax3.set_yticks([])
    
    add_moving_bound_arrows(ax3, fixed_x=ic, moving_x=rf1['max_fringe'], y_positions=y_positions, 
                            color=moving_bound_color, lw=3, mutation_scale=18, length_fraction=0.075, zorder=5)
    
    ax4.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax4.set_ylabel('Fringe $R^2$', fontsize=fs, color='black')
    ax4.set_xlim( (xlo, xhi) )
    ax4.set_ylim( (zlo, zhi) )
    
    ax4.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax4, axis_thickness)
    ax4.set_xticks([])
    ax4.set_yticks([0.9, 1.0])
    
    add_moving_bound_arrows(ax4, fixed_x=ic, moving_x=rf1['max_fringe'], y_positions=y_positions, 
                            color=moving_bound_color, lw=3, mutation_scale=18, length_fraction=0.075, zorder=5)
    
    # Update ic
    ic = max(segment_2)
    
    # Backward: 2
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        ax5.plot(x, y, '-', lw=lw, color=color, label=region)
        
    ax5.axvline(ic, color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax5.axvline(rr2['max_fringe'], color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    
    segment_3 = [ic, rr2['max_fringe']]
    segment_3_x, segment_3_y = misc_funcs.reduce_data(temperature, volume, min(segment_3), max(segment_3))
    segment_3_lr = misc_funcs.linear_regression(segment_3_x, segment_3_y)
    ax5.plot(segment_3_lr.xreg, segment_3_lr.yreg, '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression within bounds')
    
    regionized_fringe_r2 = regionize(rr2_fringe, rr2_r2, regions)
    for region in regionized_fringe_r2:
        x, y, color = regionized_fringe_r2[region]
        ax6.plot(x, y, '-', lw=lw, color=color, label=region)
        
    ax6.axvline(ic, color=fixed_bound_color, ls='--', lw=lw/2, label='Fixed bound') 
    ax6.axvline(rr2['max_fringe'], color=moving_bound_color, ls='--', lw=lw/2, label='Moving bound') 
    ax6.plot(rr2['max_fringe'], rr2['max_r2'], 'o', ms=ms, mec='black', mfc=regression_color, mew=1, label='Point of maximum $R^2$')
    
    max_r2_x = rr2['max_fringe']
    if max_r2_x == 604:# TODO: JDK override for cleaner publication plots
        max_r2_x = 600
    annotate_curve(ax6, (rr2['max_fringe']), (rr2['max_r2']), '{:.1f}'.format(max_r2_x), 1.0, (-40, 0), 
                   color=regression_color, fontsize=legend_fs_scale*fs)
    
    ax5.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax5.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    ax5.set_xlim( (xlo, xhi) )
    ax5.set_ylim( (ylo, yhi) )
    set_axis_thickness(ax5, axis_thickness)
    ax5.set_xticks([])
    ax5.set_yticks([])
    
    add_moving_bound_arrows(ax5, fixed_x=ic, moving_x=rr2['max_fringe'], y_positions=y_positions, 
                            color=moving_bound_color, lw=3, mutation_scale=18, length_fraction=0.075, zorder=5)
    
    ax6.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax6.set_ylabel('Fringe $R^2$', fontsize=fs, color='black')
    ax6.set_xlim( (xlo, xhi) )
    ax6.set_ylim( (zlo, zhi) )
    
    ax6.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax6, axis_thickness)
    ax6.set_xticks([])
    ax6.set_yticks([0.9, 1.0])
    
    add_moving_bound_arrows(ax6, fixed_x=ic, moving_x=rr2['max_fringe'], y_positions=y_positions,
                            color=moving_bound_color, lw=3, mutation_scale=18, length_fraction=0.075, zorder=5)
    
    
    # Collect legend entries from ax1 and ax2.
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    
    # Remove duplicate labels while preserving their order.
    unique_entries = {}
    for handle, label in zip(handles1 + handles2, labels1 + labels2):
        if label and not label.startswith('_'):
            unique_entries.setdefault(label, handle)
    
    row_headings = [
        (ax1, '(a) Backward Fringe R^2: First guess on lower bound'),
        (ax3, '(b) Forward Fringe R^2: Find final upper bound'),
        (ax5, '(c) Backward Fringe R^2: Find final lower bound'),
    ]
    
    # Reserve space above the top row.
    fig.tight_layout(rect=[0, 0, 1, 0.86])
    
    # Must come after tight_layout()/subplots_adjust().
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    
    for left_ax, heading in row_headings:
        position = left_ax.get_position()
        heading_y = position.y1 + 0.008
    
        # Split at the first colon and keep the colon in the bold portion.
        bold_text, colon, normal_text = heading.partition(':')
        bold_text += colon
    
        bold_artist = fig.text(
            position.x0,
            heading_y,
            bold_text,
            ha='left',
            va='bottom',
            fontsize=fs,
            fontweight='bold',
            #fontstyle='italic',
            color='black',
        )
    
        # Render and measure the bold portion.
        fig.canvas.draw()
        bold_bbox = (
            bold_artist.get_window_extent(renderer=renderer)
            .transformed(fig.transFigure.inverted())
        )
    
        # Add the remaining text immediately after the bold portion.
        fig.text(
            bold_bbox.x1,
            heading_y,
            normal_text,
            ha='left',
            va='bottom',
            fontsize=fs,
            fontweight='normal',
            fontstyle='italic',
            color='black',
        )
    
    # Legend centered above ax1 and ax2.
    fig.legend(
        unique_entries.values(),
        unique_entries.keys(),
        loc='upper center',
        bbox_to_anchor=(0.5, 1.05),
        ncol=2,
        fontsize=fs*0.90,
        frameon=True,
        edgecolor='black',
        fancybox=True,
        columnspacing=1.2,
        handletextpad=0.5,
    )
    
    fig.savefig(
        figname + '_published_4.jpeg',
        dpi=dpi,
        bbox_inches='tight',
    )
    
    
    
    #--------------------------------------#
    # IC-grid and probability density plot #
    #--------------------------------------#
    ics = np.array( csv_data['Rubbery_Guessed_ICs'] )
    tgs = np.array( csv_data['Rubbery_Possible_Tgs'] )
    lowers = np.array( csv_data['Rubbery_Possible_Lowers'] )
    uppers = np.array( csv_data['Rubbery_Possible_Uppers'] )
    slopes = np.array( csv_data['Rubbery_Possible_Slopes'] )
    
    # Deal with Tristan's annoyance
    if max(ics) == 870 and min_span == 30:
        ics = 875 - max(ics) + ics
    
    widths = uppers - lowers
    bound_penetrations = np.array( csv_data['Rubbery_Tg_penetrations'] )
    
    fs = 16
    legend_fs_scale = 0.65
    label_rel_pos   = (0.005, 0.98)
    
    # Width assigned to the colorbar relative to one plotting axis
    cbar_width = 0.055
    
    fig = plt.figure(
        figsize=((2 + cbar_width + 0.10)*dim, 2*dim + dim/4),
        layout='constrained',
    )
    
    gs = fig.add_gridspec(
        nrows=2,
        ncols=3,
        width_ratios=[1.0, 1.0, cbar_width],
        height_ratios=[1.0, 1.0],
        wspace=0.03,  # smaller value moves colorbar closer to ax2
        hspace=0.05,
    )
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    
    # Colorbar immediately to the right of ax2
    cax = fig.add_subplot(gs[0, 2])
    
    # Make the physical plotting boxes square.
    for ax in (ax1, ax2, ax3, ax4):
        ax.set_box_aspect(1)
    
    
    # Unique initial-condition values in ascending order
    ic_values = np.sort( np.unique(np.asarray(ics, dtype=float)) )
    n_colors = len(ic_values)
    
    # Place color boundaries halfway between adjacent IC values
    if n_colors == 1:
        boundaries = np.array([
            ic_values[0] - 0.5,
            ic_values[0] + 0.5,
        ])
    else:
        midpoints = 0.5*(ic_values[:-1] + ic_values[1:])
    
        boundaries = np.concatenate([
            [ic_values[0] - 0.5*(ic_values[1] - ic_values[0])],
            midpoints,
            [ic_values[-1] + 0.5*(ic_values[-1] - ic_values[-2])],
        ])
    
    # Resample coolwarm to exactly one color per IC value
    color_map = "coolwarm"
    #color_map = "hsv"
    cmap = mpl.colormaps[color_map].resampled(n_colors)
    norm = mpl.colors.BoundaryNorm(boundaries, ncolors=cmap.N, clip=True)

    # Plot possible rubbery states
    nevery = 4
    count = -nevery
    spans = []
    for ic, (_, _, guess_lr, tg_guess) in guesses.items():
        count += 1
        color = cmap(norm(ic))
        span = max(guess_lr.xreg) - min(guess_lr.xreg)
        spans.append( span )
        
        ax2.plot(guess_lr.xreg, guess_lr.yreg,  marker='|', ls='-', ms=ms/2, lw=lw/3, color=color, alpha=1.0, zorder=1/span)
        
        rfr_tg, rfr_vx, rfr_vy = find_intersection(glassy_line_x, glassy_line_y, guess_lr.xreg, guess_lr.yreg, tolerance=1e-12)
        ax2.plot(rfr_vx, rfr_vy, '--', lw=lw/6, color=regression_color, zorder=1/span)
        ax2.plot(*tg_guess, 'p', ms=ms/1.5, mec='black', mew=1.2, mfc=color, zorder=20)
        
        if count % nevery == 0:
            ax1.axvline(ic, color=color, ls='--', lw=lw/5, zorder=0.5) 

    # Normal plot
    small_zorder = 1/(2*max(spans))
    region_handles, region_labels = [], []
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        region_line_1, = ax1.plot(x, y, '-', lw=lw, color=color, zorder=small_zorder)#, label=region)
        region_line_2, = ax2.plot(x, y, '-', lw=lw, color=color, zorder=small_zorder)#, label=region)
        
        region_handles.append( region_line_1 )
        region_labels.append( str(region) )
        
    rfr_tg, rfr_vx, rfr_vy = find_intersection(glassy_line_x,  glassy_line_y, rubbery_line_x, rubbery_line_y, tolerance=1e-12)
    
    ax2.plot(glassy_line_x,  glassy_line_y, '--', lw=lw/1.5, color=regression_color, label='Glassy linear regression')
    ax1.plot(glassy_line_x,  glassy_line_y,  '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression')
    ax1.plot(rubbery_line_x, rubbery_line_y, '--', lw=lw/2, color=regression_color, zorder=10, label='Linear regression')
    ax1.plot(rfr_vx, rfr_vy, '--', lw=lw/6, color=regression_color)
    ax1.plot(*rfr_tg, 'p', ms=ms, mfc=regression_color, mec='black', mew=1.2, label="$T_g$")
    
    annotate_curve(ax1, (rfr_tg[0]), (rfr_tg[1]), "Most probable $T_g$", 0.30, (20, -20),
                   color=regression_color, fontsize=legend_fs_scale*fs*0.8)
    
    annotate_curve(ax1, glassy_line_x,  glassy_line_y, "Glassy\nlinear\nregression", 0.50, (30, -20),
                   color=regression_color, fontsize=legend_fs_scale*fs*0.8)
    
    annotate_curve(ax1, rubbery_line_x, rubbery_line_y, "Most probable\nrubbery\nlinear\nregression", 0.50, (20, -20),
                   color=regression_color, fontsize=legend_fs_scale*fs*0.8)


    ax1.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax1.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    ax1.set_xlim( (xlo, xhi) )
    ax1.set_ylim( (ylo, yhi) )
    
    ax1.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax1, axis_thickness)
    #ax1.set_xticks([])
    ax1.set_yticks([])
    
    fig.legend(
        handles=region_handles,
        labels=region_labels,
        loc='outside upper center',
        ncol=2, #len(region_handles),
        fancybox=True,
        fontsize=legend_fs_scale*fs*1.4,
        columnspacing=1.5,
        handlelength=2.5,
        #mode='expand',
    )
    
    if xpeaks is not None:
        for x in xpeaks:
            ax2.axvline(x, color=kde_color, ls='--', lw=lw/5, alpha=1, zorder=0.25) # Place a marker for labeling
    
    # Finalize color bar
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    
    # Label every fourth initial-condition value
    nevery = 5
    cbar_ticks = ic_values[::nevery]
    
    # Ensure the final value is included
    if not np.isclose(cbar_ticks[-1], ic_values[-1]):
        cbar_ticks = np.append(cbar_ticks, ic_values[-1])
    
    cbar = fig.colorbar(
        sm,
        cax=cax,
        boundaries=boundaries,
        ticks=cbar_ticks,
        spacing='uniform',
        drawedges=True,
    )
    
    cbar.set_label(
        'Initial-condition temperature (K)',
        fontsize=fs,
    )
    
    cbar.ax.tick_params(labelsize=fs)
    
    # Calculate the constrained-layout positions first
    fig.canvas.draw()
    ax2_pos = ax2.get_position()
    cax_pos = cax.get_position()
    gap = 0.125  # distance from ax1 in figure coordinates
    cax.set_in_layout(False)
    cax.set_position([
        ax2_pos.x1 + gap,  # left
        ax2_pos.y0,        # bottom
        cax_pos.width,     # width
        ax2_pos.height,    # height
    ])
    
    ax2.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax2.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    ax2.set_xlim( (xlo, xhi) )
    ax2.set_ylim( (ylo, yhi) )
    
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax2, axis_thickness)
    #ax2.set_xticks([])
    ax2.set_yticks([])
    
    ax1.text(*label_rel_pos, '(a)', transform=ax1.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
    ax2.text(*label_rel_pos, '(b)', transform=ax2.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
    
    # Representative colors from the discrete color map
    legend_colors = [cmap(0), cmap(cmap.N // 2), cmap(cmap.N - 1)]
    
    # Existing glassy-regression style
    glassy_handle = Line2D([0], [0], color=regression_color, linestyle='--', linewidth=lw/2)
    rubbery_regression_handle = Line2D([0], [0], color=regression_color, linestyle='--', linewidth=lw/1.5)
    
    # Colored dashed segments representing the IC-grid lines in ax1
    ic_grid_handles = tuple(
        Line2D(
            [0], [0],
            color=color,
            linestyle='--',
            linewidth=lw/3,
        )
        for color in legend_colors
    )
    
    ax1.legend(
        handles=[
            # glassy_handle,
            # rubbery_regression_handle,
            ic_grid_handles,
        ],
        labels=[
            # 'Glassy linear regression',
            # 'Rubbery linear regression\n(most probable)',
            'Initial-condition temperature (K)',
        ],
        handler_map={
            tuple: HandlerTuple(ndivide=None, pad=0.1),
        },
        loc='upper left',
        bbox_to_anchor=(0.075, 1.0),
        ncol=1,
        fancybox=True,
        fontsize=legend_fs_scale*fs,
    )
        
    # Several colors show that rubbery regions are colored by initial condition
    rubbery_handles = tuple(
        Line2D(
            [0], [0],
            color=color,
            linestyle='-',
            linewidth=lw/3,
            marker='|',
            markersize=ms/2,
        )
        for color in legend_colors
    )
    
    # Several colors show that Tg points use the same color mapping
    tg_handles = tuple(
        Line2D(
            [0], [0],
            linestyle='none',
            marker='p',
            markersize=ms/1.5,
            markerfacecolor=color,
            markeredgecolor='black',
            markeredgewidth=1.2,
        )
        for color in legend_colors
    )
    
    ax2.legend(
        handles=[
            glassy_handle,
            rubbery_handles,
            tg_handles,
        ],
        labels=[
            'Glassy linear regression',
            'Possible rubbery regions',
            r'Possible $T_g$ values',
        ],
        handler_map={
            tuple: HandlerTuple(ndivide=None, pad=0.1),
        },
        loc='upper left',
        bbox_to_anchor=(0.075, 1.0),
        ncol=1,
        fancybox=True,
        fontsize=legend_fs_scale*fs,
    )
    
    
    # Tgs vs IC's
    #'*', ms=ms, mfc=fill_color, mec='black', mew=1.2, 
    ax3.plot(ics, tgs, '-p', lw=lw/2, ms=ms/1.75, color=kde_color, mec='black', mew=1.2, label='Possible Tg values')
    # for ic, (_, _, guess_lr, tg_guess) in guesses.items():
    #     color = cmap(norm(ic))
    #     ax3.plot(ic, tg_guess[0], 'p', ms=ms/2, color=color)
    
    #ax3.legend(loc='upper left', bbox_to_anchor=(0, 1), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs)
    ax3.set_xlabel('Initial-condition temperature (K)', fontsize=fs, color='black')
    ax3.set_ylabel("Possible $T_g$'s", fontsize=fs, color='black')
    ax3.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax3, axis_thickness)
    #ax3.set_xticks([])
    #ax3.set_yticks([])
    ax3.text(*label_rel_pos, '(c)', transform=ax3.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
    
    ic_lo = min(ics)
    ic_hi = max(ics)
    d_ic  = ic_hi - ic_lo
    
    sf  = 0.1
    ic_lo -= sf*d_ic
    ic_hi += sf*d_ic
    ax3.set_xlim( (ic_lo, ic_hi) )
    print('ics = ', ics)
    
    tg_lo = min(tgs)
    tg_hi = max(tgs)
    dtg  = tg_hi - tg_lo
    sf  = 0.1
    tg_lo -= sf*dtg
    tg_hi += sf*dtg
    ax3.set_ylim( (tg_lo, tg_hi) )


    # Probability distribution of Tg guesses
    temp, tg_kde = csv_data['Tg_KDE']
    #csv_data['Tg_KDE_peaks'] = [xpeaks, ypeaks]
    if tg_kde is not None:   
        ax4.plot(temp, tg_kde, lw=2, color=kde_color, zorder=4, label='Convergence kernel\ndensity estimate (KDE)')
        #ax4.fill_between(temp, tg_kde, alpha=0.25, color=kde_color)
        
        ax2_twin = ax2.twinx()
        ax2_twin.set_box_aspect(1)
        ax2_twin.plot(temp, tg_kde, lw=2, color=kde_color, label='Convergence kernel\ndensity estimate (KDE)')
        ax2_twin.fill_between(temp, tg_kde, alpha=0.25, color=kde_color)
        ax2_twin.plot(xpeaks, ypeaks, '*', ms=ms/2, mfc=kde_color, mec='black', mew=1.2/3, label='KDE - Peaks')
        
        # Draw a box around this KDE
        dx = max(temp) - min(temp)
        dy = max(tg_kde) - min(tg_kde)
        sfx = 0.04
        sfy = 0.25
        
        # Box inset by 5% from the data boundaries
        ymin = -0.03
        ymax = 7.5*max(tg_kde)
        xleft   = np.min(temp)   - sfx*dx
        xright  = np.max(temp)   + sfx*dx
        ybottom = np.min(tg_kde) - sfy*dy
        ytop    = np.max(tg_kde) + sfy*dy/1.5
        ybottom = ymin - ymin/5 # OVERRIDE
        if ybottom < ymin:
            ybottom = ymin + sfy*dy/4
        
        box = np.array([
            [xleft,  ybottom],  # Bottom left
            [xright, ybottom],  # Bottom right
            [xright, ytop],     # Top right
            [xleft,  ytop],     # Top left
            [xleft,  ybottom],  # Close the box
        ])
        
        xbs, ybs = box[:, 0], box[:, 1]
        ax2_twin.plot(xbs, ybs, '-', lw=axis_thickness, color='black')        
        ax2_twin.set_ylim( ( ymin, ymax ) )
        ax2_twin.set_yticks([])
        
        connector_style = dict(
            coordsA='data',
            axesA=ax2_twin,
            coordsB='axes fraction',
            axesB=ax4,
            arrowstyle='-',
            color='black',       # Use 'black' if preferred
            linewidth=axis_thickness,
            clip_on=False,
            zorder=20,
        )
        
        # Lower-left inset corner -> upper-left corner of panel (d)
        dright = abs(xleft - xright)
        left_connector = ConnectionPatch(
            xyA=(xleft + dright*0.1, ybottom),
            xyB=(0.05, 1.0),
            **connector_style,
        )
        
        # Lower-right inset corner -> upper-right corner of panel (d)
        right_connector = ConnectionPatch(
            xyA=(xright - dright*0.0, ybottom),
            xyB=(1.0 - 0.00, 1.0),
            **connector_style,
        )
        
        # Prevent connectors from changing the constrained layout.
        left_connector.set_in_layout(False)
        right_connector.set_in_layout(False)
        
        fig.add_artist(left_connector)
        fig.add_artist(right_connector)
        
    if xpeaks is not None:
        ax4.plot(xpeaks, ypeaks, '*', ms=ms, mfc=kde_color, mec='black', mew=1.2, label='KDE - Peaks')
        
        # Annotate peaks with selection criterion labels
        ymax = max(tg_kde)
        xrange = max(temp) - min(temp)
        ytext_hi = 1.20*ymax
        ytext_lo = 0.90*ymax
        max_labeled_peak = 6
        selected_label = ''
        tg_selection = 'prob'
        for n, info in enumerate(peak_info):
            label_rank = info.get('label_rank', n)
            if label_rank > max_labeled_peak and not info['selected']:
                continue
            
            xt = info['temp']
            yd = info['density']

            # Alternate small left/right x-offsets
            dx = (-1)**n * 0.02*xrange
            dx = 0
            xtxt = xt + dx

            # Alternate between two label heights
            ytext = ytext_hi if n % 2 == 0 else ytext_lo
            #ytext = ytext_lo

            if info['selected']:
                bbox = dict(boxstyle='round,pad=0.25', fc='white', ec='black', lw=1.2)
                arrow_color = 'black'
                text = '{}\n{}'.format(selected_label, info['label'])
                text = info['label']
                weight = 'bold'
                csv_data['selected_peak'] = [[xt], [yd]]
            else:
                bbox = dict(boxstyle='round,pad=0.25', fc='white', ec='black', lw=0.8, alpha=0.9)
                arrow_color = '0.35'
                text = info['label']
                weight = 'normal'

            # Change annotation formating for publication figures
            prob = float(info['label'].split('\n')[-1])
            text = 'Probability\n{:.4f}'.format(prob)
            
            
            arrowprops = dict(arrowstyle='->', color=arrow_color, lw=1.0, shrinkA=0, shrinkB=4)
            ax4.annotate(text, xy=(xt, yd), xytext=(xtxt, ytext), textcoords='data', ha='center', va='bottom',
                         fontsize=legend_fs_scale*fs*0.9, fontweight=weight, bbox=bbox, arrowprops=arrowprops, zorder=10)
            
            
        # Plot left/right bases used for peak probability integration
        if tg_selection == 'prob':
                    
            # Resample coolwarm to exactly one color per IC value    
            kernel_effective_width = settings.get('kernel_width', ic_increment)
            color_map = "coolwarm"
            n_colors = min(max_labeled_peak, len(peak_info))
            cmap = mpl.colormaps[color_map].resampled(max(n_colors, 1))
            x_to_ys = []
            for n, info in enumerate(peak_info):
                lb = info['left_base']
                rb = info['right_base']
                xpeak = xpeaks[n]
                rounded = round(xpeak / 5) * 5
                
                # Assign colors to the labeled peaks; use gray for any remaining peaks
                if n < n_colors:
                    fill_color = cmap(n / max(n_colors - 1, 1))
                else:
                    fill_color = '0.75'
    
                xlb = temp[lb]
                xrb = temp[rb]
                ylb = tg_kde[lb]
                yrb = tg_kde[rb]
                key = 'Tg_KDE_base_{}'.format(n)
                csv_data[key] = [[xlb, xrb], [ylb, yrb]]
    
                if info['selected']:
                    base_color = 'black'
                    base_lw = 1.2
                    base_alpha = 0.8
                    shade_alpha = 0.65
                else:
                    base_color = '0.35'
                    base_lw = 0.8
                    base_alpha = 0.45
                    shade_alpha = 0.03

                    
                base_lw = 1.2
                base_alpha = 0.8
                shade_alpha = 0.75
                #ax4.axvline(xlb, color=base_color, ls=':', lw=base_lw, alpha=base_alpha, zorder=3)
                #ax4.axvline(xrb, color=base_color, ls=':', lw=base_lw, alpha=base_alpha, zorder=3)
    
                ax4.plot(xpeaks[n], ypeaks[n], '*', ms=ms, mfc=fill_color, mec='black', mew=1.2, alpha=shade_alpha, zorder=5)
                ax4.fill_between(temp[lb:rb+1], 0, tg_kde[lb:rb+1], color=fill_color, alpha=shade_alpha, zorder=2)
                
                ax4.plot(xlb, ylb, '|', ms=ms, color=base_color, alpha=base_alpha, zorder=6)
                ax4.plot(xrb, yrb, '|', ms=ms, color=base_color, alpha=base_alpha, zorder=6)
                

                
                # Change annotation formating for publication figures
                prob = float(info['label'].split('\n')[-1])
                text = 'Probability={:.4f}'.format(prob)
                label = 'KDE - Peak:{} ({})'.format(n, text)
                label = 'KDE - Peak:{} region'.format(n)
                label = 'Peak center: {} K'.format(rounded)
                print(xpeaks)
                
                # Matching horizontal peak-base region in ax3
                ax3.axhspan(xlb, xrb, facecolor=fill_color, edgecolor='none', alpha=shade_alpha, zorder=0, label=label)
                x_to_ys.append(xlb)
                x_to_ys.append(xrb)
                
                # Plot the actualy width of the kernel
                lo_peak = min(temp[lb:rb+1])
                hi_peak = max(temp[lb:rb+1])
                tgs_in_region = []
                for tg_ in tgs:
                    if lo_peak <= tg_ <= hi_peak:
                        tgs_in_region.append(tg_)
                
                if tgs_in_region:
                    label = 'KDE - Peak:{} region'.format(n)
                    kde_lo = min(tgs_in_region) - kernel_effective_width/2
                    kde_hi = max(tgs_in_region) + kernel_effective_width/2
                    #ax3.axhspan(kde_lo, kde_hi, facecolor=fill_color, edgecolor='none', alpha=shade_alpha, zorder=0, label=label)
                    #ax3.axhspan(xpeak-kernel_effective_width/2, xpeak+kernel_effective_width/2, facecolor=fill_color, edgecolor='black', alpha=1.1*shade_alpha, zorder=0)

                
                
            ax4.plot(xlb, ylb, '|', ms=ms, color='black', alpha=base_alpha, zorder=6, label=' * Left base')
            ax4.plot(xrb, yrb, '|', ms=ms, color='black', alpha=base_alpha, zorder=6, label=' * Right base')

        ax4.set_ylim(top=2.3*ymax)
        ax3.set_ylim( (min(x_to_ys), max(x_to_ys)) )
    

    ax3.legend(loc='upper left', bbox_to_anchor=(0.370, 1), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs*0.85)
    ax4.legend(loc='upper left', bbox_to_anchor=(0.075, 1), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs)
    ax4.text(*label_rel_pos, '(d)', transform=ax4.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
    ax4.set_xlabel("Possible $T_g$'s", fontsize=fs, color='black')
    ax4.set_ylabel('Estimated probablity density', fontsize=fs, color='black')
    ax4.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax4, axis_thickness)
    ax4.set_xticks( ax3.get_yticks() )
    #ax4.set_yticks([])
    
    klo = min(temp)
    khi = max(temp)
    dk  = khi - klo
    sf  = 0.1
    klo -= sf*dk
    khi += sf*dk
    ax4.set_xlim( (klo, khi) )
    
    fig.savefig(
        figname + '_published_5.jpeg',
        dpi=dpi,
        bbox_inches='tight',
    )
    
    
    #------------------------------------#
    # Probability density weighting plot #
    #------------------------------------#
    fs = 15
    legend_fs_scale = 0.65
    label_rel_pos   = (0.01, 0.98)
    
    # # Width assigned to the colorbar relative to one plotting axis
    cbar_width = 0.055
    
    fig = plt.figure(
        figsize=((2 + cbar_width + 0.50)*dim, 4*dim + dim/1.5),
        layout='constrained',
    )
    
    gs = fig.add_gridspec(
        nrows=4,
        ncols=3,
        width_ratios=[1.0, 1.0, cbar_width],
        height_ratios=[1.0, 0.5, 1.0, 0.5],
        wspace=0.08,
        hspace=0.025,
    )
    
    # Row 1: wide plot and colorbar
    ax1 = fig.add_subplot(gs[0, :2])
    cax = fig.add_subplot(gs[0, 2])
    
    # Row 2
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    
    # Row 3
    ax4 = fig.add_subplot(gs[2, 0])
    ax5 = fig.add_subplot(gs[2, 1])
    
    # Row 4: span both plotting columns
    ax6 = fig.add_subplot(gs[3, :2])
    
    # Keep panels ax4 and ax5 square and top-aligned
    for ax in (ax4, ax5):
        #ax.set_box_aspect(1)
        ax.set_anchor('N')
    
    # Allow the spanning panels to use their full width
    for ax in (ax1, ax6):
        ax.set_aspect('auto')
        
    # Zoomed in axis ax3 and ax4 (b and c)
    repeat_xlo, repeat_xhi = 320, 460
    repeat_ylo, repeat_yhi = 1.7, 2.2
    
    out_of_bound_xlo, out_of_bound_xhi = 320, 700
    out_of_bound_ylo, out_of_bound_yhi = 1.7, 4.1
    
    # Unique initial-condition values in ascending order
    ic_values = np.sort( np.unique(np.asarray(ics, dtype=float)) )
    n_colors = len(ic_values)
    
    # Place color boundaries halfway between adjacent IC values
    if n_colors == 1:
        boundaries = np.array([
            ic_values[0] - 0.5,
            ic_values[0] + 0.5,
        ])
    else:
        midpoints = 0.5*(ic_values[:-1] + ic_values[1:])
    
        boundaries = np.concatenate([
            [ic_values[0] - 0.5*(ic_values[1] - ic_values[0])],
            midpoints,
            [ic_values[-1] + 0.5*(ic_values[-1] - ic_values[-2])],
        ])
    
    # Resample coolwarm to exactly one color per IC value
    color_map = "coolwarm"
    #color_map = "hsv"
    cmap = mpl.colormaps[color_map].resampled(n_colors)
    norm = mpl.colors.BoundaryNorm(boundaries, ncolors=cmap.N, clip=True)

    # Plot possible rubbery states
    nevery = 4
    count = -nevery
    n = 0
    spans = []
    for ic, (_, _, guess_lr, tg_guess) in guesses.items():
        n += 1
        count += 1
        color = cmap(norm(ic))
        span = max(guess_lr.xreg) - min(guess_lr.xreg)
        spans.append( span )
        
        rfr_tg, rfr_vx, rfr_vy = find_intersection(glassy_line_x, glassy_line_y, guess_lr.xreg, guess_lr.yreg, tolerance=1e-12)
        ax1.plot(guess_lr.xreg, guess_lr.yreg,  marker='|', ls='-', ms=ms/2, lw=lw/3, color=color, alpha=1.0, zorder=1/span)
        ax1.plot(rfr_vx, rfr_vy, '--', lw=lw/6, color=regression_color, zorder=2)
        ax1.plot(*tg_guess, 'p', ms=ms/1.5, mec='black', mew=1.2, mfc=color, zorder=10)
    
        if n == 10:
            ax2.plot(guess_lr.xreg, guess_lr.yreg,  marker='|', ls='-', ms=ms/1.5, lw=lw/2, color=color, alpha=1.0)
            ax2.plot(rfr_vx, rfr_vy, '--', lw=lw/5, color=regression_color, zorder=2)
            ax2.plot(*tg_guess, 'p', ms=ms/1.25, mec='black', mew=1.2, mfc=color, zorder=2)
            
            bpp_x = [ min(guess_lr.xreg), tg_guess[0] ]
            bpp_y = [ min(guess_lr.yreg), min(guess_lr.yreg) ]
            bpp_y = [i - min(guess_lr.yreg)/16 for i in bpp_y]
            ax2.plot(bpp_x, bpp_y, marker='|', ls='-', ms=ms/1.25, lw=lw/1.5, color='black', alpha=1.0)
            for x in bpp_x:
                ax2.axvline(x, color='black', ls='--', lw=lw/4, alpha=1.0, zorder=0.125)
            
            bpp = max(bpp_x) - min(bpp_x)
            string = 'Rubbery bound\npenetration: {:.1f}'.format(bpp)
            annotate_curve(ax2, bpp_x, bpp_y, string, 0.25, (-90, 40), 
                           color='black', fontsize=legend_fs_scale*fs)
                
        
        if n == 11:
            ax3.plot(guess_lr.xreg, guess_lr.yreg,  marker='|', ls='-', ms=ms/1.5, lw=lw/2, color=color, alpha=1.0)
            ax3.plot(rfr_vx, rfr_vy, '--', lw=lw/5, color=regression_color, zorder=2)
            ax3.plot(*tg_guess, 'p', ms=ms/1.25, mec='black', mew=1.2, mfc=color, zorder=2)
            
            bpp_x = [ min(guess_lr.xreg), tg_guess[0] ]
            bpp_y = [ min(guess_lr.yreg), min(guess_lr.yreg) ]
            bpp_y = [i - min(guess_lr.yreg)/6 for i in bpp_y]
            ax3.plot(bpp_x, bpp_y, marker='|', ls='-', ms=ms/1.25, lw=lw/1.5, color='black', alpha=1.0)
            for x in bpp_x:
                ax3.axvline(x, color='black', ls='--', lw=lw/4, alpha=1.0, zorder=0.125)
            
            bpp = max(bpp_x) - min(bpp_x)
            string = 'Rubbery bound\npenetration: {:.1f}'.format(bpp)
            annotate_curve(ax3, bpp_x, bpp_y, string, 0.25, (-120, 0), 
                           color='black', fontsize=legend_fs_scale*fs)
            
    # Normal plot
    small_zorder = 1/(2*max(spans))
    region_handles, region_labels = [], []
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        region_line_1, = ax1.plot(x, y, '-', lw=lw, color=color, zorder=small_zorder)
        region_line_1, = ax2.plot(x, y, '-', lw=lw, color=color, zorder=small_zorder)
        region_line_2, = ax3.plot(x, y, '-', lw=lw, color=color, zorder=small_zorder)
        
        region_handles.append( region_line_1 )
        region_labels.append( str(region) )
        
    ax1.plot(glassy_line_x,  glassy_line_y, '--', lw=lw/2, color=regression_color, label='Glassy linear regression')
    ax2.plot(glassy_line_x,  glassy_line_y, '--', lw=lw/2, color=regression_color, label='Glassy linear regression')
    ax3.plot(glassy_line_x,  glassy_line_y, '--', lw=lw/2, color=regression_color, label='Glassy linear regression')
    
    # Finalize color bar
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    
    # Label every fourth initial-condition value
    nevery = 5
    cbar_ticks = ic_values[::nevery]
    
    # Ensure the final value is included
    if not np.isclose(cbar_ticks[-1], ic_values[-1]):
        cbar_ticks = np.append(cbar_ticks, ic_values[-1])
    
    cbar = fig.colorbar(
        sm,
        cax=cax,
        boundaries=boundaries,
        ticks=cbar_ticks,
        spacing='uniform',
        drawedges=True,
    )
    
    cbar.set_label(
        'Initial-condition temperature (K)',
        fontsize=fs,
    )
    
    cbar.ax.tick_params(labelsize=fs)
    
    # Calculate the constrained-layout positions first
    fig.canvas.draw()
    ax1_pos = ax1.get_position()
    cax_pos = cax.get_position()
    gap = 0.075  # distance from ax1 in figure coordinates
    cax.set_in_layout(False)
    cax.set_position([
        ax1_pos.x1 + gap,  # left
        ax1_pos.y0 - 0.05, # bottom
        cax_pos.width,     # width
        ax1_pos.height,    # height
    ])
    
    #ax1.set_xlabel('{}{}'.format(0*" ", temp_label), fontsize=fs, color='black')
    ax1.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    ax1.set_xlim( (xlo, xhi) )
    ax1.set_ylim( (ylo, yhi) )
    ax1.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax1, axis_thickness)
    ax1.set_xticks([0, 250, 500, 750, 1000])
    ax1.set_yticks([])
    
    ax2.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax2.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')    
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax2, axis_thickness)
    #ax2.set_xticks([])
    ax2.set_yticks([])
    
    ax3.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax3.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')    
    ax3.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax3, axis_thickness)
    #ax3.set_xticks([])
    ax3.set_yticks([])
    
    ax1.text(*label_rel_pos, '(a)', transform=ax1.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
    ax2.text(*label_rel_pos, '(b)', transform=ax2.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
    ax3.text(*label_rel_pos, '(c)', transform=ax3.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
    
    fig.legend(
        handles=region_handles,
        labels=region_labels,
        loc='outside upper center',
        ncol=2, #len(region_handles),
        fancybox=True,
        fontsize=legend_fs_scale*fs*1.4,
        columnspacing=1.5,
        handlelength=3.0,
        #mode='expand',
    )

    # Representative colors from the discrete color map
    legend_colors = [cmap(0), cmap(cmap.N // 2), cmap(cmap.N - 1)]
    
    # Existing glassy-regression style
    glassy_handle = Line2D([0], [0], color=regression_color, linestyle='--', linewidth=lw/2)
    
    # Several colors show that rubbery regions are colored by initial condition
    rubbery_handles = tuple(
        Line2D(
            [0], [0],
            color=color,
            linestyle='-',
            linewidth=lw/3,
            marker='|',
            markersize=ms/2,
        )
        for color in legend_colors
    )
    
    # Several colors show that Tg points use the same color mapping
    tg_handles = tuple(
        Line2D(
            [0], [0],
            linestyle='none',
            marker='p',
            markersize=ms/1.5,
            markerfacecolor=color,
            markeredgecolor='black',
            markeredgewidth=1.2,
        )
        for color in legend_colors
    )
    
    ax1.legend(
        handles=[
            glassy_handle,
            rubbery_handles,
            tg_handles,
        ],
        labels=[
            'Glassy linear regression',
            'Possible rubbery regions',
            r'Possible $T_g$ values',
        ],
        handler_map={
            tuple: HandlerTuple(ndivide=None, pad=0.1),
        },
        loc='upper left',
        bbox_to_anchor=(0.05, 1.0),
        ncol=3,
        fancybox=True,
        fontsize=legend_fs_scale*fs,
    )
    
    
    # Zoom in on repeat
    # lower_bound_color = 'tab:blue'
    # upper_bound_color = 'tab:red'
    # tgs_color         = 'tab:purple'    
    # slopes_color      = 'tab:cyan'
    # widths_color      = 'tab:gray'
    # bounds_color      = 'tab:blue'
    slopes_color      = 'tab:cyan'
    widths_color      = 'tab:gray'
    bounds_color      = 'tab:blue'
    ax2.set_xlim( (repeat_xlo, repeat_xhi) )
    ax2.set_ylim( (repeat_ylo, repeat_yhi) )
    
    # Zoom in on out_of_bound
    ax3.set_xlim( (out_of_bound_xlo, out_of_bound_xhi) )
    ax3.set_ylim( (out_of_bound_ylo, out_of_bound_yhi) )
    
    # Draw a box around repeat
    xleft   = repeat_xlo
    xright  = repeat_xhi
    ybottom = repeat_ylo
    ytop    = repeat_yhi
    box = np.array([
        [xleft,  ybottom],  # Bottom left
        [xright, ybottom],  # Bottom right
        [xright, ytop],     # Top right
        [xleft,  ytop],     # Top left
        [xleft,  ybottom],  # Close the box
    ])
    xbs, ybs = box[:, 0], box[:, 1]
    ax1.plot(xbs, ybs, '-', lw=axis_thickness, color='black')  
    
    connector_style = dict(
        coordsA='data',
        axesA=ax1,
        coordsB='axes fraction',
        axesB=ax2,
        arrowstyle='-',
        color='black',       # Use 'black' if preferred
        linewidth=axis_thickness,
        clip_on=False,
        zorder=20,
    )
    
    # Lower-left inset corner -> upper-left corner of panel (d)
    left_connector = ConnectionPatch(
        xyA=(xleft, ybottom),
        xyB=(0.00125, 1.0),
        **connector_style,
    )
    
    # Lower-right inset corner -> upper-right corner of panel (d)
    right_connector = ConnectionPatch(
        xyA=(xright, ybottom),
        xyB=(1.0 - 0.00125, 1.0),
        **connector_style,
    )
    
    # Prevent connectors from changing the constrained layout.
    left_connector.set_in_layout(False)
    right_connector.set_in_layout(False)
    fig.add_artist(left_connector)
    fig.add_artist(right_connector)
    
    # Draw a box around out_of_bound
    xleft   = out_of_bound_xlo
    xright  = out_of_bound_xhi
    ybottom = out_of_bound_ylo
    ytop    = out_of_bound_yhi
    box = np.array([
        [xleft,  ybottom],  # Bottom left
        [xright, ybottom],  # Bottom right
        [xright, ytop],     # Top right
        [xleft,  ytop],     # Top left
        [xleft,  ybottom],  # Close the box
    ])
    xbs, ybs = box[:, 0], box[:, 1]
    ax1.plot(xbs, ybs, '-', lw=axis_thickness, color='black')  
    
    connector_style = dict(
        coordsA='data',
        axesA=ax1,
        coordsB='axes fraction',
        axesB=ax3,
        arrowstyle='-',
        color='black',       # Use 'black' if preferred
        linewidth=axis_thickness,
        clip_on=False,
        zorder=20,
    )
    
    # Lower-left inset corner -> upper-left corner of panel (d)
    left_connector = ConnectionPatch(
        xyA=(xleft, ybottom),
        xyB=(0.00125, 1.0),
        **connector_style,
    )
    
    # Lower-right inset corner -> upper-right corner of panel (d)
    right_connector = ConnectionPatch(
        xyA=(xright, ybottom),
        xyB=(1.0 - 0.00125, 1.0),
        **connector_style,
    )
    
    # Prevent connectors from changing the constrained layout.
    left_connector.set_in_layout(False)
    right_connector.set_in_layout(False)
    fig.add_artist(left_connector)
    fig.add_artist(right_connector)
    
    # Possibilities
    #ax3.plot(ics, lowers, '-v', lw=lw/2, ms=ms/2, color=lower_bound_color,  label='Lower bound')
    #ax3.plot(ics, uppers, '-^', lw=lw/2, ms=ms/2, color=upper_bound_color,  label='Upper bound')
    #ax3.plot(ics, tgs,    '-p', lw=lw/2, ms=ms/2, color=tgs_color,     label="$T_g$'s")
    ax4.plot(ics, widths, '-s', lw=lw/3, ms=ms/3, color=widths_color,  label="Rubbery widths")
    ax4.plot(ics, bound_penetrations, '-o', lw=lw/3, ms=ms/3, color=bounds_color,  label="Bound penetrations")
    

    ax4_twin = ax4.twinx()
    #ax4_twin.set_box_aspect(1)
    ax4_twin.set_anchor('N')
    ax4_twin.plot(ics, slopes, '-d', lw=lw/3, ms=ms/3, color=slopes_color, label="Rubbery slopes")
    
    ax4_twin.set_xlim( (ic_lo, ic_hi) )
    ax4_twin.set_ylabel('Possible rubbery slopes ($Å^3$/K)', fontsize=fs, color=slopes_color)
    #ax4_twin.set_yticks([])
    ax4_twin.tick_params(axis='y', labelsize=fs, colors=slopes_color)
    
    ax4_twin.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax4_twin, axis_thickness)

    
    handles_1, labels_1 = ax4.get_legend_handles_labels()
    handles_2, labels_2 = ax4_twin.get_legend_handles_labels()
    ax4.legend(handles_2 + handles_1, labels_2 + labels_1,
        loc='upper left', bbox_to_anchor=(0.1, 1), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs*0.95)

    ax4.set_xlabel('Initial-condition temperature (K)', fontsize=fs, color='black')
    ax4.set_ylabel('Possible rubbery states (K)', fontsize=fs, color='black')
    ax4.tick_params(axis='x', labelsize=fs, colors='black')
    ax4.tick_params(axis='y', labelsize=fs, colors='black')
    
    nlo = min([min(widths),  min(bound_penetrations) ])
    nhi = max([max(widths),  max(bound_penetrations) ])
    dn  = nhi - nlo
    
    sf  = 0.1
    nlo -= sf*dn
    nhi += 5*sf*dn
    ax4.set_xlim( (ic_lo, ic_hi) )
    ax4.set_ylim( (nlo, nhi) )
    
    
    # Weighting
    glassy_line = outputs['glassy_line'] # [glassy_lr.b0,  glassy_lr.b1,  glassy_lr.r2]
    weighted_rubbery_widths = widths/max(widths)    # Normalized by rubbery-max-length
    weighted_rubbery_slopes = slopes/glassy_line[1] # Normalized by glassy-slope
    
    tau = min_span/4
    weighted_bound_penatration = np.exp(-bound_penetrations / tau)
    #print(weighted_bound_penatration)
    
    weights = weighted_bound_penatration * weighted_rubbery_widths *  weighted_rubbery_slopes    
    ax5.plot(ics, weighted_rubbery_widths,    '-s', lw=lw/3, ms=ms/3, color=widths_color,  label='Rubbery width weight (RWW)', zorder=1)
    ax5.plot(ics, weighted_bound_penatration, '-o', lw=lw/3, ms=ms/3, color=bounds_color,  label='Bound penetration weight (BPW)', zorder=1)
    ax5.plot(ics, weights,                    '-p', lw=lw/3, ms=ms/2, color=kde_color,     label='Weights = RSW*RWW*BPW', zorder=100)
    
    ax5_twin = ax5.twinx()
    #ax5_twin.set_box_aspect(1)
    ax5_twin.set_anchor('N')
    ax5.set_zorder(ax5_twin.get_zorder() + 1) # Draw ax5 above ax5_twin
    ax5.patch.set_visible(False) # Prevent ax5's background from hiding the twin-axis curve
    
    ax5_twin.set_xlim( (ic_lo, ic_hi) )
    ax5_twin.set_ylabel('Rubbery slopes weight', fontsize=fs, color=slopes_color)
    ax5_twin.tick_params(axis='y', labelsize=fs, colors=slopes_color)
    ax5_twin.plot(ics, weighted_rubbery_slopes,    '-d', lw=lw/3, ms=ms/3, color=slopes_color,  label='Rubbery slope weight (RSW)', zorder=1)

    
    handles_1, labels_1 = ax5.get_legend_handles_labels()
    handles_2, labels_2 = ax5_twin.get_legend_handles_labels()
    ax5.legend(handles_2 + handles_1, labels_2 + labels_1,
        loc='upper left', bbox_to_anchor=(0.1, 1), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs*0.95)
    
    ax5.set_xlabel('Initial-condition temperature (K)', fontsize=fs, color='black')
    ax5.set_ylabel('KDE - weights', fontsize=fs, color='black')
    ax5.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax5, axis_thickness)
    
    nlo = min([min(weighted_rubbery_widths),  min(weighted_bound_penatration), min(weights) ])
    nhi = max([max(weighted_rubbery_widths),  max(weighted_bound_penatration), max(weights) ])
    dn  = nhi - nlo
    
    sf  = 0.1
    nlo -= sf*dn
    nhi += 10*sf*dn
    
    ax5.set_xlim( (ic_lo, ic_hi) )
    ax5.set_ylim( (nlo, nhi) )
    ax4.text(*label_rel_pos, '(d)', transform=ax4.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
    ax5.text(*label_rel_pos, '(e)', transform=ax5.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
    
    # Bound penetration equation
    tau = min_span/4
    bound_penetration = np.linspace(0, 100, 10_000)
    bound_weight = np.exp(-bound_penetration / tau)
    
    label = r'$y = \exp\left(-\frac{x}{\tau}\right)$'
    ax6.plot(bound_penetration, bound_weight, '-o', lw=lw/3, ms=ms/3, color=bounds_color,  label=label)
    
    ax6.set_xlabel("Bound penetrations", fontsize=fs, color='black')
    ax6.set_ylabel('Bound\npenetration\nweight', fontsize=fs, color='black')
    ax6.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax6, axis_thickness)
    ax6.text(*label_rel_pos, '(f)', transform=ax6.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
    ax6.legend(loc='upper left', bbox_to_anchor=(0.075, 1), fancybox=True, ncol=1, fontsize=fs)
    
    fig.savefig(
        figname + '_published_6.jpeg',
        dpi=dpi,
        bbox_inches='tight',
    )
    
    #---------------------------------------------------------#
    # Alternative methods - Patrone and Segemented regression #
    #---------------------------------------------------------#
    regression_color = '#0097A7' #regression_color
    patrone_color    = '#E69F00' #'tab:green'
    pr_color         = '#009E73' #'tab:cyan'
    fs = 16
    legend_fs_scale = 0.9
    label_rel_pos   = (0.02, 0.98)
    
    fig = plt.figure(figsize=(2.5*dim, 2.5*dim))
    outer_gs = fig.add_gridspec(
        nrows=2,
        ncols=1,
        left=0.10,      # room for ax1 and ax3 y-labels
        right=0.98,
        bottom=0.10,
        top=0.76,
        hspace=0.40,
    )
    
    # Independent spacing for each row
    top_gs    = outer_gs[0].subgridspec(1, 2, wspace=0.35)
    bottom_gs = outer_gs[1].subgridspec(1, 3, wspace=0.55)
    
    ax1 = fig.add_subplot(top_gs[0, 0])
    ax2 = fig.add_subplot(top_gs[0, 1])
    
    ax3 = fig.add_subplot(bottom_gs[0, 0])
    ax4 = fig.add_subplot(bottom_gs[0, 1])
    ax5 = fig.add_subplot(bottom_gs[0, 2])
    # for ax in (ax1, ax2):#, ax3, ax4, ax5):
    #     ax.set_box_aspect(1)
    
    # RFR Tg and CTE's
    rfr_tg, rfr_vx, rfr_vy = find_intersection(glassy_line_x,  glassy_line_y, rubbery_line_x, rubbery_line_y, tolerance=1e-12)
    
    # Compute some CTE's
    def compute_CTE(xpoints, ypoints):
        x1, x2 = xpoints
        y1, y2 = ypoints
        return (y2 - y1) / (x2 - x1)
    rfr_glassy_cte  = compute_CTE(glassy_line_x,  glassy_line_y)
    rfr_rubbery_cte = compute_CTE(rubbery_line_x, rubbery_line_y)

    
    # Reduce temperature and volume to rfr_bounds
    rfr_xlo = min(glassy_line_x)
    rfr_xhi = max(rubbery_line_x)
    reduced_temp, reduced_vol = misc_funcs.reduce_data(temperature, volume, rfr_xlo, rfr_xhi)
    
    # Piecewise Regression (pr) Tg and CTE's
    pr_x_unbounded, pr_y_unbounded, xbreaks, ybreaks, slopes = misc_funcs.piecewise_regression(temperature, volume, rfr_xlo, rfr_xhi, 1)
    pr_unbounded_tg = (xbreaks[1], ybreaks[1])
    pr_unbounded_glassy_cte  = compute_CTE(xbreaks[0:2], ybreaks[0:2])
    pr_unbounded_rubbery_cte = compute_CTE(xbreaks[1:3], ybreaks[1:3])

    pr_x_bounded, pr_y_bounded, xbreaks, ybreaks, slopes = misc_funcs.piecewise_regression(reduced_temp, reduced_vol, rfr_xlo, rfr_xhi, 1)
    pr_bounded_tg = (xbreaks[1], ybreaks[1])
    pr_bounded_glassy_cte  = compute_CTE(xbreaks[0:2], ybreaks[0:2])
    pr_bounded_rubbery_cte = compute_CTE(xbreaks[1:3], ybreaks[1:3])
    
    # Patrone (h for hyperbola) Tg and CTE's
    slope_method = 'tangent'
    slope_method = 'asymptote'
    out = misc_funcs.fit_hyperbola(temperature, volume, rfr_xlo, rfr_xhi, minimum_convergence=None, initial_guess=True, maxiter=10**4, slope_method=slope_method)
    h_x_unbounded, h_y_unbounded, params, center, slopes, transition, h_unbounded_tg, h_unbounded_tangents = out
    h_unbounded_glassy_cte  = compute_CTE(h_unbounded_tangents[0][0:2], h_unbounded_tangents[1][0:2])
    h_unbounded_rubbery_cte = compute_CTE(h_unbounded_tangents[0][1:3], h_unbounded_tangents[1][1:3])
    
    out = misc_funcs.fit_hyperbola(reduced_temp, reduced_vol, rfr_xlo, rfr_xhi, minimum_convergence=None, initial_guess=True, maxiter=10**4, slope_method=slope_method)
    h_x_bounded, h_y_bounded, params, center, slopes, transition, h_bounded_tg, h_bounded_tangents = out
    h_bounded_glassy_cte  = compute_CTE(h_bounded_tangents[0][0:2], h_bounded_tangents[1][0:2])
    h_bounded_rubbery_cte = compute_CTE(h_bounded_tangents[0][1:3], h_bounded_tangents[1][1:3])
    
    
    # Normal plot
    for region in regionized_temp_volume:
        x, y, color = regionized_temp_volume[region]
        ax1.plot(x, y, '-', lw=lw*1.2, color=color)#, label=region)
        ax2.plot(x, y, '-', lw=lw*1.2, color=color)#, label=region)
        
        
    
    # Segmented regression - plotting
    ax1.plot(pr_x_unbounded, pr_y_unbounded, '--', lw=lw/1.2, color=pr_color, label='Piecewise regression')
    ax2.plot(pr_x_bounded,   pr_y_bounded,   '--', lw=lw/1.2, color=pr_color, label='Piecewise regression')
    #ax2.plot(reduced_temp, reduced_vol, '-', lw=lw, color='black')
    
    ax1.plot(*pr_unbounded_tg, 'p', ms=ms, mfc=pr_color, mec='black', mew=1.2, label="Piecewise regression $T_g$")
    ax2.plot(*pr_bounded_tg, 'p', ms=ms, mfc=pr_color, mec='black', mew=1.2, label="Piecewise regression $T_g$")
        
    # Patrone plotting
    ax1.plot(h_x_unbounded, h_y_unbounded, '-', lw=lw/1.7, color=patrone_color, label='Hyperbola fit')
    ax1.plot(*h_unbounded_tangents, '--', lw=lw/3, color=patrone_color, label='Hyperbola tangent intersection')
    ax1.plot(*h_unbounded_tg, 'p', ms=ms, mfc=patrone_color, mec='black', mew=1.2, label="Hyperbola $T_g$")
    
    ax2.plot(h_x_bounded, h_y_bounded, '-', lw=lw/1.7, color=patrone_color, label='Hyperbola fit')
    ax2.plot(*h_bounded_tangents, '--', lw=lw/3, color=patrone_color, label='Hyperbola tangent intersection')
    ax2.plot(*h_bounded_tg, 'p', ms=ms, mfc=patrone_color, mec='black', mew=1.2, label="Hyperbola $T_g$")

        
    # RFR - plotting
    ax1.plot(glassy_line_x,  glassy_line_y, '--', lw=lw/2, color=regression_color, label='Glassy linear regression (RFR-thermal)')
    ax2.plot(glassy_line_x,  glassy_line_y, '--', lw=lw/2, color=regression_color, label='Glassy linear regression (RFR-thermal)')
    
    ax1.plot(rubbery_line_x, rubbery_line_y, '--', lw=lw/2, color=regression_color, label='Rubbery linear regression (RFR-thermal)')
    ax2.plot(rubbery_line_x, rubbery_line_y, '--', lw=lw/2, color=regression_color, label='Rubbery linear regression (RFR-thermal)')
    
    ax1.plot(*rfr_tg, 'p', ms=ms, mfc=regression_color, mec='black', mew=1.2, label="$T_g$ (RFR-thermal)")
    ax2.plot(*rfr_tg, 'p', ms=ms, mfc=regression_color, mec='black', mew=1.2, label="$T_g$ (RFR-thermal)")
    
    ax1.plot(rfr_vx, rfr_vy, '--', lw=lw/6, color=regression_color)
    ax2.plot(rfr_vx, rfr_vy, '--', lw=lw/6, color=regression_color)

    
    ax1.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax1.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    ax1.set_xlim( (xlo, xhi) )
    ax1.set_ylim( (ylo, yhi) )
    
    ax1.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax1, axis_thickness)
    #ax1.set_xticks([])
    ax1.set_yticks([])
    
    ax2.set_xlabel('{}'.format(temp_label), fontsize=fs, color='black')
    ax2.set_ylabel('{}'.format(volume_label), fontsize=fs, color='black')
    ax2.set_xlim( (xlo, xhi) )
    ax2.set_ylim( (ylo, yhi) )
    
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    set_axis_thickness(ax2, axis_thickness)
    #ax2.set_xticks([])
    ax2.set_yticks([])
    
    ax1.text(*label_rel_pos, '(a) - Unbounded', transform=ax1.transAxes, fontsize=fs*legend_fs_scale, fontweight='bold', va='top', ha='left')
    ax2.text(*label_rel_pos, '(b) - Bounded to RFR-regions', transform=ax2.transAxes, fontsize=fs*legend_fs_scale, fontweight='bold', va='top', ha='left')
    
    # After plotting everything:
    handles, labels = [], []
    for ax in (ax1, ax2):#, ax3, ax4, ax5):
        ax_handles, ax_labels = ax.get_legend_handles_labels()
        for handle, label in zip(ax_handles, ax_labels):
            if label not in labels:  # avoid duplicate legend entries
                handles.append(handle)
                labels.append(label)
    
    # Combined horizontal extent of ax1 and ax2
    pos1 = ax1.get_position()
    pos2 = ax2.get_position()
    
    legend_left = pos1.x0
    legend_width = pos2.x1 - pos1.x0
    
    fig.legend(
        handles,
        labels,
        loc='lower left',
        bbox_to_anchor=(
            legend_left,     # left edge
            0.79,            # bottom edge
            legend_width,    # total width
            0.10,            # available height
        ),
        bbox_transform=fig.transFigure,
        mode='expand',
        ncol=2,
        borderaxespad=0,
        fancybox=True,
        fontsize=legend_fs_scale * fs,
    )
    
    # Plot Tg results
    labels    = np.array(['RFR\nthermal', 'Hyperbola',       'Piecewise\nregression'])
    bounded   = np.array([np.nan,         h_bounded_tg[0],   pr_bounded_tg[0]]) # RFR index-0 does not have bounded test
    unbounded = np.array([rfr_tg[0],      h_unbounded_tg[0], pr_unbounded_tg[0]])
    colors = np.array([regression_color,  patrone_color,     pr_color])
    
    x = np.arange(len(labels))
    width = 0.35
    has_bounded = ~np.isnan(bounded)
    
    # Paired bars are offset; a lone unbounded bar is centered
    bounded_x = x[has_bounded] - width / 2
    unbounded_x = np.where(has_bounded, x + width / 2, x)
    
    ax3.bar(bounded_x, bounded[has_bounded], width, color=colors[has_bounded], edgecolor='black', linewidth=0.8, hatch='///')
    ax3.bar(unbounded_x, unbounded, width, color=colors, edgecolor='black', linewidth=0.8, hatch='\\\\\\')
    ax3.axhline(unbounded[0], color=regression_color, ls='--', lw=lw/5, alpha=1)
    
    ax3.set_xticks(x)
    ax3.set_xticklabels(labels)
    ax3.set_ylabel('$T_g$ (K)', fontsize=fs, color='black')
    
    tlo = 0
    thi = 1.5*max([ np.nanmax(bounded), np.nanmax(unbounded)  ])
    ax3.set_ylim( (tlo, thi) )
    
    method_handles = [
        Patch(facecolor='white', edgecolor='black', hatch='///', label='Bounded'),
        Patch(facecolor='white', edgecolor='black', hatch='\\\\\\', label='Unbounded')]
    ax3.legend(handles=method_handles, fontsize=legend_fs_scale*fs*0.95)
    
    ax3.tick_params(axis='both', which='major', labelsize=fs/1.35)
    set_axis_thickness(ax3, axis_thickness)
    ax3.text(*label_rel_pos, '(c)', transform=ax3.transAxes, fontsize=fs*legend_fs_scale, fontweight='bold', va='top', ha='left')
    
    
    # Glassy CTE's
    labels    = np.array(['RFR\nthermal', 'Hyperbola',           'Piecewise\nregression'])
    bounded   = np.array([np.nan,         h_bounded_glassy_cte,   pr_bounded_glassy_cte]) # RFR index-0 does not have bounded test
    unbounded = np.array([rfr_glassy_cte, h_unbounded_glassy_cte, pr_unbounded_glassy_cte])
    colors = np.array([regression_color,  patrone_color,          pr_color])
    # print(bounded)
    # print(unbounded)
    
    x = np.arange(len(labels))
    width = 0.35
    has_bounded = ~np.isnan(bounded)
    
    # Paired bars are offset; a lone unbounded bar is centered
    bounded_x = x[has_bounded] - width / 2
    unbounded_x = np.where(has_bounded, x + width / 2, x)
    
    ax4.bar(bounded_x, bounded[has_bounded], width, color=colors[has_bounded], edgecolor='black', linewidth=0.8, hatch='///')
    ax4.bar(unbounded_x, unbounded, width, color=colors, edgecolor='black', linewidth=0.8, hatch='\\\\\\')
    ax4.axhline(unbounded[0], color=regression_color, ls='--', lw=lw/5, alpha=1)
    
    ax4.set_xticks(x)
    ax4.set_xticklabels(labels)
    ax4.set_ylabel('Glassy CTE ($K^{-1}$)', fontsize=fs, color='black')
    
    tlo = 0
    thi = 1.5*max([ np.nanmax(bounded), np.nanmax(unbounded)  ])
    ax4.set_ylim( (tlo, thi) )
    
    method_handles = [
        Patch(facecolor='white', edgecolor='black', hatch='///', label='Bounded'),
        Patch(facecolor='white', edgecolor='black', hatch='\\\\\\', label='Unbounded')]
    ax4.legend(handles=method_handles, fontsize=legend_fs_scale*fs*0.95)
    
    ax4.tick_params(axis='both', which='major', labelsize=fs/1.35)
    set_axis_thickness(ax4, axis_thickness)
    ax4.text(*label_rel_pos, '(d)', transform=ax4.transAxes, fontsize=fs*legend_fs_scale, fontweight='bold', va='top', ha='left')
    
    # Rubbery CTE's
    labels    = np.array(['RFR\nthermal', 'Hyperbola',             'Piecewise\nregression'])
    bounded   = np.array([np.nan,          h_bounded_rubbery_cte,   pr_bounded_rubbery_cte]) # RFR index-0 does not have bounded test
    unbounded = np.array([rfr_rubbery_cte, h_unbounded_rubbery_cte, pr_unbounded_rubbery_cte])
    colors = np.array([regression_color,   patrone_color,           pr_color])
    
    x = np.arange(len(labels))
    width = 0.35
    has_bounded = ~np.isnan(bounded)
    
    # Paired bars are offset; a lone unbounded bar is centered
    bounded_x = x[has_bounded] - width / 2
    unbounded_x = np.where(has_bounded, x + width / 2, x)
    
    ax5.bar(bounded_x, bounded[has_bounded], width, color=colors[has_bounded], edgecolor='black', linewidth=0.8, hatch='///')
    ax5.bar(unbounded_x, unbounded, width, color=colors, edgecolor='black', linewidth=0.8, hatch='\\\\\\')
    ax5.axhline(unbounded[0], color=regression_color, ls='--', lw=lw/5, alpha=1)
    
    ax5.set_xticks(x)
    ax5.set_xticklabels(labels)
    ax5.set_ylabel('Rubbery CTE ($K^{-1}$)', fontsize=fs, color='black')
    
    tlo = 0
    thi = 1.5*max([ np.nanmax(bounded), np.nanmax(unbounded)  ])
    ax5.set_ylim( (tlo, thi) )
    
    method_handles = [
        Patch(facecolor='white', edgecolor='black', hatch='///', label='Bounded'),
        Patch(facecolor='white', edgecolor='black', hatch='\\\\\\', label='Unbounded')]
    ax5.legend(handles=method_handles, fontsize=legend_fs_scale*fs*0.95)
    
    ax5.tick_params(axis='both', which='major', labelsize=fs/1.35)
    set_axis_thickness(ax5, axis_thickness)
    ax5.text(*label_rel_pos, '(e)', transform=ax5.transAxes, fontsize=fs*legend_fs_scale, fontweight='bold', va='top', ha='left')



    fig.savefig(
        figname + '_published_7.jpeg',
        dpi=dpi,
        bbox_inches='tight',
    )
    
    return None

##########################################################################
# Function for reducing dense data to sparser data for plotting purposes #
##########################################################################
def reduce_xy_for_plot(xdata, ydata, plot_dx):
    """
    Reduce discrete x/y data for plotting by keeping the point nearest each
    plot_dx-spaced x location.

    Parameters
    ----------
    xdata, ydata : array-like
        Original discrete data.
    plot_dx : float
        Desired spacing of plotted x-values.

    Returns
    -------
    xplot, yplot : ndarray
        Reduced plotting arrays.
    """
    xdata = np.asarray(xdata)
    ydata = np.asarray(ydata)
    
    dx = np.median(np.diff(xdata))
    if dx < plot_dx:
        stride = max(1, int(round(plot_dx / dx)))
        xplot = xdata[::stride]
        yplot = ydata[::stride]
        print('Reducing plot resolution: {} -> {} points (median dx = {:.5g}, plotting every {} points)'.format(
              len(xdata), len(xplot), dx, stride))
    else:
        xplot = xdata
        yplot = ydata
    return xplot, yplot