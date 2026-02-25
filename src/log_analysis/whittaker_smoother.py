# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
October 19th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    ****************************************************************************
    * Requirements:                                                            *
    *   python whittaker-eilers module (for testing only within this script):  *
    *    - pip3 install whittaker-eilers (if pip manager is installed)         *
    *                                                                          *
    *   python numpy module (for using this script as a module):               *
    *    - pip3 install numpy (if pip manager is installed)                    *
    *                                                                          *
    *   python scipy module (for using this script as a module):               *
    *    - pip3 install scipy (if pip manager is installed)                    *
    *                                                                          *
    *   python matplotlib module (for using this script as a module):          *
    *    - pip3 install matplotlib (if pip manager is installed)               *
    ****************************************************************************
    
Please see for more details:
    https://pubs.acs.org/doi/10.1021/ac034173t
    https://numpy.org/devdocs/user/numpy-for-matlab-users.html
    https://scipy.github.io/old-wiki/pages/NumPy_for_Matlab_Users.html
    https://towardsdatascience.com/how-to-tune-the-perfect-smoother-bcc5a67660b1
    https://towardsdatascience.com/the-perfect-way-to-smooth-your-noisy-data-4f3fe6b44440
"""
##############################
# Import Necessary Libraries #
##############################
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import math



############################################################################################
# Function implementating base Whittaker Smoother to call from PyPi (for testing purpose)  #
#   ydata = Y-data to smoooth                                                              #
#   order = order of smoother (number of data points smoother uses - 2 is a good default)  #
#   lmbda = controls "smoothness" can be an integer or 'op' for optimize                   #
############################################################################################
def PyPi_smooth(ydata, order, lmbda):
    from whittaker_eilers import WhittakerSmoother
    # Automatic lambda optimization
    if str(lmbda).startswith('op'):
        smoother = WhittakerSmoother(lmbda=1, order=order, data_length=len(ydata))
        results = smoother.smooth_optimal(ydata, break_serial_correlation=False)
        optimal_lambda = results.get_optimal().get_lambda()
        smoothed = results.get_optimal().get_smoothed()
        #cve = results.get_optimal().get_cross_validation_error()
    else: # Manual control
        smoother = WhittakerSmoother(lmbda=lmbda, order=order, data_length=len(ydata))
        smoothed = smoother.smooth(ydata) 
        optimal_lambda = None
        #cve = None    
    return smoothed, optimal_lambda


###########################################################
# Function to build D = (m-d)*N sparse difference matrix  #
# (drop-in replacement for MATLab's diff() funciton)      #
###########################################################
def sparse_eye_diff(m, d, format='csc'):
    diagonals = np.zeros(2*d + 1)
    diagonals[d] = 1
    for i in range(d):
        diagonals = diagonals[:-1] - diagonals[1:]
    offsets = np.arange(d+1)
    return sp.sparse.diags(diagonals, offsets, shape=(m-d, m), format=format)


################################################################################################################
# Function implementing Whittaker-Eilers smoothing function w/o interopolation. The parameters are as follows: #
#    y = ydata to smooth in numpy array                                                                        #
#    d = integer order of the smoother                                                                         #
#    lmbda = lmbda float smoothing constant                                                                    #
#    compute_cve = True or False Boolean to set whether or not to compute the Cross-validation error (CVE)     #
#    cve_mode = a string of 'numpy' 'scipy' or 'fast' to set how CVE is being computed. The meaning of each:   #
#      'numpy' will compute the full h-matrix by converting the sparse C-matrix to a fully matrix and then     #
#              take the inverse.                                                                               #
#      'scipy' will use the sparse C-matrix to inverse. The 'numpy' mode is quicker for small data sets, but   #
#              uses a lot more memory as fully matrices are used. I would recommend 'scipy' as the default due #
#              to both time and space constraints.                                                             #
#      'fast'  will renormalize the y-data to be multiples of length n (set in the function) to compute the    #
#              partial h-matrix and scale that up to get the full h-matrix. Eilers default is to use fast      #
#              when the length of y-data is larger then 1000.                                                  #
################################################################################################################
def Whittaker_Eilers_without_interpolation(x, y, d, lmbda, compute_cve=False, cve_mode='scipy', nevery=1):
    # If nevery > 1 and compute_cve we are trying to break serial correlated data (when data is correlated
    # with the lagged version of itself). Eiliers suggests sampling nevery data points in such a case, by
    # weighting different nevery points
    if compute_cve and nevery > 1:
        w = np.zeros_like(y)
        w[::nevery] = 1.0
        y = w*y
        #x = x[::nevery]
        #y = y[::nevery]
    
    # If we want to compute the CVE in 'fast' mode, we need to re-normalize y to be multiples of n
    if compute_cve and cve_mode == 'fast':
        n = 50 # set the number of n-points for "smaller" h-matrix. Eiler's suggests 100.
        mn = np.floor(len(y)/n)
        mn = int(mn*n)
        y = y[:mn]
        x = x[:mn]
    
    # Base Whittaker Eilers smoother
    m = len(y)
    E = sp.sparse.eye(m, dtype=int, format='csc')
    D = sparse_eye_diff(m, d, format='csc')
    C = E + lmbda*D.conj().T.dot(D)
    z = sp.sparse.linalg.spsolve(C, y)
    
    # Computation of hat diagonal and cross-validation.
    if compute_cve:
        if cve_mode == 'numpy':
            C = C.todense().T
            H = np.linalg.inv(C)
            h = np.diagonal(H)
        if cve_mode == 'scipy':
            H = sp.sparse.linalg.inv(C)
            h = sp.sparse.csr_matrix.diagonal(H)
        if cve_mode == 'fast':
            E1 = sp.sparse.eye(n, dtype='int32', format='csc')
            D1 = sparse_eye_diff(n, d, format='csc')
            lambda1 = lmbda*(n/m)**(2*d)
            H1 = sp.sparse.linalg.inv(E1 + lambda1*D1.conj().T.dot(D1))
            h1 = sp.sparse.csr_matrix.diagonal(H1)
            u = np.zeros(m)
            k = int(np.floor(m/2) - 1)
            k1 = int(np.floor(n/2) - 1)
            u[k] = 1
            v = sp.sparse.linalg.spsolve(C, u)
            f = np.round( (np.arange(0,m) - 1)*(n - 1)/(m - 1)).astype(int)
            h = h1[f]*v[k]/h1[k1]
        r = (y - z)/(1 - h)
        s = np.sum(r**2)
        cve = np.sqrt(s/m)
    else: cve = None
    return x, z, cve


##########################################################################################
# Function to parse out MinLambda, MaxLambda, and NumLambda from string formatted like:  #
# lmbda = 'op<MinLambda,MaxLambda,NumLambda>' or 'op<MinLambda,MaxLambda,NumLambda>-p'   #
# where MinLambda is a float, MaxLambda is a float and NumLambda is an integer.          #
##########################################################################################
def parse_lambda_settings(lmbda):
    # Set defaults
    min_lambda = ''; max_lambda = ''; num_lambda = ''
    
    # Start parsing
    starting_char = '<'; ending_char = '>'
    starting_flag = False; ending_flag = False
    string = ''
    for i in lmbda:
        if i == ' ': continue
        if i == starting_char:
            starting_flag = True
            continue
        if i == ending_char: 
            ending_flag = True
            continue
        if starting_flag and not ending_flag:
            string += i
    
    # Split string to get values
    values = string.split(',')
    if len(values) >= 3:
        min_lambda = float(values[0])
        max_lambda = float(values[1])
        num_lambda = int(float(values[2]))
        successful_parse = True
    else: successful_parse = False
    return min_lambda, max_lambda, num_lambda, successful_parse, values


#################################################################
# Function to increment lambda's, compute cve, and minimize cve #
#################################################################
def Whittaker_Eilers_optimize_lambda(x, y, d, lmbda_method, xlabel='', ylabel='', basename='', colors=None):
    
    #--------------------------------------------------------------------------------------#
    # Check for user values based on lmbda_method = 'op<MinLambda,MaxLambda,NumLambda>' or #
    # 'op<MinLambda,MaxLambda,NumLambda>-p' and let user know. If not set defaults and let #
    # user know what defaults are being used for lambda spacing.                           #
    #--------------------------------------------------------------------------------------#
    successful_parse = False
    if '<' in lmbda_method and '>' in lmbda_method and lmbda_method.count(',') >= 2:
        min_lambda, max_lambda, num_lambda, successful_parse, values = parse_lambda_settings(lmbda_method)
        nevery = 1
        print(f'User supplied inputs to set lambda spacing and range via {lmbda_method}. Spacing parameters:')
        print('{:>12}: {}'.format('MinLambda', min_lambda))
        print('{:>12}: {}'.format('MaxLambda', max_lambda))
        print('{:>12}: {}'.format('NumLambda', num_lambda))  
        if len(values) == 4:
            try:
                nevery = int(values[3])
                print('{:>12}: {}'.format('Nevery',    nevery)) 
            except: pass
        if not successful_parse:
            print(f'WARNING could not parse {lmbda_method}. Please ensure format is "op<MinLambda,MaxLambda,NumLambda>" or "op<MinLambda,MaxLambda,NumLambda>-p"')
    if not successful_parse:
        min_lambda = 1e-2
        max_lambda = 1e12
        num_lambda = 50
        nevery = 1
        print('Using defaults to set lambda spacing. Spacing parameters:')
        print('{:>12}: {}'.format('MinLambda', min_lambda))
        print('{:>12}: {}'.format('MaxLambda', max_lambda))
        print('{:>12}: {}'.format('NumLambda', num_lambda)) 
        print('{:>12}: {}'.format('Nevery',    nevery)) 
        
    #-----------------------------------------------------------#
    # Check for errors in min_lamba, max_lambda, and multiplier #
    #-----------------------------------------------------------#
    if min_lambda >= max_lambda:
        raise Exception(f'ERROR MinLambda >= MaxLambda. {min_lambda} >= {max_lambda}.')
    if min_lambda < 0:
        raise Exception(f'ERROR MinLambda < 0. MinLambda = {min_lambda} and must be greater than 0.')
    if num_lambda < 2:
        raise Exception(f'ERROR NumLambda < 2. NumLambda = {num_lambda} and must be greater than 2.')
    
    #----------------------------------------------------------------------------------------#
    # Set the mode in which to compute the Cross-validation error (CVE) based on length of   #
    # y-data. Eilers recommends anything larger then 1000, be computed with the 'fast' mode. #
    #----------------------------------------------------------------------------------------#
    if len(y) > 500:
        cve_mode = 'fast'
        print(f'\n\nCross-validation mode is set to "{cve_mode}" because length of y-data is > 500. Length of data {len(y)}.')
    else:
        cve_mode = 'scipy'
        print(f'\n\nCross-validation mode is set to "{cve_mode}" because length of y-data is <= 500. Length of data {len(y)}.')      
        
    #-------------------#
    # Plotting settings #
    #-------------------#
    fs = 12
    legend_fs_scale = 0.8
    label_rel_pos = (0.005, 0.99)
    label_rel_pos = ()
    
    #------------------------------------------------------------------#
    # Check for nevery == 0, if so attempt to break serial-correlation #
    #------------------------------------------------------------------#
    if nevery == 0:
        n_min, n_max, n_inc = 1, 10, 1
        serial = {'x':   [],
                  'z':   [],
                  'n':   [],
                  'cve': [],
                  'lambda': [],
                  'min_diff': []}
        break_count = 0 # Count of break out conditions
        nover       = 2 # Number of times we went over the break out conditions
        for i, n in enumerate(range(n_min, n_max+1, n_inc)):
            lmbda_method_loop = 'op<{}, {}, {}, {}>'.format(min_lambda, max_lambda, num_lambda, n)
            lambda_loop, fig_loop = Whittaker_Eilers_optimize_lambda(x, y, d, lmbda_method_loop)
            x_loop, z_loop, cve_loop = Whittaker_Eilers_without_interpolation(x, y, d, lambda_loop, compute_cve=True, cve_mode=cve_mode, nevery=n)
            min_diff = lambda_loop - min_lambda
            serial['x'].append( x_loop )
            serial['z'].append( z_loop )
            serial['n'].append( n )
            serial['cve'].append( cve_loop )
            serial['lambda'].append( lambda_loop )
            serial['min_diff'].append( lambda_loop - min_lambda )
            
            # Setup break out conditions
            if i > 0:
                if min_diff > min(serial['min_diff']):
                    break_count += 1
            if break_count > nover: break
            
            
        # Find the nevery, that bumps the computed lambda
        # off of the lower bound for the first time.
        min_diff = np.array(serial['min_diff'])
        possible_indexes = np.where(min_diff > 0)[0]
        if len(possible_indexes) == 0:
            index = 0
        else: index = min(possible_indexes)
        opt_nevery  = serial['n'][index]
        opt_mindiff = serial['min_diff'][index]
        
        # Re-run with the newly choosen nevery
        lmbda_method_final = 'op<{}, {}, {}, {}>'.format(min_lambda, max_lambda, num_lambda, opt_nevery)
        if lmbda_method.endswith('-p'): lmbda_method_final += '-p'
        opt_lambda, fig = Whittaker_Eilers_optimize_lambda(x, y, d, lmbda_method_final, xlabel=xlabel, ylabel=ylabel, basename=basename, colors=colors)
        
        # Add a row to the bottom of the plot for this analysis
        old_axes     = fig.axes
        ncols, nrows = 2, 2
        
        # Resize figure to make room
        w, h = fig.get_size_inches()
        fig.set_size_inches(w, h*1.5)
        gs = gridspec.GridSpec(nrows + 1, ncols, figure=fig)
        
        # Reposition old axes into new grid
        for i, ax in enumerate(old_axes):
            row = i // ncols
            col = i % ncols
            ax.set_position(gs[row, col].get_position(fig))
            ax.set_subplotspec(gs[row, col])
        
        # Add new spanning subplot
        ax5 = fig.add_subplot(gs[nrows, :])        
        # ax5.set_yscale('log')
        # ax5.yaxis.set_major_formatter(ScalarFormatter())
        ax5.set_title('Automatically breaking serial correlation')
        
        ax5.plot(serial['n'], serial['min_diff'], '-o', lw=5, ms=10, color='tab:blue', label='')
        ax5.axhline(y=0, ls='--', lw=3, zorder=1, color='black', alpha=1.0, label='Zero-line')
        ax5.plot(opt_nevery, opt_mindiff, '*', ms=15, color='tab:orange', label='Selected nevery={}'.format(opt_nevery))

        
        ax5.set_xlabel('nevery', fontsize=fs)
        ax5.set_ylabel(r'$\lambda_{CVE, optimized}$ - $\lambda_{Min, bound}$', fontsize=fs)
        ax5.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs)
        if label_rel_pos: ax5.text(*label_rel_pos, '(f)', transform=ax5.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')

        return opt_lambda, fig
    
    #-------------------------------------------------------#
    # Normal optimization that does not do an nevery search #
    #-------------------------------------------------------#
    else:
        #------------------------------#
        # Compute course CVE "spectra" #
        #------------------------------#
        # Compute the course Cross-validation error (CVE) for each lamba
        if min_lambda > 0:
            course_lambdas = np.geomspace(min_lambda, max_lambda, num=num_lambda, endpoint=True, dtype='float64', axis=0)
        else:
            course_lambdas = np.geomspace(1e-16, max_lambda, num=num_lambda, endpoint=True, dtype='float64', axis=0)
            course_lambdas = np.insert(course_lambdas, 0, 0)
        print(f'\nComputing Course Cross-validation error for {len(course_lambdas)} different lambdas ...')
        course_cves, course_zs = [], []
        for n, lmbda in enumerate(course_lambdas):
            x_out, z, cve = Whittaker_Eilers_without_interpolation(x, y, d, lmbda, compute_cve=True, cve_mode=cve_mode, nevery=nevery)
            course_cves.append(cve)
            course_zs.append(z)
            
        # Find minimum cve and use as optimized lambda
        course_minimum_index = course_cves.index(min(course_cves))
        course_optimized_cve = course_cves[course_minimum_index]
        course_optimized_lambda = course_lambdas[course_minimum_index]
        
        #----------------------------#
        # Compute fine CVE "spectra" #
        #----------------------------#
        # Find left/right indexes and reset to bounds if needed
        lo_index = course_minimum_index - 1
        hi_index = course_minimum_index + 1
        if lo_index <= 0: 
            lo_index = 0
            print('  WARNING course minimum is near the start of lambda testing, decrease MinLambda.')     
        if hi_index >= len(course_lambdas)-1:
            hi_index = len(course_lambdas)-1
            print('  WARNING course minimum is near the end of lambda testing, increase MaxLambda.') 
        
        # Reset the min_lambda and max_lamdba values and generate lo/hi fine_lambdas based on
        # len(course_lambdas) to perform a "moving from course minimum outward search"
        min_lambda_fine = course_lambdas[lo_index] 
        max_lambda_fine = course_lambdas[hi_index] 
        lo_fine_lambdas = np.flip(np.linspace(min_lambda_fine, course_optimized_lambda, len(course_lambdas))[:-1])
        hi_fine_lambdas = np.linspace(max_lambda_fine, course_optimized_lambda, len(course_lambdas))[:-1]
        
        # Compute the fine Cross-validation error (CVE) for each lambda on
        lo_pass = False; hi_pass = False; npoints = len(course_lambdas)-1; fine_lambdas_cves = {} # {lambda:cve}
        print(f'\nComputing Fine Cross-validation error for {2*npoints} different lambdas moving from course minimum outward ...');
        for n in range(npoints):
            # Compute CVE on both sides of the minimum
            lo_lambda = lo_fine_lambdas[n]; hi_lambda = hi_fine_lambdas[n]
            lo_x, lo_z, lo_cve = Whittaker_Eilers_without_interpolation(x, y, d, lo_lambda, compute_cve=True, cve_mode=cve_mode, nevery=nevery)
            lo_y, hi_z, hi_cve = Whittaker_Eilers_without_interpolation(x, y, d, hi_lambda, compute_cve=True, cve_mode=cve_mode, nevery=nevery)
            fine_lambdas_cves[lo_lambda] = lo_cve
            fine_lambdas_cves[hi_lambda] = hi_cve
            
            # Update when lo_cve and hi_cve get larger then course_optimized_cve
            if lo_cve >= course_optimized_cve: lo_pass = True
            if hi_cve >= course_optimized_cve: hi_pass = True
            
            # If both lo_pass and hi_pass we have started increasing on 
            # both sides of the course minimum and break out of this loop
            #if lo_pass and hi_pass: break 
        
        # Sort fine_lambdas_cves in ascending order based on lambda for cleaner plotting
        fine_lambdas_cves = dict(sorted(fine_lambdas_cves.items()))
        fine_lambdas = list(fine_lambdas_cves.keys())
        fine_cves = list(fine_lambdas_cves.values())
            
        # Find minimum cve and use as optimized lambda
        fine_minimum_index = fine_cves.index(min(fine_cves))
        fine_optimized_cve = fine_cves[fine_minimum_index]
        fine_optimized_lambda = fine_lambdas[fine_minimum_index]
            
        # Plot the CVE vs lambda plot
        fig = None
        if str(lmbda_method).endswith('-p'):
            if isinstance(colors, (tuple , list)):
                colors = colors
            else:
                # Color wheel defined by matplotlib:
                #   colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
                # However, we can construct our own color wheel to prioritize the colors we want first and we can 
                # have way more colors defined than what matplotlib defines. Below is Josh's preferred color wheel
                colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple', 'tab:red', 'tab:gray','tab:olive', 'tab:cyan', 'tab:pink', 'teal',
                          'crimson', 'lime', 'tomato',  'blue', 'orange', 'green', 'purple', 'red', 'gray', 'olive', 'cyan', 'pink', 'tab:brown']
            
            # Function to walk around the color wheel defined by colors and get the next color,
            # if the color index is exceeds the color wheel, it will reset the color index to 0
            def walk_colors(color_index, colors):
                color = colors[color_index]
                color_index += 1
                if color_index + 1 > len(colors): color_index = 0
                return color, color_index
            
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(2*5, 2*4))
            
            ax1.semilogx(course_lambdas, course_cves, ls='-', lw=4, color='tab:blue', zorder=0, label=r"CVEs based on course $\lambda$'s")
            ax1.semilogx(fine_lambdas, fine_cves, ls='-', lw=4, color='tab:green', zorder=0, label=r"CVEs based on fine $\lambda$'s  ")
            ax1.plot(course_optimized_lambda, course_optimized_cve, 'o', mfc='tab:cyan', mec='black', ms=12, lw=3, zorder=2, label='Course Minimum CVE (x={:.2f}, y={:.2f}; nevery={})'.format(course_optimized_lambda, course_optimized_cve, nevery))
            ax1.plot(fine_optimized_lambda, fine_optimized_cve, 'o', mfc='lime', mec='black', ms=12, lw=3, zorder=2, label='Fine Minimum CVE (x={:.2f}, y={:.2f}; nevery={})'.format(fine_optimized_lambda, fine_optimized_cve, nevery))
            ax1.set_xlabel(r'$\lambda$', fontsize=fs)
            ax1.set_ylabel('Cross-validation error (CVE)', fontsize=fs)
            if label_rel_pos: ax1.text(*label_rel_pos, '(a)', transform=ax1.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
            
            #        (center, height, span, height)
            anchor = (0.0,    1.1,    2.2,  0.1)
            ax1.legend(bbox_to_anchor=anchor, loc=2, ncol=2, mode='expand', fontsize=legend_fs_scale*fs)
            #ax1.xaxis.set_major_formatter(ScalarFormatter())
            
            if lo_index > 1: lo_zoomedex = lo_index - 1
            else: lo_zoomedex = 0
            if hi_index < len(course_lambdas)-1:  hi_zoomedex = hi_index + 1
            else: hi_zoomedex = len(course_lambdas)-1
            ax2.semilogx(course_lambdas[lo_zoomedex:hi_zoomedex+1], course_cves[lo_zoomedex:hi_zoomedex+1], ls='-', lw=4, color='tab:blue', label=r"CVEs based on course $\lambda$'s")
            ax2.semilogx(fine_lambdas, fine_cves, ls='-', lw=4, color='tab:green', label=r"CVEs based on fine $\lambda$'s  ")
            ax2.plot(course_optimized_lambda, course_optimized_cve, 'o', mfc='tab:cyan', mec='black', ms=12, lw=3, label='Course Minimum CVE (x={:.2f}; y={:.2f}; nevery={})'.format(course_optimized_lambda, course_optimized_cve, nevery))
            ax2.plot(fine_optimized_lambda, fine_optimized_cve, 'o', mfc='lime', mec='black', ms=12, lw=3, label='Fine Minimum CVE (x={:.2f}; y={:.2f}; nevery={})'.format(fine_optimized_lambda, fine_optimized_cve, nevery))
            ax2.set_xlabel(r'$\lambda$', fontsize=fs)
            ax2.set_ylabel('Cross-validation error (CVE)', fontsize=fs)
            if label_rel_pos: ax2.text(*label_rel_pos, '(b)', transform=ax2.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
            #ax2.xaxis.set_major_formatter(ScalarFormatter())
            
            nplotted = 8 # maximum number of plotted lambdas
            modulus = math.floor(len(course_lambdas)/nplotted)
            color_index = 0
            ax3.plot(x, y, '-', lw=5, color='tab:blue', label='Raw data')
            for n, (lamda, z, cve) in enumerate(zip(course_lambdas, course_zs, course_cves), 1):
                if n%modulus == 0:
                    course_x, course_z, course_cve = Whittaker_Eilers_without_interpolation(x, y, d, lamda, compute_cve=False, cve_mode=cve_mode)
                    label = '{}={:.1e}'.format(r'$\lambda$', lamda)
                    color, color_index = walk_colors(color_index, colors)
                    ax1.plot(lamda, cve, 's', ms=8, mec='black', lw=3, color=color, zorder=1)
                    ax3.plot(x, course_z, '-', lw=2, color=color, label=label)
            ax3.set_xlabel('{}'.format(xlabel), fontsize=fs)
            ax3.set_ylabel('{}'.format(ylabel), fontsize=fs)
            ax3.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs)
            if label_rel_pos: ax3.text(*label_rel_pos, '(c)', transform=ax3.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
            
            fine_x, fine_z, fine_cve = Whittaker_Eilers_without_interpolation(x, y, d, fine_optimized_lambda, compute_cve=False, cve_mode=cve_mode)
            ax4.plot(x, y, '-', lw=5, color='tab:blue', label='Raw data')
            ax4.plot(x, fine_z, '-', lw=2, color='lime', label='Fine Minimum CVE Z-series')
            ax4.set_xlabel('{}'.format(xlabel), fontsize=fs)
            ax4.set_ylabel('{}'.format(ylabel), fontsize=fs)
            ax4.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0), fancybox=True, ncol=1, fontsize=legend_fs_scale*fs)
            if label_rel_pos: ax4.text(*label_rel_pos, '(d)', transform=ax4.transAxes, fontsize=fs, fontweight='bold', va='top', ha='left')
                    
    
            fig.tight_layout()
            if basename:
                figname = basename+'_WE_LOO-CVE.jpeg'
                print(f'Rendering {figname}')
                fig.savefig(figname, dpi=300)
        return fine_optimized_lambda, fig


############################################################################################
# Function implementating base Whittaker Smoother to call from LUNAR/log_analysis          #
#   ydata = Y-data to smoooth                                                              #
#   order = order of smoother (number of data points smoother uses - 2 is a good default)  #
#   lmbda = controls "smoothness" can be an integer or 'op' for optimize                   #
############################################################################################
def smooth(xdata, ydata, order, lmbda, xlabel='', ylabel='', basename='', colors=None):
    x = np.array(xdata.copy())
    y = np.array(ydata.copy())
    # Automatic lambda optimization
    if str(lmbda).startswith('op'):
        optimal_lambda, fig = Whittaker_Eilers_optimize_lambda(x, y, order, lmbda, xlabel=xlabel, ylabel=ylabel, basename=basename, colors=colors)
        lmbda = optimal_lambda 
    else: optimal_lambda, fig = None, None
    x_out, smoothed, cve = Whittaker_Eilers_without_interpolation(x, y, order, lmbda, compute_cve=False)
    return smoothed, optimal_lambda, fig




########################
# Testing the smoother #
########################
if __name__ == "__main__": 
    # Seed for reproducability
    seed = 1
    
    # Generate "real" data
    x = np.linspace(0, 0.1, 1200)
    y_real = (1/(1+np.exp(-x*50))-0.5)*2000
    
    # Generate "noisy" data with an underlying sine wave
    if seed > 0: np.random.seed(seed)
    sine_wave = (np.mean(y_real)/25)*np.sin(250*(2*np.pi*x)) # Sine wave with amplitude 1/25th the mean of y_real at 250 hz
    y_noisy = y_real + np.random.normal(0, 100, len(x)) + sine_wave
    
    # Smooth "noisy" data using Josh's version. lmbda can be used to set the
    # MinLambda, MaxLambda, and NumLambda value for the range of Lambda's to
    # test. For example:
    #    'op<MinLambda,MaxLambda,NumLambda>'
    # where the MinLambda will be is the minimum lamba to test, MaxLambda is the
    # maximum lambda to test and NumLambda sets the number of equally spaced 
    # lambda's to test on a log scale. The '-p' flag can be added to plot the
    # CVE "spectra", for example: 'op<MinLambda,MaxLambda,NumLambda>-p'.
    order = 2
    in1_lmbda = 1_000_000 # Manually take control of lambda
    in1_lmbda = 'op' # use default spacing and do not plot
    in1_lmbda = 'op-p' # use default spacing and plot
    in1_lmbda = 'op<1e-2, 1e10, 50>-p' # Adjust spacing and plot
    y_smooth, out1_lmbda = smooth(y_noisy, order, in1_lmbda)
    if out1_lmbda is None: out1_lmbda = in1_lmbda
    
    # Smooth "noisy" data using PyPi version
    order = 2
    in2_lmbda = 'op'
    y_PyPi, out2_lmbda = PyPi_smooth(y_noisy, order, in2_lmbda)
    if out2_lmbda is None: out2_lmbda = in2_lmbda
    
    # Plot the data
    fig, ax = plt.subplots()
    ax.plot(x, y_noisy,'tab:blue', lw=2, label='Noisy y')
    ax.plot(x, y_real, 'tab:orange', lw=4, label='Real y')
    ax.plot(x, y_smooth, 'tab:cyan', ls='--', lw=3, label='Joshs implementation (order={}; lambda={})'.format(order, out1_lmbda))
    ax.plot(x, y_PyPi, 'tab:green', ls='-.', lw=2, label='PyPi implementation (order={}; lambda={})'.format(order, out2_lmbda))
    ax.legend()