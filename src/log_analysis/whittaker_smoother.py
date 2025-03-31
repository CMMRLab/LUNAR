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
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp



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
#              benchmarking different ways of computing the h-matrix.                                          #
#      'fast'  will renormalize the y-data to be multiples of length n (set in the function) to compute the    #
#              partial h-matrix and scale that up to get the full h-matrix. Eilers default is to use fast      #
#              when the length of y-data is larger then 1000.                                                  #
################################################################################################################
def Whittaker_Eilers_without_interpolation(y, d, lmbda, compute_cve=False, cve_mode='scipy'):
    # If we want to compute the CVE in 'fast' mode, we need to re-normalize y to be multiples of n
    if compute_cve and cve_mode == 'fast':
        n = 50 # set the number of n-points for "smaller" h-matrix. Eiler's suggests 100.
        mn = np.floor(len(y)/n)
        mn = int(mn*n)
        y = y[:mn]
    
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
    return z, cve


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
    if len(values) == 3:
        min_lambda = float(values[0])
        max_lambda = float(values[1])
        num_lambda = int(float(values[2]))
        successful_parse = True
    else: successful_parse = False
    return min_lambda, max_lambda, num_lambda, successful_parse


#################################################################
# Function to increment lambda's, compute cve, and minimize cve #
#################################################################
def Whittaker_Eilers_optimize_lambda(y, d, lmbda_method):
    from tqdm import tqdm
    
    #--------------------------------------------------------------------------------------#
    # Check for user values based on lmbda_method = 'op<MinLambda,MaxLambda,NumLambda>' or #
    # 'op<MinLambda,MaxLambda,NumLambda>-p' and let user know. If not set defaults and let #
    # user know what defaults are being used for lambda spacing.                           #
    #--------------------------------------------------------------------------------------#
    successful_parse = False
    if '<' in lmbda_method and '>' in lmbda_method and lmbda_method.count(',') == 2:
        min_lambda, max_lambda, num_lambda, successful_parse = parse_lambda_settings(lmbda_method)
        print(f'User supplied inputs to set lambda spacing and range via {lmbda_method}. Spacing parameters:')
        print('{:>12}: {}'.format('MinLambda', min_lambda))
        print('{:>12}: {}'.format('MaxLambda', max_lambda))
        print('{:>12}: {}'.format('NumLambda', num_lambda))  
        if not successful_parse:
            print(f'WARNING could not parse {lmbda_method}. Please ensure format is "op<MinLambda,MaxLambda,Multiplier>" or "op<MinLambda,MaxLambda,Multiplier>-p"')
    if not successful_parse:
        min_lambda = 1e-2
        max_lambda = 1e12
        num_lambda = 50
        print('Using defaults to set lambda spacing. Spacing parameters:')
        print('{:>12}: {}'.format('MinLambda', min_lambda))
        print('{:>12}: {}'.format('MaxLambda', max_lambda))
        print('{:>12}: {}'.format('NumLambda', num_lambda))        
    
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
    
    #------------------------------#
    # Compute course CVE "spectra" #
    #------------------------------#
    # Compute the course Cross-validation error (CVE) for each lamba
    course_lambdas = np.geomspace(min_lambda, max_lambda, num=num_lambda, endpoint=True, dtype='float64', axis=0)
    print(f'\nComputing Course Cross-validation error for {len(course_lambdas)} different lambdas ...'); course_cves = []
    for n, lmbda in enumerate(tqdm(course_lambdas)):
        z, cve = Whittaker_Eilers_without_interpolation(y, d, lmbda, compute_cve=True, cve_mode=cve_mode)
        course_cves.append(cve)
        
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
    min_lambda = course_lambdas[lo_index] 
    max_lambda = course_lambdas[hi_index] 
    lo_fine_lambdas = np.flip(np.linspace(min_lambda, course_optimized_lambda, len(course_lambdas))[:-1])
    hi_fine_lambdas = np.linspace(max_lambda, course_optimized_lambda, len(course_lambdas))[:-1]
    
    # Compute the fine Cross-validation error (CVE) for each lambda on
    lo_pass = False; hi_pass = False; npoints = len(course_lambdas)-1; fine_lambdas_cves = {} # {lambda:cve}
    print(f'\nComputing Fine Cross-validation error for {2*npoints} different lambdas moving from course minimum outward ...');
    for n in tqdm(range(npoints)):
        # Compute CVE on both sides of the minimum
        lo_lambda = lo_fine_lambdas[n]; hi_lambda = hi_fine_lambdas[n]
        lo_z, lo_cve = Whittaker_Eilers_without_interpolation(y, d, lo_lambda, compute_cve=True, cve_mode=cve_mode)
        hi_z, hi_cve = Whittaker_Eilers_without_interpolation(y, d, hi_lambda, compute_cve=True, cve_mode=cve_mode)
        fine_lambdas_cves[lo_lambda] = lo_cve
        fine_lambdas_cves[hi_lambda] = hi_cve
        
        # Update when lo_cve and hi_cve get larger then course_optimized_cve
        if lo_cve >= course_optimized_cve: lo_pass = True
        if hi_cve >= course_optimized_cve: hi_pass = True
        
        # If both lo_pass and hi_pass we have started increasing on 
        # both sides of the course minimum and break out of this loop
        if lo_pass and hi_pass: break 
    
    # Sort fine_lambdas_cves in ascending order based on lambda for cleaner plotting
    fine_lambdas_cves = dict(sorted(fine_lambdas_cves.items()))
    fine_lambdas = list(fine_lambdas_cves.keys())
    fine_cves = list(fine_lambdas_cves.values())
        
    # Find minimum cve and use as optimized lambda
    fine_minimum_index = fine_cves.index(min(fine_cves))
    fine_optimized_cve = fine_cves[fine_minimum_index]
    fine_optimized_lambda = fine_lambdas[fine_minimum_index]
        
    # Plot the CVE vs lambda plot
    if str(lmbda_method).endswith('-p'):
        fig, ax = plt.subplots()
        ax.semilogx(course_lambdas, course_cves, ls='-', lw=4, color='tab:blue', label="CVEs based on course $\lambda$'s")
        ax.semilogx(fine_lambdas, fine_cves, ls='-', lw=4, color='tab:green', label="CVEs based on fine $\lambda$'s  ")
        ax.plot(course_optimized_lambda, course_optimized_cve, 'o', mfc='tab:cyan', mec='black', ms=12, lw=3, label='Course Minimum CVE (x={:.4f}; y={:.4f})'.format(course_optimized_lambda, course_optimized_cve))
        ax.plot(fine_optimized_lambda, fine_optimized_cve, 'o', mfc='lime', mec='black', ms=12, lw=3, label='Fine Minimum CVE (x={:.4f}; y={:.4f})'.format(fine_optimized_lambda, fine_optimized_cve))
        ax.set_xlabel(r'$\lambda$', fontsize=12)
        ax.set_ylabel('Cross-validation error (CVE)', fontsize=12)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, ncol=2, fontsize=10)
    return fine_optimized_lambda


############################################################################################
# Function implementating base Whittaker Smoother to call from LUNAR/log_analysis          #
#   ydata = Y-data to smoooth                                                              #
#   order = order of smoother (number of data points smoother uses - 2 is a good default)  #
#   lmbda = controls "smoothness" can be an integer or 'op' for optimize                   #
############################################################################################
def smooth(ydata, order, lmbda):
    ydata = np.array(ydata.copy())
    # Automatic lambda optimization
    if str(lmbda).startswith('op'):
        optimal_lambda = Whittaker_Eilers_optimize_lambda(ydata, order, lmbda)
        lmbda = optimal_lambda 
    else: optimal_lambda = None
    smoothed, cve = Whittaker_Eilers_without_interpolation(ydata, order, lmbda, False)
    return smoothed, optimal_lambda




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
    in1_lmbda = 'op<1e-2, 1e10, 100>-p' # Adjust spacing and plot
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