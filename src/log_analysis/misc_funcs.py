# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
July 28, 2026
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
import numpy as np
import math
import csv


#############################################
# Josh's hand built linear regression model #
#############################################
class linear_regression:
    def __init__(self, xdata, ydata):
        # Check if xdata and ydate exists. If it does not generate some defaults
        if not xdata and not ydata:
            print('WARNING no data for linear regression. Setting xdata=[0.0, 0.5, 1.0]; ydata=[0.0, 0.0, 0.0] to avoid div by zero errors - BUT DONT TRUST OUTPUT.')
            xdata = [0.0, 0.5, 1.0]
            ydata = [0.0, 0.0, 0.0]
        
        # Compute cummulative parameters
        sum_xi = 0; sum_yi = 0; sum_xi_2 = 0; sum_yi_2 = 0; sum_xi_yi = 0; n = 0
        for x, y in zip(xdata, ydata):
            sum_xi += x; sum_yi += y; n += 1
            sum_xi_2 += x*x; sum_yi_2 += y*y; 
            sum_xi_yi += x*y
        if n == 0: n = 1 # Avoid div by zero errors
                
        # SSxy and Sxx
        self.SSxy = (sum_xi_yi) - ((1/n)*(sum_xi*sum_yi))
        self.SSxx = sum_xi_2 - (1/n)*sum_xi*sum_xi
        
        # y bar and x bar
        self.x_bar = sum_xi/n
        self.y_bar = sum_yi/n 

        # b0 and b1
        self.b1 = self.SSxy/self.SSxx
        self.b0 = self.y_bar - self.b1*self.x_bar
        self.n = n
    
        # Errors in residuals and regressions
        self.SStotal = sum_yi_2 - (sum_yi*sum_yi/n)
        self.SSreg = self.b1*self.b1*self.SSxx
        self.SSres = self.SStotal - (self.b1*self.b1*self.SSxx)
        self.MSres = self.SSres/(n - 2)

        # variances and r^2
        if self.SSxx != 0:
            self.b0_variance = self.MSres*( (1/n) + (self.x_bar*self.x_bar/self.SSxx) )
            self.b1_variance = self.MSres/self.SSxx
            self.r2 = self.SSreg/self.SStotal
        else:
            self.b0_variance = 0
            self.b1_variance = 0
            self.r2 = 0
        
        # linear regression data
        self.xreg = [min(xdata), max(xdata)]
        self.yreg = [i*self.b1 + self.b0 for i in self.xreg]
        
    # Method for computing confidence interval for b0
    def confidence_interval_b0(self, alpha):
        from scipy import stats    
        df_res = self.n - 2
        t_statistic = stats.t.ppf(alpha, df_res)
        std_dev = math.sqrt(self.b0_variance)
        plus_minus = t_statistic*std_dev
        return sorted([self.b0-plus_minus, self.b0+plus_minus])
    
    # Method for computing confidence interval for b1
    def confidence_interval_b1(self, alpha):
        from scipy import stats
        df_res = self.n - 2
        t_statistic = stats.t.ppf(alpha, df_res)
        std_dev = math.sqrt(self.b1_variance)
        plus_minus = t_statistic*std_dev
        return sorted([self.b1-plus_minus, self.b1+plus_minus])
    
    
#########################################
# Function implementing LOWESS smoother #
#########################################
def lowess(xdata, ydata, fraction=0.2, max_iter=10):    
    import statsmodels.api as sm
    out = sm.nonparametric.lowess(np.array(ydata), np.array(xdata), frac=fraction, it=max_iter, is_sorted=False)   
    xout = out[:, 0]
    yout = out[:, 1]
    return xout, yout


###################################
# Josh's avg/moving avg functions #
###################################
def moving_average(x, w):
  return np.convolve(x, np.ones(w)/w, 'valid')

def avg(lst):
    return sum(lst)/len(lst)


###############################
# Josh's reduce data function #
###############################
def reduce_data(xdata, ydata, xlo, xhi):
    xlo, xhi = sorted([xlo, xhi])
    if xlo > xhi:
        print(f'ERROR (reduce_data) no data between {xlo} - {xhi}')
    reducedx = []; reducedy = [];  
    for x, y in zip(xdata, ydata):
        if x >= xlo and x <= xhi:
            reducedx.append(x); reducedy.append(y);
    return reducedx, reducedy


####################################################################
# Function to sort independant variable X-list in ascending order  #
# and to re-construct Y-list based on the newly sorted X-list. For #
# example:                                                         #
#   x = [0, 1, 2, 2, 4, 3]                                         #
#   y = [0, 1, 2, 3, 4, 5]                                         #
#                                                                  #
# Perform increasing independant sort:                             #
#   newx, newy, index_backmap = increasing_independant_sort(x, y)  #
#   newx = [0, 1, 2, 2, 3, 4]                                      #
#   newy = [0, 1, 2, 3, 5, 4]                                      #
#                                                                  #
# Reconstruct x and y data from  newx, newy, and index_backmap:    #
#   reconx = []; recony = []                                       #
#   for i in index_backmap:                                        #
#       reconx.append(newx[i])                                     #
#       recony.append(newy[i])                                     #
#                                                                  #
#   reconx = [0, 1, 2, 2, 4, 3]                                    #
#   reconx = [0, 1, 2, 3, 4, 5]                                    #
####################################################################
def increasing_independant_sort(xdata, ydata):
    paired = {(n, i):j for n, (i, j) in enumerate(zip(xdata, ydata))}
    paired = dict(sorted(paired.items(), key=lambda x:abs(x[0][1]) )) # [0=keys;1=values][index position in tuple-key]
    newx = []; newy = []; index_backmap = [] # {}
    for n, key in enumerate(paired):
        index_backmap.append(key[0])
        value = paired[key]
        newx.append(key[1])
        newy.append(value)
    return newx, newy, index_backmap


###################################################################################
# Function to walk along x- and y-data and remove any x-data that backtracks. For #
# example:                                                                        #
#   x = [0, 1, 2, 2, 4, 3]                                                        #
#   y = [0, 1, 2, 3, 4, 5]                                                        #
#                                                                                 #
# Perform cleaning:                                                               #
#   newx, newy, removed_indexes = increasing_independant_cleaning(x, y)           #
#   newx = [0, 1, 2, 4]                                                           #
#   newy = [0, 1, 2, 4]                                                           #
###################################################################################
def increasing_independant_cleaning(xdata, ydata):
    newx = [xdata[0]]; newy = [ydata[0]]; 
    removed_indexes = []; passed_x = xdata[0]
    for n, (x, y) in enumerate(zip(xdata, ydata)):
        if x > passed_x:
            newx.append(x)
            newy.append(y)
            passed_x = x
        else: removed_indexes.append(n)
    return newx, newy, removed_indexes
    

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
            ix.append( x[i] )
            iy.append(  summation )
    else: print('ERROR (compute_integral) inconsistent number of data points between X and Y arrays')
    return ix, iy


##################################################################################################################################
# Function to compute a first and second derivative using central finite difference.                                             #
# https://www.mathworks.com/matlabcentral/answers/494553-first-and-second-order-central-difference                               #
# https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter20.02-Finite-Difference-Approximating-Derivatives.html #
##################################################################################################################################
def compute_derivative(xdata, ydata):
    dxn = xdata[1:-1] # ignore first and last point
    dy1 = [0]*len(dxn)
    dy2 = [0]*len(dxn)
    if len(xdata) == len(ydata):
        for i in range(1, len(xdata)-1):
            dx = (xdata[i+1] - xdata[i-1])/2
            if dx == 0:
                print('WARNING finite difference dx was zero at x={}. Derivative was set to zero to avoid infinite derivative.'.format(xdata[i]))
            else:
                dy1[i-1] = (ydata[i+1] - ydata[i-1])/(2*dx)
                dy2[i-1] = (ydata[i+1] - 2*ydata[i] + ydata[i-1])/(dx*dx) 
    else: print('ERROR (compute_derivative) inconsistent number of data points between X and Y arrays')
    return dxn, dy1, dy2


############################################################
# Function to compute derivatives based on polynomial fits #
############################################################
def poly_derivative(xdata, ydata, span=2, deg=3):
    xdata = np.array(xdata); ydata = np.array(ydata)
    npoints = min([len(xdata), len(ydata)]) - span
    drxn = []; dry1 = []; dry2 = []
    for i in range(span, npoints):
        px = xdata[i-span:i+span+1]
        py = ydata[i-span:i+span+1]
        
        coefficients = np.polynomial.polynomial.polyfit(px, py, deg)
        x0 = xdata[i]; y0 = 0
        yp1 = 0; yp2 = 0
        for n in range(len(coefficients)):
            y0 += coefficients[n]*x0**n
            if n >= 1: 
                yp1 += coefficients[n]*n*x0**(n-1)
            if n >= 2: 
                yp2 += coefficients[n]*n*(n-1)*x0**(n-2)
    
        drxn.append(x0)
        dry1.append(yp1)
        dry2.append(yp2)           
    return drxn, dry1, dry2


#############################################################
# This was is Josh's method for computing euler's curvature #
#############################################################
def compute_euler_curvature(xdata, ydata):
    npoints = min([len(xdata), len(ydata)]) - 3
    xcurvature = []; ycurvature = []
    for i in range(1, npoints):
        x0, x1, x2, x3 = xdata[i:i+4]
        y0, y1, y2, y3 = ydata[i:i+4]
        
        # Compute the tangent vectors using central difference
        t1 = [x2-x0, y2-y0]
        t2 = [x3-x1, y3-y1]
        
        # Normalize the tangent vectors
        t1_magnitude = math.sqrt(t1[0]*t1[0] + t1[1]*t1[1])
        t2_magnitude = math.sqrt(t2[0]*t2[0] + t2[1]*t2[1])
        h1_hat = [t1[0]/t1_magnitude, t1[1]/t1_magnitude]
        h2_hat = [t2[0]/t2_magnitude, t2[1]/t2_magnitude]
        
        # Compute the change in tangent vectors "t3"
        t3 = [h1_hat[0]-h2_hat[0], h1_hat[1]-h2_hat[1]]
        t3_magnitude = math.sqrt(t3[0]*t3[0] + t3[1]*t3[1])
        
        # Set the direction of the of the t3 vector
        if t1[1] > t2[1]:
            direction = -1
        else:
            direction = 1
        
        # Compute the arc length
        ds = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        
        # Finally compute kappa
        kappa = direction*(t3_magnitude/ds)
        
        # Log values
        xcurvature.append(x2)
        ycurvature.append(kappa)
    return xcurvature, ycurvature


##########################################
# Function to compute Newton's curvature #
##########################################
def compute_newton_curvature(xdata, ydata, method='signed', span=1, deg=1):
    if span == 1 or deg == 1:
        print('Computing derivatives using central difference')
        dxn, dy1, dy2 = compute_derivative(xdata, ydata)
    else:
        print(f'Computing derivatives polyfit with span={span} and deg={deg}')
        dxn, dy1, dy2 = poly_derivative(xdata, ydata, span, deg)
    yp1 = np.array(dy1); yp2 = np.array(dy2)
    if method == 'signed':
        kappa = yp2/((1 + yp1*yp1)**(3/2))
    elif method == 'unsigned':
        kappa = np.abs(yp2)/((1 + yp1*yp1)**(3/2))
    else: raise Exception(f'ERROR unsupported method {method}. Currently supported methods: "signed" or "unsigned"')
    return dxn, kappa


##########################################################
# Function to find zero crossing (intial version came    #
# from Tristan Muzzy - This version is fully vectorized) #
##########################################################
def value_crossing(xdata, ydata, yvalue=0, style='low'):
    # Convert to numpy arrays
    x = np.array(xdata)
    y = np.array(ydata) - yvalue
        
    # Use Tristan's method of signs to find a zero crossing
    s = np.sign(y)
    flips = s[1:]*s[:-1]
    lo = np.where(flips == -1)[0]
    hi = lo[lo <= np.size(x)-2] + 1
    if style == 'low':
        x_cross = x[lo]
        y_cross = y[lo] + yvalue
    elif style == 'high':
        x_cross = x[hi]
        y_cross = y[hi] + yvalue
    elif style == 'interp':
        x_cross = -y[lo]*( (x[hi] - x[lo])/(y[hi] - y[lo]) ) + x[lo]
        y_cross = np.zeros_like(x_cross) + yvalue
    else:
        raise Exception(f'ERROR style {style} not recognized. Currently supported styles are "low" or "high" or "interp".')
    return list(x_cross), list(y_cross)


######################################################################
# Function to find peaks and valleys in a data set using derivatives #
######################################################################
def find_peaks_and_valleys(xdata, ydata):
    # Compute derivatives and find critical points
    dxn, dy1, dy2 = compute_derivative(xdata, ydata)
    critical_x, critical_y = value_crossing(dxn, dy1, yvalue=0, style='low')
    
    # Use 1st and 2nd derivative test to determine if the point is a peak or a valley
    xpeaks = []; ypeaks = []; xvalleys = []; yvalleys = []
    for x, y in zip(critical_x, critical_y):
        i = xdata.index(x)
        if dy2[dxn.index(x)] <= 0:
            xpeaks.append(xdata[i])
            ypeaks.append(ydata[i])
        else:
            xvalleys.append(xdata[i])
            yvalleys.append(ydata[i])
    return xpeaks, ypeaks, xvalleys, yvalleys


#########################################################################
# Function to fit y = c0 + c1*x + ... + c_n*x^n to x and y numpy arrays #
#########################################################################
def np_curve_fit(x, y, degree=2, domain=None):
    # Define polynomial function
    def poly(x, coef):
        yout = np.zeros_like(x)
        for n in range(len(coef)):
            yout += coef[n]*(x**n)
        return yout
    
    # Ensure X and Y are numpy arrays and find domain
    x = np.array(x); y = np.array(y)
    if domain is not None:
        indices = np.where( np.logical_and(x >= domain[0], x <= domain[1]) )[0]
        x = x[indices]; y = y[indices];
    else:
        x = x.copy(); y = y.copy();
    
    # Find best fit
    param = np.polynomial.polynomial.polyfit(x, y, degree)
    xfit = np.linspace(min(x), max(x), num=len(x))
    yfit = poly(xfit, param)
    return x, y, xfit, yfit, param


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
def edge_tan_intersect(x, y, N=1):
    '''
    Discrete derivative tangent lines at edge(N steps from)
    Returns analytical intersections between tangent lines
    '''
    yp = np.gradient(y, x)
    xi = (-y[N] + y[-N] + x[N] * yp[N] - x[-N] * yp[-N]) / (yp[N] - yp[-N])
    yi = yp[N] * (xi - x[N]) + y[N]
    return xi, yi

def fit_hyperbola(x, y, xlo, xhi, minimum_convergence=None, initial_guess=False, maxiter=10**4, slope_method='asymptote'):
    from scipy import optimize
      
    # define default outputs (in-case something goes wrong)
    xout = list(np.linspace(min(x), max(x), num=len(x)))
    yout = len(y)*[0];
    params = [0, 0, 0, 0, 0]; center = [0, 0];
    slopes = [0, 0]; transition = [0, 0];
    tangent_intersection = [0, 0]
    tangents = [(0, 0, 0), (0, 0, 0)] # {(x-points), (ypoints)}
    
    # Convert to float64 and then lists to numpy arrays
    xx = np.array(x, dtype=np.float64); yy = np.array(y, dtype=np.float64);
    
    # Setup intial guess
    p0 = None
    if initial_guess:
        slopeguess = (yy[-1]-yy[0])/(xx[-1]-xx[0])
        p0 = (np.mean(xx), np.mean(yy), slopeguess, slopeguess, np.log((xx[-1]-xx[0])**2/100))

    # Define the hyperbola equation (eqn 1 in paper)
    def hyberbola(t, t0, v0, a, b, c):
        h0 = 0.5*(t-t0) + np.sqrt( ((t-t0)*(t-t0))/4 + np.exp(c, dtype=np.float64) )
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
    
    # Find tangent intersection   
    if slope_method == 'asymptote':
        m1, m2 = slopes # slopes at lo-end and hi-end
    
    # Tristan said derivitive method is better?
    elif slope_method == 'tangent':
        dxn, dy1, dy2 = compute_derivative(xout, yout)
        m1 = dy1[0]  # numerical derivative at lo-end
        m2 = dy1[-1] # numerical derivative at hi-end
    else:
        raise Exception('ERROR slope_method needs to be "{}" or "{}", but is "{}"'.format('asymptote', 'tangent', slope_method))
    
    b1 = yout[0] - m1*xout[0] # Y-intercept at lo-end
    b2 = yout[-1] - m2*xout[-1] # Y-intercept at hi-end
    if m1 != m2:
        x_intersection = (b2 - b1) / (m1 - m2)
        y_intersection = m1*x_intersection + b1
        tangent_intersection = [x_intersection, y_intersection]
        
        # Upate tangent points for plotting
        tangents[0] = (xout[0], x_intersection, xout[-1])
        tangents[1] = (yout[0], y_intersection, yout[-1])
        
        
    
    # Find where function is transitioning (eqn 5 in paper)
    if minimum_convergence is not None:
        dtemp = (np.exp(c/2)*(2*minimum_convergence-1))/(np.sqrt(minimum_convergence*(1-minimum_convergence)))
        transition = [t0-dtemp, t0+dtemp]
    return xout, yout, params, center, slopes, transition, tangent_intersection, tangents


##########################################
# Function to find closest value in list #
##########################################
def closest(lst, value):
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-value))]  


######################################################################
# Statistic's functions to use why you dont want numpy as dependancy #
######################################################################
def compute_mean(data):
    if data:
        return sum(data)/len(data)
    else:
        return 0
 
def compute_variance(data):
    if data:
      mean = compute_mean(data)
      deviations = [(x - mean)**2 for x in data]
      variance = sum(deviations)/len(data)
      return variance
    else:
        return 0
 
def compute_standard_deviation(data):
  variance = compute_variance(data)
  return math.sqrt(variance)  


#########################################################################
# Function to apply a bin-average to effectively "resample" LAMMPS data #
#########################################################################
def bin_average(x, y, dx, start=None, stop=None):
    """Average y in bins of width dx along x.

    Returns bin centers and mean y values. Empty bins are NaN.
    The final bin may be narrower than dx.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError("x and y must be one-dimensional and the same length")
    if dx <= 0:
        raise ValueError("dx must be greater than zero")

    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) == 0:
        return np.array([]), np.array([])

    start = np.min(x) if start is None else start
    stop = np.max(x) if stop is None else stop
    if stop <= start:
        raise ValueError("stop must be greater than start")

    n_bins = int(np.ceil((stop - start) / dx))
    edges = start + np.arange(n_bins + 1) * dx
    edges[-1] = stop  # The final bin may be narrower than dx

    centers = (edges[:-1] + edges[1:]) / 2
    in_range = (x >= start) & (x < stop)

    counts, _ = np.histogram(x[in_range], bins=edges)
    sums, _ = np.histogram(x[in_range], bins=edges, weights=y[in_range])

    averages = np.full(n_bins, np.nan)
    occupied = counts > 0
    averages[occupied] = sums[occupied] / counts[occupied]

    return centers, averages


###################################################
# Function to interpolate data onto constant grid #
###################################################
def interp_constant_dx(xdata, ydata, dx, n_points):
    # Look for duplicate x-values
    unique = sorted(set(xdata))
    if len(xdata) != len(unique):
        print('WARNING interp-constant-dx found duplicate X-values. Grouping and averaging Y-values')
        groups = {i:[] for i in unique}
        for x, y in zip(xdata, ydata):
            groups[x].append(y)
            
        averages = {x:sum(ys)/len(ys) for x, ys in groups.items()}
        xdata    = [x for x in sorted(averages.keys())]
        ydata    = [averages[x] for x in xdata]        
        
        
    
    # Covnert to numpy arrays
    xdata = np.asarray(xdata, dtype=float)
    ydata = np.asarray(ydata, dtype=float)
    xmin  = min(xdata)
    xmax  = max(xdata)
    #print(dx, n_points)
    
    # Error check
    if isinstance(n_points, int) and isinstance(dx, (int, float)):
        raise Exception('ERROR np and dx can not be used simulationly')
    
    # Compute the desired number of points
    if dx is None and n_points is None:
        n_points = len(xdata)
    elif isinstance(dx, (int, float)):
        xspan    = xmax - xmin
        n_points = int(math.ceil(xspan/dx))
    elif isinstance(n_points, int):
        n_points = n_points
    
    # Dont let number of points dip below the orignal data
    if n_points < len(xdata):
        n_points = len(xdata)
        
    xinterp = np.linspace(xmin, xmax, n_points)    
    yinterp = np.interp(xinterp, xdata, ydata)
    return xinterp, yinterp
    


#######################################################
# Function to write Kemppainen-Muzzy data to csv file #
#######################################################
def savedata_to_csv(csv_data, csv_name):
    """
    Write data using either format:

        csv_data[label] = [xdata, ydata]
        csv_data[label] = series
        csv_data[label] = [series]       # Also supported

    Paired data receives "(x)" and "(y)" column labels.
    Standalone series use only the provided label.
    """

    def is_series(value):
        """Return True for non-string iterable objects."""
        if isinstance(value, (str, bytes)):
            return False

        try:
            iter(value)
            return True
        except TypeError:
            return False

    data_groups = []
    for label, data in csv_data.items():
        label = str(label)[:75]

        if not is_series(data):
            raise TypeError(f'CSV data associated with "{label}" must be iterable.')

        data = list(data)

        # Paired x/y data: [xdata, ydata]
        if len(data) == 2 and all(is_series(values) for values in data):
            xdata = list(data[0])
            ydata = list(data[1])

            columns = [(f'{label} (x)', xdata), (f'{label} (y)', ydata)]

        # Wrapped standalone series: [series]
        elif len(data) == 1 and is_series(data[0]):
            columns = [(label, list(data[0]))]

        # Direct standalone series: [value1, value2, ...]
        else:
            columns = [(label, data)]

        group_length = max(len(values) for _, values in columns)
        data_groups.append((group_length, columns))

    # Put the longest data groups first, matching the original behavior
    data_groups.sort(key=lambda group: group[0], reverse=True)

    # Flatten the groups into individual CSV columns
    columns = [
        column
        for _, group_columns in data_groups
        for column in group_columns
    ]

    nrows = max((len(values) for _, values in columns), default=0)
    with open(csv_name, 'w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)

        # Write column headings
        writer.writerow([title for title, _ in columns])

        # Write data, leaving blank cells for shorter columns
        for row_index in range(nrows):
            row = [
                values[row_index] if row_index < len(values) else ''
                for _, values in columns
            ]
            writer.writerow(row)
    return