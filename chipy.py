# -*- coding: utf-8 -*-
"""
Siddharth Sule
April 2020
--------------------------------------------------
ChiPy - A Module for Undergraduate Laboratory work
--------------------------------------------------

This module is designed to do all the required work for undergraduate lab :
    
    - Apply fits to data and judge its quality #INCOMPLETE
    - Plot data in laboratory style #NOT STARTED
    - Propagate errors #NOT STARTED
"""
#------------------------------------------------------------------------------
#IMPORTED MODULES
import numpy as np
#------------------------------------------------------------------------------
def linear_fit(x, y, yerr):
    """
    Parameters
    ----------
    x : float 1D array, x values of the experiment.
    y : float 1d array, y values of the experiment.
    yerr : float 1D array, error in y values of the experiment.

    Returns
    -------
    p : float 1D array, linear fit y = p[0]*x + p[1].
    perr : float 1D array, error on elements of p.
    reduced_chi_sqr : float, reduced chi squared value of fit.
    
    Description
    -----------
    This function applies a linear fit to the data. It returns the data for the
    linear fit and the errors on the parameters, as well as the reduced chi
    squared value of the fit.
    """
    p, c = np.polyfit(x, y, 1, w = 1/yerr, cov = True)
    perr = np.array([np.sqrt(c[0,0]), np.sqrt(c[1,1])])
    
    f = lambda x : (p[0] * x) + p[1]
    
    chi_sqr = 0
    
    for i in range(len(x)):
        
        chi_sqr += ((y[i] - f(x[i]))/(yerr[i]))**2
    
    reduced_chi_sqr = chi_sqr / (len(x) - 2)
    
    return p, perr, reduced_chi_sqr
#------------------------------------------------------------------------------
def quadratic_fit(x, y, yerr):
    """
    Parameters
    ----------
    x : float 1D array, x values of the experiment.
    y : float 1d array, y values of the experiment.
    yerr : float 1D array, error in y values of the experiment.

    Returns
    -------
    p : float 1D array, linear fit y = p[0]*x**2 + p[1]*x + p[2].
    perr : float 1D array, error on elements of p.
    reduced_chi_sqr : float, reduced chi squared value of fit.
    
    Description
    -----------
    This function applies a quadratic fit to the data. It returns the data for 
    the quadratic fit and the errors on the parameters, as well as the reduced 
    chi squared value of the fit.
    """
    p, c = np.polyfit(x, y, 2, w = 1/yerr, cov = True)
    perr = np.array([np.sqrt(c[0,0]), np.sqrt(c[1,1]), np.sqrt(c[2,2])])
    
    f = lambda x : (p[0] * x**2) + (p[1] * x) + p[2]
    
    chi_sqr = 0
    
    for i in range(len(x)):
        
        chi_sqr += ((y[i] - f(x[i]))/(yerr[i]))**2
    
    reduced_chi_sqr = chi_sqr / (len(x) - 3)
    
    return p, perr, reduced_chi_sqr
#------------------------------------------------------------------------------
def cubic_fit(x, y, yerr):
    """
    Parameters
    ----------
    x : float 1D array, x values of the experiment.
    y : float 1d array, y values of the experiment.
    yerr : float 1D array, error in y values of the experiment.

    Returns
    -------
    p : float 1D array, linear fit y = p[0]*x**3 + p[1]*x**2 + p[2]*x + p[3].
    perr : float 1D array, error on elements of p.
    reduced_chi_sqr : float, reduced chi squared value of fit.
    
    Description
    -----------
    This function applies a cubic fit to the data. It returns the data for 
    the cubic fit and the errors on the parameters, as well as the reduced chi
    squared value of the fit.
    """
    p, c = np.polyfit(x, y, 3, w = 1/yerr, cov = True)
    perr = np.array([np.sqrt(c[0,0]), np.sqrt(c[1,1]), 
                     np.sqrt(c[2,2]), np.sqrt(c[3,3])])
    
    f = lambda x : (p[0] * x**3) + (p[1] * x**2) + (p[2] * x) + p[3]
    
    chi_sqr = 0
    
    for i in range(len(x)):
        
        chi_sqr += ((y[i] - f(x[i]))/(yerr[i]))**2
    
    reduced_chi_sqr = chi_sqr / (len(x) - 4)
    
    return p, perr, reduced_chi_sqr
#------------------------------------------------------------------------------
def exponential_fit(x, y, yerr):
    """
    Parameters
    ----------
    x : float 1D array, x values of the experiment.
    y : float 1d array, y values of the experiment.
    yerr : float 1D array, error in y values of the experiment.

    Returns
    -------
    e : float 1D array, exponential fit y = A*exp(B*x)
    eerr : float 1D array, error on elements of e.
    reduced_chi_sqr : float, reduced chi squared value of fit.
    
    Description
    -----------
    This function applies np.log to the data, and then applies a linear fit. 
    It returns the data for the exponential fit and the errors on the 
    parameters, as well as the reduced chi squared value of the fit.
    """
    #FIT : y = Ae^(Bx) => NATURAL LOG : ln(y) = ln(A) + Bx
    ln_y = np.log(y)
    
    # ERROR IN LN(Y) = YERR / Y
    ln_yerr = yerr / y
    
    p, c = np.polyfit(x, ln_y, 1, w = 1/ln_yerr, cov = True)
    perr = np.array([np.sqrt(c[0,0]), np.sqrt(c[1,1])])
    
    f = lambda x : (p[0] * x) + p[1]
    
    chi_sqr = 0
    
    for i in range(len(x)):
        
        chi_sqr += ((y[i] - f(x[i]))/(yerr[i]))**2
    
    reduced_chi_sqr = chi_sqr / (len(x) - 2)
    
    #RETURN A AND B
    e = np.array([np.exp(p[1]), p[0]])
    eerr = np.array([(perr[1] * np.exp(p[1])), perr[0]])
    
    return e, eerr, reduced_chi_sqr
#------------------------------------------------------------------------------    