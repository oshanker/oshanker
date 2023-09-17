#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 22:37:20 2023

@author: uorugant
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit


gauss_norm = 1/(math.sqrt(2*math.pi))
def gauss(x, sigma ):
    """

    Parameters
    ----------
    x : arg
    
    sigma : std dev

    Returns
    -------
    normal distribution

    """
    A = gauss_norm/(sigma)
    return A*np.exp(-(x)**2/(2*sigma**2))

def exp_pdf_raw(x, lam):
    """
    

    Parameters
    ----------
    x : arg
    
    lam : constant for exponential distribution

    Returns
    -------
    exponential distribution

    """
    xx = np.abs(np.copy(x))
    return  np.exp(-lam * xx) * lam/2
    
def exp_gauss(x, p, lam, sigma):
    """    

    Parameters
    ----------
    x : arg
    
    p : coefficient for exponential distribution
    
    lam : constant for exponential distribution
    
    sigma :  std dev

    Returns
    -------
    ret : linear combination of exponential distribution and normal distribution.

    """
    ret = p * exp_pdf_raw(x, lam) + (1-p) * gauss(x, sigma )
    return ret



def exp_pdf(x, lam, corr):
    """
    

    Parameters
    ----------
    x : arg
    
    lam : constant for exponential distribution

    Returns
    -------
    exponential distribution

    """
    xx = np.abs(np.copy(x))
    return  corr*np.exp(-lam * xx) * lam/2

def do_fit(func, xdata, ydata, bounds=([0.7, 0.8], [3.3, 2.0])):
    print(np.sum(ydata))
    popt, pcov = curve_fit(func, xdata, ydata, bounds=bounds)
    print('param', popt)
    cond = np.linalg.cond(pcov)
    print("cond", cond)
    print("diag cov", np.diag(pcov))
    plt.plot(xdata, ydata, 'b-', label='data')
    plt.plot(xdata, func(xdata, *popt), 'g--',
         label='fit' )
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    plt.show()
    return popt

def do_plot_func(func, popt, xdata, label='func', style='g--'):
    plt.plot(xdata, func(xdata, *popt), style,
         label=label )
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid()
    plt.legend()
    #plt.show()
    return popt

def do_plot_func_with_coeff(func, coeff, popt, xdata, label='func', style='g--'):
    plt.plot(xdata, coeff*func(xdata, *popt), style,
         label=label )
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid()
    plt.legend()
    #plt.show()
    return popt

def main():
    bins_in = []
    for i in np.arange(-16.125, 16.5, 0.05):
        bins_in.append(i)
    
    xaxis_zero = 322
    xaxis = []
    for i in range(0, len(bins_in)-1):
        xaxis.append((bins_in[i]+bins_in[i+1])/2)
    print(len(xaxis))
    size = 300000
    lam = 1.0
    data1 = np.random.exponential(scale=1.0/lam, size=size)
    data2 = np.random.exponential(scale=1.0/lam, size=size)
    sigma = 2.0
    data3 = np.random.normal(loc=0.0, scale=sigma, size=size)
    data = np.concatenate([data1, -data2, data3])
    #data = np.concatenate([data1, -data2])
    hist, bins_range = np.histogram(data, bins=bins_in, density=True)
    histnorm = np.sum(hist)
    print('np.sum(hist)',  histnorm)
    print('np.var(data, ddof=0)',  np.var(data, ddof=0))
    xdata = np.array(xaxis)
    p = 2/3;
    print('p*2/(lam*lam) ', p*2/(lam*lam) , " + (1-p)*sigma*sigma ", (1-p)*sigma*sigma
          , ' = ', p*2/(lam*lam) + (1-p)*sigma*sigma)
    print('<x**2> (from hist)',  np.sum(hist*xdata*xdata)/histnorm)
    #popt = do_fit(exp_pdf, xdata, hist)
    # bounds=([0.5, 0.5, 1.3], [1.0, 4.7, 3.0])
    # popt = do_fit(exp_gauss, xdata, hist, bounds=bounds) #0.9
    #do_plot_func(der_gauss, popt, x, 'der')
    xx = []
    for i in range(-8, 9):
        xx.append(i)
    xx_array = np.array(xx)
    
    # popt = [0.22968019,  0.69744049,  7.73814023]
    # do_plot_func(exp_gauss, popt, xx_array, label='exp_gauss')
    # # popt_exp = [1.00132733, 0.666666667]
    # # do_plot_func(exp_pdf, popt_exp, xx_array, label='exp_pdf',style='r--')
    # popt_exp = [0.69744049 ]
    # do_plot_func_with_coeff(exp_pdf_raw, 0.22968019, popt_exp, xx_array, label='exp_pdf_with_coeff',style='r--')
    # popt_exp = [7.73814023 ]
    # do_plot_func_with_coeff(gauss, 1-0.22968019, popt_exp, xx_array, label='gauss_with_coeff',style='b--')

    popt = [0.61859396,  1.42590918,  9.4 ]
    do_plot_func(exp_gauss, popt, xx_array, label='exp_gauss')
    popt_exp = [1.42590918 ]
    do_plot_func_with_coeff(exp_pdf_raw, 0.61859396, popt_exp, xx_array, label='exp_pdf_with_coeff',style='r--')
    popt_exp = [9.4 ]
    do_plot_func_with_coeff(gauss, 1-0.61859396, popt_exp, xx_array, label='gauss_with_coeff',style='b--')
    
    plt.show()
    # simul = gauss(x, 2.0)
    # plt.plot(xaxis, simul, 'b-', label='data')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.legend()
    # plt.show()

main()
