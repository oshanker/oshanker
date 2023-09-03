#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 10:24:20 2023

    The simulation output from characteristic.py ('../out/fitcue.txt') is processed 
    by smoothed.py to extract the distribution parameters p, lambda and sigma.

@author: uorugant
"""


import flip_and_switch
import my_functions
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math

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

def exp_pdf(x, lam):
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
    ret = p * exp_pdf(x, lam) + (1-p) * gauss(x, sigma )
    return ret

def do_fit(func, xdata, ydata):
    """
    fit function to data

    Parameters
    ----------
    func : function to fit
    
    xdata : independent variables
    
    ydata : dependent variable (distribution to be fitted)

    Returns
    -------
    popt : fitted parameters 

    """
    print(np.sum(ydata))
    #popt, pcov = curve_fit(func, xdata, ydata, bounds=([0.01], [2.5]))
    
    #p, lam, sigma
    popt, pcov = curve_fit(func, xdata, ydata, bounds=([0.28, 3.6, 1.5], [0.31, 3.7, 1.52]))
    # [0.3        3.55578876 1.52809399]
    # [0.29117862 3.64662775 1.50870849]
    # [0.29117866 3.64662728 1.50870858]
    
    print('param', popt)
    cond = np.linalg.cond(pcov)
    print("cond", cond)
    print("diag cov", np.diag(pcov))
    plt.plot(xdata, ydata, 'b-', label='data')
    plt.plot(xdata, func(xdata, *popt), 'g--',
         label='fit: ' )
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()
    return popt


def main():
   """
    The simulation output from characteristic.py ('../out/fitcue.txt') is processed 
    by smoothed.py to extract the distribution parameters p, lambda and sigma.

    Returns
    -------
    None.

   """
   file_name =  '../out/rawcue.txt'
   data = np.loadtxt(file_name)
   n = int(data[-1][0]) 
   sample_size = int(data[-1][1]) 
   #-2.0 to 2.0, values only
   data = data[0:-1,:]
   rows = data.shape[0]
   cols = data.shape[1]
   col0 = np.copy(data[:,0:1])
   zdf = []
   for i in np.arange(-1.0, 1.05, 0.05):
        zdf.append(np.around(i, decimals=2))
   
   smooth = True
   if smooth:
       data = np.concatenate([data,col0], axis=1)
       rev = np.flip(data, axis=1)
       print(data[0,:])
       print(rev[0,:])
       data = (data+rev)/2
       data = flip_and_switch.flip_and_switch(data)
   
   print(data.shape) 
   testfit = []
   phi_values = []
   for j in range(0, cols+1):
        phi_values.append(int(j*360/(cols)))
        
   print(phi_values, len(phi_values))

   file_name = '../out/smoothedfitcue'+ str(n) + '.txt'
   header = "n " + str(n) + " sample " + str(sample_size)    

   for j in range(0, rows):
        out = my_functions.fit(zdf[j], data[j], phi_values)
        testfit.append(out)
   output_array = np.asarray(testfit) 
   np.savetxt( file_name, 
               output_array, 
               fmt='%4.2f %7.3f %7.3f %8.4f %7.3f', 
               header=header,
               delimiter=',')
   ydata = output_array[:,1] 
   #do_fit(exp_pdf, np.asarray(zdf), ydata) 1.17
   #do_fit(gauss, np.asarray(zdf), ydata) #0.9
   popt = do_fit(exp_gauss, np.asarray(zdf), ydata) #0.9
   print('param', popt)


main()