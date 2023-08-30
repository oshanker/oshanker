#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 10:24:20 2023

@author: uorugant
"""


import flip_and_switch
import my_functions
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def exp_pdf(x, lam):
    xx = np.copy(x)
    for i in range(0, xx.shape[0]):
        if(xx[i] < 0):
            xx[i] = -xx[i]
    return  np.exp(-lam * xx) * lam/2
    
def do_exp_fit(xdata, ydata):
    print(np.sum(ydata))
    popt, pcov = curve_fit(exp_pdf, xdata, ydata, bounds=([1], [3]))
    print(popt)
    cond = np.linalg.cond(pcov)
    print(cond)
    print(np.diag(pcov))
    plt.plot(xdata, ydata, 'b-', label='data')
    plt.plot(xdata, exp_pdf(xdata, *popt), 'g--',
         label='fit: lam=%5.3f' % tuple(popt))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()


def main():
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
   do_exp_fit(zdf, output_array[:,1]) 


main()