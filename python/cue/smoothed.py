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
   bounds=([0.3, 2.0, 1.5], [0.4, 2.5, 1.8])
   popt = my_functions.do_fit(my_functions.exp_gauss, np.asarray(zdf), ydata, bounds=bounds) #0.9
   print('param', popt)


main()