#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 10:24:20 2023

@author: uorugant
"""

import pandas as pd
import my_functions
import numpy as np

def main():
   file_name =  '../out/rawcue.txt'
   data = np.loadtxt(file_name)
   #-2.0 to 2.0, values only
   data = data[0:-2,:]
   col0 = np.copy(data[:,0:1])
   
   smooth = True
   if smooth:
       data = np.concatenate([data,col0], axis=1)
       rev = np.flip(data, axis=1)
       print(data[0,:])
       print(rev[0,:])
       data = (data+rev)/2
   
   print(data.shape) 
   testfit = []
   rows = data.shape[0]
   cols = data.shape[1]
   phi_values = []
   for j in range(0, cols):
        phi_values.append(int(j*360/(cols-1)))
        
   print(phi_values)

   
   for j in range(0, rows):
        out = my_functions.fit(data[j], phi_values)
        testfit.append(out)
   np.savetxt('../out/smoothedfitcue.txt', 
               np.asarray(testfit), 
               fmt='%7.3f&%7.3f&%8.4f&%7.3f\\\\', 
               delimiter=',')

main()