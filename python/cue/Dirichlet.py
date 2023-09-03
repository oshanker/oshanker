#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 15:01:10 2023

@author: uorugant
"""

import my_functions
import numpy as np

def main():
   """
    The Dirichlet data is processed 
    to extract the distribution parameters p, lambda and sigma.
    T = 1.6E9

    Returns
    -------
    None.

   """

   input_data = [
[-4, 0.0245, -0.021, 0.001],
[-3, 0.041, -0.0335, 0.001],
[-2, 0.0725, -0.056, 0.002],
[-1, 0.1475, -0.0995, 0.001],
[0, 0.325,  0.00, -0.024],
[1, 0.1475, 0.0995, 0.001],
[2, 0.0725, 0.056, 0.002],
[3, 0.041, 0.0335, 0.001],
[4, 0.0245, 0.021, 0.001],
            ]
   output_array = np.asarray(input_data) 
   
   print(output_array.shape)
   xdata = output_array[:,0] 
   ydata = output_array[:,1] 
   popt = my_functions.do_fit(my_functions.exp_gauss, xdata, ydata,
                              bounds=([0.65, 0.8, 5.8], [0.75, 1.1, 6.1])) #0.9
   print('param', popt)
   #[0.59999999 1.01975752 4.99994054]

main()