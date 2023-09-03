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
    The T12 data is processed 
    to extract the distribution parameters p, lambda and sigma.

    Returns
    -------
    None.

   """

   input_data = [
[-3,0.039,-0.0285,0.001,0.99982],
[-2.5,0.051,-0.036,0.001,0.99973],
[-2,0.068,-0.047,0.001,0.99988],
[-1.5,0.095,-0.0625,0.001,0.99981],
[-1,0.14,-0.086,0,0.99991],
[-0.5,0.223,-0.1145,-0.004,0.99993],
[0,0.332,0,-0.014,0.99403],
[0.5,0.223,0.1145,-0.004,0.99994],
[1,0.14,0.086,0,0.99992],
[1.5,0.095,0.0625,0.001,0.99989],
[2,0.068,0.047,0.001,0.99989],
[2.5,0.051,0.036,0.001,0.9997],
[3,0.039,0.0285,0.001,0.99982]
            ]
   output_array = np.asarray(input_data) 
   
   print(output_array.shape)
   xdata = output_array[:,0] 
   ydata = output_array[:,1] 
   popt = my_functions.do_fit(my_functions.exp_gauss, xdata, ydata,
                              bounds=([0.65, 0.8, 5.9], [0.70, 1.1, 6.2])) #0.9
   print('param', popt)
   #[0.59999999 1.01975752 4.99994054]

main()