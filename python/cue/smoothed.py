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
   data = np.concatenate([data,col0], axis=1)
   print(data.shape) 
   rev = np.flip(data, axis=1)
   print(data[0,:])
   print(rev[0,:])

main()