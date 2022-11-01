#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 17:37:27 2022

@author: uorugant
"""
import sys
import numpy as np 
from tensorflow import keras
from tensorflow.keras import layers

import pandas
import matplotlib.pyplot as plt

def main():
    print("-------")
    print(sys.argv)

    def getdata():
        print("check alignment")
        dataset = pandas.read_csv('../../oldriemann/data/zetaE12.csv', header=0)
        dataset.drop(dataset.columns[0], axis=1, inplace=True)
        
        raw_data = dataset.values[1:]
        print('raw_data.shape', raw_data.shape)
        print('raw_data[0]', raw_data[0])
        print(raw_data[:5])
        
        dataset1 = pandas.read_csv('../../oldriemann/data/zerosE12.csv.max', header=None)
    
        temperature = dataset1.values
        print('temperature[0]', temperature[0])
        return raw_data, temperature
    raw_data, temperature = getdata()
    print('temperature shape', temperature.shape)

    
if __name__   == '__main__':
     main()