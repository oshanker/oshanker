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
        print(dataset.head())

        dataset.drop(dataset.columns[0], axis=1, inplace=True)
        
        raw_data = dataset.values[1:]
        print('raw_data.shape', raw_data.shape)
        print('raw_data[0]', raw_data[0])
        print(raw_data[:5])
        
       # dataset1 = pandas.read_csv('../../oldriemann/data/zerosE12.csv.max', header=None)
    
        dataset1 = pandas.read_csv('../../oldriemann/data/zetaE12.csv', header=0)
        dataset1.drop(dataset1.columns[1], axis=1, inplace=True)
        temperature = dataset1.values[1:]
        print('temperature[0]', temperature[0])
        print(temperature[:5])
        return raw_data, temperature
    
    def example1():
        limit = 15
        time_sequence = np.arange(limit)  
        sequence_length = 3
        sampling_rate = 1
        delay = sampling_rate * (sequence_length + 1 - 1)
        batch_size = 3 
        
        raw_data = np.zeros((time_sequence.shape[0] ), dtype=int)  
        y_sequence = np.zeros((time_sequence.shape[0] ), dtype=int)  
        for i in time_sequence:
            raw_data[i] =  i
            if i >= sequence_length - 1:
                y_sequence[i] = np.max(raw_data[i+1-sequence_length:i+1])
        print(raw_data)
        

        print(y_sequence)
        dummy_dataset = keras.utils.timeseries_dataset_from_array(
            data=raw_data,                                 
            targets=y_sequence[delay-1:],                               
            sequence_length=sequence_length,                                      
            batch_size=batch_size,                                           
            sampling_rate=sampling_rate,
            shuffle=True,
       )
     
        print('===========')
        print('dummy_dataset, batch_size', batch_size)
        count = 0
        for inputs, targets in dummy_dataset:
            for i in range(inputs.shape[0]):
                count = count + 1
                print([int(x) for x in inputs[i]], int(targets[i]))
            print('count', count)
        print('===========')
        print('sampling_rate != 1 messes up the order')
        print('===========')
       
        for iterval in np.arange(1):
            print('===========')
            print('iterval', iterval)
            for samples, targets in dummy_dataset:
                print("samples shape:", samples.shape)
                print("targets shape:", targets.shape)
                for i in np.arange( samples.shape[0]):
                    print("samples ", [int(x) for x in samples[i]],
                          end =" ")
                    print("targets ", int(targets[i]))
                
            print('===========')
           

    example1()
    
    # raw_data, temperature = getdata()
    # print('temperature shape', temperature.shape)
    
    
if __name__   == '__main__':
     main()