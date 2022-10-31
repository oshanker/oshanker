#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 20:31:03 2022

@author: uorugant
"""
import sys
import numpy as np 
from tensorflow import keras
import pandas

# https://www.manning.com/books/deep-learning-with-python-second-edition
# https://github.com/fchollet/deep-learning-with-python-notebooks
# https://s3.amazonaws.com/keras-datasets/jena_climate_2009_2016.csv.zip

def main():
    print("-------")
    print(sys.argv)
    time_sequence = np.arange(10)                                
    y_sequence = np.arange(3,18,3)                                
    dummy_dataset = keras.utils.timeseries_dataset_from_array(
        data=time_sequence[:-3],                                 
        targets=y_sequence,                               
        sequence_length=3,                                      
        batch_size=3,                                           
    )
 
    for inputs, targets in dummy_dataset:
        for i in range(inputs.shape[0]):
            print([int(x) for x in inputs[i]], int(targets[i]))
    
    dataset = pandas.read_csv('../../../jena_climate_2009_2016.csv', header=0)
    dataset.drop(dataset.columns[0], axis=1, inplace=True)
    
    print(dataset.head)
    values = dataset.values



if __name__   == '__main__':
     main()