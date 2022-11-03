#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 17:37:27 2022

@author: uorugant
"""
import sys
import numpy as np 
from tensorflow import keras

import pandas
import matplotlib.pyplot as plt

def plotzeta(raw_data):
    plt.plot(range(1,15), raw_data[:14])
    plt.grid(True)
    plt.xlabel('offset gram index')
    plt.ylabel('zeta')

    plt.title('zeta at gram points')

def timeseries_dataset(raw_data, temperature, delay, sequence_length, 
                       batch_size, start_index, end_index):
    train_dataset = keras.utils.timeseries_dataset_from_array(
        raw_data[:-delay],
        targets=temperature[delay-1:],
        sampling_rate=1,
        sequence_length=sequence_length,
        shuffle=True,
        batch_size=batch_size,
        start_index=start_index,
        end_index=end_index)
    return train_dataset 
    
def getZetadata(upper, sequence_length ):
    dataset = pandas.read_csv('../../oldriemann/data/zetaE12.csv', header=0)
    print(dataset.head())
    
    #drop extra column
    dataset.drop(dataset.columns[0], axis=1, inplace=True)
    
    raw_data = dataset.values[1:upper]
    #raw_data = dataset.values
    print('raw_data.shape', raw_data.shape)
    print('raw_data[0]', raw_data[0])
    print(raw_data[:5])
    
    # dataset1 = pandas.read_csv('../../oldriemann/data/zerosE12.csv.max', header=None)

    dataset1 = pandas.read_csv('../../oldriemann/data/zetaE12.csv', header=0)
    #drop extra column
    dataset1.drop(dataset1.columns[1], axis=1, inplace=True)
    
    #temperature = dataset1.values[1:upper]
    #temperature = dataset1.values
    
    temperature = np.zeros((raw_data.shape[0] ))  
    for i in np.arange(raw_data.shape[0]) :
        if i >= sequence_length - 1:
            temperature[i] = np.max(
                np.abs(raw_data[i+1-sequence_length:i+1]) )
    
    print('temperature[sequence_length - 1]', temperature[sequence_length - 1])
    print(temperature[sequence_length - 1:sequence_length + 5])
    
    return raw_data, temperature

def evaluate_naive_method(dataset, do_print = False):
    total_abs_err = 0. 
    samples_seen = 0 
    for samples, targets in dataset:
        total = 0
        for i in np.arange( samples.shape[0]):
            pred  = np.max( np.abs(samples[i]) )
            if do_print:
                print("pred ", pred,  end =" ")
                print("targets ", float(targets[i]))
            total  += np.abs( pred -  float(targets[i]))

        total_abs_err += total
        samples_seen += samples.shape[0]
    return total_abs_err / samples_seen


def main():
    print("-------")
    print(sys.argv)

    def example2():
        limit = 15
        sequence_length = 3
        time_sequence = np.arange(limit)  
        raw_data = np.zeros((time_sequence.shape[0] ), dtype=int)  
        y_sequence = np.zeros((raw_data.shape[0] ), dtype=int)  
        sign = 1
        for i in np.arange(raw_data.shape[0]) :
            raw_data[i] =  i * sign
            sign = -sign
            if i >= sequence_length - 1:
                y_sequence[i] = np.max(
                    np.abs(raw_data[i+1-sequence_length:i+1]) )
        print(raw_data)
        print(y_sequence)

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
           

    #example1()
    #example2()
    sequence_length = 3
    raw_data, temperature = getZetadata(500001, sequence_length)
    
    delay = sequence_length 
       
    batch_size = 3 

    train_dataset = timeseries_dataset(raw_data, temperature, 
                               delay, sequence_length, 
                               batch_size, 0, 15)

    print('===========')
    mae = evaluate_naive_method(train_dataset)
    print('mae', mae)
    
    print('===========')
    print('train_dataset, batch_size', batch_size)
    count = 0
    for inputs, targets in train_dataset:
        for i in range(inputs.shape[0]):
            count = count + 1
            print([float(x) for x in inputs[i]], float(targets[i]))
        print('count', count)
    print('===========')
    
    
if __name__   == '__main__':
     main()
     
     """
                    learning_rate=0.001, momentum=0.9), 

Epoch 20/20
977/977 [==============================] - 28s 29ms/step - loss: 0.9157 - mae: 0.5804 - val_loss: 0.6761 - val_mae: 0.4750
489/489 [==============================] - 8s 15ms/step - loss: 0.6510 - mae: 0.4762
Test MAE: 0.48

                    learning_rate=0.001, momentum=0.905), 

Epoch 20/20
977/977 [==============================] - 21s 21ms/step - loss: 1.0820 - mae: 0.6391 - val_loss: 0.8594 - val_mae: 0.6069
489/489 [==============================] - 6s 11ms/step - loss: 0.8213 - mae: 0.6074
Test MAE: 0.61

smooth
                    learning_rate=0.001, momentum=0.895), 
Epoch 20/20
977/977 [==============================] - 21s 22ms/step - loss: 0.9322 - mae: 0.5925 - val_loss: 0.5451 - val_mae: 0.4496
489/489 [==============================] - 7s 13ms/step - loss: 0.5637 - mae: 0.4519
Test MAE: 0.45

                    learning_rate=0.001, momentum=0.8975), 
Epoch 20/20
977/977 [==============================] - 21s 22ms/step - loss: 0.8620 - mae: 0.5647 - val_loss: 0.8698 - val_mae: 0.4838
489/489 [==============================] - 6s 11ms/step - loss: 0.8443 - mae: 0.5605
Test MAE: 0.56

    
     """     
     
     
     
     
     
     