#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 20:31:03 2022

@author: uorugant
"""
import sys
import numpy as np 
from tensorflow import keras
from tensorflow.keras import layers

import pandas
import matplotlib.pyplot as plt

def plot(history, prefix):
    
    loss = history.history["mae"]
    val_loss = history.history["val_mae"]
    epochs = range(2, len(loss) + 1)
    plt.figure()
    plt.plot(epochs, loss[1:], "bo", label=prefix+" Training MAE")
    plt.plot(epochs, val_loss[1:], "b", label=prefix+" Validation MAE")
    plt.title("Training and validation MAE")
    plt.legend()
    plt.show()
    
# https://www.manning.com/books/deep-learning-with-python-second-edition
# https://github.com/fchollet/deep-learning-with-python-notebooks
# https://s3.amazonaws.com/keras-datasets/jena_climate_2009_2016.csv.zip

def main():
    print("-------")
    print(sys.argv)
    
    def example1():
        time_sequence = np.arange(10)  
        
        raw_data = np.zeros((7, 2 ))  
        for i in time_sequence[:-3]:
            raw_data[i,:] = [i, i*i]
        print(raw_data)
        y_sequence = np.arange(3,21,3)                                
        dummy_dataset = keras.utils.timeseries_dataset_from_array(
            data=raw_data,                                 
            targets=y_sequence,                               
            sequence_length=3,                                      
            batch_size=2,                                           
        )
     
        for inputs, targets in dummy_dataset:
            for i in range(inputs.shape[0]):
                print(inputs[i], int(targets[i]))
    
    def getdata():
        dataset = pandas.read_csv('../../../jena_climate_2009_2016.csv', header=0)
        dataset.drop(dataset.columns[0], axis=1, inplace=True)
        
        print(dataset.head())
        raw_data = dataset.values
        print('raw_data.shape', raw_data.shape)
        print('raw_data[0]', raw_data[0])
    
        temperature = np.array(raw_data[:,1], copy=True)
        temperature[0] = float("nan")
        print('raw_data[0]', raw_data[0])
        print('temperature[0]', temperature[0])
        return raw_data, temperature
    
    raw_data, temperature = getdata()
    print('temperature shape', temperature.shape)
    
    num_train_samples = int(0.5 * len(raw_data))
    num_val_samples = int(0.25 * len(raw_data))
    num_test_samples = len(raw_data) - num_train_samples - num_val_samples
    print("num_train_samples:", num_train_samples)
    print("num_val_samples:", num_val_samples)
    print("num_test_samples:", num_test_samples)

    mean = raw_data[:num_train_samples].mean(axis=0)
    raw_data -= mean
    std = raw_data[:num_train_samples].std(axis=0)
    raw_data /= std
    
    sampling_rate = 1 
    sequence_length = 120 
    #delay = sampling_rate * (sequence_length + 24 - 1)
    delay = (sequence_length + 6*24 - 1)
    batch_size = 256 
      
    train_dataset = keras.utils.timeseries_dataset_from_array(
        raw_data[:-delay],
        targets=temperature[delay:],
        sampling_rate=sampling_rate,
        sequence_length=sequence_length,
        shuffle=True,
        batch_size=batch_size,
        start_index=0,
        end_index=num_train_samples)
      
    print(
        """
        To preserve order,  sampling_rate has to be 1! 
        model is quickly overfitting, despite only having very 
          few units: the training and validation losses start 
          to diverge considerably after a few epochs. 
          You’re already familiar with a classic technique 
          for fighting this phenomenon: dropout, which randomly zeros out 
          input units of a layer in order to break happenstance correlations 
          in the training data 
        """
          )
    print('type(train_dataset)', type(train_dataset))
    
    val_dataset = keras.utils.timeseries_dataset_from_array(
        raw_data[:-delay],
        targets=temperature[delay:],
        sampling_rate=sampling_rate,
        sequence_length=sequence_length,
        shuffle=True,
        batch_size=batch_size,
        start_index=num_train_samples,
        end_index=num_train_samples + num_val_samples)
      
    test_dataset = keras.utils.timeseries_dataset_from_array(
        raw_data[:-delay],
        targets=temperature[delay:],
        sampling_rate=sampling_rate,
        sequence_length=sequence_length,
        shuffle=True,
        batch_size=batch_size,
        start_index=num_train_samples + num_val_samples)
    
    for samples, targets in train_dataset:
        print("samples shape:", samples.shape)
        print("targets shape:", targets.shape)
        break

    def evaluate_naive_method(dataset):
        total_abs_err = 0. 
        samples_seen = 0 
        for samples, targets in dataset:
            preds = samples[:, -1, 1] * std[1] + mean[1]         
            total_abs_err += np.sum(np.abs(preds - targets))
            samples_seen += samples.shape[0]
        return total_abs_err / samples_seen
      
    print(f"evaluate_naive_method Validation MAE: {evaluate_naive_method(val_dataset):.2f}") 
    print(f"evaluate_naive_method Test MAE: {evaluate_naive_method(test_dataset):.2f}")

    inputs = keras.Input(shape=(sequence_length, raw_data.shape[-1]))
    optimizer=['adam',"rmsprop"]
    
    #######################################
    
    file_name = "../out/jena_lstm.keras"
    def train_keras():
        x = layers.LSTM(16)(inputs)
        outputs = layers.Dense(1)(x)
        model = keras.Model(inputs, outputs, name="LSTM")
          
        callbacks = [
            keras.callbacks.ModelCheckpoint(file_name,
                                            save_best_only=True)
        ]
        
        epochs=10
        
        model.compile(optimizer=optimizer[1], loss="mse", metrics=["mae"])
        
        history = model.fit(train_dataset,
                            epochs=epochs,
                            validation_data=val_dataset,
                            callbacks=callbacks)
        plot(history, "LSTM")
    
    train_keras()
    model = keras.models.load_model(file_name) 
    print(f"LSTM Test MAE: {model.evaluate(test_dataset)[1]:.2f}")
    print(model.summary())
    
    #######################################
    
    
    
if __name__   == '__main__':
     main()