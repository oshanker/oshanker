#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 20:31:03 2022

@author: uorugant
"""
import sys
import numpy as np 
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from test001 import getZetadata
import os.path as path

import pandas
import matplotlib.pyplot as plt

def myplot(history, prefix, minidx = 2):
    
    loss = history.history["mae"]
    val_loss = history.history["val_mae"]
    epochs = range(int(minidx), len(loss) + 1)
    plt.figure()
    plt.plot(epochs, loss[1:], "bo", label=prefix+" Training MAE")
    plt.plot(epochs, val_loss[1:], "b", label=prefix+" Validation MAE")
    plt.title(prefix+" Training and validation MAE")
    plt.legend()
    plt.show()
    

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



# https://www.manning.com/books/deep-learning-with-python-second-edition
# https://github.com/fchollet/deep-learning-with-python-notebooks
# https://s3.amazonaws.com/keras-datasets/jena_climate_2009_2016.csv.zip

def main():
    print("-------")
    print(sys.argv)
    
    reload = False
    
    
    #sequence_length = 120 
    sequence_length = 25 
    
    #raw_data, temperature = getdata()
    raw_data, temperature = getZetadata(500001, sequence_length)
    print('temperature shape', temperature.shape)
    
    num_train_samples = int(0.5 * len(raw_data))
    num_val_samples = int(0.25 * len(raw_data))
    num_test_samples = len(raw_data) - num_train_samples - num_val_samples
    print("num_train_samples:", num_train_samples)
    print("num_val_samples:", num_val_samples)
    print("num_test_samples:", num_test_samples)

    mean = raw_data[:num_train_samples].mean(axis=0)
    #raw_data -= mean
    print('raw_data mean', mean)
    std = raw_data[:num_train_samples].std(axis=0)
    #raw_data /= std
    print('raw_data std', std)
    
    print('temperature mean', temperature[:num_train_samples].mean(axis=0))
    print('temperature std', temperature[:num_train_samples].std(axis=0))
    
    sampling_rate = 1 
    
    
    #delay = sampling_rate * (sequence_length + 24 - 1)
    #delay = (sequence_length + 6*24 - 1)
    
    delay = sampling_rate * (sequence_length + 1 - 1)
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
      
    #print(f"evaluate_naive_method Validation MAE: {evaluate_naive_method(val_dataset):.2f}") 
    #print(f"evaluate_naive_method Test MAE: {evaluate_naive_method(test_dataset):.2f}")

    inputs = keras.Input(shape=(sequence_length, raw_data.shape[-1]))
    optimizer=['adam',"rmsprop"]
    
    #######################################
    id = "LSTM"
    file_name = "../out/jena_xxx.keras"
    def train_keras():
        if reload:
            model = keras.models.load_model(file_name) 
        else:
            x = layers.LSTM(16)(inputs)
            outputs = layers.Dense(1, name="output_layer")(x)
            model = keras.Model(inputs, outputs, name=id)
            model.compile(optimizer=optimizer[1], loss="mse", metrics=["mae"])
        
          
        callbacks = [
            keras.callbacks.ModelCheckpoint(file_name,
                                            save_best_only=True)
        ]
        
        epochs=3
        #epochs=3
        
        history = model.fit(train_dataset,
                            epochs=epochs,
                            validation_data=val_dataset,
                            callbacks=callbacks)
        myplot(history, id)
        
    def train_dropout():
        x = layers.LSTM(32, recurrent_dropout=0.25)(inputs)
        x = layers.Dropout(0.5)(x)                             

        outputs = layers.Dense(1, name="output_layer")(x)
        
        model = keras.Model(inputs, outputs, name=id)
          
        callbacks = [
            keras.callbacks.ModelCheckpoint(file_name,
                                            save_best_only=True)
        ]
        
        epochs=8
        #epochs=3
        
        model.compile(optimizer=optimizer[1], loss="mse", metrics=["mae"])
        
        history = model.fit(train_dataset,
                            epochs=epochs,
                            validation_data=val_dataset,
                            callbacks=callbacks)
        myplot(history, id)
    
    def train_bidirectional():
        if reload and path.exists(file_name):
            model = keras.models.load_model(file_name) 
        else:
            # https://towardsdatascience.com/a-look-at-gradient-descent-and-rmsprop-optimizers-f77d483ef08b#:~:text=The%20difference%20between%20RMSprop%20and,is%20usually%20set%20to%200.9.
            x = layers.Bidirectional(layers.LSTM(16))(inputs)
            outputs = layers.Dense(1, name="output_layer")(x)
            model = keras.Model(inputs, outputs, name=id)
            model.compile(
                optimizer=tf.keras.optimizers.RMSprop(
                    learning_rate=0.0015, momentum=0.85), 
                loss="mse", metrics=["mae"])
          
        callbacks = [
            keras.callbacks.ModelCheckpoint(file_name,
                                            save_best_only=True)
        ]
        
        epochs=20
        #epochs=3
        
        history = model.fit(train_dataset,
                            epochs=epochs,
                            validation_data=val_dataset,
                            callbacks=callbacks)
        myplot(history, id)
    
    
    def train_expt():
        x = layers.Bidirectional(layers.LSTM(16))(inputs)

        outputs = layers.Dense(1, name="output_layer")(x)
        
        model = keras.Model(inputs, outputs, name=id)
          
        callbacks = [
            keras.callbacks.ModelCheckpoint(file_name,
                                            save_best_only=True)
        ]
        
        epochs=8
        #epochs=3
        
        model.compile(optimizer=optimizer[1], loss="mse", metrics=["mae"])
        
        history = model.fit(train_dataset,
                            epochs=epochs,
                            validation_data=val_dataset,
                            callbacks=callbacks)
        myplot(history, id)
    
   
    #train_keras()
    #train_dropout()
    train_bidirectional()
    
    #######################################
    
    model = keras.models.load_model(file_name) 
    print(f"Test MAE: {model.evaluate(test_dataset)[1]:.2f}")
    print(model.summary())
   
    
if __name__   == '__main__':
     main()