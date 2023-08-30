#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 10:24:21 2023

@author: uorugant
"""
import numpy as np


def flip_and_switch(arr):
    shape = arr.shape
    #col90_idx = int((arr.shape[1]-1)/4)
    #col270_idx = int(3*(arr.shape[1]-1)/4)
    col180_idx = int((arr.shape[1]-1)/2)
    col360_idx = int((arr.shape[1]-1))
    flipped = np.flip(arr, axis=0)
    ret = np.zeros(shape)
    for i in range(0, col180_idx+1):
       ret[:,i] = flipped[:,col180_idx-i ]
    for i in range(col180_idx+1, col360_idx+1):
       ret[:,i] = ret[:,col360_idx-i ]
    return (ret+arr)/2

def test_flip_and_switch():
    A = [
          [ 1, 1, 5, 1, 1],
          [ 2, 4, 4, 4, 2],
          [ 3, 5, 3, 5, 3],
          [ 4, 4, 2, 4, 4],
          [ 5, 1, 1, 1, 5]
        ]
    inputA = [
          [ 1, 1, 5, 1, 1],
          [ 3, 4, 4, 4, 3],
          [ 3, 5, 3, 5, 3],
          [ 4, 4, 1, 4, 4],
          [ 5, 1, 1, 1, 5]
        ]
    arr = np.array(A)    
    in_arr = np.array(inputA)  
    flipped = flip_and_switch (in_arr)
    i = 0
    for col in in_arr.T:
        print("Column ",i,":", col)
        i+=1
    print('flipped')
    i = 0
    for col in flipped.T:
        print("Column ",i,":", col)
        i+=1
    print((flipped-arr))

