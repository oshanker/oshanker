#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 21:24:57 2022

@author: uorugant
"""

import matplotlib.pyplot as plt
import numpy as np
import my_functions
import time

def updateDiffs(n, diff):
    angles = np.degrees(my_functions.sample(n))
    for i in range(1, n):
        diff.append(angles[i] - angles[i-1])
    diff.append( angles[0] + 360 - angles[n-1] )

def main():
    print("-------")
    n = 76
    np.random.seed(2021)
    #print(angles)
    diff = []
    start = time.time()
    for i in range(0, 10000):
        updateDiffs(n, diff)
    
    elapsed = time.time()-start
    print('elapsed', elapsed)
    bins_in = []
    for i in np.arange(0,27,0.25):
        bins_in.append(i)
    
    #print (bins_in)
    
    print(np.median(diff), np.mean(diff), np.std(diff))
    #hist, bins_range = np.histogram(diff, bins=bins_in)
    #print(hist)
    #print(hist.sum())
    plt.hist(diff, bins=bins_in)
    plt.title("Histogram Diffs")
    plt.grid()
    plt.show()
    # printing eigen vectors
    # print("Printing Right eigenvectors of the given square array:\n",
    #       v)
    # print(q.dot(v))
    
main()