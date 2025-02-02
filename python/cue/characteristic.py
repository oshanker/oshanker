#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 17:51:50 2022

@author: uorugant
"""
import math
import numpy as np
import pandas as pd
from numpy.linalg import qr
import time

import cmath
import matplotlib.pyplot as plt
import my_functions



def sampleeigvals(n):
    """
    eigenvalues of cue sample matrix.
    
    uses qr factorization     

    Parameters
    ----------
    n : size of cue matrix

    Returns
    -------
    angles : phases of eigenvalues

    """
    z = (np.random.randn(n,n) + 1j*np.random.randn(n,n))/ math.sqrt(2.0)
    q,r = qr(z)
    
    d = np.diagonal(r)
    ph = d/np.absolute(d)
    q = np.multiply(q, ph, q)
    w = np.linalg.eigvals(q)
    angles = []
    for i in range(0, n):
        angles.append(cmath.phase(w[i]))
        
    return angles



def generate_cue_distribution():
    """
    The real enchilada. RMT eigenfunction value distribution.
    The simulation output ('../out/fitcue.txt') is processed by smoothed.py
    to extract the distribution parameters p, lambda and sigma.

    Returns
    -------
    None.

    """
    print("-------")
    n = 80
    start = time.time_ns()
    np.random.seed(int(start/1.0E9)-1693197846)
    
    cos_incr = math.cos(-2*math.pi/n)
    sin_incr = math.sin(-2*math.pi/n)

    bins_in = []
    for i in np.arange(-16.125, 16.15, 0.05):
        bins_in.append(i)
    xaxis = []
    for i in range(0, len(bins_in)-1):
        xaxis.append((bins_in[i]+bins_in[i+1])/2)
    
    xaxis_zero = 322
    print('xaxis[', xaxis_zero , '] = ', xaxis[xaxis_zero])
    
    sample_size = 60000
    
    sums = [
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0
            ]
    vars = np.zeros(len(sums))
    print("-- cue --")
    abs_sum = np.zeros(len(sums))
    grams = []
    for j in range(0, len(sums)):
        grams.append(np.zeros((len(bins_in)-1,), dtype=int))
    phi_incr =(4/len(sums)) *math.pi/n
    
    offsets = np.arange(0, len(sums))*phi_incr
    #horz_rad = horz*2*math.pi/n - math.pi 
    
    for i in range(0, sample_size):
        angles = sampleeigvals(n)
        sum1n = np.sum(angles)/n
        #xbase = horz_rad + sum1n
        xbase0 = - math.pi + sum1n
        
        for j in range(0, len(sums)):
            #x = xbase + offsets[j]
            x0 = xbase0 + offsets[j]
            arg = (angles-x0)/2
            sin1 = 2*np.sin(arg)
            cos1 = 2*np.cos(arg)
            # if i==0 and j==0:
            #     print(len(x), int(n/2))
            #     print(np.diff(x), 2*2*math.pi/n)
            y1list = []
            z = np.prod(sin1)
           
            y1list.append(z)
             
            for k in range(1, int((n+1)/2)):
                s1_temp = sin1*cos_incr + cos1*sin_incr
                cos1 = cos1*cos_incr - sin1*sin_incr
                sin1 = s1_temp
                z = np.prod(sin1)
               
                y1list.append(z)
                
            values_array = np.array(y1list)
            
            hist_incr, bins_range = np.histogram(values_array, bins=bins_in, density=True)
            grams[j] = grams[j] + hist_incr
            sums[j] = sums[j] + np.mean(values_array)
            vars[j] = vars[j] + np.var(values_array, ddof=0)
            abs_sum[j] = abs_sum[j] + np.mean(values_array*values_array)
        
        
    for i in range(0, len(sums)):
        grams[i] = grams[i]/sample_size
        # plt.plot(xaxis, grams[i], '-x', color = 'black')
        # plt.grid()
        # plt.show()
        # input("waiting")
        # plt.clf()    
    
    means = np.array(sums)/sample_size
    vars = vars/sample_size
    abs_sum = abs_sum/sample_size
    print('vars', vars)
    hist = grams[0]
    histnorm = np.sum(hist)
    print('np.sum(hist)',  histnorm)
    xdata = np.array(xaxis)
    print('np.var(from hist)?? ',  np.sum(hist*xdata*xdata)/histnorm)
    print('abs_sum ',  abs_sum)
    
    index90 = int(len(sums)/4)
    print(grams[index90][xaxis_zero - 1], grams[index90][xaxis_zero], 
          grams[index90][xaxis_zero + 1])
    # plt.plot(xaxis, grams[0], '-x', color = 'black')
    # plt.grid()
    
    all_data = np.asarray(grams).T
    print('all_data.shape', all_data.shape)
    header = "n " + str(n) + " sample " + str(sample_size) + \
        " " + str(xaxis[0]) + " " + str(xaxis[-1])

    np.savetxt('../out/all_data' + str(n) + '.gz', 
               all_data,
               header=header)

    
    data_symm = []
    zdf = []
    for zindex in range(xaxis_zero - 160, xaxis_zero + 161):
        row = []
        for j in range(0, len(sums)):
            row.append(grams[j][zindex])
        data_symm.append(row)
        zdf.append(np.around(xaxis[zindex], decimals=2))

    data_symm.append(means)
    zdf.append(-100)
    
    phi_values = []
    for j in range(0, len(sums)):
        phi_values.append(int(j*360/len(sums)))
        
    x2array = np.zeros((3,len(phi_values)))
    x2array[0] = phi_values
    x2array[1] = abs_sum
    x2array[2] = means
    np.savetxt('../out/x2_' + str(n) + '.txt', 
               x2array, fmt='%3.15f')
        
    df = pd.DataFrame(data_symm, columns=phi_values, index=zdf)
    file_name =  '../out/cue.txt'
    df.to_csv(file_name, sep=' ')
    
    testfit = []
    for j in range(0, len(data_symm), 5):
        out = my_functions.fit(zdf[j], df.iloc[j].values, phi_values)
        testfit.append(out)
    np.savetxt('../out/fitcue.txt', 
               np.asarray(testfit), 
               fmt='%.2f %7.3f %7.3f %8.4f %7.3f', 
               delimiter=',')
    # changing data, put at very end
    data_symm[-1][0] = n
    data_symm[-1][1] = sample_size
    data_symm[-1][2] = zdf[0]
    header = "n " + str(n) + " sample " + str(sample_size)    

    np.savetxt('../out/rawcue.txt', 
               np.asarray(data_symm),
               header=header)
    
    elapsed = (time.time_ns()-start)/1.0E9
    print('elapsed', elapsed)

    
    
def plotValues(horz, y1, horz1, y2):
    
    plt.plot(horz, y1, '-x', color = 'red')
    plt.grid()
    ty1 = (np.max(y1)+np.min(y1))/2
    plt.text(0.25, ty1, "0", fontsize=12, color = 'red')
    plt.show()
    input("waiting")
    plt.clf()
    plt.plot(horz1, y2, '-.o', color = 'green')
    plt.grid()
    ty2 = (np.max(y2)+np.min(y2))/2
    plt.text(0.25, ty2, "1", fontsize=12, color = 'green')
    

   
generate_cue_distribution()

