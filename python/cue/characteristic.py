#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 17:51:50 2022

@author: uorugant
"""
import math
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from numpy.linalg import qr
import time

import cmath
import matplotlib.pyplot as plt
import my_functions

q = my_functions.sampleQ()
w = np.linalg.eigvals(q)
eigenlist = []
for e in w:
    eigenlist.append( cmath.phase(e))
eigen = np.array(eigenlist)

epsilon = 0.000001
conv=180/math.pi


def sampleeigvals(n):
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



def test4():
    print("-------")
    n = 63
    np.random.seed(2020)
    
    cos_incr = math.cos(-2*math.pi/n)
    sin_incr = math.sin(-2*math.pi/n)
    

    bins_in = []
    for i in np.arange(-3.15, 3.5, 0.1):
        bins_in.append(i)
    xaxis = []
    for i in range(0, len(bins_in)-1):
        xaxis.append((bins_in[i]+bins_in[i+1])/2)
    
    xaxis_zero = 31
    print('index', xaxis[xaxis_zero])
    sample_size = 20000
    
    sums = [
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0
            ]
    print("-- cue --")
    grams = []
    for j in range(0, len(sums)):
        grams.append(np.zeros((len(bins_in)-1,), dtype=int))
    phi_incr =(4/len(sums)) *math.pi/n
    
    offsets = np.arange(0, len(sums))*phi_incr
    #horz_rad = horz*2*math.pi/n - math.pi 
    
    start = time.time()
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
                
            y_1 = np.array(y1list)
            
            hist, bins_range = np.histogram(y_1, bins=bins_in)
            grams[j] = grams[j] + hist
            sums[j] = sums[j] + np.mean(y_1)
        
        
    for i in range(0, len(sums)):
        converted = grams[i].sum()*(xaxis[1]-xaxis[0])
        grams[i] = grams[i]/converted
        # plt.plot(xaxis, grams[i], '-x', color = 'black')
        # plt.grid()
        # plt.show()
        # input("waiting")
        # plt.clf()    
    
    means = np.array(sums)/sample_size
    index90 = int(len(sums)/4)
    print(grams[index90][xaxis_zero - 1], grams[index90][xaxis_zero], 
          grams[index90][xaxis_zero + 1])
    # plt.plot(xaxis, grams[0], '-x', color = 'black')
    # plt.grid()
    
    data = []
    zdf = []
    for zindex in range(xaxis_zero - 20, xaxis_zero + 22, 5):
        row = []
        for j in range(0, len(sums)):
            row.append(grams[j][zindex])
        data.append(row)
        zdf.append(np.around(xaxis[zindex], decimals=1))

    data.append(means)
    zdf.append(-100)
    
    phi_values = []
    for j in range(0, len(sums)):
        phi_values.append(int(j*360/len(sums)))
        
    df = pd.DataFrame(data,columns=phi_values,index=zdf)
    file_name =  '../out/cue.txt'
    df.to_csv(file_name, sep=' ')

    testfit = []
    for j in range(0, len(data)):
        out = fit(df.iloc[j].values, phi_values)
        testfit.append(out)
    np.savetxt('../out/fitcue.txt', 
               np.asarray(testfit), 
               fmt='%7.3f&%7.3f&%8.4f&%7.3f\\\\', 
               delimiter=',')
    elapsed = time.time()-start
    print('elapsed', elapsed)

def fit(values, phi_values):
    x = []
    for phi in phi_values:
        x.append([math.cos(math.pi*phi/180), math.cos(math.pi*phi/90)])
        #x.append([math.cos(math.pi*phi/180)])
    x = np.array(x)
    reg = LinearRegression().fit(x, values)
    r2 = reg.score(x, values)
    #
    coeff = reg.coef_
    intercept = reg.intercept_
    return [intercept, coeff[0], coeff[1], r2]
    
    
    # sum_hist = grams[0]+grams[2]
    # print(sum_hist)
    # print(sum_hist[15], sum_hist[16])
    
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
    

def test3():
    sum1 = np.sum(eigen)
    x = np.array([0, 1, 2])*2*math.pi/3 + sum1/3-math.pi
    print(x*conv)
    for theta in x:
       f = my_functions.eval_factor(theta, eigen)
       print('f', f)
       z = my_functions.Z(theta, eigen)
       print('z', z)
    
def test2():
    theta =  1.5* 2*math.pi/6
    f = my_functions.eval_direct(theta, eigen)
    print(f)
    f = my_functions.eval_factor(theta, eigen)
    print(f)
    
def test1():
    base = 2*math.pi/30
    x = np.array([-1,  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
                  10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                  20, 21, 22])*base
    y1list = []
    y2list = []
    y3list = []
    for theta in x:
        # print("eval at", math.degrees(theta))
        # print(eval(theta, eigen))
        # print(complex(1-math.cos(3*theta), math.sin(3*theta)))
        f = my_functions.eval_direct(theta, eigen)
        y1list.append(f.real)
        y2list.append(f.imag)
        if abs(f)<epsilon:
            y3list.append(math.pi/2)
        else:
            y3list.append(cmath.phase(f))
        
    y1 = np.array(y1list)
    plt.plot(x, y1, '-x', color = 'red')
    
    y2 = np.array(y2list)
    plt.plot(x, y2, '-.o', color = 'green')
    
    y3 = np.array(y3list)
    plt.plot(x, y3, '--+', color = 'black')
    
    plt.text(0.25, 1.3, "phase", fontsize=12)
    plt.text(0.25, 0, "real", fontsize=12, color = 'red')
    plt.text(0.92, 0.4, "imag", fontsize=12, color = 'green')
    plt.grid()
    
    plt.show()
    print(x*conv)
    print(y3*conv)
    
test4()

