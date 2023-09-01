#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 22:37:20 2023

@author: uorugant
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

gauss_norm = 1/(math.sqrt(2*math.pi))
def gauss(x, sigma ):
    A = gauss_norm/(sigma)
    return A*np.exp(-(x)**2/(2*sigma**2))
    
def do_fit(func, xdata, ydata):
    print(np.sum(ydata))
    popt, pcov = curve_fit(func, xdata, ydata, bounds=([1.5], [2.5]))
    print('param', popt)
    cond = np.linalg.cond(pcov)
    print("cond", cond)
    print("diag cov", np.diag(pcov))
    plt.plot(xdata, ydata, 'b-', label='data')
    plt.plot(xdata, func(xdata, *popt), 'g--',
         label='fit: lam=%5.3f' % tuple(popt))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()

def main():
    bins_in = []
    for i in np.arange(-3.125, 3.5, 0.05):
        bins_in.append(i)
    xaxis = []
    for i in range(0, len(bins_in)-1):
        xaxis.append((bins_in[i]+bins_in[i+1])/2)
    size = 50000
    data = np.random.normal(loc=0.0, scale=2.0, size=size)
    hist, bins_range = np.histogram(data, bins=bins_in, density=True)
    x = np.array(xaxis)
    do_fit(gauss, x, hist)
    
    # simul = gauss(x, 2.0)
    # plt.plot(xaxis, simul, 'b-', label='data')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.legend()
    # plt.show()

main()