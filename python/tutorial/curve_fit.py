#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 13:40:09 2023

@author: uorugant
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

def exp_pdf(x, lam):
    xx = np.copy(x)
    for i in range(0, xx.shape[0]):
        if(xx[i] < 0):
            xx[i] = -xx[i]
    return  np.exp(-lam * xx) * lam/2

def curve_fit_test():
    xdata = np.linspace(-2, 2, 50)
    y = func(xdata, 2.5, 1.3, 0.5)
    rng = np.random.default_rng()
    y_noise = 0.1 * rng.normal(size=xdata.size)
    ydata = y + y_noise
    plt.plot(xdata, ydata, 'b-', label='data')
    popt, pcov = curve_fit(func, xdata, ydata)
    print(popt)
    plt.plot(xdata, func(xdata, *popt), 'r-',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    popt, pcov = curve_fit(func, xdata, ydata, bounds=([2.2, 1., 0.4], [2.6, 1.4, 0.6]))
    print(popt)
    plt.plot(xdata, func(xdata, *popt), 'g--',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()
    
    cond = np.linalg.cond(pcov)
    print(cond)
    print(np.diag(pcov))
    
def exp_fit_test():
    xdata = np.linspace(-1, 1, 41)
    ydata = exp_pdf(xdata, 1.45)
    do_exp_fit(xdata, ydata)
    
def do_exp_fit(xdata, ydata):
    print(np.sum(ydata))
    popt, pcov = curve_fit(exp_pdf, xdata, ydata, bounds=([1], [3]))
    print(popt)
    cond = np.linalg.cond(pcov)
    print(cond)
    print(np.diag(pcov))
    plt.plot(xdata, ydata, 'b-', label='data')
    plt.plot(xdata, exp_pdf(xdata, *popt), 'g--',
         label='fit: lam=%5.3f' % tuple(popt))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()


exp_fit_test()