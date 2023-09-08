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

def exp_pdf(x, lam, corr):
    """
    

    Parameters
    ----------
    x : arg
    
    lam : constant for exponential distribution

    Returns
    -------
    exponential distribution

    """
    xx = np.abs(np.copy(x))
    return  corr*np.exp(-lam * xx) * lam/2

def do_fit(func, xdata, ydata):
    print(np.sum(ydata))
    popt, pcov = curve_fit(func, xdata, ydata, bounds=([1.7, 0.8], [3.3, 2.0]))
    print('param', popt)
    cond = np.linalg.cond(pcov)
    print("cond", cond)
    print("diag cov", np.diag(pcov))
    plt.plot(xdata, ydata, 'b-', label='data')
    plt.plot(xdata, func(xdata, *popt), 'g--',
         label='fit' )
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid()
    plt.show()
    return popt

def do_plot_func(func, popt, xdata, label):
    plt.plot(xdata, func(xdata, *popt), 'g--',
         label=label )
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid()
    plt.legend()
    plt.show()
    return popt

def main():
    bins_in = []
    for i in np.arange(-3.125, 3.5, 0.05):
        bins_in.append(i)
    xaxis = []
    for i in range(0, len(bins_in)-1):
        xaxis.append((bins_in[i]+bins_in[i+1])/2)
    size = 300000
    lam = 3.0
    data1 = np.random.exponential(scale=1.0/lam, size=size)
    data2 = np.random.exponential(scale=1.0/lam, size=size)
    data = np.concatenate([data1, -data2])
    hist, bins_range = np.histogram(data, bins=bins_in, density=True)
    histnorm = np.sum(hist)
    print('np.sum(hist)',  histnorm)
    print('np.var(data)',  np.var(data))
    xdata = np.array(xaxis)
    print('np.var(from hist)?? 2/lam**2',  np.sum(hist*xdata*xdata)/histnorm)
    popt = do_fit(exp_pdf, xdata, hist)
    #do_plot_func(der_gauss, popt, x, 'der')
    
    # simul = gauss(x, 2.0)
    # plt.plot(xaxis, simul, 'b-', label='data')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.legend()
    # plt.show()

main()