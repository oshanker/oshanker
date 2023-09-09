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
def gauss(x, sigma, corr ):
    A = gauss_norm/(sigma)
    return corr*A*np.exp(-(x*x)/(2*sigma*sigma))

def der_gauss(x, sigma ):
    A = gauss_norm/(sigma**3)
    return -A*x*np.exp(-(x)**2/(2*sigma**2))
    
def do_fit(func, xdata, ydata):
    print(np.sum(ydata))
    popt, pcov = curve_fit(func, xdata, ydata, bounds=([0.5, 0.8], [2.5, 2.0]))
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
    for i in np.arange(-16.125, 16.5, 0.05):
        bins_in.append(i)
    xaxis = []
    for i in range(0, len(bins_in)-1):
        xaxis.append((bins_in[i]+bins_in[i+1])/2)
    
    print(len(xaxis))
    xaxis_zero = 322
    print('xaxis[', xaxis_zero , '] = ', xaxis[xaxis_zero])
    size = 300000
    sigma = 2.0
    data = np.random.normal(loc=0.0, scale=sigma, size=size)
    hist, bins_range = np.histogram(data, bins=bins_in, density=True)
    histnorm = np.sum(hist)
    print('sigma*hist[', xaxis_zero , '] = ', sigma*hist[xaxis_zero])
    halfmax = int(1.1775*sigma/0.05)
    print('sigma*hist[', xaxis_zero+halfmax , '] = ', sigma*hist[xaxis_zero+halfmax ])
    print('sigma*hist[', xaxis_zero-halfmax , '] = ', sigma*hist[xaxis_zero-halfmax ])
    print('np.sum(hist)',  histnorm)
    print('np.var(data)',  np.var(data))
    print('abs(sigma - np.std(s, ddof=1))', abs(sigma - np.std(data, ddof=1)))
    xdata = np.array(xaxis)
    x2 = xdata*xdata
    print('truncated np.var(from hist)??',  np.sum(hist*x2)/histnorm)
    popt = do_fit(gauss, xdata, hist)
    #do_plot_func(der_gauss, popt, x, 'der')
    
    # simul = gauss(x, 2.0)
    # plt.plot(xaxis, simul, 'b-', label='data')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.legend()
    # plt.show()

main()