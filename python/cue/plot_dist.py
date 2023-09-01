#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 10:09:04 2023

@author: uorugant
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math


gauss_norm = 1/(math.sqrt(2*math.pi))
def gauss(x, sigma ):
    A = gauss_norm/(sigma)
    return A*np.exp(-(x)**2/(2*sigma**2))

def exp_pdf(x, lam):
    xx = np.abs(np.copy(x))
    return  np.exp(-lam * xx) * lam/2
    
def exp_gauss(x, p, lam, sigma):
    ret = p * exp_pdf(x, lam) + (1-p) * gauss(x, sigma )
    return ret

def der_gauss(x, sigma ):
    A = gauss_norm/(sigma**3)
    return -A*x*np.exp(-(x)**2/(2*sigma**2))

def der_exp_pdf(x, lam):
    ret = []
    for arg in x:
        if arg == 0:
            ret.append(0)
        elif arg > 0:
            ret.append(-math.exp(-lam * arg) * lam*lam/2)
        else:
            ret.append(math.exp(lam * arg) * lam*lam/2)
    return  np.asarray(ret)
    
def der_exp_gauss(x, p, lam, sigma):
    ret = p * der_exp_pdf(x, lam) + (1-p) * der_gauss(x, sigma )
    return ret

def do_plot_with_data(func, popt, xdata, ydata, label):
    plt.plot(xdata, ydata, 'b-', label=label)
    plt.plot(xdata, func(xdata, *popt), 'g--',
         label='func' )
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
    plt.legend()
    plt.show()
    return popt

def do_plot_data(xdata, ydata, label):
    plt.plot(xdata, ydata, 'g--',
         label=label )
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()


def main():
    n = 95
    file_name = '../out/smoothedfitcue'+ str(n) + '.txt'
    data = np.loadtxt(file_name)
    popt = [0.29117866, 3.64662728, 1.50870858]
    output_array = np.asarray(data) 
    xdata = output_array[:,0] 
    ydata = output_array[:,1] 
    #do_plot_with_data(exp_gauss, popt, xdata, ydata, 'fit to A')
    #do_plot_with_data(der_exp_gauss, popt, xdata, ydata, 'der  A')
    b_data = output_array[:,2] 
    #do_plot_data(xdata, b_data, 'B')
    do_plot_with_data(der_exp_gauss, popt, xdata, -10*b_data, 'eval der  A')
    
main()