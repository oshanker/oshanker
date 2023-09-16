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
import my_functions

    
def b_der_exp_gauss(x, b, p, lam, sigma):
    ret = b*my_functions.der_exp_gauss(x, p, lam, sigma)
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
    n = 80
    file_name = '../out/smoothedfitcue'+ str(n) + '.txt'
    data = np.loadtxt(file_name)
    output_array = np.asarray(data) 
    xdata = output_array[:,0] 
    ydata = output_array[:,1] 
    #do_plot_with_data(exp_gauss, popt, xdata, ydata, 'fit to A')
    #do_plot_with_data(der_exp_gauss, popt, xdata, ydata, 'der  A')
    b_data = output_array[:,2] 
    #do_plot_data(xdata, b_data, 'B')
    
    #popt = [0.83159422, 2.05637868, 0.56744309]
    #do_plot_with_data(my_functions.der_exp_gauss, popt, xdata, -13.0*b_data, 'ydata')
    
    # popt = [-0.076923076923077, 0.83159422, 2.05637868, 0.56744309]
    # do_plot_with_data(b_der_exp_gauss, popt, xdata, b_data, 'ydata')
    
    bounds=([-2.005, 0.2, 0.6, 7.5], [-1.999, 0.3, 0.8, 8.5])
    popt = my_functions.do_fit(b_der_exp_gauss, xdata, b_data, bounds=bounds) #0.9
    print('param', popt)

    #param [       0.43159422, 2.05637868, 3.56744309]
    #param [-1.997       0.22968019  0.69744049  7.73814023]
    
main()