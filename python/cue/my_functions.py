#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 10:21:08 2022

@author: uorugant
"""
import math
import numpy as np
from numpy.linalg import qr
import cmath
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

gauss_norm = 1/(math.sqrt(2*math.pi))
def gauss(x, sigma ):
    """

    Parameters
    ----------
    x : arg
    
    sigma : std dev

    Returns
    -------
    normal distribution

    """
    A = gauss_norm/(sigma)
    return A*np.exp(-(x)**2/(2*sigma**2))

def exp_pdf(x, lam):
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
    return  np.exp(-lam * xx) * lam/2
    
def exp_gauss(x, p, lam, sigma):
    """    

    Parameters
    ----------
    x : arg
    
    p : coefficient for exponential distribution
    
    lam : constant for exponential distribution
    
    sigma :  std dev

    Returns
    -------
    ret : linear combination of exponential distribution and normal distribution.

    """
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

def do_fit(func, xdata, ydata, bounds=([0.28, 3.6, 1.5], [0.31, 3.7, 1.52])):
    """
    fit function to data

    Parameters
    ----------
    func : function to fit
    
    xdata : independent variables
    
    ydata : dependent variable (distribution to be fitted)

    Returns
    -------
    popt : fitted parameters 

    """
    print(np.sum(ydata))
    #popt, pcov = curve_fit(func, xdata, ydata, bounds=([0.01], [2.5]))
    
    #p, lam, sigma
    popt, pcov = curve_fit(func, xdata, ydata, bounds=bounds)
    # [0.3        3.55578876 1.52809399]
    # [0.29117862 3.64662775 1.50870849]
    # [0.29117866 3.64662728 1.50870858]
    
    print('param', popt)
    cond = np.linalg.cond(pcov)
    print("cond", cond)
    print("diag cov", np.diag(pcov))
    plt.plot(xdata, ydata, 'b-', label='data')
    plt.plot(xdata, func(xdata, *popt), 'g--',
         label='fit: ' )
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()
    return popt


def sample(n):
    """
    sorted eigenvalues of cue sample matrix.
    
    uses qr factorization 
    """
    z = (np.random.randn(n,n) + 1j*np.random.randn(n,n))/ np.lib.scimath.sqrt(2.0)
    q,r = qr(z)
    #print(z)
    #print(q.dot(r))
    #print("-------")
    #print(q)
    #print("unitarity")
    #print(q.dot(q.T.conjugate()))
    
    d = np.diagonal(r)
    ph = d/np.absolute(d)
    q = np.multiply(q,ph,q)
    #print(q)
    w, eigenvect = np.linalg.eig(q)
    # printing eigen values
    angles = []
    for i in range(0, n):
        angles.append(cmath.phase(w[i]))
        
    angles = np.sort( angles)
    return angles

def fit(z, values, phi_values):
    """
    fit A B C

    Parameters
    ----------
    values : array
        DESCRIPTION.
    phi_values : array
        angles phi.

    Returns
    -------
    list
        DESCRIPTION.

    """
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
    return [z, intercept, coeff[0], coeff[1], r2]


def eval_direct(theta, root_angles):
    #n = eigen.size
    arg = root_angles-theta
    ret = 1+0j
    for x in arg:
        ret = ret * complex(1-math.cos(x), -math.sin(x))
    return ret

def eval_factor(theta, root_angles):
    #n = eigen.size
    arg = (root_angles-theta)/2
    sum1 = np.sum(arg + math.pi/2)
    sin1 = -2*np.sin(arg)
    prod = np.prod(sin1)
    ret = prod*cmath.exp((sum1)*1j)    
    return ret

def Z(theta, root_angles):
    arg = (root_angles-theta)/2
    sin1 = 2*np.sin(arg)
    prod = np.prod(sin1)
    return prod

def sampleQ():
    omega =  cmath.exp((2*cmath.pi/3)*1j)
    root2 = (1/math.sqrt(2))
    root3 = (1/math.sqrt(3))
    root6 = (1/math.sqrt(6))
    omegasq = omega*omega
    a = np.array([[root3, root3,root3], [root6, root6, -2*root6], 
                  [root2, -root2, 0]])
    diag = np.diagflat([1, omega, omegasq])
    
    q = a.dot(diag).dot(a.T)
    # print("q")
    # print(q)
    # print("unitarity")
    # print(q.dot(q.T.conjugate()))
    return q

