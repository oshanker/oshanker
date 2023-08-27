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

def sample(n):
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

