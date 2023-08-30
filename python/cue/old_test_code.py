#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 22:15:25 2023

@author: uorugant
"""
import math
import numpy as np

import cmath
import matplotlib.pyplot as plt
import my_functions

conv=180/math.pi


def used():
    q = my_functions.sampleQ()
    w = np.linalg.eigvals(q)
    eigenlist = []
    for e in w:
        eigenlist.append( cmath.phase(e))
    eigen = np.array(eigenlist)
    return eigen

def test3():
    eigen = used()
    sum1 = np.sum(eigen)
    x = np.array([0, 1, 2])*2*math.pi/3 + sum1/3-math.pi
    print(x*conv)
    for theta in x:
       f = my_functions.eval_factor(theta, eigen)
       print('f', f)
       z = my_functions.Z(theta, eigen)
       print('z', z)
    
def test2():
    eigen = used()
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
    eigen = used()
    for theta in x:
        # print("eval at", math.degrees(theta))
        # print(eval(theta, eigen))
        # print(complex(1-math.cos(3*theta), math.sin(3*theta)))
        f = my_functions.eval_direct(theta, eigen)
        y1list.append(f.real)
        y2list.append(f.imag)
        epsilon = 0.000001
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
 