#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 11:02:57 2020

@author: shankero
"""


import numpy as np
from sklearn.linear_model import LinearRegression
X = np.array([[0,1],[1, 1], [1, 2], [2, 2], [2, 3]])
# y = 1 * x_0 + 2 * x_1 + 3
y = np.dot(X, np.array([1, 2])) + 3
reg = LinearRegression().fit(X, y)
print('Score ', reg.score(X, y))
#1.0
coeff = reg.coef_
intercept = reg.intercept_
print('X*', coeff, '+', intercept)
print('y:', y)
print('check:', np.dot(X, coeff) + intercept)
#array([1., 2.])
#3.0000...
print('predict',reg.predict(np.array([[3, 5]])))
#array([16.])
