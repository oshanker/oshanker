# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# https://towardsdatascience.com/simple-and-multiple-linear-regression-in-python-c928425168f9 
import math

import numpy as np
import pandas as pd
import os


from sklearn.datasets import fetch_california_housing
from sklearn import linear_model 
from sklearn import datasets ## imports datasets from scikit-learn

def starter():
    data = fetch_california_housing()    
    df = pd.DataFrame(data.data, columns=data.feature_names)
    
    target = pd.DataFrame(data.target, columns=["MEDV"])
    
    print( target)
    
    
    X = df
    y = target["MEDV"]
    
    lm = linear_model.LinearRegression()
    model = lm.fit(X,y)
    predictions = lm.predict(X)
    print(predictions[0:5])
    
    print(lm.score(X,y))
    print(lm.coef_)
    print(lm.intercept_)

def main():
    print(os.getcwd())
    a = 161
    b = 28
    inverse = inv(a,b)
    print(inverse)

def inv(a, b):
    q = a//b
    r = a%b
    row = np.array([a,	b,	q,	r,	1,	0,	1,	0,	1,	-q])
    print(row)
    for i in range(100):
        row = bezout(row)
        print(i, row)
        if row[3] == 0:
            break
    gcd = row[1]
    print("gcd", gcd)
    
    t = row[8]
    s = row[5]
    print(s, t, "s × a + t × b", s * a + t * b)
    if gcd == 1:
        inv = row[8]%a
        print('inv', inv, a, b, (inv*b)%a)
        return inv
    else:
        return math.nan

def bezout(row):
    # see https://extendedeuclideanalgorithm.com/xea.php
    if row[3] == 0:
        print("divide by zero", row)
    newrow = np.zeros(10, dtype=int)
    newrow[0] = row[1]
    newrow[1] = row[3]
    newrow[2] = newrow[0]//newrow[1]
    newrow[3] = newrow[0]%newrow[1]
    newrow[4] = row[5]
    newrow[5] = row[6]
    newrow[6] = newrow[4] - newrow[2]*newrow[5]
    
    newrow[7] = row[8]
    newrow[8] = row[9]
    newrow[9] = newrow[7] - newrow[2]*newrow[8]
    
    return newrow

if __name__   == '__main__':
     main()