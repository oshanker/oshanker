#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 12:06:35 2025

@author: uorugant
"""
import math
from temp import inv

primes = [2, 3, 5, 7, 11, 13, 17]

def init(p, q, e):
    n = p*q
    phi = (p-1)*(q-1)
    d = inv(phi,e)
    if math.isnan(d):
        print ("invalid e", e)
    return(n, phi, d)

def message(m, n, e, d):
    print('-------')
    print('message', m)
    
    c = 1
    for i in range(e):
        c = (c*m)%n
    print('c', c)
    decoded = 1
    for i in range(d):
        decoded = (decoded*c)%n
    print('decoded', decoded)
    
def main():
    p = 11
    q = 13
    e = 37
    init(p, q, 3)
    n, phi, d = init(p, q, e)
    print(n, phi, d)
    message(5, n, e, d)
    message(17, n, e, d)

if __name__   == '__main__':
     main()