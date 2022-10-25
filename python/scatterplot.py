#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 21:56:47 2022

@author: uorugant
"""
import sys
from pandas import read_csv
import matplotlib.pyplot as plt

def main():
    print("-------")
    print(sys.argv)
    # ~/junk/oshanker/oshanker/oldriemann/data/gzetaE12/zerosE12.csv
    file_name =  sys.argv[1]     
    dataset = read_csv(file_name, skip_blank_lines=True, header=None)
    print(dataset)
    plot = dataset.plot(x=0, y=1)
    plot.grid()


if __name__   == '__main__':
     main()