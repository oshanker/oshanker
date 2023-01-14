#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 21:56:47 2022

@author: uorugant
"""
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    print("-------")
    # ~/junk/oshanker/oshanker/oldriemann/data/gzetaE12/zerosE12.csv
    data =  [
            [2 , 2.268576, 1.504897, 0.999090, 0.663992, 0.441400 , ],
            [3 , 3.624520, 1.896329, 0.998864, 0.526210, 0.274796 , ],
            [4 , 5.588850, 2.367189, 1.001220, 0.426549, 0.178314 , ],
            [5 , 8.849518, 2.923269, 1.011907, 0.343228, 0.115100 , ],
            [6 , 14.373004, 3.728921, 1.008597, 0.266224, 0.070801 , ],
            [7 , 23.961623, 4.721631, 0.974761, 0.205256, 0.041497 , ],
            [8 , 51.790514, 6.714706, 0.983726, 0.149338, 0.018210 , ],
            [9 , 116.632353, 10.129730, 0.975052, 0.103749, 0.008631 , ],
            [10 , 361.615385, 13.306122, 0.949264, 0.057506, 0.004515 , ],
            [11 , 1397.000000, 34.533333, 1.027888, 0.041000, 0.001084 , ],
        ]
    x = np.arange(-0.2, 0.3, 0.1)
    df = pd.DataFrame(data)
    #df = df.T
    df.drop(df.columns[0], axis=1, inplace=True)
    print(df)
    plt.figure()
    # https://rowannicholls.github.io/python/graphs/symbols_linestyles_colours.html
    plt.plot(x, df.iloc[0], 
             linestyle='-',
             color='darkblue',
             #marker='o',
             label="2")
    plt.plot(x, df.iloc[1], 
             linestyle='-',
             color='springgreen',
             #marker='v',
             label="3")
    plt.plot(x, df.iloc[2], 
             linestyle='-',
             color='black',
             #marker='o',
             label="4")
    plt.plot(x, df.iloc[3], 
             linestyle='-',
             color='red',
             #marker='o',
             label="5")
    plt.plot(x, df.iloc[4], 
             linestyle='-',
             color='springgreen',
             #marker='v',
             label="6")
    plt.plot(x, df.iloc[5], 
             linestyle='-',
             color='black',
             #marker='o',
             label="7")
    plt.plot(x, df.iloc[6], 
             linestyle='-',
             color='darkblue',
             #marker='^',
             label="8")
    plt.plot(x, df.iloc[7], 
             linestyle='-',
             color='darkgoldenrod',
             #marker='+',
             label="9")
    plt.plot(x, df.iloc[8], 
             linestyle='-',
             color='red',
             #marker='^',
             label="10")
    plt.plot(x, df.iloc[9], 
             linestyle='-',
             color='darkgoldenrod',
             #marker='+',
             label="11")
    plt.grid(True)
    plt.xlabel('displacement')
    plt.ylabel('ratio')
    plt.ylim([0, 10])

    plt.title('sharp transition')
    plt.legend()
    plt.show()


if __name__   == '__main__':
     main()