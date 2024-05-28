#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 13:26:02 2024

@author: uorugant
"""

import matplotlib.pyplot as plt

def plot_list(data, xlabel='Index', ylabel='Value'):
    """
    Plot a list of values.

    Parameters:
        data (list): The list of values to be plotted.
    """
    plt.plot(data)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('Plot of List')
    plt.grid(True)
    plt.show()

def runperson():
# Example usage:
    my_list = [1, 4, 3, 4, 1]
    plot_list(my_list)
	
if __name__ == "__main__":
    runperson()
    
# python3 base_transformer_shanker/functions.py
