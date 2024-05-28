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
    
def plot_multiple_lists(*lists, xlabel='Index', ylabel='Value', labels=None):
    """
    Plot multiple lists on the same graph.

    Parameters:
        *lists (list): Variable length list of lists to be plotted.
        labels (list): Optional list of labels for each list.
    """
    for i, data in enumerate(lists):
        if labels:
            label = labels[i]
        else:
            label = f'Data {i+1}'
        plt.plot(data, label=label)
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('Plot of Multiple Lists')
    plt.legend()
    plt.grid(True)
    plt.show()

#

def runperson():
# Example usage:
    # my_list = [1, 4, 3, 4, 1]
    # plot_list(my_list)
    list1 = [1, 2, 3, 4, 5]
    list2 = [2, 3, 4, 5, 6]
    list3 = [3, 4, 5, 6, 7]
    
    plot_multiple_lists(list1, list2, list3, labels=['List 1', 'List 2', 'List 3'])

	
if __name__ == "__main__":
    runperson()
    
# python3 base_transformer_shanker/functions.py
