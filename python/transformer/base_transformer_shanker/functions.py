#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 13:26:02 2024


@author: uorugant
"""

import matplotlib.pyplot as plt
import random
from torch.nn.utils.rnn import pad_sequence
from base_transformer_shanker.constants import PAD_IDX


def iterate_rows(tensor):
    """
    Iterate over the rows of a PyTorch tensor.

    Args:
    - tensor (torch.Tensor): The input tensor

    Yields:
    - torch.Tensor: Each row of the input tensor
    """
    for row in tensor:
        yield row

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

def tokens_to_str(tokens):
	return "".join([chr(x+94) for x in tokens])


def runplots():
    # my_list = [1, 4, 3, 4, 1]
    # plot_list(my_list)
    list1 = [1, 2, 3, 4, 5]
    list2 = [2, 3, 4, 5, 6]
    list3 = [3, 4, 5, 6, 7]
    plot_multiple_lists(list1, list2, list3, labels=['List 1', 'List 2', 'List 3'])
    
def choose_symbol(symbols, probabilities):
    """
    Chooses a symbol from a set of symbols according to a given probability distribution.
    
    Args:
        symbols (list): List of symbols to choose from.
        probabilities (list): List of probabilities corresponding to each symbol.
        
    Returns:
        symbol: The chosen symbol.
    """
    return random.choices(symbols, weights=probabilities, k=1)[0]

def gd_collate_fn(batch):
    """ 
    This function pads inputs with PAD_IDX to have batches of equal length
    """
    src_batch, tgt_batch = [], []
    for src_sample, tgt_sample in batch:
        src_batch.append(src_sample)
        tgt_batch.append(tgt_sample)

    src_batch = pad_sequence(src_batch, padding_value=PAD_IDX, batch_first=True)
    tgt_batch = pad_sequence(tgt_batch, padding_value=PAD_IDX, batch_first=True)
    return src_batch, tgt_batch



def runperson():
    symbols = ['A', 'B', 'C']
    probabilities = [0.3, 0.4, 0.3]
    
    chosen_symbol = choose_symbol(symbols, probabilities)
    print("Chosen symbol:", chosen_symbol)
	
if __name__ == "__main__":
    runperson()
    
# python3 base_transformer_shanker/functions.py
