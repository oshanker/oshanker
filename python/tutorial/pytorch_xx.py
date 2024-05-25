#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 16:40:54 2023

@author: uorugant
"""
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.data as data
import math
import copy


x = torch.rand(5, 3)
print(x)


x = np.array([1, 2, 3])
y = np.array([2, 3, 4])

print(np.dot(x, y))

# https://statistics.berkeley.edu/computing/faqs/git-auth
# https://towardsdatascience.com/a-complete-guide-to-write-your-own-transformers-29e23f371ddd
# https://towardsdatascience.com/build-your-own-transformer-from-scratch-using-pytorch-84c850470dcb