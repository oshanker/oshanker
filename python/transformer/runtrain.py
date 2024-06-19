#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 14:13:24 2024

@author: uorugant
"""
import torch

from train import Train
from base_transformer_shanker.data.multipleIntervalsData import  MultipleIntervalsDataset 
from base_transformer_shanker.constants import args 
from transformer import Transformer
from base_transformer_shanker.data.intervalsData import  IntervalsDataset 
from torch.utils.data import DataLoader
import numpy as np



class RunTrain(Train):
    def __init__(self, model, model_loaded = False):
        super(RunTrain, self).__init__(model, model_loaded)
        print("init RunTrain")


def runThreesIntervalTrain(train, train_iter, train_iter_1, eval_iter, eval_iter_1, 
                     model_path, ignore_index: int = -100, epochsToRun = 2):
    dataloader_train = DataLoader(train_iter, batch_size=256)
    dataloader_train_1 = DataLoader(train_iter_1, batch_size=256)
    dataloader_val = DataLoader(eval_iter, batch_size=256)
    dataloader_val_1 = DataLoader(eval_iter_1, batch_size=256)

    print("=== ROUND 1 ===")
    train.run( dataloader_train, dataloader_val, 
              ignore_index=ignore_index,
        title='First round of training', epochsToRun = epochsToRun)
    
    print("=== ROUND 2 ===")
    train.run( dataloader_train_1, dataloader_val, 
              ignore_index=ignore_index,
        title='Second round of training', filename = '../out/errorsthree.csv', 
        epochsToRun = 2)
    
    print("=== TEST ===")
    train.evaluate2( dataloader_val_1)
    torch.save(train.model.state_dict(), model_path)


def my_runIntervalTrain(train, train_iter, train_iter_1, eval_iter, eval_iter_1, 
                     path, ignore_index: int = -100, epochsToRun = 2):
    dataloader_train = DataLoader(train_iter, batch_size=256)
    dataloader_train_1 = DataLoader(train_iter_1, batch_size=256)
    dataloader_val = DataLoader(eval_iter, batch_size=256)
    dataloader_val_1 = DataLoader(eval_iter_1, batch_size=256)

    print("=== ROUND 1 ===")
    train.run( dataloader_train, dataloader_val, 
              ignore_index=ignore_index,
        title='First round of training', epochsToRun = epochsToRun)
    
    print("=== ROUND 2 ===")
    train.run( dataloader_train_1, dataloader_val, 
              ignore_index=ignore_index,
        title='Second round of training', filename = '../out/errors.csv', 
        epochsToRun = 2)
    
    print("=== TEST ===")
    train.evaluate2( dataloader_val_1)     
    torch.save(train.model.state_dict(), path)


    
def runperson(train):
    S=10
    L = 10
    path = "../data/intervalsTestE28.csv"
    train_iter = IntervalsDataset(100, path, 0,  S = S, L = L)
    train_iter_1 = IntervalsDataset(600, path, 100, S = S, L = L)
    eval_iter = IntervalsDataset(200000, path, 1000, S = S, L = L)
    eval_iter_1 = IntervalsDataset(4, path, 585, S = S, L = L)
    
    model_path = "../out/intervalsE28.pt" 
    my_runIntervalTrain(train, train_iter, train_iter_1, eval_iter, 
                        eval_iter_1, model_path)
    print("=====================================")
    print(path)
    print("=====================================")
        
    
def runThrees(train):
    S=10
    L = 10
    path = "../data/intervalsTestE28Threes.csv"
    train_iter = MultipleIntervalsDataset(3800, path, 0,  S = S, L = L)
    train_iter_1 = MultipleIntervalsDataset(2600, path, 7000, S = S, L = L)
    eval_iter = MultipleIntervalsDataset(10000, path, 9600, S = S, L = L)
    eval_iter_1 = MultipleIntervalsDataset(4, path, 466, S = S, L = L)
    
    model_path = "../out/intervalsThree.pt" 
    runThreesIntervalTrain(train, train_iter, train_iter_1, eval_iter, 
                     eval_iter_1, model_path, epochsToRun = 4)
    print("=====================================")
    print(path)
    print("=====================================")
    
if __name__ == "__main__":
    np.random.seed(0)
    my_model = Transformer(**args)
    model_loaded = False
    if model_loaded:
        # we are saying model loaded, so load model!
        # this actually overtrains
        model_path = "../data/intervalsThree.pt" 
        my_model.load_state_dict(torch.load(model_path))
    #runperson()
    train = RunTrain(my_model, model_loaded)
    runThrees(train)
    #runperson(train)
    print("runtrain was run!")
