#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 14:13:24 2024

@author: uorugant
"""
from train import Train, runIntervalTrain
from base_transformer_shanker.data.multipleIntervalsData import  MultipleIntervalsDataset 
from base_transformer_shanker.constants import args 
from transformer import Transformer
from base_transformer_shanker.data.intervalsData import  IntervalsDataset 
from torch.utils.data import DataLoader



class RunTrain(Train):
    def __init__(self, model):
        super(RunTrain, self).__init__(model)
        print("init RunTrain")


def runThreesIntervalTrain(train_iter, train_iter_1, eval_iter, eval_iter_1, 
                     path, ignore_index: int = -100, epochsToRun = 2):
    my_model = Transformer(**args)
    train = RunTrain(my_model)
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
#     torch.save(model.state_dict(), path)

    
def runperson():
    S=10
    L = 10
    path = "../out/intervalsTestE28.csv"
    train_iter = IntervalsDataset(6400, path, 0,  S = S, L = L)
    train_iter_1 = IntervalsDataset(8000, path, 6400, S = S, L = L)
    eval_iter = IntervalsDataset(10000, path, 14400+75, S = S, L = L)
    eval_iter_1 = IntervalsDataset(3, path, 14400+75, S = S, L = L)
    
    path = "../out/intervals.pt" 
    runIntervalTrain(train_iter, train_iter_1, eval_iter, eval_iter_1, path)
        
    
def runThrees():
    S=10
    L = 10
    path = "../out/intervalsTestE28Threes.csv"
    train_iter = MultipleIntervalsDataset(1133, path, 0,  S = S, L = L)
    train_iter_1 = MultipleIntervalsDataset(333, path, 733, S = S, L = L)
    eval_iter = MultipleIntervalsDataset(500, path, 1466, S = S, L = L)
    eval_iter_1 = MultipleIntervalsDataset(3, path, 1470, S = S, L = L)
    
    path = "../out/intervals.pt" 
    runThreesIntervalTrain(train_iter, train_iter_1, eval_iter, 
                     eval_iter_1, path, epochsToRun = 4)
    
if __name__ == "__main__":
    #runperson()
    #runThrees()
    runperson()
