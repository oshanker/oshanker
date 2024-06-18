import torch
import time
import torch.nn as nn
import numpy as np
import os

from tqdm import tqdm
from torch.utils.data import DataLoader
from transformer import Transformer
import base_transformer_shanker.functions as functions
from base_transformer_shanker.data.intervalsData import  IntervalsDataset 
from base_transformer_shanker.data.multipleIntervalsData import  MultipleIntervalsDataset 
from base_transformer_shanker.data.stringdata1 import GenerateNoMarkerDataset
from base_transformer_shanker.data.fixedData import  FixedDataset 
from base_transformer_shanker.constants import args 

# Code is based on
# https://towardsdatascience.com/a-complete-guide-to-write-your-own-transformers-29e23f371ddd 
# https://pytorch.org/tutorials/beginner/introyt/trainingyt.html
"""
model.train() tells your model that you are training the model. This helps inform layers such as Dropout and BatchNorm, which are designed to behave differently during training and evaluation. For instance, in training mode, BatchNorm updates a moving average on each new batch; whereas, for evaluation mode, these updates are frozen.

More details: model.train() sets the mode to train (see source code). You can call either model.eval() or model.train(mode=False) to tell that you are testing. It is somewhat intuitive to expect train function to train model but it does not do that. It just sets the mode.

https://stackoverflow.com/questions/51433378/what-does-model-train-do-in-pytorch
"""

# Define model here
my_model = Transformer(**args)

class Train():
    def __init__(self, model):
        self.model = model

        # Initialize model parameters
        for p in self.model.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)
    
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=0.001, 
                                          betas=(0.9, 0.98), eps=1e-9)

    def run(self, dataloader_train, dataloader_val, 
            ignore_index, title, filename=None, epochsToRun = 2):
        loss_fn = torch.nn.CrossEntropyLoss()
        # Save history to dictionnary
        history = {
            'train_loss': [],
            'eval_loss': [],
            'train_acc': [],
            'eval_acc': []
        }
    
        for epoch in range(1, epochsToRun):
            start_time = time.time()
            train_loss, train_acc, hist_loss, hist_acc = self.train( self.optimizer, 
                                    dataloader_train, loss_fn, epoch, ignore_index)
            history['train_loss'] += hist_loss
            history['train_acc'] += hist_acc
            end_time = time.time()
            filename1 = filename
            if filename is not None:
                directory = os.path.dirname(filename)
                basefilename = os.path.basename(filename)
                filename1 = directory + '/' + str(epoch) + basefilename
            val_loss, val_acc, hist_loss, hist_acc = self.evaluate(self.model, dataloader_val, 
                                    loss_fn, ignore_index, filename1)
            history['eval_loss'] += hist_loss
            history['eval_acc'] += hist_acc
            print((f"Epoch: {epoch}, Train loss: {train_loss:.3f}, Train acc: {train_acc:.3f}, Val loss: {val_loss:.3f}, Val acc: {val_acc:.6f} "f"Epoch time = {(end_time - start_time):.3f}s"))
        
        functions.plot_multiple_lists(history['train_acc'][5:], history['eval_acc'][5:],
                                      xlabel='Batch Iteration', ylabel='Accuracy', 
                                    labels=['train_acc','eval_acc'], title=title)
    

    def train(self, optimizer, dataloader, loss_fn, epoch, ignore_index):
        self.model.train()
        losses = 0
        acc = 0
        history_loss = []
        history_acc = [] 
    
        with tqdm(dataloader, position=0, leave=True) as tepoch:
            for x, y in tepoch:
                tepoch.set_description(f"Epoch {epoch}")
    
                optimizer.zero_grad()
                logits = self.model(x, y)
                loss = loss_fn(logits.contiguous().view(-1, self.model.vocab_size), 
                               y.contiguous().view(-1))
                loss.backward()
                optimizer.step()
                losses += loss.item()
                
                preds = logits.argmax(dim=-1)
                masked_pred = preds * (y != ignore_index) if ignore_index > -50 else preds
                accuracy = (masked_pred == y).float().mean()
                acc += accuracy.item()
                
                history_loss.append(loss.item())
                history_acc.append(accuracy.item())
                tepoch.set_postfix(loss=loss.item(), accuracy=100. * accuracy.item())
    
        length = len(list(dataloader)) 
        if ( length == 0): length = 1
        return losses / length, acc / length, history_loss, history_acc
    
    
    def evaluate(self, model, dataloader, loss_fn, ignore_index, 
                 filename=None, printerrors=False):
        model.eval()
        losses = 0
        acc = 0
        history_loss = []
        history_acc = [] 
        rowcount = 0
        error_count = 0
        file = None if filename is None else open(filename, 'w')
        for x, y in tqdm(dataloader, position=0, leave=True):
    
            logits = model(x, y)
            loss = loss_fn(logits.contiguous().view(-1, model.vocab_size), 
                           y.contiguous().view(-1))
            losses += loss.item()
            
            preds = logits.argmax(dim=-1)
            masked_pred = preds * (y != ignore_index) if ignore_index > -50 else preds
            accuracy = (masked_pred == y)
            
            L = y.size()[1]
            
            if filename is not None:
                empty_tensor = np.empty([0, L+2], dtype=int) 
        
                for rowidx in range(0, y.size()[0]):
                    rowcount = rowcount + 1
                    row_error_count = L - accuracy[rowidx,:].int().detach().numpy().sum()
                    if row_error_count > 0:
                        error_count = error_count + row_error_count
                        yy = np.concatenate(
                            (y[rowidx,:].int().detach().numpy(),
                                             np.array([100, row_error_count]) ) )
                        tocat = [yy]
                        zz = np.concatenate((masked_pred[rowidx,:].int().detach().numpy(),
                                             np.array([200, 0]) ))
                        tocat_1 = [zz]
                        empty_tensor = np.concatenate((empty_tensor, tocat,tocat_1), 
                                                      axis=0)
                functions.write_integers_to_open_file(empty_tensor, file)
            accuracy = accuracy.float().mean()
            acc += accuracy.item()
            
            history_loss.append(loss.item())
            history_acc.append(accuracy.item())
        
        if filename is not None:
            print("=====================================")
            print("error_count, rowcount", error_count, rowcount)
            print("=====================================")
            file.close()
        length = len(list(dataloader)) 
        return losses / length, acc / length, history_loss, history_acc
    
    def evaluate1(self, eval_iter, ignore_index, collate_fn):
        self.model.eval()
        loss_fn = torch.nn.CrossEntropyLoss()
        losses = 0
        acc = 0
        dataloader = DataLoader(eval_iter, batch_size=256, collate_fn=collate_fn)
    
        for x, y in tqdm(dataloader, position=0, leave=True):
            print("x --> ", x.size())
            for row in functions.iterate_rows(x):
                print(functions.tokens_to_str(row))
            print("y (expected) --> ", y.size())
            for row in functions.iterate_rows(y):
                print(functions.tokens_to_str(row))
                
            
            logits = self.model(x, y)
            loss = loss_fn(logits.contiguous().view(-1, self.model.vocab_size), 
                           y.contiguous().view(-1))
            losses += loss.item()
            
            preds = logits.argmax(dim=-1)
            masked_pred = preds * (y != ignore_index) if ignore_index > -50 else preds
            print("validate ", (masked_pred == y).float())
            accuracy = (masked_pred == y).float().mean()
            print("masked_pred --> ", masked_pred.size())
            for row in functions.iterate_rows(masked_pred):
                print(functions.tokens_to_str(row))
            acc += accuracy.item()
            print("accuracy --> ", accuracy.detach().numpy())
            
        
        length = len(list(dataloader)) 
        return losses / length, acc / length
    
    def evaluate2(self, dataloader, filename=None):
        self.model.eval()
    
        for x, y in dataloader:
            logits = self.model(x, y)
            #print(logits[0,0,:].detach().numpy())
            preds = logits.argmax(dim=-1)
            masked_pred = preds  
            accuracy = (masked_pred == y)
            
            L = y.size()[1]
            
            if filename is not None:
                empty_tensor = np.empty([0, L+2], dtype=int) 
        
                for rowidx in range(0, y.size()[0]):
                    print("--", rowidx, accuracy[rowidx,:].int().detach().numpy().sum())
                    if rowidx < 2:
                        yy = np.concatenate((y[rowidx,:].int().detach().numpy(),
                                             np.array([100, rowidx*100])))
                        tocat = [yy]
                        zz = np.concatenate((masked_pred[rowidx,:].int().detach().numpy(),
                                             np.array([200, rowidx*200])))
                        tocat_1 = [zz]
                        empty_tensor = np.concatenate((empty_tensor, tocat,tocat_1), 
                                                      axis=0)
                functions.write_integers_to_file(empty_tensor, 
                                                  filename)
            accuracy = accuracy.float().mean()
            print('accuracy', accuracy.detach().numpy() )
            for step in range(0, y.size()[0]):
                actual = y[step, :].detach().numpy()
                mypred = masked_pred[step, :].detach().numpy() 
                print("actual    ", actual)
                print("prediction", mypred)

def runTrain(train_iter, train_iter_1, eval_iter, eval_iter_1,
             path, ignore_index: int = -100):
    """
    B = Batch size
    S = Source sequence length
    L = Target sequence length
    E = Model dimension

    encode Input
        x: (B, S) with elements in (0, C) where C is num_classes
    encode Output
        (B, S, E) embedding
    
    decode Input
        encoded_x: (B, S, E)
        y: (B, L) with elements in (0, C) where C is num_classes
    decode Output
        (B, L, C) logits
    """
    train = Train(my_model)

    # Instantiate datasets
    collate_fn=functions.gd_collate_fn if ignore_index > -50 else None
    dataloader_train = DataLoader(train_iter, batch_size=256, collate_fn=collate_fn)
    dataloader_train_1 = DataLoader(train_iter_1, batch_size=256, collate_fn=collate_fn)
    dataloader_val = DataLoader(eval_iter, batch_size=256, collate_fn=collate_fn)

    print("=== ROUND 1 ===")
    train.run( dataloader_train, dataloader_val, ignore_index=ignore_index,
        title='First round of training')
    print("=== TEST 1 ===")
    train.evaluate1( eval_iter_1, ignore_index, collate_fn)
    print("=== ROUND 2 ===")
    train.run( dataloader_train_1, dataloader_val, ignore_index=ignore_index,
        title='Second round of training')
    print("=== TEST ===")
    train.evaluate1( eval_iter_1, ignore_index, collate_fn)
#     torch.save(model.state_dict(), path)


def runIntervalTrain(train_iter, train_iter_1, eval_iter, eval_iter_1, 
                     path, ignore_index: int = -100, epochsToRun = 2):
    train = Train(my_model)
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


def runstring():
    reverseString=True
    path = "../out/stringreverse" + str(reverseString) + ".pt" 
    train_iter = GenerateNoMarkerDataset(6400, reverseString=reverseString)
    train_iter_1 = GenerateNoMarkerDataset(6400, reverseString=reverseString)
    eval_iter = GenerateNoMarkerDataset(12800, reverseString=reverseString)
    eval_iter_1 = FixedDataset(3, reverseString=reverseString, drop = 0)
    runTrain(train_iter, train_iter_1, eval_iter, eval_iter_1, path)    
    
def runperson():
    S=10
    L = 10
    path = "../out/intervals.csv"
    train_iter = IntervalsDataset(6400, path, 0,  S = S, L = L)
    train_iter_1 = IntervalsDataset(8000, path, 6400, S = S, L = L)
    eval_iter = IntervalsDataset(10000, path, 14400+75, S = S, L = L)
    eval_iter_1 = IntervalsDataset(3, path, 14400+75, S = S, L = L)
    
    path = "../out/intervals.pt" 
    runIntervalTrain(train_iter, train_iter_1, eval_iter, eval_iter_1, path)
    
def runThrees():
    S=10
    L = 10
    path = "../out/intervalsTest.csv"
    train_iter = MultipleIntervalsDataset(1133, path, 0,  S = S, L = L)
    train_iter_1 = MultipleIntervalsDataset(333, path, 733, S = S, L = L)
    eval_iter = MultipleIntervalsDataset(500, path, 1466, S = S, L = L)
    eval_iter_1 = MultipleIntervalsDataset(3, path, 1470, S = S, L = L)
    
    path = "../out/intervals.pt" 
    runIntervalTrain(train_iter, train_iter_1, eval_iter, 
                     eval_iter_1, path, epochsToRun = 4)
    
if __name__ == "__main__":
    #runperson()
    runThrees()
    runperson()
    print("- train was run!")
