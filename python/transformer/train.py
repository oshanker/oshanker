import torch
import time
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from torch.utils.data import DataLoader
from torch.nn.utils.rnn import pad_sequence
from mpl_toolkits.axes_grid1 import ImageGrid
from base_transformer_shanker.data import GenerateDataset as GD
from transformer import Transformer
from constants import args
import base_transformer_shanker.functions as functions
from base_transformer_shanker.data.stringdata1 import GenerateNoMarkerDataset 

# Code is based on
# https://towardsdatascience.com/a-complete-guide-to-write-your-own-transformers-29e23f371ddd 
# https://pytorch.org/tutorials/beginner/introyt/trainingyt.html
"""
model.train() tells your model that you are training the model. This helps inform layers such as Dropout and BatchNorm, which are designed to behave differently during training and evaluation. For instance, in training mode, BatchNorm updates a moving average on each new batch; whereas, for evaluation mode, these updates are frozen.

More details: model.train() sets the mode to train (see source code). You can call either model.eval() or model.train(mode=False) to tell that you are testing. It is somewhat intuitive to expect train function to train model but it does not do that. It just sets the mode.

https://stackoverflow.com/questions/51433378/what-does-model-train-do-in-pytorch
"""

PAD_IDX = GD.PAD_IDX
SOS_IDX = GD.SOS_IDX
EOS_IDX = GD.EOS_IDX

def collate_fn(batch):
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

def train(model, optimizer, dataloader, loss_fn, epoch):
    model.train()
    losses = 0
    acc = 0
    history_loss = []
    history_acc = [] 

    with tqdm(dataloader, position=0, leave=True) as tepoch:
        for x, y in tepoch:
            tepoch.set_description(f"Epoch {epoch}")

            optimizer.zero_grad()
            logits = model(x, y[:, :-1])
            loss = loss_fn(logits.contiguous().view(-1, model.vocab_size), y[:, 1:].contiguous().view(-1))
            loss.backward()
            optimizer.step()
            losses += loss.item()
            
            preds = logits.argmax(dim=-1)
            masked_pred = preds * (y[:, 1:]!=PAD_IDX)
            accuracy = (masked_pred == y[:, 1:]).float().mean()
            acc += accuracy.item()
            
            history_loss.append(loss.item())
            history_acc.append(accuracy.item())
            tepoch.set_postfix(loss=loss.item(), accuracy=100. * accuracy.item())

    length = len(list(dataloader)) 
    
    print("train history_loss",  len(history_loss))
#     print("history_acc", history_acc)
     

    return losses / length, acc / length, history_loss, history_acc


def evaluate(model, dataloader, loss_fn):
    model.eval()
    losses = 0
    acc = 0
    history_loss = []
    history_acc = [] 

    for x, y in tqdm(dataloader, position=0, leave=True):

        logits = model(x, y[:, :-1])
        loss = loss_fn(logits.contiguous().view(-1, model.vocab_size), y[:, 1:].contiguous().view(-1))
        losses += loss.item()
        
        preds = logits.argmax(dim=-1)
        masked_pred = preds * (y[:, 1:]!=PAD_IDX)
        accuracy = (masked_pred == y[:, 1:]).float().mean()
        acc += accuracy.item()
        
        history_loss.append(loss.item())
        history_acc.append(accuracy.item())
    
    length = len(list(dataloader)) 
    
    print("eval history_loss",  len(history_loss))
#     print("history_acc", history_acc)
     

    return losses / length, acc / length, history_loss, history_acc


def runTrain(train_iter, eval_iter, path):
    # Define model here
    model = Transformer(**args)

    # Instantiate datasets
    dataloader_train = DataLoader(train_iter, batch_size=256, collate_fn=collate_fn)
    dataloader_val = DataLoader(eval_iter, batch_size=256, collate_fn=collate_fn)


    # Initialize model parameters
    for p in model.parameters():
        if p.dim() > 1:
            nn.init.xavier_uniform_(p)

    # Define loss function : we ignore logits which are padding tokens
    loss_fn = torch.nn.CrossEntropyLoss(ignore_index=PAD_IDX)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001, betas=(0.9, 0.98), eps=1e-9)
    
    # Save history to dictionnary
    history = {
        'train_loss': [],
        'eval_loss': [],
        'train_acc': [],
        'eval_acc': []
    }

    # Main loop
    epochsToRun = 2
    for epoch in range(1, epochsToRun):
        start_time = time.time()
        train_loss, train_acc, hist_loss, hist_acc = train(model, optimizer, dataloader_train, loss_fn, epoch)
        history['train_loss'] += hist_loss
        history['train_acc'] += hist_acc
        end_time = time.time()
        val_loss, val_acc, hist_loss, hist_acc = evaluate(model, dataloader_val, loss_fn)
        history['eval_loss'] += hist_loss
        history['eval_acc'] += hist_acc
        print((f"Epoch: {epoch}, Train loss: {train_loss:.3f}, Train acc: {train_acc:.3f}, Val loss: {val_loss:.3f}, Val acc: {val_acc:.3f} "f"Epoch time = {(end_time - start_time):.3f}s"))
    
    # history of loss and accuracy for each batch.
    for key in history:
        print(key, len(history[key]))
    functions.plot_list(history['train_acc'][75:],ylabel='train_acc')
#     torch.save(model.state_dict(), path)


def runperson():
    # train_iter = GD(50000, reverseString=False)
    # eval_iter = GD(20000, reverseString=False)
    train_iter = GenerateNoMarkerDataset(50000, reverseString=False)
    eval_iter = GenerateNoMarkerDataset(20000, reverseString=False)
    path = "../out/forward.pt"
    runTrain(train_iter, eval_iter, path)
    
if __name__ == "__main__":
    runperson()
