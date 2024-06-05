import torch
import time
import torch.nn as nn

from tqdm import tqdm
from torch.utils.data import DataLoader
from transformer import Transformer
import base_transformer_shanker.functions as functions
#from base_transformer_shanker.data.intervalsData import  IntervalsDataset 
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



def train(model, optimizer, dataloader, loss_fn, epoch, ignore_index):
    model.train()
    losses = 0
    acc = 0
    history_loss = []
    history_acc = [] 

    with tqdm(dataloader, position=0, leave=True) as tepoch:
        for x, y in tepoch:
            tepoch.set_description(f"Epoch {epoch}")

            optimizer.zero_grad()
            logits = model(x, y)
            loss = loss_fn(logits.contiguous().view(-1, model.vocab_size), 
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
    
     

    return losses / length, acc / length, history_loss, history_acc


def evaluate(model, dataloader, loss_fn, ignore_index):
    model.eval()
    losses = 0
    acc = 0
    history_loss = []
    history_acc = [] 

    for x, y in tqdm(dataloader, position=0, leave=True):

# =============================================================================
        # the marker disease needs to unravel from here.
        # same should be done for train and evaluate.
        # damn the markers!
# =============================================================================
        logits = model(x, y)
        loss = loss_fn(logits.contiguous().view(-1, model.vocab_size), 
                       y.contiguous().view(-1))
        losses += loss.item()
        
        preds = logits.argmax(dim=-1)
        masked_pred = preds * (y != ignore_index) if ignore_index > -50 else preds
        accuracy = (masked_pred == y).float().mean()
        acc += accuracy.item()
        
        history_loss.append(loss.item())
        history_acc.append(accuracy.item())
    
    length = len(list(dataloader)) 
    
     

    return losses / length, acc / length, history_loss, history_acc

def evaluate1(model, dataloader, loss_fn, ignore_index):
    model.eval()
    losses = 0
    acc = 0

    for x, y in tqdm(dataloader, position=0, leave=True):
        print("x --> ", x.size())
        for row in functions.iterate_rows(x):
            print(functions.tokens_to_str(row))
        print("y (expected) --> ", y.size())
        for row in functions.iterate_rows(y):
            print(functions.tokens_to_str(row))
            
        
        logits = model(x, y)
        loss = loss_fn(logits.contiguous().view(-1, model.vocab_size), 
                       y.contiguous().view(-1))
        losses += loss.item()
        
        preds = logits.argmax(dim=-1)
        masked_pred = preds * (y != ignore_index) if ignore_index > -50 else preds
        accuracy = (masked_pred == y).float().mean()
        print("masked_pred --> ", masked_pred.size())
        for row in functions.iterate_rows(masked_pred):
            print(functions.tokens_to_str(row))
        acc += accuracy.item()
        print("accuracy --> ", accuracy.detach().numpy())
        
    
    length = len(list(dataloader)) 
    

    return losses / length, acc / length

def evaluate2(model, dataloader, loss_fn, ignore_index):
    model.eval()
    losses = 0
    acc = 0

    for x, y in tqdm(dataloader, position=0, leave=True):
        print("x --> ", x.size())
        for row in functions.iterate_rows(x):
            print((row))
        print("y (expected) --> ", y.size())
            
        
        logits = model(x, y)
        functions.plot_list(logits[0,0,:].detach().numpy())
        #print(logits[0,0,:].detach().numpy())
        loss = loss_fn(logits.contiguous().view(-1, model.vocab_size), 
                       y.contiguous().view(-1))
        losses += loss.item()
        
        preds = logits.argmax(dim=-1)
        masked_pred = preds * (y != ignore_index) if ignore_index > -50 else preds
        accuracy = (masked_pred == y).float().mean()
        
        for step in range(0, y.size()[0]):
            actual = y[step, :].detach().numpy()-1
            mypred = masked_pred[step, :].detach().numpy() - 1
            print("actual    ", actual)
            print("prediction", mypred)

        acc += accuracy.item()
    
    length = len(list(dataloader)) 
    

    return losses / length, acc / length


def runTrain(train_iter, eval_iter, path, ignore_index: int = -100):
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
    # Define model here
    model = Transformer(**args)

    # Instantiate datasets
    collate_fn=functions.gd_collate_fn if ignore_index > -50 else None
    dataloader_train = DataLoader(train_iter, batch_size=256, collate_fn=collate_fn)
    dataloader_val = DataLoader(eval_iter, batch_size=256, collate_fn=collate_fn)


    # Initialize model parameters
    for p in model.parameters():
        if p.dim() > 1:
            nn.init.xavier_uniform_(p)

    # Define loss function : we ignore logits which are padding tokens
    loss_fn = torch.nn.CrossEntropyLoss()
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
        train_loss, train_acc, hist_loss, hist_acc = train(model, optimizer, 
                                dataloader_train, loss_fn, epoch, ignore_index)
        history['train_loss'] += hist_loss
        history['train_acc'] += hist_acc
        end_time = time.time()
        
        val_loss, val_acc, hist_loss, hist_acc = evaluate(model, dataloader_val, loss_fn, ignore_index)
        history['eval_loss'] += hist_loss
        history['eval_acc'] += hist_acc
        print((f"Epoch: {epoch}, Train loss: {train_loss:.3f}, Train acc: {train_acc:.3f}, Val loss: {val_loss:.3f}, Val acc: {val_acc:.3f} "f"Epoch time = {(end_time - start_time):.3f}s"))
        
        # evaluate2(model, dataloader_val, loss_fn, ignore_index)
        # print((f"Epoch: {epoch}, Train loss: {train_loss:.3f}, Train acc: {train_acc:.3f},  "f"Epoch time = {(end_time - start_time):.3f}s"))
    
    functions.plot_multiple_lists(history['train_acc'][5:], history['eval_acc'][5:],
                                labels=['train_acc','eval_acc'])
#     torch.save(model.state_dict(), path)


def runIntervalTrain(train_iter, train_iter_1, eval_iter, path, ignore_index: int = -100):
    # Define model here
    model = Transformer(**args)

    # Instantiate datasets
    dataloader_train = DataLoader(train_iter, batch_size=256)
    dataloader_train_1 = DataLoader(train_iter_1, batch_size=256)
    dataloader_val = DataLoader(eval_iter, batch_size=256)


    # Initialize model parameters
    for p in model.parameters():
        if p.dim() > 1:
            nn.init.xavier_uniform_(p)

    # Define loss function : we ignore logits which are padding tokens
    loss_fn = torch.nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001, betas=(0.9, 0.98), eps=1e-9)
    
    # Save history to dictionnary
    history = {
        'train_loss': [],
        'eval_loss': [],
        'train_acc': [],
        'eval_acc': []
    }

    # Main loop
    print("=== ROUND 1 ===")
    epochsToRun = 2
    for epoch in range(1, epochsToRun):
        start_time = time.time()
        train_loss, train_acc, hist_loss, hist_acc = train(model, optimizer, 
                                dataloader_train, loss_fn, epoch, ignore_index)
        history['train_loss'] += hist_loss
        history['train_acc'] += hist_acc
        end_time = time.time()
        
        val_loss, val_acc, hist_loss, hist_acc = evaluate(model, dataloader_val, 
                                                          loss_fn, ignore_index)
        history['eval_loss'] += hist_loss
        history['eval_acc'] += hist_acc
        print((f"Epoch: {epoch}, Train loss: {train_loss:.3f}, Train acc: {train_acc:.3f}, Val loss: {val_loss:.3f}, Val acc: {val_acc:.3f} "f"Epoch time = {(end_time - start_time):.3f}s"))
        

    functions.plot_multiple_lists(history['train_acc'][5:], history['eval_acc'][5:],
                                  xlabel='Batch Iteration', ylabel='Accuracy', 
                                labels=['train_acc','eval_acc'], title='First round of training')

    # Save history to dictionnary
    history = {
        'train_loss': [],
        'eval_loss': [],
        'train_acc': [],
        'eval_acc': []
    }

    # ROUND 2
    print("=== ROUND 2 ===")
    epochsToRun = 2
    for epoch in range(1, epochsToRun):
        start_time = time.time()
        train_loss, train_acc, hist_loss, hist_acc = train(model, optimizer, 
                                dataloader_train_1, loss_fn, epoch, ignore_index)
        history['train_loss'] += hist_loss
        history['train_acc'] += hist_acc
        end_time = time.time()
        
        val_loss, val_acc, hist_loss, hist_acc = evaluate(model, dataloader_val, 
                                                          loss_fn, ignore_index)
        history['eval_loss'] += hist_loss
        history['eval_acc'] += hist_acc
        print((f"Epoch: {epoch}, Train loss: {train_loss:.3f}, Train acc: {train_acc:.3f}, Val loss: {val_loss:.3f}, Val acc: {val_acc:.3f} "f"Epoch time = {(end_time - start_time):.3f}s"))
        
    print("dataloader_val dataset", dataloader_val.dataset)
    print("=== DONE ===")
    #evaluate2(model, dataloader_val, loss_fn, ignore_index)

    functions.plot_multiple_lists(history['train_acc'][5:], history['eval_acc'][5:],
                                  xlabel='Batch Iteration', ylabel='Accuracy', 
                                labels=['train_acc','eval_acc'], title='Second round of training')
#     torch.save(model.state_dict(), path)


def runperson():
    reverseString=True
    path = "../out/stringreverse" + str(reverseString) + ".pt" 
    #train_iter = GenerateNoMarkerDataset(6400, reverseString=reverseString)
    train_iter = GenerateNoMarkerDataset(12800, reverseString=reverseString)
    eval_iter = GenerateNoMarkerDataset(12800, reverseString=reverseString)
    #eval_iter = FixedDataset(3, reverseString=reverseString, drop = 0)
    runTrain(train_iter, eval_iter, path)    
    
    # S=30
    # L = 20
    # path = "../out/intervals.csv"
    # train_iter = IntervalsDataset(6400, path, 0,  S = S, L = L)
    # train_iter_1 = IntervalsDataset(8400, path, 6400, S = S, L = L)
    # eval_iter = IntervalsDataset(10006, path, 14800+100, S = S, L = L)
    
    # path = "../out/intervals.pt" 
    # runIntervalTrain(train_iter, train_iter_1, eval_iter, path)
    # print('train_iter', train_iter)
    # print('train_iter_1', train_iter_1)
    # print('eval_iter', eval_iter)
    
if __name__ == "__main__":
    runperson()
