import torch
import time
import torch.nn as nn

from tqdm import tqdm
from torch.utils.data import DataLoader
from transformer import Transformer
import base_transformer_shanker.functions as functions
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
# =============================================================================
        # the marker disease needs to unravel from here.
        # same should be done for train and evaluate.
# Epoch: 1, Train loss: 2.016, Train acc: 0.475,  Epoch time = 4.564s
# =============================================================================
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
    
    print("train history_loss",  len(history_loss))
#     print("history_acc", history_acc)
     

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
    
    print("eval history_loss",  len(history_loss))
#     print("history_acc", history_acc)
     

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
# =============================================================================
#    y -->  torch.Size([3, 10])
#    masked_pred -->  torch.Size([3, 9])
#    here are we seeing the implicit EOS/SOS remnant? 
# =============================================================================
        masked_pred = preds * (y != ignore_index) if ignore_index > -50 else preds
        accuracy = (masked_pred == y).float().mean()
        print("masked_pred --> ", masked_pred.size())
        for row in functions.iterate_rows(masked_pred):
            print(functions.tokens_to_str(row))
        acc += accuracy.item()
        print("accuracy --> ", accuracy)
        
    
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
    loss_fn = torch.nn.CrossEntropyLoss(ignore_index=ignore_index)
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
        
        # val_loss, val_acc, hist_loss, hist_acc = evaluate(model, dataloader_val, loss_fn, ignore_index)
        # history['eval_loss'] += hist_loss
        # history['eval_acc'] += hist_acc
        # print((f"Epoch: {epoch}, Train loss: {train_loss:.3f}, Train acc: {train_acc:.3f}, Val loss: {val_loss:.3f}, Val acc: {val_acc:.3f} "f"Epoch time = {(end_time - start_time):.3f}s"))
        
        evaluate1(model, dataloader_val, loss_fn, ignore_index)
        print((f"Epoch: {epoch}, Train loss: {train_loss:.3f}, Train acc: {train_acc:.3f},  "f"Epoch time = {(end_time - start_time):.3f}s"))
    
    # history of loss and accuracy for each batch.
    functions.plot_list(history['train_acc'][5:],ylabel='train_acc')
#     torch.save(model.state_dict(), path)


def runperson():
    # train_iter = GD(50000, reverseString=False)
    # eval_iter = GD(20000, reverseString=False)
    reverseString=True
    train_iter = GenerateNoMarkerDataset(12800, reverseString=reverseString)
    eval_iter = FixedDataset(3, reverseString=reverseString, drop = 0)
    path = "../out/reverse.pt" if reverseString else "../out/forward.pt"
    runTrain(train_iter, eval_iter, path)
    
if __name__ == "__main__":
    runperson()
