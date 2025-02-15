import torch
import base_transformer_shanker.functions as functions
from transformer import Transformer
from base_transformer_shanker.constants import args
from base_transformer_shanker.data.intervalsData import  IntervalsDataset 
from torch.utils.data import DataLoader
from base_transformer_shanker.data.multipleIntervalsData import  MultipleIntervalsDataset 
import numpy as np


class Predict():
    """
    https://pytorch.org/tutorials/beginner/saving_loading_models.html
    """

    def __init__(self, transformer):
        self.model = transformer
    
    
    def evaluate(self, dataloader, ignore_index: int = -100, 
                 filename=None, printerrors=False):
        self.model.eval()
        loss_fn = torch.nn.CrossEntropyLoss()
        losses = 0
        acc = 0
        history_loss = []
        history_acc = [] 
        rowcount = 0
        error_count = 0
        file = None if filename is None else open(filename, 'w')
        for x, y in dataloader:
    
            logits = self.model(x, y)
            loss = loss_fn(logits.contiguous().view(-1, self.model.vocab_size), 
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
    
    def evaluate1(self, dataloader, ignore_index: int = -100):
        self.model.eval()
        acc = 0
    
        for x, y in dataloader:
            print("x --> ", x.size())
                
            
            logits = self.model(x, y)
            
            preds = logits.argmax(dim=-1)
            masked_pred = preds * (y != ignore_index) if ignore_index > -50 else preds
            accuracy = (masked_pred == y).float().mean()
            print("masked_pred --> ", masked_pred.size())
            acc += accuracy.item()
            print("accuracy --> ", accuracy.detach().numpy())
            test = np.array([1,1,1,1,0,2,1,1,1,1])
            for step in range(0, y.size()[0]):
                input = x[step, :].detach().numpy()
                #do_print = np.array_equal(test, input)
                do_print = True
                if (do_print):
                    actual = y[step, :].detach().numpy()
                    mypred = masked_pred[step, :].detach().numpy() 
                    
                    print("----------", step)
                    print("input     ", np.array2string(input, separator=","))
                    print("actual    ", np.array2string(actual, separator=","))
                    print("prediction", np.array2string(mypred, separator=","))
            
        
    
    
    def evaluate2(self, dataloader, filename=None):
        self.model.eval()
    
        for x, y in dataloader:
            #  def _call_impl(self, *input, **kwargs) in nn.module
            # https://stephencowchau.medium.com/pytorch-module-call-vs-forward-c4df3ff304b1
            # If we don't have any hooks, we want to skip the rest of the logic in
            # this function, and just call forward.
            logits = self.model(x, y)
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
                print("----------")
                print("actual    ", actual)
                print("prediction", mypred)
    

def runperson():
    # Define model here
    np.random.seed(0)
    model = Transformer(**args)
    model_path = "../data/intervals4500_500.pt" 
    model.load_state_dict(torch.load(model_path))
    model.eval()
    L = 10
    train = Predict(model)
    # path = "../data/intervalsTestE28Threes.csv"
    # eval_iter_1 = MultipleIntervalsDataset(4, path, 466, S = S, L = L)
    # dataloader_val_1 = DataLoader(eval_iter_1, batch_size=256)
    
    # print("=== TEST ===")
    # train.evaluate2( dataloader_val_1)    
    # eval_iter = MultipleIntervalsDataset(10000, path, 9600, S = S, L = L)
    # dataloader_val = DataLoader(eval_iter, batch_size=256)
    # val_loss, val_acc, hist_loss, hist_acc = train.evaluate(dataloader_val, 
    #                                     filename = '../out/errorsthree.csv')
    # print((f"  Val loss: {val_loss:.3f}, Val acc: {val_acc:.6f} "))
    
    #path = "../data/intervalsTestE28.csv"
    path = "../data/intervalsE12.csv"
    print("=====================================")
    print(f"regular intervals {path}")
    print("=====================================")
    # eval_iter = IntervalsDataset(200000, path, 1000, S = S, L = L)
    # dataloader_val = DataLoader(eval_iter, batch_size=256)
    # val_loss, val_acc, hist_loss, hist_acc = train.evaluate(dataloader_val, 
    #                                     filename = '../out/errorsthree.csv')
    # print((f"  Val loss: {val_loss:.3f}, Val acc: {val_acc:.6f} "))

    eval_iter_1 = IntervalsDataset(10, path, 0, S = 0, L = L)
    dataloader_val_1 = DataLoader(eval_iter_1, batch_size=256)
    train.evaluate1( dataloader_val_1)    
    #train.evaluate2( dataloader_val_1)    
       
    
    
if __name__ == "__main__":
    runperson()
    print("=====================================")
    print("- predict was run!")
    print("=====================================")
