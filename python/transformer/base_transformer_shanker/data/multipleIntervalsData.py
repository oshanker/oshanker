
import numpy as np
import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader

np.random.seed(0)
def tokens_to_str(tokens):
    return "".join([chr(x+94) for x in tokens])



class MultipleIntervalsDataset(Dataset):
    S = 20
    L = 10
    def __init__(self, n_samples, path, offset, S =20, L = 10):
        super(MultipleIntervalsDataset, self).__init__()
        self.arr = np.loadtxt(path,
                     delimiter=",", dtype=np.int_) 
        # max_index = len(self.arr) - S - L
        
        self.n_samples = n_samples
        self.offset = offset
        self.S = S
        self.L = L
        self.row_len = len(self.arr[0,:])
        if self.row_len < self.S+self.L:
            raise IndexError(f"Size out of bounds: row len {self.row_len}")

    def getRowCol(self, index):
        divisor = self.row_len-self.S-self.L+1
        
        quotient = index // divisor
        remainder = index % divisor
        ret = np.array([quotient, remainder])
        return ret


    def __len__(self):
        return self.n_samples  # number of samples in the dataset

    def __getitem__(self, index):
        my_index = index + self.offset
        max_index = my_index+self.S+self.L
        if max_index >= len(self.arr):
            raise IndexError(f"Index out of bounds: {index}")
        return torch.tensor(self.arr[my_index :my_index+self.S], dtype=torch.int64), \
            torch.tensor(self.arr[my_index+self.S :max_index], dtype=torch.int64)

    def __str__(self):
        return f"MultipleIntervalsDataset:[ S: {self.S}, L: {self.L}, \
              length: {self.n_samples}, offset: {self.offset}] "

def runpersonMain():
    path = "../../../out/integers.csv"
    train_iter = MultipleIntervalsDataset(4, path, 0, S = 3, L =5 )
    dataloader_train = DataLoader(train_iter, batch_size=3)
    print("train_iter.size", len(train_iter))
    print (train_iter)
    i = 1
    for index in range(0,5):
        print (train_iter.getRowCol(index))
    # for s, t in iter(dataloader_train):
    #     print(f"=== {i} ===")
    #     i = i + 1
    #     for step in range(0, s.size()[0]):
    #         print("in", s[step, :].detach().numpy() )
    #         print("out", t[step, :].detach().numpy() )

if __name__ == "__main__":
    runpersonMain()

