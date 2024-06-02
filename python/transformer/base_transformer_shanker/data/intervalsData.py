
import numpy as np
import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
import os

np.random.seed(0)
def tokens_to_str(tokens):
    return "".join([chr(x+94) for x in tokens])

def generate_random_string1(x):
    #len = np.random.randint(10, 20)
    store = [
        "helloworld",
        "yelloworld",
        "hellospain",
        ]
    return store[x%len(store)]


class IntervalsDataset(Dataset):
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
# =============================================================================
# generate data for model    
# =============================================================================
    S = 10
    L = 10
    def __init__(self, n_samples, path, offset):
        super(IntervalsDataset, self).__init__()
        self.arr = np.loadtxt(path,
                     delimiter=",", dtype=np.int_) + 1
        print("arr ", self.arr.size)
        self.n_samples = n_samples
        self.offset = offset

    def __len__(self):
        #return len(self.arr)-(self.S+self.L)  # number of samples in the dataset
        return self.n_samples  # number of samples in the dataset

    def __getitem__(self, index):
        #return torch.zeros(self.S, dtype=torch.int64),  torch.ones(self.L, dtype=torch.int64)
        my_index = index + self.offset
        return torch.tensor(self.arr[my_index :my_index+self.S], dtype=torch.int64), \
            torch.tensor(self.arr[my_index+self.S :my_index+self.S+self.L], dtype=torch.int64)

def runpersonMain():
    path = "../../../out/intervals.csv"
    train_iter = IntervalsDataset(6, path, 0 )
    dataloader_train = DataLoader(train_iter, batch_size=3)
    s, t = next(iter(dataloader_train))
    print(s[:, ...])
    print(t[:, ...])
    print("s.size", s.size())
    
    i = 1
    for s, t in iter(dataloader_train):
        for step in range(0, s.size()[0]):
            print("in", s[step, :])
            print("out", t[step, :])

if __name__ == "__main__":
    runpersonMain()

