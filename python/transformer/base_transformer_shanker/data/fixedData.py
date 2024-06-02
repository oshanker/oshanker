
import numpy as np
import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from base_transformer_shanker.constants import PAD_IDX

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

def drop_last_n_chars(strings, n):
    return [string[:-n] for string in strings]


class FixedDataset(Dataset):
# =============================================================================
# generate data for model    
# =============================================================================
    def __init__(self, n_samples, pad_idx=PAD_IDX,  reverseString=True, drop=0):
        print('drop', drop)
        super(FixedDataset, self).__init__()
        self.pad_idx = pad_idx
        self.values = [generate_random_string1(x) for x in range(n_samples)]
        self.labels = [x[::-1] for x in self.values] if reverseString else self.values
        if  drop > 0:
            self.labels = drop_last_n_chars(self.labels, drop)

    def __len__(self):
        return len(self.values)  # number of samples in the dataset

    def __getitem__(self, index):
        return self.text_transform(self.values[index].rstrip("\n")), \
            self.text_transform(self.labels[index].rstrip("\n"))
        
    def text_transform(self, x):
        return torch.tensor([ord(z)-97+3 for z in x] )

def runpersonMain():
    train_iter = FixedDataset(6, drop = 2 )
    dataloader_train = DataLoader(train_iter, batch_size=3)
    s, t = next(iter(dataloader_train))
    print(s[:, ...])
    print(t[:, ...])
    print("s.size", s.size())
    
    for s, t in iter(dataloader_train):
        for step in range(0, s.size()[0]):
            print("?", tokens_to_str(s[step, :]))
            print("?", tokens_to_str(t[step, :]))

if __name__ == "__main__":
    runpersonMain()

