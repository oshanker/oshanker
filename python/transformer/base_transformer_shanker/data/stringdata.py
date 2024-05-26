import numpy as np
import torch
from torch.utils.data import Dataset


np.random.seed(0)

def generate_random_string():
    #len = np.random.randint(10, 20)
    len = 12
    return "".join([chr(x) for x in np.random.randint(97, 97+26, len)])

class ReverseDataset(Dataset):

    PAD_IDX = 0
    SOS_IDX = 1
    EOS_IDX = 2
    
    def __init__(self, n_samples, pad_idx=PAD_IDX, sos_idx=SOS_IDX, eos_idx=EOS_IDX):
        super(ReverseDataset, self).__init__()
        self.pad_idx = pad_idx
        self.sos_idx = sos_idx
        self.eos_idx = eos_idx
        self.values = [generate_random_string() for _ in range(n_samples)]
        self.labels = [x[::-1] for x in self.values]

    def __len__(self):
        return len(self.values)  # number of samples in the dataset

    def __getitem__(self, index):
        return self.text_transform(self.values[index].rstrip("\n")), \
            self.text_transform(self.labels[index].rstrip("\n"))
        
    def text_transform(self, x):
        return torch.tensor([self.sos_idx] + [ord(z)-97+3 for z in x] + [self.eos_idx])
