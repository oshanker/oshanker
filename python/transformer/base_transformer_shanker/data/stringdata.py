import numpy as np
import torch
from torch.utils.data import Dataset
from base_transformer_shanker.constants import PAD_IDX, SOS_IDX, EOS_IDX
#import base_transformer_shanker.constants 
from torch.utils.data import DataLoader
from base_transformer_shanker.functions import gd_collate_fn, tokens_to_str

np.random.seed(0)

def generate_random_string():
    len = np.random.randint(10, 20)
    #len = 12
    return "".join([chr(x) for x in np.random.randint(97, 97+26, len)])

class GenerateDataset(Dataset):
    
    def __init__(self, n_samples, pad_idx=PAD_IDX, sos_idx=SOS_IDX, eos_idx=EOS_IDX, reverseString=True):
        super(GenerateDataset, self).__init__()
        self.pad_idx = pad_idx
        self.sos_idx = sos_idx
        self.eos_idx = eos_idx
        self.values = [generate_random_string() for _ in range(n_samples)]
        self.labels = [x[::-1] for x in self.values] if reverseString else self.values

    def __len__(self):
        return len(self.values)  # number of samples in the dataset

    def __getitem__(self, index):
        return self.text_transform(self.values[index].rstrip("\n")), \
            self.text_transform(self.labels[index].rstrip("\n"))
        
    def text_transform(self, x):
        return torch.tensor([self.sos_idx] + [ord(z)-97+3 for z in x] + [self.eos_idx])

def runpersonMain():
	train_iter = GenerateDataset(6, pad_idx=PAD_IDX, sos_idx=SOS_IDX, eos_idx=EOS_IDX)
	dataloader_train = DataLoader(train_iter, batch_size=3, collate_fn=gd_collate_fn)
	s, t = next(iter(dataloader_train))
	print(s[:, ...])
	print(t[:, ...])
	print("s.size", s.size())
	print("?", tokens_to_str(s[0, :]))
	print("?", tokens_to_str(t[0, :]))
	
	for s, t in iter(dataloader_train):
	     print(s.size())

if __name__ == "__main__":
    runpersonMain()

