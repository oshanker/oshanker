import numpy as np
import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from base_transformer_shanker.data.stringdata import GenerateDataset as GD


np.random.seed(0)
def tokens_to_str(tokens):
	return "".join([chr(x+94) for x in tokens])

def generate_random_string1():
    #len = np.random.randint(10, 20)
    len = 12
    return "".join([chr(x) for x in np.random.randint(97, 97+26, len)])

class GenerateNoMarkerDataset(Dataset):

    
    def __init__(self, n_samples, pad_idx=GD.PAD_IDX,  reverseString=True):
        super(GenerateNoMarkerDataset, self).__init__()
        self.pad_idx = pad_idx
        self.values = [generate_random_string1() for _ in range(n_samples)]
        self.labels = [x[::-1] for x in self.values] if reverseString else self.values

    def __len__(self):
        return len(self.values)  # number of samples in the dataset

    def __getitem__(self, index):
        return self.text_transform(self.values[index].rstrip("\n")), \
            self.text_transform(self.labels[index].rstrip("\n"))
        
    def text_transform(self, x):
        return torch.tensor([ord(z)-97+3 for z in x] )

def runpersonMain():
	train_iter = GenerateNoMarkerDataset(6 )
	dataloader_train = DataLoader(train_iter, batch_size=3)
	s, t = next(iter(dataloader_train))
	print(s[:, ...])
	print(t[:, ...])
	print("s.size", s.size())
	print("?", tokens_to_str(s[0, :]))
	print("?", tokens_to_str(t[0, :]))
	
	for s, t in iter(dataloader_train):
	     print(s.size())



runpersonMain()

