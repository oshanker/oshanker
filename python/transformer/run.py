import torch
import numpy as np

def runpersonMain():
    arr = np.loadtxt("../out/intervals.csv",
                 delimiter=",", dtype=np.int32)
    print(arr)    
    tensor = torch.tensor(arr, dtype=torch.int32)
    print('size', tensor.size(), 'sum', torch.sum(tensor))
    idx = 0
    tensor = torch.tensor(arr[idx :idx+10], dtype=torch.int32)
    print('size', tensor.size(), 'sum', torch.sum(tensor))
    tensor = torch.tensor(arr[idx+10 :idx+15], dtype=torch.int32)
    print('size', tensor.size(), 'sum', torch.sum(tensor))
    



if __name__ == "__main__":
    runpersonMain()


# export PYTHONPATH="$PWD"
# python3 run.py