import torch
import numpy as np

def runpersonMain():
    arr = np.loadtxt("../out/intervals.csv",
                 delimiter=",", dtype=np.int32)
    print(arr)    
    tensor = torch.tensor(arr, dtype=torch.int32)
    print(tensor)
    print(torch.sum(tensor))



if __name__ == "__main__":
    runpersonMain()


# export PYTHONPATH="$PWD"
# python3 run.py