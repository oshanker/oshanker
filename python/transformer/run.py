import torch
import numpy as np
import base_transformer_shanker.functions as functions

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
    #runpersonMain()
    integer_lists = [
           [1, 2, 3, 4, 5, 6, 7, 8, 9],
           [ 10, 11, 12, 13, 14, 15, 16, 17, 18]
        ]
    functions.write_integers_to_file(integer_lists, '../out/integers.csv')
    #functions.write_integers_to_file(integer_lists)
    print("Integers have been written to 'integers.csv'")

# export PYTHONPATH="$PWD"
# python3 run.py