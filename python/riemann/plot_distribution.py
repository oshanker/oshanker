import sys
import numpy as np
from pandas import read_csv
from matplotlib import pyplot

def inspect(colIndex):
    print('col header', dataset.columns[colIndex] )
    col = values[:, colIndex]
    print('col', colIndex, type(col), col.shape)
    print(col)
    print('sum', np.sum(col))
    colMean = np.mean(col)
    colsd = np.std(col)
    sys.stdout.write("Mean = " + '\t' + str(colMean) + '\t\t' +
                "Standard Deviation = " + '\t ' + str(colsd) + "\n")
    
def getColPhi():
    # specify columns to plot
    groups = np.arange(1, dataset.shape[1], 1)
    print('groups', type(groups))
    print(groups)
    i = 1
    
    inspect(0)
    inspect(1)
    inspect(groups[-1]-1)

    


# python -i plot_distribution.py 
# load dataset
dataset = read_csv('../../oldriemann/out/gzetaE12/calcHist12.csv', header=0)
dataset.drop('Unnamed: 25', axis = 1, inplace = True)
print('dataset.columns', dataset.columns, dataset.columns.shape)
#dataset.dtypes
print('dataset', type(dataset))
print('values for: ', dataset.at[0,'Unnamed: 0'])
row0 = dataset.iloc[0][1:].to_frame() #.reset_index(drop=True)
print(row0.head())
print(row0.tail())
summary = row0.describe()
print(summary)

# values = dataset.values
# print('values', type(values))


# plot each column
'''
pyplot.figure()
for group in groups:
	pyplot.subplot(len(groups), 1, i)
	pyplot.plot(values[:, group])
	pyplot.title(dataset.columns[group], y=0.5, loc='right')
	i += 1
pyplot.show()
'''

