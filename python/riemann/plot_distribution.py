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
    


# python -i plot_distribution.py 
# load dataset
dataset = read_csv('../../oldriemann/out/gzetaE12/calcHist12.csv', header=0)
print('dataset.columns', dataset.columns, dataset.columns.shape)
print('dataset', type(dataset))
values = dataset.values
print('values', type(values))

# specify columns to plot
groups = np.arange(1, dataset.shape[1], 1)
print('groups', type(groups))
print(groups)
i = 1

inspect(0)
inspect(1)
inspect(groups[-1]-1)

# col = values[:, groups[0]]
# print('col', type(col), col.shape)
# print(col)
# colMean = np.mean(col)
# colsd = np.std(col)
# sys.stdout.write("Mean = " + '\t' + str(colMean) + '\t\t' +
#             "Standard Deviation = " + '\t ' + str(colsd) + "\n")



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

