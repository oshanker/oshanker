import numpy as np
from pandas import read_csv
from matplotlib import pyplot
# load dataset
dataset = read_csv('../../oldriemann/data/gzetaE12/percentile_calc.csv', header=0)
values = dataset.values
# specify columns to plot
groups = np.arange(1, dataset.shape[1], 1)
i = 1
# plot each column
pyplot.figure()
for group in groups:
	pyplot.subplot(len(groups), 1, i)
	pyplot.plot(values[:, group])
	pyplot.title(dataset.columns[group], y=0.5, loc='right')
	i += 1
pyplot.show()
