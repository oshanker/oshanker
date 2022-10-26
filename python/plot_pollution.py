import sys
from pandas import read_csv
from matplotlib import pyplot


def main():
    # load dataset
    print("-------")
    print(sys.argv)
    # ~/junk/oshanker/oshanker/python/data/pollution.csv
    file_name =  sys.argv[1]     
    dataset = read_csv(file_name, header=0, index_col=0)
    print(dataset)
    values = dataset.values
    # specify columns to plot
    groups = [0, 1, 2, 3, 5, 6, 7]
    i = 1
    # plot each column
    pyplot.figure()
    for group in groups:
    	pyplot.subplot(len(groups), 1, i)
    	pyplot.plot(values[:, group])
    	pyplot.title(dataset.columns[group], y=0.5, loc='right')
    	i += 1
    pyplot.show()
    

if __name__   == '__main__':
     main()