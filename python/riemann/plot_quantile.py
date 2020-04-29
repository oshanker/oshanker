import sys
import numpy as np
import math
from pandas import read_csv
from matplotlib import pyplot
from sklearn.linear_model import LinearRegression

def inspect(colIndex):
    print('col header', dataset.columns[colIndex] )
    quantile = values[0, colIndex]
    print('quantile', quantile)
    col = values[1:, colIndex]
    angle = values[1:, 0]*360/col.shape[0]
    print('angle',angle)
    print('col index', colIndex, type(col), col.shape)
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
    
    
    inspect(0)
    inspect(2)
    inspect(groups[-1]-1)

def inspectRow(index):
    z = values[0, index]
    print('values for: ', z)
    debug = True
    
    row0 = dataset.iloc[index][1:].to_frame() #.reset_index(drop=True)
    x = []
    col = values[1:, index]
    angle = values[1:, 0]*360/col.shape[0]
    for phi in angle:
        x.append([math.cos(math.pi*phi/180), math.cos(math.pi*phi/90)])
    x = np.array(x)
    #x = x.reshape(-1,1)
    y = col
    reg = LinearRegression().fit(x, y)
    #print(reg.)
    r2 = reg.score(x, y)
    #
    coeff = reg.coef_
    intercept = reg.intercept_
    
    pred = np.dot(x, coeff) + intercept
    if debug:
        print(angle[3], x[3],y[3])
        print(coeff, intercept)
        print(x[3,0]*coeff[0] + x[3,1]*coeff[1] + intercept, pred[3])
    savedCoeff.append([z, coeff[0], coeff[1], intercept, r2])
    
    if doAction == 'plot':
        option = 'save'
        pyplot.figure()
        pyplot.subplot(1, 1, 1)
        pyplot.scatter(angle, pred,color='black')
        pyplot.plot(angle, y)
        title = 'f=' + str(z)
        pyplot.legend(['fit','actual'])
        pyplot.xlabel('phi')
        pyplot.ylabel('quantile')
        pyplot.title(title, y=0.75, loc='right')
        if(option == 'save'):
            fname = 'f' + str(z).replace('.','') + '.eps'
            pyplot.savefig('../../oldriemann/out/gzetaE12/' + fname)
        else:
            pyplot.show()
    elif doAction == 'testfit':
        out = [z]
        for i in [0, 2, 3, 4, 6]:
            out.append(y[i])
            out.append(pred[i])
        print(len(out))
        testfit.append(out)

    
def plotDistribution(groups):
    
    values = dataset.values
    # print('values', type(values))
    
    i = 1
    # plot each column
    pyplot.figure()
    for group in groups:
    	pyplot.subplot(len(groups), 1, i)
    	pyplot.plot(values[:, 0], values[:, group])
    	pyplot.title(dataset.columns[group], y=0.75, loc='right')
        #pyplot.
    	i += 1
    pyplot.show()
    


# python -i plot_distribution.py 
# load dataset
savedCoeff = []
dataset = read_csv('../../oldriemann/out/gzetaE12/percentile_calc.csv', header=0)
#dataset.drop('Unnamed: 25', axis = 1, inplace = True)
print('dataset.columns', dataset.columns, dataset.columns.shape)
#dataset.dtypes
print('dataset', type(dataset))
print(dataset)
doAction = 'saveCoeff'
values = dataset.values
if doAction == 'inspect':
    getColPhi()
elif doAction == 'plot':
    inspectRow(dataset.shape[1]-2)
elif doAction == 'saveCoeff':
    for index in range(2,dataset.shape[1]-1):
        inspectRow(index)
    
    np.savetxt('../../oldriemann/out/gzetaE12/fitQuantileCoeff.csv', np.asarray(savedCoeff), 
               fmt='%7.2f, %7.3f, %7.3f, %7.3f, %8.5f', delimiter=',')
elif doAction == 'testfit':
    testfit = []
    for index in range(13,26):
        inspectRow(index)
    np.savetxt('../../oldriemann/out/gzetaE12/testQuantilefit.csv', np.asarray(testfit), 
               fmt='%7.2f&%7.3f&%7.3f&%7.3f&%7.3f&%7.3f&%7.3f&%7.3f&%7.3f&%7.3f&%7.3f \\\\', 
               delimiter=',')

groups = [1]
#plotDistribution(groups)


