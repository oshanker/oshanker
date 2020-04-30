import sys
import numpy as np
import math
from pandas import read_csv
from matplotlib import pyplot
from sklearn.linear_model import LinearRegression

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
    
    
    inspect(0)
    inspect(1)
    inspect(groups[-1]-1)

def inspectRow(index):
    z = dataset.at[index,'Unnamed: 0']
    print('values for: ', z)
    
    row0 = dataset.iloc[index][1:].to_frame() #.reset_index(drop=True)
    # print(row0.head())
    # print(row0.tail())
    # summary = row0.describe()
    # print(summary)
    values = row0.values
    # print (values)
    x = []
    angle = []
    for phi in range(0,360,15):
        x.append([math.cos(math.pi*phi/180), math.cos(math.pi*phi/90)])
        angle.append(phi)
    x = np.array(x)
    angle = np.array(angle)
    #x = x.reshape(-1,1)
    y = values[:, 0]
    reg = LinearRegression().fit(x, y)
    r2 = reg.score(x, y)
    #
    coeff = reg.coef_
    intercept = reg.intercept_
    
    
    savedCoeff.append([z, coeff[0], coeff[1], intercept, r2])
    
    pred = np.dot(x, coeff) + intercept
    if doAction == 'plot':
        option = 'save'
        pyplot.figure()
        pyplot.subplot(1, 1, 1)
        pyplot.scatter(angle, pred,color='black')
        pyplot.plot(angle, y)
        title = 'z=' + str(z)
        pyplot.legend(['fit from \nuniversality','actual'])
        pyplot.xlabel('phi')
        pyplot.ylabel('probability density')
        pyplot.title(title, y=0.75, loc='right')
        if(option == 'save'):
            fname = 'z' + str(z).replace('.','') + '.eps'
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
# input generated in math.GSeriesTest.testWriteZetaPhiE12()
# or in math.MoreGSeriesTest.testInterpolate()
savedCoeff = []
runFor = 1
fitCoeff_in = [
    '../../oldriemann/out/gzetaE12/calcHist12.csv',
    '../../oldriemann/out/gzetaE12/calcHist_fine12.csv',
    '../../oldriemann/out/gzetaE28/calcHist12.csv']
outFit = [
    '../../oldriemann/out/gzetaE12/fitCoeff.csv',
    '../../oldriemann/out/gzetaE12/fitCoeff_fine.csv',
    '../../oldriemann/out/gzetaE28/fitCoeff.csv']
dataset = read_csv(fitCoeff_in[runFor], header=0)
dataset.drop('Unnamed: 25', axis = 1, inplace = True)
print('dataset.columns', dataset.columns, dataset.columns.shape)
#dataset.dtypes
print('dataset', type(dataset))
doAction = 'saveCoeff'
if doAction == 'plot':
    inspectRow(19)
elif doAction == 'saveCoeff':
    for index in range(0,39):
        inspectRow(index)
    
    np.savetxt(outFit[runFor], np.asarray(savedCoeff), 
               fmt='%7.2f, %7.3f, %7.3f, %7.3f, %8.5f', delimiter=',')
elif doAction == 'testfit':
    testfit = []
    for index in range(13,26):
        inspectRow(index)
    np.savetxt('../../oldriemann/out/gzetaE12/testfit.csv', np.asarray(testfit), 
               fmt='%7.2f&%7.3f&%7.3f&%7.3f&%7.3f&%7.3f&%7.3f&%7.3f&%7.3f&%7.3f&%7.3f \\\\', 
               delimiter=',')

groups = [1]
#plotDistribution(groups)


