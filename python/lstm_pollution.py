import sys
from math import sqrt
from numpy import concatenate
from matplotlib import pyplot
from pandas import read_csv
from pandas import DataFrame
from pandas import concat
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import mean_squared_error
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.utils.vis_utils import plot_model

import pandas as pd
import numpy as np
 
# https://machinelearningmastery.com/multivariate-time-series-forecasting-lstms-keras/
# convert series to supervised learning
def series_to_supervised(df, n_in=1, n_out=1, dropnan=True):
	n_vars = df.shape[1]
	cols, names = list(), list()
	# input sequence (t-n, ... t-1)
	for i in range(n_in, 0, -1):
		cols.append(df.shift(i))
		names += [('var%d(t-%d)' % (j+1, i)) for j in range(n_vars)]
	# forecast sequence (t, t+1, ... t+n)
	for i in range(0, n_out):
		cols.append(df.shift(-i))
		if i == 0:
			names += [('var%d(t)' % (j+1)) for j in range(n_vars)]
		else:
			names += [('var%d(t+%d)' % (j+1, i)) for j in range(n_vars)]
	# put it all together
	agg = concat(cols, axis=1)
	agg.columns = names
	# drop rows with NaN values
	if dropnan:
		agg.dropna(inplace=True)
	return agg
 
def main():
    print("-------")
    print(sys.argv)
    # ~/junk/oshanker/oshanker/python/data/
    file_name =  sys.argv[1]     
    # load dataset
    dataset = read_csv(file_name+'pollution.csv', header=0, index_col=0)
    values = dataset.values
    # integer encode direction
    encoder = LabelEncoder()
    values[:,4] = encoder.fit_transform(values[:,4])
    # ensure all data is float
    values = values.astype('float32')
    dataset = DataFrame(values)
    
    
    # normalize features
    scaler_x = MinMaxScaler(feature_range=(0, 1))
    scaler_y = MinMaxScaler(feature_range=(0, 1))
    steps = 1
    # frame as supervised learning
    reframed = series_to_supervised( dataset, steps, 1)
    # drop columns we don't want to predict
    reframed.drop(reframed.columns[[9,10,11,12,13,14,15]], axis=1, inplace=True)
    print(reframed.head())
     
    # split into train and test sets
    values = reframed.values
    n_train_hours = 365 * 24
    train = values[:n_train_hours, :]
    test = values[n_train_hours:, :]
    # split into input and outputs
    #drop pollution as predictor if start = 1
    start = 0
    train_X = scaler_x.fit_transform(train[:, start:-1])
    train_y = scaler_y.fit_transform(train[:, -1].reshape(train.shape[0], 1))
    test_X = scaler_x.fit_transform(test[:, start:-1])
    test_y = scaler_y.fit_transform(test[:, -1].reshape(test.shape[0], 1))
    # reshape input to be 3D [samples, timesteps, features]
    train_X = train_X.reshape((train_X.shape[0], 1, train_X.shape[1]))
    test_X = test_X.reshape((test_X.shape[0], 1, test_X.shape[1]))
    print('train_X.shape', train_X.shape, train_y.shape, 
          'test_X.shape', test_X.shape, test_y.shape)
    print(train_X[0, :], train_y[0])
     
    # design network
    model = Sequential()
    model.add(LSTM(50, input_shape=(train_X.shape[1], train_X.shape[2]),activation='relu'))
    model.add(Dense(32,activation='relu'))
    model.add(Dense(1))
    print(model.summary())
    plot_model(model, to_file='out/model_lstm.png', show_shapes=True, show_layer_names=True)

    model.compile(loss='mae', optimizer='adam')
    # fit network
    epochs = 20
    history = model.fit(train_X, train_y, epochs=epochs, batch_size=72, validation_data=(test_X, test_y), verbose=2, shuffle=False)
    # plot history
    def plot_history():
    	pyplot.plot(history.history['loss'], label='train')
    	pyplot.plot(history.history['val_loss'], label='test')
    	pyplot.legend()
    	pyplot.show()
    plot_history()
    val = input("Enter your value: ")
    if (val == "q"):
        return;
 
    # make a prediction
    plotsize = 10000
    yhat = model.predict(test_X)
    # invert scaling for forecast
    inv_yhat = scaler_y.inverse_transform(yhat)
    inv_yhat = inv_yhat[:,0]
    forecast = inv_yhat[0:plotsize]
    #print("forecast ", inv_yhat[0:5])
    # invert scaling for actual
    test_y = test_y.reshape((len(test_y), 1))
    #inv_y = concatenate((test_y, test_X[:, (1-start):]), axis=1)
    inv_y = scaler_y.inverse_transform(test_y)
    inv_y = inv_y[:,0]
    actual = inv_y[0:plotsize]
    #print("actual ", actual)
    # calculate RMSE
    rmse = sqrt(mean_squared_error(inv_y, inv_yhat))
    print('Test RMSE: %.3f' % rmse, ' actual mean %.3f' % inv_y.mean(), 
       ' prediction mean %.3f' % inv_yhat.mean())
    joined = np.array([forecast, actual])
    joined = joined.transpose()
    df = pd.DataFrame(joined,  columns =["forecast", "actual"])
    # save to file
    df.to_csv(file_name + 'pollution_out.csv')

    values=df.values
    groups = [0, 1]
    # plot each column
    pyplot.figure()
    for group in groups:
    	pyplot.plot(values[:, group], label=df.columns[group])
    
    pyplot.legend()
    pyplot.show()
    

if __name__   == '__main__':
     main()