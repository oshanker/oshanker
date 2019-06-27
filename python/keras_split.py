import numpy as np
import pandas as pd
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasRegressor
# python3 -m pip install --user scikit-learn
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

from sklearn.model_selection import train_test_split

# https://machinelearningmastery.com/multivariate-time-series-forecasting-lstms-keras/
# https://machinelearningmastery.com/regression-tutorial-keras-deep-learning-library-python/
# load dataset
dataframe = pd.read_csv("data/housing.data", delim_whitespace=True, header=None)
dataset = dataframe.values
# split into input (X) and output (Y) variables
X = dataset[:,0:13]
Y = dataset[:,13]

# fix random seed for reproducibility
seed = 7
np.random.seed(seed)


X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=.2)

#The efficient ADAM optimization algorithm is used and a mean squared error loss function is optimized.
# define base model
def baseline_model():
	# create model
	model = Sequential()
	model.add(Dense(13, input_dim=13, kernel_initializer='normal', activation='relu'))
	model.add(Dense(1, kernel_initializer='normal'))
	# Compile model
	model.compile(loss='mean_squared_error', optimizer='adam')
	return model

# evaluate model with standardized dataset
estimator = KerasRegressor(build_fn=baseline_model, epochs=100, batch_size=5, verbose=0)

estimator.fit(X_train, Y_train)
prediction = estimator.predict(X_test)

test_error =  np.abs(Y_test - prediction)
mean_error = np.mean(test_error)
min_error = np.min(test_error)
max_error = np.max(test_error)
std_error = np.std(test_error)
print(mean_error, ",  np.mean(test_error)")
print(min_error, ",  np.min(test_error)")
print(max_error, ",  np.max(test_error)")
print(std_error, ",  np.std(test_error)")

prediction = estimator.predict(X_train)
train_error =  np.abs(Y_train - prediction)
mean_error = np.mean(train_error)
min_error = np.min(train_error)
max_error = np.max(train_error)
std_error = np.std(train_error)
print()
print(mean_error, ",  np.mean(train_error)")
print(min_error, ",  np.min(train_error)")
print(max_error, ",  np.max(train_error)")
print(std_error, ",  np.std(train_error)")





