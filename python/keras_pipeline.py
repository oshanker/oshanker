import numpy as np
import pandas
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasRegressor
# python3 -m pip install --user scikit-learn
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

# https://towardsdatascience.com/recurrent-neural-networks-by-example-in-python-ffd204f99470
# https://machinelearningmastery.com/regression-tutorial-keras-deep-learning-library-python/
# load dataset
dataframe = pandas.read_csv("data/housing.data", delim_whitespace=True, header=None)
dataset = dataframe.values
# split into input (X) and output (Y) variables
X = dataset[:,0:13]
print(X.shape)
Y = dataset[:,13]

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

# fix random seed for reproducibility
seed = 7
# evaluate model with standardized dataset
np.random.seed(seed)
estimators = []
estimators.append(('standardize', StandardScaler()))
estimators.append(('mlp', KerasRegressor(build_fn=baseline_model, epochs=50, batch_size=5, verbose=0)))
# a pipeline is also a sklearn Estimator. 
#  the sample above does *not* provide a trained model as output.
# you have to train the pipeline: pipeline.fit(X,Y)
# We can then call .predict()  (same function name, different API),
# http://scikit-learn.org/stable/modules/generated/sklearn.pipeline.Pipeline.html#sklearn.pipeline.Pipeline.predict
pipeline = Pipeline(estimators)
kfold = KFold(n_splits=10, random_state=seed)
results = cross_val_score(pipeline, X, Y, cv=kfold)
print("Standardized: %.2f (%.2f) MSE" % (results.mean(), results.std()))

estimator = estimators[1][1]
X_test = X[0:100,:]
Y_test = Y[0:100]

#you have to fit the estimator again after cross_val_score to evaluate on the new data
#estimator.fit(X, Y)
#prediction = estimator.predict(X_test)
pipeline.fit(X, Y)
prediction = pipeline.predict(X_test)

# https://stackoverflow.com/questions/44132652/keras-how-to-perform-a-prediction-using-kerasregressor
test_error =  np.abs(Y_test - prediction)
mean_error = np.mean(test_error)
min_error = np.min(test_error)
max_error = np.max(test_error)
std_error = np.std(test_error)
print(mean_error, ",  np.mean(test_error)")
print(min_error, ",  np.min(test_error)")
print(max_error, ",  np.max(test_error)")
print(std_error, ",  np.std(test_error)")

prediction = pipeline.predict(X)
train_error =  np.abs(Y - prediction)
mean_error = np.mean(train_error)
min_error = np.min(train_error)
max_error = np.max(train_error)
std_error = np.std(train_error)
print()
print(mean_error, ",  np.mean(train_error)")
print(min_error, ",  np.min(train_error)")
print(max_error, ",  np.max(train_error)")
print(std_error, ",  np.std(train_error)")

t1 = pipeline.predict(X[0:2,:])
t2 = estimator.predict(X[0:2,:])
print(t1, t2, " actual ", Y[0:2])







