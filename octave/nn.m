%% Machine Learning  - Neural Network Learning
% run by calling nn (it is a script)
%set up a simple neural net, train it and verify.
% normalize y to be between 0 and 1

%% Initialization
clear ; close all; clc

%% Setup the parameters you will use for this exercise
input_layer_size  = 3;  % 
hidden_layer_size = 5;   % 
num_labels = 2;          %  
                          % 
% Load Training Data
fprintf('Loading and Visualizing Data ...\n')
m = 30;
% X - function arguments (rows: number of training examples)
X = 15*(rand(m, input_layer_size)-0.5) ;
% y - function values
y = zeros(m, num_labels);


W1 = weights(input_layer_size, hidden_layer_size);
W2 = weights(hidden_layer_size, num_labels);
[rms, y] = predict(W1, W2, X, y);

% Weight regularization parameter (we set this to 0 here).
lambda = 0.00;

%% ================ Part 6: Initializing Pameters ================
%  In this part of the exercise, you will be starting to implment a two
%  layer neural network.  Initialize the weights of the neural network
%  (randInitializeWeights.m)

fprintf('\nInitializing Neural Network Parameters ...\n')

initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
initial_Theta2 = randInitializeWeights(hidden_layer_size, num_labels);
[rms, pred] = predict(initial_Theta1, initial_Theta2, X, y);
fprintf('initial rms %12.3f\n', rms);
% Unroll parameters
initial_nn_params = [initial_Theta1(:) ; initial_Theta2(:)];

nnCostFunction(initial_nn_params,  input_layer_size, ...
       hidden_layer_size, X, y);

%% =============== Part 7: Implement Backpropagation ===============
%  backpropagation algorithm for the neural network.  nnCostFunction.m returns the partial
%  derivatives of the parameters.
%
fprintf('Checking Backpropagation... \n');
%%{
%  Check gradients by running checkNNGradients

checkNNGradients(lambda);

fprintf('\nProgram paused (Backpropagation). Press enter to continue.\n');
pause;
%%}

%% =================== Part 8: Training NN ===================
%  You have now implemented all the code necessary to train a neural 
%  network. To train your neural network, we will now use "fmincg", which
%  is a function which works similarly to "fminunc". Recall that these
%  advanced optimizers are able to train our cost functions efficiently as
%  long as we provide them with the gradient computations.
%
fprintf('\nTraining Neural Network... \n')

        [Theta1 Theta2] = trainNeuralNet(input_layer_size, ...
             hidden_layer_size,  400, initial_nn_params, X, y, lambda);
                     
%% ================= Part 10: Implement Predict =================
%  After training the neural network, we would like to use it to predict
%  the labels. You will now implement the "predict" function to use the
%  neural network to predict the labels of the training set. This lets
%  you compute the training set accuracy.

[rms, pred] = predict(Theta1, Theta2, X, y);
fprintf('lambda %12.3f final rms %12.3f\n', lambda, rms);
diff = y-pred;
disp('train pred diff')
disp([y pred diff])






%{
Xval = 15*(rand(m, input_layer_size)-0.5) ;
yval = zeros(m, num_labels);
[rms, yval] = predict(W1, W2, Xval, yval)
lambda_vec = [ 0.003 0.01 0.03 0.1 ]';
[error_train, error_val] = ...
    validationCurve(input_layer_size, ...
	     hidden_layer_size, ...
         X, y, Xval, yval, lambda_vec);
close all;
label = 'lambda';
	plot(lambda_vec, error_train, 'r', lambda_vec, error_val);
	legend('Train', 'Cross Validation');
	xlabel(label);
	ylabel('Error');
%}