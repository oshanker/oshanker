function [error_train, error_val, Theta1, Theta2] = ...
    learningCurve( ...
    initial_nn_params, hidden_layer_size, X, y, Xval, yval, sizes, lambda, iter)
    
%LEARNINGCURVE Generates the train and cross validation set errors needed 
%to plot a learning curve
%   [error_train, error_val] = ...
%       LEARNINGCURVE(X, y, Xval, yval, lambda) returns the train and
%       cross validation set errors for a learning curve. In particular, 
%       it returns two vectors of the same length - error_train and 
%       error_val. Then, error_train(i) contains the training error for
%       i examples (and similarly for error_val(i)).
%
%   In this function, you will compute the train and test errors for
%   dataset sizes from 1 up to m. 

% Number of training examples
m = size(sizes, 2);
num_labels = size(y, 2);            
input_layer_size = (size(initial_nn_params, 1) - num_labels)...
     /hidden_layer_size - num_labels -1;

error_train = zeros(m, 1);
error_val   = zeros(m, 1);
%maxIter = 500;
		if exist('iter', 'var')   
           maxIter = iter;
		else
           maxIter = 100;
        end
idx = 1;
   for i = sizes

		if exist('lambda', 'var') && lambda > 0 
		    disp('lambda')
			[Theta1 Theta2] = trainNeuralNet( 
					input_layer_size, ...
					hidden_layer_size, ...
					maxIter, initial_nn_params,	X(1:i, :), y(1:i, :), lambda);
		else
			[Theta1 Theta2] = trainNeuralNet(
				 input_layer_size, ...
			  	 hidden_layer_size, ...
				 maxIter, initial_nn_params, X(1:i, :), y(1:i, :));
		end
		error_train(idx) =  predict(Theta1, Theta2, ...
						X(1:i, :), y(1:i, :));           
		error_val(idx) =  predict(Theta1, Theta2, ...
						 Xval, yval); 
		idx++;          
    end


%


end
