function [Theta1 Theta2] = trainNeuralNet(
          input_layer_size, ...
          hidden_layer_size, ...
          maxIter, initial_nn_params, X, y, lambda, wt)
% lambda is Weight regularization parameter
% wt is a set of weights to be used in calculating the cost function.

% Initialize Theta
		num_labels = size(y, 2);            
if exist('wt', 'var')  
   wt1 = wt;
else
   wt1 = [];
end 

% Create "short hand" for the cost function to be minimized
    if exist('lambda', 'var')  
        costFunction = @(p) nnCostFunction(p, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                    X, y, lambda, wt1);
    else
        costFunction = @(p) nnCostFunction(p, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                    X, y, 0, wt1);
    end

train_options = optimset('MaxIter', maxIter, 'GradObj', 'on');

% Minimize using fmincg
[nn_params, fx, it] = fmincg(costFunction, initial_nn_params, train_options);
%[nn_params, fx, it] = fminunc(costFunction, initial_nn_params, train_options);


% Obtain Theta1 and Theta2 back from nn_params
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));


end
