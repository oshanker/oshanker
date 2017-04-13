function [error_train, error_val] = ...
    validationCurve(
	   input_layer_size,  hidden_layer_size, ...
       X, y, Xval, yval, lambda_vec)
%VALIDATIONCURVE Generate the train and validation errors needed to
%plot a validation curve that we can use to select lambda
%   [lambda_vec, error_train, error_val] = ...
%       VALIDATIONCURVE(X, y, Xval, yval) returns the train
%       and validation errors (in error_train, error_val)
%       for different values of lambda. You are given the training set (X,
%       y) and validation set (Xval, yval).
%

error_train = zeros(length(lambda_vec), 1);
error_val = zeros(length(lambda_vec), 1);
m = size(X,1);      

% ====================== YOUR CODE HERE ======================
%  return training errors in 
%               error_train and the validation errors in error_val. The 
%               vector lambda_vec contains the different lambda parameters 
%               to use for each calculation of the errors, i.e, 
%               error_train(i), and error_val(i) should give 
%               you the errors obtained after training with 
%               lambda = lambda_vec(i)
%
% Note: You can loop over lambda_vec with the following:
%
   for i = 1:length(lambda_vec)
        lambda = lambda_vec(i);
		num_labels = size(y, 2);            
		initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
		initial_Theta2 = randInitializeWeights(hidden_layer_size, num_labels);
		initial_nn_params = [initial_Theta1(:) ; initial_Theta2(:)];
        fprintf(['Lambda: %g\n'], lambda);

		[Theta1 Theta2] = trainNeuralNet( 
								 input_layer_size, ...
								 hidden_layer_size, ...
						600, initial_nn_params, X, y, lambda);
		error_train(i) =  predict(Theta1, Theta2, ...
						X, y);           
		error_val(i) =  predict(Theta1, Theta2, ...
						 Xval, yval); 
   end
%
%

% =========================================================================

end
