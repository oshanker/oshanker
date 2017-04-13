function checkNNGradients(lambda)
%CHECKNNGRADIENTS Creates a small neural network to check the
%backpropagation gradients
%   CHECKNNGRADIENTS(lambda) Creates a small neural network to check the
%   backpropagation gradients, it will output the analytical gradients
%   produced by your backprop code and the numerical gradients (computed
%   using computeNumericalGradient). These two gradient computations should
%   result in very similar values.
%    Need computeNumericalGradient

input_layer_size = 3;
hidden_layer_size = 2;
num_labels = 1;
m = 5;

% We generate some 'random' test data
Theta1 = randInitializeWeights( input_layer_size, hidden_layer_size);
Theta2 = randInitializeWeights( hidden_layer_size, num_labels);

X = zeros(m, input_layer_size);
y = zeros(m, 1);

for i = 1:m
yy = 1/11;
for j = 1:input_layer_size
X(i, j) = (i+j)/55;
yy += X(i, j);
end
y(i) = yy;
end

% Unroll parameters
nn_params = [Theta1(:) ; Theta2(:)];

% Short hand for cost function
if ~exist('lambda', 'var') || isempty(lambda)
costFunc = @(p) nnCostFunction(p, input_layer_size, hidden_layer_size, ...
                                X, y);
else
costFunc = @(p) nnCostFunction(p, input_layer_size, hidden_layer_size, ...
                                X, y, lambda);
end
%nnCostFunction(nn_params, input_layer_size, hidden_layer_size, X, y);
[cost, grad] = costFunc(nn_params);
numgrad = computeNumericalGradient(costFunc, nn_params);

% Visually examine the two gradient computations.  The two columns
% you get should be very similar. 
disp([numgrad grad]);
fprintf(['The above two columns you get should be very similar.\n' ...
         '(Left-Your Numerical Gradient, Right-Analytical Gradient)\n']);

% Evaluate the norm of the difference between two solutions.  
% If you have a correct implementation, and assuming you used EPSILON = 0.0001 
% in computeNumericalGradient.m, then diff below should be less than 1e-9
diff = norm(numgrad-grad)/norm(numgrad+grad);

fprintf(['If your backpropagation implementation is correct, then \n' ...
         'the relative difference will be small (less than 1e-9). ' ...
         'Relative Difference: %g\n'], diff);

end
