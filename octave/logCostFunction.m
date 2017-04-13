function [J grad] = logCostFunction(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, ...
                                   X, y, lambda)
%NNCOSTFUNCTION Implements the neural network cost function for a two layer
%neural network which performs classification
%   [J grad] = NNCOSTFUNCTON(nn_params, hidden_layer_size, num_labels, ...
%   X, y, lambda) computes the cost and gradient of the neural network. The
%   parameters for the neural network are "unrolled" into the vector
%   nn_params and need to be converted back into the weight matrices. 

Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));

% Setup some useful variables
m = size(X, 1);
         
% You need to return the following variables correctly 
J = 0;

% ====================== YOUR CODE HERE ======================
% Part 1: Feedforward the neural network and return the cost in the
%         variable J. 
%         Note: The vector y passed into the function is a vector of labels
%               containing values from 1..K. You need to map this vector into a 
%               binary vector of 1's and 0's to be used with the neural network
%               cost function.
%

a1 = [ones(m, 1) X];

z2 = a1*Theta1';

sz2 = sigmoid(z2);

a2 = [ones(m, 1) sz2];

z3 = [a2*Theta2'];

htheta = sigmoid(z3);

% size(htheta) m X num_labels

[zz, p] = max(htheta, [], 2);

lh = -log(htheta);
l1h = -log(1-htheta);
for i = 1:m
	for k = 1:num_labels
		if(y(i)==k)
		   J = J + lh(i, k);
		else
		   J = J + l1h(i, k);
		endif
	endfor
endfor
J=J/m;

if exist('lambda', 'var') && lambda > 0  
	reg = 0;
	for j =1:size(Theta1,1)
		for k = 2:size(Theta1,2)
		reg = reg + Theta1(j,k)*Theta1(j,k);
		endfor
	endfor
	for j =1:size(Theta2,1)
		for k = 2:size(Theta2,2)
		reg = reg + Theta2(j,k)*Theta2(j,k);
		endfor
	endfor

   J = J + lambda*reg/(2*m);
endif

%
% Part 2: Implement the backpropagation algorithm to compute the gradients
%         Theta1_grad and Theta2_grad. 

d3 = htheta;

for t = 1:num_labels
   d3(:,t) = d3(:,t) - (y==t);
endfor
z2grad = [zeros(m, 1) sz2.*(1-sz2)];
d2 = (d3*Theta2).*z2grad;

if exist('lambda', 'var') && lambda > 0  
	reg2 = lambda*([zeros(size(Theta2, 1), 1) Theta2(:, 2:end)])/m;
	reg1 = lambda*([zeros(size(Theta1, 1), 1) Theta1(:, 2:end)])/m;
	Theta2_grad = d3'*a2/m + reg2;
	Theta1_grad = d2(:, 2:end)'*a1/m + reg1;
else
	Theta2_grad = d3'*a2/m;
	Theta1_grad = d2(:, 2:end)'*a1/m;
endif

grad = [Theta1_grad(:) ; Theta2_grad(:)];

end
