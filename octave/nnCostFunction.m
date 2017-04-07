function [J grad] = nnCostFunction(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   X, y, lambda, wt)
%NNCOSTFUNCTION Implements the neural network cost function for a two layer
%neural network which performs classification
%   [J grad] = NNCOSTFUNCTON(nn_params, hidden_layer_size, num_labels, ...
%   X, y, lambda) computes the cost and gradient of the neural network. The
%   parameters for the neural network are "unrolled" into the vector
%   nn_params and need to be converted back into the weight matrices. 
% 
%   The returned parameter grad should be a "unrolled" vector of the
%   partial derivatives of the neural network.
%
% lambda is Weight regularization parameter
% wt is a set of weights to be used in calculating the cost function.
% The gradient needs to be checked for correctness in handling wt
% (numerical grad check seems to show a problem when both wt and lambda
% are present. We seem to be ok if only of the two is present,
% or when both are absent)

% Reshape nn_params back into the parameters Theta1 and Theta2, the weight matrices
% for our 2 layer neural network
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));
xxx = nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end);
Theta2 = reshape(xxx, size(y, 2), (hidden_layer_size + 1));

% Setup some useful variables
m = size(X, 1);
         
% You need to return the following variables correctly 
J = 0;
Theta1_grad = zeros(size(Theta1));
Theta2_grad = zeros(size(Theta2));

% Part 1: Feedforward the neural network and return the cost in the
%         variable J. 

a1 = [ones(m, 1) X];

z2 = a1*Theta1';

sz2 = sigmoid(z2);

a2 = [ones(m, 1) sz2];

z3 = [a2*Theta2'];

htheta = sigmoid(z3);


%{
for i = 1:m
J += -y(i)*log(htheta(i))-(1-y(i))*log(1-htheta(i));
end
J=J/m;
%}

% size(htheta) m X num_labels
if exist('wt', 'var')  && isequal(size(wt), size(y)) 
	d3 = wt.*(htheta-y);
	norm = sum(sum(wt));
	J = sum(sumsq(d3./sqrt(wt)))/(2*norm);
else
	d3 = htheta-y;
	norm = m;
	J = sum(sumsq(d3))/(2*norm);
end
d3 = d3.*(htheta.*(1-htheta));

if exist('lambda', 'var')  && lambda > 0  
	reg = 0;
	reg +=  sum(sum(Theta1(:,2:end).*Theta1(:,2:end)));
	reg +=  sum(sum(Theta2(:,2:end).*Theta2(:,2:end)));
	
	J = J + lambda*reg/(2*m);
end

z2grad = [zeros(m, 1) sz2.*(1-sz2)];

d2 = (d3*Theta2).*z2grad;

if exist('lambda', 'var') && lambda > 0  
	reg2 = lambda*([zeros(size(Theta2, 1), 1) Theta2(:, 2:end)])/m;
	reg1 = lambda*([zeros(size(Theta1, 1), 1) Theta1(:, 2:end)])/m;
	Theta2_grad = d3'*a2/m + reg2;
	Theta1_grad = d2(:, 2:end)'*a1/m + reg1;
else
	Theta2_grad = d3'*a2/norm;
	Theta1_grad = d2(:, 2:end)'*a1/norm;
endif

grad = [Theta1_grad(:) ; Theta2_grad(:)];

end
