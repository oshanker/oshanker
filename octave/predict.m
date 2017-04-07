function [rms, p] = predict(Theta1, Theta2, X, y, excite)
%PREDICT Predict the label of an input given a trained neural network
%   p = PREDICT(Theta1, Theta2, X) outputs the predicted label of X given the
%   trained weights of a neural network (Theta1, Theta2)

% Useful values
m = size(X, 1);
num_labels = size(Theta2, 1);

if !exist('excite', 'var') 
   excite = @(p) sigmoid(p);
end   

in1 = [ones(m, 1) X] * Theta1';
h1 = excite(in1);
h2 = excite([ones(m, 1) h1] * Theta2');
p = h2;
if !isempty(y)
	diff = p-y;
	xx = sumsq(diff)/m;
	rms = sqrt(sum(xx)/size(xx,1));
else
   rms = 0;
endif

% =========================================================================


end
