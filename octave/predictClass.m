function p = predictClass(Theta1, Theta2, X, options)
%PREDICT Predict the label of an input given a trained neural network
%   p = PREDICT(Theta1, Theta2, X) outputs the predicted label of X given the
%   trained weights of a neural network (Theta1, Theta2)

% Useful values
m = size(X, 1);
num_labels = size(Theta2, 1);

% You need to return the following variables correctly 
p = zeros(size(X, 1), 1);

h1 = sigmoid([ones(m, 1) X] * Theta1');
h2 = sigmoid([ones(m, 1) h1] * Theta2');
#h2 m x num_labels

if exist('options', 'var') && ~isempty(options) && isfield(options, 'classEvalRule')
switch (options.classEvalRule)
case 'evencheck'
   for idxr = 1:m
      begin = 1;
      if sign(X(idxr,1)) == sign(X(idxr,21))
         begin = 2;
      endif
      array = begin:2:num_labels;
      for idxclass = 1:num_labels
         for j = array
            h2(idxr, j) = 0;
         endfor
      endfor
   endfor
otherwise
endswitch
endif
[dummy, p] = max(h2, [], 2);

% =========================================================================


end
