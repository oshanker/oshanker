function [zeroDiffValues] = scrip
# Knowing some diffs, predict next diff.
# http://www.opengardensblog.futuretext.com/wp-content/uploads/2014/03/OctaveTutorialAndrewNg.pdf 
#http://www.holehouse.org/mlclass/

#need data, randInit, learningcurve, trainNeuralNet, fmincg
#nnCostFunction, sigmoid, predict

warning('off', 'Octave:possible-matlab-short-circuit-operator');
id = fopen('data/zeros12.txt');
[rzeros] = fscanf(id, '%f');
offset = 267653395647;

norm = 2*pi*3;
factor = log(offset+rzeros(1)/(2*pi))/(norm);
#find zero diffs, normalized
j = 0;
max = 0.015/factor;
for i = 1:(size(rzeros, 1)-1)
   diff = rzeros(i+1) - rzeros(i);
   if(diff > max)
      continue;
   endif
   j++;
   zeroDiffValues(j, 1:3) = [(diff)*factor, rzeros(i), rzeros(i+1)];
end
# minimum, first quartile, median, third
# quartile, maximum, mean, standard deviation, skewness, and kurtosis
# of the columns.
statistics(zeroDiffValues)'
save(['data/yt' int2str(offset) '.txt'], 'zeroDiffValues');
return;

input_layer_size = 100;
hidden_layer_size = 300;   %  hidden units
train = 2000;
val = 3000;
[X, y, Xval, yval, valbegin] = calcXY(zeroDiffValues, input_layer_size, train, val);
disp(['valbegin', ' diff']);
disp([valbegin, zeroDiffValues(valbegin+1+input_layer_size), yval(1,1)]);
sz = 40
Xtest = Xval(1:sz, :);
ytest = yval(1:sz, :);

%sizes = linspace((train-400), train, 2);
sizes  = train;

num_labels = size(y, 2);   
%{         
		load weights.mat;
%}
Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
Theta2 = randInitializeWeights(hidden_layer_size, num_labels);


initial_nn_params = [Theta1(:) ; Theta2(:)];
%disp([size(initial_Theta1), size(initial_Theta2), size(initial_nn_params)])

[error_train, error_val, Theta1, Theta2] = ...
    learningCurve(initial_nn_params, hidden_layer_size, X, y, Xval, yval, sizes);
disp('sizes, error_train, error_val');
disp([sizes', error_train, error_val]);

%{
return;
old_error_val = error_val;
save weights.mat Theta1 Theta2 old_error_val
%}

%{
label = ('Number of training examples');
ttl = ' Learning Curve (lambda = 0)';
show(sizes, error_train, error_val, label,ttl);

fprintf('Program paused (Number of training examples). Press enter to continue.\n');
pause;
%}  
  

[pred ] = predictLocal(Xtest, Theta1, Theta2);
disp('pred(1), pred(2), y(1) y(2)');
disp([pred(1:10,:), ytest(1:10,:)])
#plot(1:size(pred, 1), pred(:,1), 'r', 1:size(pred, 1), ytest(:,1))
#legend('pred', 'actual');

fprintf('Program paused (predict). Press enter to continue.\n');
pause;
#close all;

%{
disp("forward")
[pred ] = forward(Xtest(1, :), Theta1, Theta2, ytest);
%disp([pred, ytest])
plot(1:size(pred, 1), pred(:,1), 'r', 1:size(pred, 1), ytest(:,1))
legend('pred', 'actual');

fprintf('Program paused (forward). Press enter to continue.\n');
pause;
%}



%% =========== Part 8: Validation for Selecting Lambda =============
%  You will now implement validationCurve to test various values of 
%  lambda on a validation set. You will then use this to select the
%  "best" lambda value.
%

% Selected values of lambda 
%{
lambda_vec = [ 0.003 0.01 0.03 0.1 ]';
[error_train, error_val] = ...
    validationCurve(input_layer_size, ...
	     hidden_layer_size, ...
         X, y, Xval, yval, lambda_vec);
close all;
label = 'lambda';
show(lambda_vec, error_train, error_val, label);

fprintf('lambda\t\tTrain Error\tValidation Error\n');
for i = 1:length(lambda_vec)
	fprintf(' %f\t%f\t%f\n', ...
            lambda_vec(i), error_train(i), error_val(i));
end

fprintf('Program paused. Press enter to continue.\n');
pause;
%}
fprintf('done.\n');

end;

##########################################################################
function [pred ] = predictLocal(X, Theta1, Theta2)
	[junk, pred] = predict(Theta1, Theta2, X, []);
end

##########################################################################
function [pred ] = forward(X, Theta1, Theta2, y)
	sz = size(y, 1);
	rowx = X(1,:);
	pred = zeros(size(y));
	for i = 1:sz
		[junk, val] = predict(Theta1, Theta2, rowx, []);
		rowx = [rowx(2:end) val(1,1)];
		pred(i,:) = val;
	end
	return
end

##########################################################################
# L = input features count, train = m, 
function [X, y, Xval, yval, valbegin] = calcXY(zeroDiffValues, L, train, val)
	begin = 1;
	
	[X, y] = shape(zeroDiffValues, L, train, begin);
	
	valbegin = (size(zeroDiffValues, 1) - val - L - 1);
	[Xval, yval] = shape(zeroDiffValues, L, val, valbegin);
end

##########################################################################
function show(lambda_vec, error_train, error_val, label,ttl)
	plot(lambda_vec, error_train, 'r', lambda_vec, error_val);
	if exist('ttl', 'var') && ~isempty(ttl) 
	   title(ttl);
	end    
	legend('Train', 'Cross Validation');
	xlabel(label);
	ylabel('Error');
end

##########################################################################
function [X, y] = shape(rzeros, L, train, begin)
	X = zeros(train, L);
	y = zeros(train, 2);
	for i = 1:train
		idxr = i + begin;
		X(i, 1:(L)) = rzeros(idxr:(idxr+L-1));
		y(i, 1) = rzeros(idxr+L);
		y(i, 2) = rzeros(idxr+L+1);
	end
	return;
end

##########################################################################
function [X, y] = easyshape(rzeros, L, train, begin)
	X = reshape(rzeros(begin:(begin + L*train - 1)), L, train);
	X = X';
	y = [X(2:end,1); rzeros(begin+numel(X))];
	return;
end



