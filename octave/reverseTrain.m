
 %% Initialization
warning('off', 'Octave:possible-matlab-short-circuit-operator');
clear ; close all;

options = struct('mindiff', 1, 'plotx', 1, 'divnorm', 3, ...
     'savex', 'data/weights-reverseRegressiondiff.mat', 'loadx', 'data/weights-reverseRegressiondiff.mat',  ...
     'offset', 267653395647,  'train', 9500, 'val', 1495, 'readRows', 11000, ...
     'maxIter', 400, 'classEvalRule', 0, 'ignoreCutoff', 1); 

    fid = fopen('data/gram12.txt');
    [grampts] = fscanf(fid, '%f',[22, options.readRows])';
    fclose(fid);
    fid = fopen('data/zeros12.txt');
    [zerovals] = fscanf(fid, '%f', options.readRows);
    fclose(fid);
    rev_gram = flipud(grampts);
    rev_zero = flipud(zerovals);
 disp([size(grampts), size(zerovals)]);
 hidden_layer_size = 200;   % 
 num_labels = 1;          %  
 lambda = 0;
 sample = zeros( size(zerovals, 1)-1, 1);
 rev_sample = zeros( size(zerovals, 1)-1, 1);
    disp([size(rev_gram), size(rev_zero), size(rev_sample)]);
 %return;
 norm = 2*pi*options.divnorm;
 for i = 1:size(sample, 1)
   %sample(i) = (zerovals(i+1) - zerovals(i));
   value = (zerovals(i+1) - zerovals(i));
   sample(i) = value;
   rev_sample(i) = (rev_zero(i+1) - rev_zero(i));
 end
    
    izmatch = 1;
    
    [Xval, yval, izmatch] = XyVals(grampts, zerovals, sample, 1, options.val, options, izmatch);

    %{
    [X, y, test] = XyVals(rev_gram, rev_zero, rev_sample, options.val + 1, ...
                 options.train+options.val, options, izmatch);
   fprintf('... izmatch: %f \n', izmatch);
   %}             
   [X, y,test] = XyVals(grampts, zerovals, sample, options.val+1, ...
                 options.train+options.val, options, options.val - 1);
%    4055      21    4055       1   10996
% 1. gramidx is doing double duty : fix
% 2. zIdx is using zerovals, sensitive to forward bakward

     disp([size(X), size(y), test]);
    clear grampts zerovals fid;

    
 fprintf('Initializing Neural Network Parameters ...\n')
 input_layer_size = size(X, 2);

if exist('options', 'var') && ~isempty(options) && isfield(options, 'classEvalRule') ...
        && (options.classEvalRule == 1)
    num_labels = max(y)          %  
	initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
	initial_Theta2 = randInitializeWeights(hidden_layer_size, num_labels);
	% Unroll parameters
	initial_nn_params = [initial_Theta1(:) ; initial_Theta2(:)];
	J = logCostFunction(initial_nn_params, input_layer_size, hidden_layer_size, ...
					   num_labels, X, y, lambda);
	
	fprintf('Cost at init(log): %f \n', J);
	train_options = optimset('MaxIter', options.maxIter);
    		
	costFunction = @(p) logCostFunction(p, ...
									   input_layer_size, ...
									   hidden_layer_size, ...
									   num_labels, X, y, lambda);
	
	nn_params = fmincg(costFunction, initial_nn_params, train_options);
	
	Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
					 hidden_layer_size, (input_layer_size + 1));
	
	Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
					 num_labels, (hidden_layer_size + 1));
	pred = predictClass(Theta1, Theta2, X);
	error_train = (1-mean(double(pred == y))) * 100;
	fprintf('train percent error %g : confusion: \n',  error_train);
	confusion = zeros(num_labels);
	for i = 1:num_labels
	   for j = 1:num_labels
	      confusion(i,j) = sum(y==i & pred==j);
	   endfor
	endfor    
	disp(confusion)  
	predval = predictClass(Theta1, Theta2, Xval);
	error_val = (1-mean(double(predval == yval))) * 100;
	fprintf('val percent error %g : confusion: \n',  error_val);
	confusion = zeros(num_labels);
	for i = 1:num_labels
	   for j = 1:num_labels
	      confusion(i,j) = sum(yval==i & predval==j);
	   endfor
	endfor    
	disp(confusion)  
   return;
endif

%regression
    if isfield(options, 'load')
		load(options.load);
        initial_Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                         hidden_layer_size, (input_layer_size + 1));

        initial_Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                         num_labels, (hidden_layer_size + 1));
	else
        initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
        initial_Theta2 = randInitializeWeights(hidden_layer_size, num_labels);
        iterations = 0;
    endif
    % Unroll parameters
    initial_nn_params = [initial_Theta1(:) ; initial_Theta2(:)];
    wt = 1 + 3*(y<0.15) + 3*(y<0.1) + 3*(y<0.05) + 5*(y<0.007);
    excite = @(p) sigmoid(p);
    costFunction = @(p) nnCostFunction(p, ...
                                       input_layer_size, ...
                                       hidden_layer_size, ...
                                        X, y, lambda, wt);

rms = predict(initial_Theta1, initial_Theta2, X, y, excite);
initCost = costFunction(initial_nn_params);
fprintf('initial rms %12.3f initial cost %f ysize %d\n', rms, initCost, size(y,1));

train_options = optimset('MaxIter', options.maxIter, 'GradObj', 'on');

% Minimize using fmincg 400-0.000864102 800-0.000642887
[nn_params, fx, it] = fmincg(costFunction, initial_nn_params, train_options);
iterations = iterations + it;

% Obtain Theta1 and Theta2 back from nn_params
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));

[rms, pred] = predict(Theta1, Theta2, X, y, excite);
% Unroll parameters
nn_params = [Theta1(:) ; Theta2(:)];
finalCost = costFunction(nn_params);
fprintf('final rms %g final cost %g\n', rms, finalCost);
[rms, predval] = predict(Theta1, Theta2, Xval, yval, excite);
fprintf('val rms %g count %d\n', rms, size(yval,1));
if isfield(options, 'save')
    save ( options.save, 'nn_params', 'iterations');
endif
disp('plot  (+ actual, * validate');
figure;
hold on;
yvalactual = options.divnorm*(yval);
predvalactual = options.divnorm*(predval);
yactual = options.divnorm*(y);
predactual = options.divnorm*(pred);
plot(yactual, predactual,'k+', yvalactual, predvalactual,'r*')

