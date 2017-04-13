
 %% Initialization
warning('off', 'Octave:possible-matlab-short-circuit-operator');
clear ; close all;

options = struct('mindiff', 1, 'plotx', 1, 'divnorm', 3, ...
     'savex', '',  'offset', 267653395647, ...
      'train', 9500, 'val', 1495, 'readRows', 11000, ...
     'maxIter', 400, 'classEvalRule', 1); 

    fid = fopen('data/gram12.txt');
    [grampts] = fscanf(fid, '%f',[22, options.readRows])';
    %disp(grampts(:,1))
    fclose(fid);
%history 10 -q
    fid = fopen('data/zeros12.txt');
    [zerovals] = fscanf(fid, '%f', options.readRows);
    %disp(zerovals)
    fclose(fid);
 disp([size(grampts), size(zerovals)]);
 hidden_layer_size = 200;   % 
 num_labels = 1;          %  
 lambda = 0;
 sample = zeros( size(zerovals, 1)-1, 1);
 norm = 2*pi*options.divnorm;
 for i = 1:size(sample, 1)
   %sample(i) = (zerovals(i+1) - zerovals(i));
   value = (((zerovals(i+1) - zerovals(i))*log((options.offset+zerovals(i))/(2*pi))/(norm)));
   sample(i) = value;
 end
    
    izmatch = 1;
    
    [Xval, yval, izmatch] = XyVals(grampts, zerovals, sample, 1, options.val, options, izmatch);
    [X, y] = XyVals(grampts, zerovals, sample, options.val+1, ...
                 options.train+options.val, options, izmatch);
    
fprintf('Initializing Neural Network Parameters ...\n')
 input_layer_size = size(X, 2);
%wt = 1./(1.0+y);


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
initial_Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
initial_Theta2 = randInitializeWeights(hidden_layer_size, num_labels);
% Unroll parameters
initial_nn_params = [initial_Theta1(:) ; initial_Theta2(:)];
wt = 1 + 3*(y<0.15) + 3*(y<0.1) + 3*(y<0.05) + 5*(y<0.007);
excite = @(p) sigmoid(p);
costFunction = @(p) nnCostFunction(p, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                    X, y, lambda, wt);
%{
wt = [];
excite = @(p) tanh(p);
costFunction = @(p) tanhnnCostFunction(p, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                    X, y, lambda, wt);
%}

rms = predict(initial_Theta1, initial_Theta2, X, y, excite);
initCost = costFunction(initial_nn_params);
fprintf('initial rms %12.3f initial cost %f ysize %d\n', rms, initCost, size(y,1));

train_options = optimset('MaxIter', options.maxIter, 'GradObj', 'on');

% Minimize using fmincg
[nn_params, fx, it] = fmincg(costFunction, initial_nn_params, train_options);


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
%a= [yval(1:10), predval(1:10)]
figure;
hold on;
yvalactual = options.divnorm*(yval);
predvalactual = options.divnorm*(predval);
plot(yvalactual, predvalactual,'r*')
disp('Press return to continue');
pause;

yactual = options.divnorm*(y);
predactual = options.divnorm*(pred);

plot( yactual, predactual,'k+')
validx = find(yvalactual<0.2);
yidx = find(yactual<0.2);
%%{
figure;
plot( yactual(yidx), predactual(yidx),'k+', ...
   yvalactual(validx), predvalactual(validx),'r*')
%%}