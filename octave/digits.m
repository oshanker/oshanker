%% Neural Network Learning (classification)

%  https://sites.google.com/site/riemannzetazeros/Home/programs
%  if digits, need digits.mat digitweights.mat logCostFunction
%   predictClass

%% Initialization
clear ; close all;

 # zcount digits classdiff
 options = struct('digits', 1, 'plot', 1, 'divnorm', 3, 'offset', 267653395647, ...
     'savex', 'weights-classdiff.mat', 'load', 'weights-classdiff.mat',  ...
     'zeros', 'data/zeros12.txt', 'gram', 'data/gram12.txt', 'train', 9500, ...
     'maxIter', 1000, 'lambda', 0);
%{
    options.zeros = 'data/zeros12.txt';
    options.gram = 'data/gram12.txt';
    options.offset = 267653395647;
    options.train = 0; 
%}        

cutoff = 0.13;     


%% ================ Part 2: Loading Parameters ================

if isfield(options, 'digits')  
	input_layer_size  = 400;  % 20x20 Input Images of Digits
	hidden_layer_size = 25;   % 25 hidden units
	num_labels = 10;          % 10 labels, from 1 to 10   

    load('data/digits.mat');
    disp([size(X), size(y)]);	
	% Load the weights into variables Theta1 and Theta2
	load('data/digitweights.mat');
    %Theta1 = randInitializeWeights(400, 25);
    %Theta2 = randInitializeWeights(25, 10);
    disp([size(Theta1), size(Theta2)]);	
	nn_params = [Theta1(:) ; Theta2(:)];
	
	% Weight regularization parameter (we set this to 0 here).
	lambda = 0;
	
	J = logCostFunction(nn_params, input_layer_size, hidden_layer_size, ...
					   num_labels, X, y, lambda);
	
	fprintf(['Cost at parameters (loaded from digitweights): %f '...
			 '\n(this value should be about 0.103068) lambda %f \n'], J, lambda);
	
	%% =============== Part 4: Implement Regularization ===============
	%  Once your cost function implementation is correct, you should now
	%  continue to implement the regularization with the cost.
	%
	
	fprintf('\nChecking Cost Function (w/ Regularization) ... \n')
	
	% Weight regularization parameter (we set this to 1 here).
	lambda = 1;
	
	J = logCostFunction(nn_params, input_layer_size, hidden_layer_size, ...
					   num_labels, X, y, lambda);
	
	fprintf(['Cost at parameters (loaded from weights) (w/ lambda = 1): %f '...
			 '\n(this value should be about 9.717137 )\n'], J);
	
	% Also output the costFunction debugging values
	lambda = 3;
	debug_J  = logCostFunction(nn_params, input_layer_size, ...
							  hidden_layer_size, num_labels, X, y, lambda);
	val = -1;
	fprintf(['\n\nCost at (fixed) debugging parameters (w/ lambda = 3): %f ' ...
			 '\n(this value should be about 28.945276 )\n\n'], debug_J);
    pred = predictClass(Theta1, Theta2, X);
    fprintf('init pred %g\n', mean(double(pred == y)) * 100);

else
	 [grampts, rzeros, sample] = init(options); 
	 if isfield(options, 'load')
		load(options.load);
		input_layer_size = size(Theta1,2)-1;
		hidden_layer_size = size(Theta1,1);
		num_labels = size(Theta2, 1);
	 else
		 input_layer_size = 34;
		 hidden_layer_size = 250;	
		 num_labels = 2;           
		Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
		Theta2 = randInitializeWeights(hidden_layer_size, num_labels);
	 endif
	
	 val = 0;
	 [X, yy] = calcXy(grampts, rzeros, sample, 1, options.train, ...
				   input_layer_size, num_labels, options, Theta1, Theta2);
	 y =  2 - (yy<cutoff);
endif

%% ================ Part 6: Initializing Pameters ================
 if (options.maxIter > 0 && options.train > 1)
    %{
     # C 18580341990000 : 81958236867985 + .922933955444789782293          
     nextInput = {'13', 'A', 'B', 'M', 'D'};
     nextOffset = [1034741742900, 7954022502370, 323393653040, 2124447368580, ...
         18523741991600];
     for i = 1:(size(nextInput, 2))
		 nopt = rmfield(options, {'train'});
		 nopt.zeros = ['data/zeros' nextInput{1,i} '.txt'];
		 nopt.gram = ['data/gram' nextInput{1,i} '.txt'];
		 nopt.offset = nextOffset(i);
		 [Xadd, yyadd] = calcXy([], [], [], 1, -1, ...
			input_layer_size, num_labels, nopt);
	     yadd =  2 - (yyadd<cutoff);
		 X = [X; Xadd];
		 y = [y; yadd];
		 yy = [yy; yyadd];
     endfor 
     %}

	fprintf('\nInitializing Neural Network Parameters ...\n')
    % val = floor((size(rzeros, 1) - options.train - num_labels - 2)/2)
	initial_nn_params = [Theta1(:) ; Theta2(:)];

	train_options = optimset('MaxIter', options.maxIter);
    		
	costFunction = @(p) logCostFunction(p, ...
									   input_layer_size, ...
									   hidden_layer_size, ...
									   num_labels, X, y, options.lambda);
	
	[nn_params, error_val] = fmincg(costFunction, initial_nn_params, train_options);
	
	Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
					 hidden_layer_size, (input_layer_size + 1));
	
	Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
					 num_labels, (hidden_layer_size + 1));
	
	 if isfield(options, 'save')
		old_error_val = error_val;
		save ( options.save, 'Theta1', 'Theta2', 'old_error_val');
	 endif
endif	
	
pred = predictClass(Theta1, Theta2, X);
error_train = mean(double(pred == y)) * 100;
error_val = -1;
error_test = -1;

if val > 1    
    [Xval, yyval] = calcXy(grampts, rzeros, sample, 1+options.train, val+options.train, ...
                       input_layer_size, num_labels, options);
	yval =  2 - (yyval<cutoff);
	pred_val = predictClass(Theta1, Theta2, Xval);
	error_val = mean(double(pred_val == yval)) * 100;
    
    [Xtest, yytest] = calcXy(grampts, rzeros, sample, 1+val+options.train, ...
                     2*val+options.train, input_layer_size, num_labels, options);
	ytest =  2 - (yytest<cutoff);
	pred_test = predictClass(Theta1, Theta2, Xtest);
	error_test = mean(double(pred_test == ytest)) * 100;
 
    disp([sum(y),sum(yval), sum(ytest), sum([sum(y),sum(yval), sum(ytest)])])	
    disp([sum(pred),sum(pred_val), sum(pred_test), ... 
         sum([sum(pred),sum(pred_val), sum(pred_test)])])
endif	

#plotData(X, y>1)
disp(' error_train, error_val error_test');
disp([ error_train, error_val,  error_test]);

  fprintf('TP %d FN %d TN %d FP %d\n', sum((pred==1) & (y == 1)), ...
  sum((pred==2) & (y == 1)),  sum((pred==2) & (y == 2)), ...
  sum((pred==1) & (y == 2)))

if isfield(options, 'plotyyyyyyyy')
    maxpts =  size(y, 1);
	plot(1:maxpts, pred(1:maxpts,1), '-*k', ...
	        1:maxpts, y(1:maxpts,1), '-+k')
	legend('predicted', 'actual');
	title(ttl);
	grid();
end


function comp(pred, y)
  fprintf('TP %d FN %d TN %d FP %d\n', sum((pred==1) & (y == 1)), ...
  sum((pred==2) & (y == 1)),  sum((pred==2) & (y == 2)), ...
  sum((pred==1) & (y == 2)))
endfunction