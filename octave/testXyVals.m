
 %% Initialization
warning('off', 'Octave:possible-matlab-short-circuit-operator');
clear ; close all;

options = struct('mindiff', 1,  'divnorm', 3, ...
      'val', 95, 'readRows', 100, ...
     'maxIter', 400, 'classEvalRule', 0, 'ignoreCutoff', 1); 

    fid = fopen('data/gram12.txt');
    [grampts] = fscanf(fid, '%f',[22, options.readRows])';
    %disp(grampts(:,1))
    fclose(fid);
%history 10 -q
    fid = fopen('data/zeros12.txt');
    [zerovals] = fscanf(fid, '%f', options.readRows);
    %disp(zerovals)
    fclose(fid);
 hidden_layer_size = 200;   % 
 num_labels = 1;          %  
 lambda = 0;
 sample = zeros( size(zerovals, 1)-1, 1);
 norm = 2*pi*options.divnorm;
 for i = 1:size(sample, 1)
   sample(i) = (zerovals(i+1) - zerovals(i));
   %{
   value = (((zerovals(i+1) - zerovals(i))*log((options.offset+zerovals(i))/(2*pi))/(norm)));
   sample(i) = value;
   %}
 end
    
    izmatch = 1;
    
    [Xval, yval, imatch] = XyVals(grampts, zerovals, sample, 1, options.val, options, izmatch);
    disp([size(Xval), size(yval)])
    pts = 2;
    beginIdx = 1;
    if(beginIdx <= pts)
       beginIdx = pts+1;
    endif   
    if(mod(beginIdx,2 ) == 1)
       beginIdx++;
    endif   
    idxToGram = beginIdx; 
    while idxToGram <= options.val;
		[izmatch, rangez] = zIdx(grampts, zerovals, idxToGram, izmatch);
		if isempty(rangez)
		   disp(izmatch);
    else
		   rangez = [ zerovals(rangez)'];
       rangez = [grampts(idxToGram,1), izmatch,  rangez];
       disp( rangez);
		endif
		idxToGram++;
		idxToGram++;
    endwhile
    clear grampts zerovals fid;
save out/XyVals.mat;