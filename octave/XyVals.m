function [Xval, yval, izmatch] = XyVals(grampts, zerovals, sample,  
  beginIdx, endIdx, options, izmatch)
    cutoff = 1.0/options.divnorm;
    pts = 2;
    if(beginIdx <= pts)
       beginIdx = pts+1;
    endif   
    if(mod(beginIdx,2 ) == 0)
       fprintf('even %ld, \n', beginIdx );
    else   
       fprintf('odd %ld, changing it \n', beginIdx );
       beginIdx++;
    endif   
    input_layer_size = 3*(3+2*pts);
    idxX = 1;
    idxToGram = beginIdx; 
    while idxToGram <= endIdx;
        Xin = [grampts((idxToGram-pts):(idxToGram+2+pts), 2) ...
               ; grampts((idxToGram-pts):(idxToGram+2+pts), 4) ... 
               ; grampts((idxToGram-pts):(idxToGram+2+pts), 5)
           ];
		[izmatch, rangez] = zIdx(grampts, zerovals, idxToGram, izmatch);	
		idxToGram++;
		idxToGram++;
		%{
		[izmatch, rangez2] = zIdx(grampts, zerovals, idxToGram, izmatch);	
		rangez = [rangez, rangez2];
		%}
		if isempty(rangez)
		   continue
		endif
		value = min(sample(rangez));
		%%{
		if (value > cutoff)
		   continue
		endif
		%%}
		Xval(idxX, 1:input_layer_size) = Xin; 
        if exist('options', 'var') && ~isempty(options) && isfield(options, 'classEvalRule') ...
        && (options.classEvalRule == 1)
			if(value < 0.25/options.divnorm)
			   yval(idxX,1) = 1;
			elseif(value < 0.35/options.divnorm)
			   yval(idxX,1) = 2;
			elseif(value < 0.45/options.divnorm)
			   yval(idxX,1) = 3;
			elseif(value < 0.55/options.divnorm)
			   yval(idxX,1) = 4;
			else    
			   yval(idxX,1) = 5;
			endif 
		else   
			yval(idxX,1) = value;
		endif   
		idxX++;
	    %fprintf('mindiff %g\n', value);
    endwhile
endfunction