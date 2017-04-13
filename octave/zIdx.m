function [izmatch, range] = zIdx(grampts, rzeros, gbegin, izcand)
% izmatch  - closest zero idx after grampt(gbegin) in forward direction
%             could be outside the gram interval 
% range - zeros in gram interval beginning at gbegin (empty if closest zero
%              is outside the gram interval.
% izcand - initial guess for zero idx

    izmatch = idx(grampts, rzeros, gbegin, izcand);
    if (nargout > 1)
       izmax = size(rzeros, 1);
       range = [];
       found = 0;
	   if rzeros(izmatch) < grampts(gbegin, 1)
		  return;
       endif  
       next = izmatch;
       while (rzeros(next) < grampts(gbegin+1, 1))
          found++;
          range(found) = next;
          next++;
          if(next > izmax)
             break;
          endif   
       endwhile
    endif
endfunction


##########################################################################
function izmatch = idx(grampts, rzeros, gbegin, izcand)

    if izcand > size(rzeros, 1)-1
       izcand = size(rzeros, 1)-1;
    endif
    
	izmatch = izcand;
	
	izmax = size(rzeros, 1)-1;
	
	if rzeros(izmatch) == grampts(gbegin, 1)
	   return;
	end   
	
	if rzeros(izcand) >= grampts(gbegin, 1)
	   while (true)
		   previdx = izmatch - 1;
		   if izmatch == 1
			  return;
		   end  
		   if rzeros(previdx) < grampts(gbegin, 1)
			  return;
		   end  
		   izmatch = previdx;
	   endwhile
	else
	   while (true)
		   previdx = izmatch + 1;
		   if izmatch == izmax
			  return;
		   endif  
		   izmatch = previdx;
		   if rzeros(izmatch) >= grampts(gbegin, 1)
			  return;
		   endif  
	   endwhile
	endif

endfunction
