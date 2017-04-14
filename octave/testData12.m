    % helps to show what is in data/gram12.txt
    fid = fopen('data/gram12.txt');
    [grampts] = fscanf(fid, '%f',[22, 2]);
    disp([ grampts(1,2) - grampts(1,1)]);
    for i = 4:12
       disp([int32(i), (i-2)*(grampts(i,1)^2 + grampts(i+9,1)^2), ... 
          (i-2)*(grampts(i,2)^2 + grampts(i+9,2)^2)])
    endfor
