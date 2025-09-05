function y = readTestSig(fname,type)
% Function to read signal into MATLAB from file.
fileID = fopen(fname,'r');
formatSpec = '%f';
x = fscanf(fileID,formatSpec);

if (strcmp(type,'real'))
    y = x; 
else 
    j = 1;
    y = complex(zeros(numel(x)/2,1));
    for i = 1:numel(x)/2
        y(i) = x(j) + 1i*x(j+1);
       j = j + 2;
   end
end 

end
