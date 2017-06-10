function [] = print2file(variable,folder,numformat,separator,fileformat)
% This function print a variable into a single .txt with the name of the variable itself 
% example print2file(variable,'folder\subfolder','%d2.1','\n','dat')

% Argument management:
%default value for
if nargin < 5
   fileformat = 'txt'
end
if nargin < 4
   separator = '\n'
end
if nargin < 3
   numformat = '%f4.3'
end

filename = inputname(1); %get the name of the variable
fileID = fopen( ['report\result\' filename '.txt'],'wt');
    fprintf(fileID,[numformat separator],variable);
fclose(fileID);

end

