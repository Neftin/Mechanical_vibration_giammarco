function [] = print2file(variable,folder,numformat,separator,fileformat,filename)
% This function print a variable into a single .txt with the name of the variable itself 
% example print2file(variable,'folder\subfolder','%2.1f','\n','dat')

% Argument management:
%default value for
if nargin < 6
   filename = inputname(1);
end
if nargin < 5
   fileformat = 'txt';
end
if nargin < 4
   separator = '\n';
end
if nargin < 3
   numformat = '%f4.3';
end

%get the name of the variable
fileID = fopen( ['report\result\' filename '.' fileformat],'wt');
    fprintf(fileID,[numformat separator],variable);
fclose(fileID);

end

