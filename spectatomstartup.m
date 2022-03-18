function spectatomstartup()
% spectatomstartup  Startup file for spectatom
% MAKE adds paths of the spectatom to Matlab.

file_path = mfilename('fullpath');
tmp = strfind(file_path,'spectatomstartup');
file_path = file_path(1:(tmp(end)-1));

% Folder for all test functions
addpath([file_path 'test']);

% Folder for all source files recursively
addpath(genpath([file_path 'src']));

end
