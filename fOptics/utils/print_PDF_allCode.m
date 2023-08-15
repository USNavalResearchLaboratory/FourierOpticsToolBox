% Create pdf of every file in folder 

clc; clear all; 

%%
pth=fullfile('.','*.m');
files=filefun(pth,Inf);
%%
options = struct('format','pdf','outputDir','./PDFsOfCode/', 'evalCode', false);
for i = 1:length(files)
publish(files{i}, options)
end