function [ X, Y ] = genGrids(M, N, dx, dy)
%genGrids - generate numerical meshgird arrays
%
% This function creates the meshgrids X and Y with M rows and N cols. If no
% dx or dy is given, then the default spacing between points is one.
% Otherwise the spacing between points is dx along the horizontal
% (x-direction) and dy along the vertical (y-dierection. If only dx is
% given, then the dy is equal to dx. The center of the grid (zero) is at:
% (row, col) = (M/2+1, N/2+1).
%
% Syntax: [ X, Y ] = genGrids(M, N, dx, dy)
%   [X, Y] = genGrids(128, 128); 
%               or
%   [X, Y] = genGrids(128, 128, 1E-3); 
%               or
%   [X, Y] = genGrids(128, 64, 1E-3, 2E-2); 
%
% Author: Dennis Gardner
%
% Initial Creation: 9/2/2016
% Modified on 9/6/2016 - one can fed in M and N or an array
% Modified on 8/07/2017 - more helper file
% modified 10/03/2017

if ~exist('dx', 'var') && ~exist('dy', 'var')
    dx = 1;
    dy = 1; 
    
elseif exist('dx', 'var') && ~exist('dy', 'var')
    dy = dx;
end

% grid vectors
xgv = ((-N/2:1:N/2-1)-0.5)*dx;
ygv = ((-M/2:1:M/2-1)-0.5)*dy;

% grid arrays
[X, Y] =  meshgrid(xgv, ygv);





end