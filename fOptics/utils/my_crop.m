function [ Icrop ] = my_crop( I, newArraySize )
%crop - crop the image
%   I can be a 2D or 3D array

% original size of image
[oldRows, oldCols, ~] = size(I); 


% newRows and newCols
newRows = newArraySize(1); 
newCols = newArraySize(2);

% center position
cenRow = oldRows/2+1; 
cenCol = oldCols/2+1;

Icrop = I(cenRow - newRows/2 : cenRow + newRows/2-1, ...
          cenCol - newCols/2 : cenCol + newCols/2-1,:);

end

