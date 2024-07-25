function [ circle ] = createCirclePx(M, N, radius)
%creatCirclePx - create circle with given radius
%   Detailed explanation goes here



x = -N/2:N/2-1;
y = -M/2:M/2-1;

[X, Y] = meshgrid(x,y);


circle = sqrt(X.^2 + Y.^2); 

circle(circle>radius) = 0; 
circle(circle>0) = 1; 
if radius >= 1
    
    circle(round(M/2+1), round(N/2+1)) = 1; 
end


end

