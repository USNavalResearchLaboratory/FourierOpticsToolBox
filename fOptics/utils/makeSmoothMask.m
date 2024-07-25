function [mask_smooth] = makeSmoothMask(efield, diameter)
%makeSmoothMask Make smooth mask based on Field dimensions

arguments
    efield
    diameter = 10;
end
radius = (efield.N/2);
xx = -(radius-0.5):(radius-0.5);
yy = -(radius-0.5):(radius-0.5);
[XX,YY] = meshgrid(xx,yy);
[~,RR] = cart2pol(XX,YY);
mask_smooth = 1 - exp(-2*(RR./(0.5*diameter)).^8);

end