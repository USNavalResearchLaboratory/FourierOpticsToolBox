%RANDGUASS  Generate a random number of guassians 
%    grayImage = randguass(gsize,r,numberOfGaussians,windowSize)

%               gsize               Size of the output, should be a 1x2 vector
%               radius              radius of guassian 
%               numberOfGaussians   Number of guassians generated
%               windowSize          Window of potential gaussian locations
%               factor              The scaled peak of the gaussians




function grayImage = randguass(gsize,radius,numberOfGaussians,windowSize,factor)

% windowSize = 1000; % Could be random if you want.
% r = 5;
sigma = sqrt(((radius+1)^2)/(2*log(255)));
% sigma = 10; % Could be random if you want.
% numberOfGaussians = 3;

% Set up some parameters.
fontSize = 20;
backgroundGrayLevel = 0;


% Create one Gaussian.
g = fspecial('gaussian', windowSize, sigma);
grayImage = backgroundGrayLevel * ones(gsize(1), gsize(2));
% Create random signs so that the Gaussians are
% randomly brighter or darker than the background.
%s = 2*randi(2, [1 numberOfGaussians])-3;
s = factor * ones(1,numberOfGaussians);
% Note: g and grayImage are floating point images, not uint8,
% though you could modify the program to have them be uint8 if you wanted.
% Get a list of random locations.
randomRow = randi(gsize(2)-windowSize+1, [1 numberOfGaussians]);
randomCol = randi(gsize(1)-windowSize+1, [1 numberOfGaussians]);
% Place the Gaussians on the image at those random locations.
for k = 1 : numberOfGaussians
  grayImage(randomRow(k):randomRow(k)+windowSize-1, randomCol(k):randomCol(k)+windowSize-1) = ...
    grayImage(randomRow(k):randomRow(k)+windowSize-1, randomCol(k):randomCol(k)+windowSize-1) + ...
    s(k) * g;
end

% % Display the final image.
% imshow(grayImage, []);
% caption = sprintf('%d Gaussians, Randomly Placed', numberOfGaussians);
% title(caption, 'FontSize', fontSize);
% axis on;
% % Enlarge figure to full screen.
% set(gcf, 'units','normalized','outerposition',[0 0 1 1])
% set(gcf,'name','Demo by ImageAnalyst','numbertitle','off')