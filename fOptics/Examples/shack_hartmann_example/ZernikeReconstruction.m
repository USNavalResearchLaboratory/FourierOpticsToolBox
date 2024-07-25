% This function reconstructs wavefront from slope matrices dZx and dZy
% with spacing dx up to ’terms’ Zernike polynomials.
%
function A = ZernikeReconstruction(SH, terms)

dZx = SH.slopes(:,:,1); 
dZy= SH.slopes(:,:,2); 
dx = SH.dx; 


[n, m] = size(dZx);
nn = n;
xx = reshape(dZx, nn^2, 1);
yy = reshape(dZy, nn^2, 1);
ss = [xx yy]';
S = reshape(ss, 2*nn^2, 1);
clear xx; clear yy; clear ss;
X = calcMatrix(nn, terms, dx);
[U, W, V] = svd(X, 0);
clear X;
W = pinv(W);
A = V*W*U'*S;
clear U; clear W; clear V; clear S;
end
% This function calculate the matrix for Zernike derivatives
function Z = calcMatrix(nn, terms, dx);
R = (nn-1)*dx/2;
[X,Y] = meshgrid(-R:dx:R);
r = sqrt(X.^2+Y.^2);
Z = zeros(2*nn^2, terms);

for i = 1:terms
z = zeros(1, terms+1);
z(i+1) = 1;
S = ZernikePolynomials(nn, z);
[dZx, dZy] = gradient(S, dx);
dZx(r>R) = 0;
dZy(r>R) = 0;
xx = reshape(dZx, nn^2, 1);
yy = reshape(dZy, nn^2, 1);
ss = [xx yy]';
Z(:,i) = reshape(ss, 2*nn^2, 1);
end
% This function calculates the wavefront based on a set of Zernike
% coefficients z to get frame size (nn+1)x(nn+1).
%
end
function S = ZernikePolynomials(nn, z)
terms = length(z)-1;
[X, Y] = meshgrid(-1:2/nn:1);
r = sqrt(X.^2+Y.^2);
r(X.^2+Y.^2>1) = 0;
Theta = atan2(Y, X);
S = zeros(nn+1);
for i = 0:terms
[n, m] = single2doubleZ(i);
if (m == 0)
pa = sqrt(n+1);
else
pa = sqrt(2*(n+1));
end
coef = pa;
Surf = zeros(nn+1);
for s = 0:(n-abs(m))/2
c1 = n-s;
c2 = (n+m)/2-s;
c3 = (n-m)/2-s;
Surf = Surf + (-1)^s*factorial(c1)/factorial(s)...
/factorial(c2)/factorial(c3)*power(r, n-2*s);
end
if (m < 0)
Surf = Surf.*sin(abs(m)*Theta);
elseif (m > 0)
    Surf = Surf.*cos(m*Theta);
end
S = S + z(i+1)*coef*Surf;
end
S(r > 1) = 0;
% This function converts single->double index in Zernike polynomials
%
end
function [n, m] = single2doubleZ(jj)
n = floor(sqrt(2*jj+1)+0.5)-1;
m = 2*jj-n*(n+2);
end
