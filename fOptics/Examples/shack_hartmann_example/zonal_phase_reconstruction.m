function W = zonal_phase_reconstruction(SH)

Sx = SH.slopes(:,:,1); 
Sy= SH.slopes(:,:,2); 
ds = SH.dx; 

[n, n] = size(Sx);
S = [reshape(Sx', 1, n*n) reshape(Sy', 1, n*n)]';
E = getE(n);
[U, D, V] = svd(E, 0);
D = pinv(D);
C = getC(n);
W = V*D*U'*C*S;
W = reshape(W', n, n)./ds;
%W = exp(1i.*W); 
% This function obtains the matrix E for zonal reconstruction
%
function E = getE(n)
E = zeros(2*n*(n-1),n*n);
for i = 1:n
for j = 1:(n-1)
E((i-1)*(n-1)+j,(i-1)*n+j) = -1;
E((i-1)*(n-1)+j,(i-1)*n+j+1) = 1;
E((n+i-1)*(n-1)+j,i+(j-1)*n) = -1;
E((n+i-1)*(n-1)+j,i+j*n) = 1;
end
end
% This function obtains the matrix C for zonal reconstruction
%
function C = getC(n)

C = zeros(2*n*(n-1),2*n*n);
for i = 1:n
for j = 1:(n-1)
C((i-1)*(n-1)+j,(i-1)*n+j) = 0.5;
C((i-1)*(n-1)+j,(i-1)*n+j+1) = 0.5;
C((n+i-1)*(n-1)+j,n*(n+j-1)+i) = 0.5;
C((n+i-1)*(n-1)+j,n*(n+j)+i) = 0.5;
end
end