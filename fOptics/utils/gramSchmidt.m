function U = gramSchmidt(V)
        % code copied from Wiki on 02/15/2018
        % https://en.wikipedia.org/wiki/Gram-Schmidt_process
        k = size(V,1);
        l = size(V,2);
        U = zeros(k,l);
        U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
        for ii = 2:l
          U(:,ii) = V(:,ii);
          for jj = 1:ii-1
            U(:,ii) = U(:,ii) - ( U(:,ii)'*U(:,jj) )/( U(:,jj)'*U(:,jj) )*U(:,jj);
          end
          U(:,ii) = U(:,ii)/sqrt(U(:,ii)'*U(:,ii));
        end
end