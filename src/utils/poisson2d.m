function [p] = poisson2d(n,m)
% poisson2d(n,m)
% k superblocks, m blocks, order n
% construct Block tridiagonal matrix for the poisson operator
a=sparse(toeplitz([-2 1 zeros(1,n-2)]));
b=sparse(toeplitz([-2 1 zeros(1,m-2)]));
idm = sparse([1:m],[1:m],1);
idn = sparse([1:n],[1:n],1);
p = -kron(b,idn)-kron(idm,a);

