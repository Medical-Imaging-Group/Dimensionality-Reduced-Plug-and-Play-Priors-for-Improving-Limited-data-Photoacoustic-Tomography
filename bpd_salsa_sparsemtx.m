
function [x, cost] = bpd_salsa_sparsemtx(y, A, lambda, mu, Nit)

% x = bpd_salsa_sparsemtx(y, A, lambda, mu, Nit)
%
% BASIS PURSUIT DENOISING
% Minimize ||y - A x||_2^2 + lambda * || x ||_1
% where A is a sparse matrix
%
% INPUT
%   y      : data
%   A      : sparse matrix
%   lambda : regularization parameter
%   mu     : ADMM parameter
%   Nit    : Number of iterations
%
% OUTPUT
%   x      : solution to BPD problem
%
% [x, cost] = bpd_salsa_sparsemtx(...) returns cost function history

% Ivan Selesnick
% NYU-Poly
% selesi@poly.edu
% March 2012

% The program implements SALSA (Afonso, Bioucas-Dias, Figueiredo,
% IEEE Trans Image Proc, 2010, p. 2345)


if nargout > 1
    ComputeCost = true;
    cost = zeros(1,Nit);
else
    ComputeCost = false;
end    

ATy = A'*y;
x = ATy;
d = zeros(size(x));

[M, N] = size(A);

F = A'*A + mu*speye(N);
issparse(F)

for i = 1:Nit
    
%     u = soft(x + d, 0.5*lambda/mu);
%     x = F \ (ATy + mu*(u-d));
%     d = d - u + x;
    
    u = soft(x + d, 0.5*lambda/mu) - d;
    x = F \ (ATy + mu*u);
    d = x - u;
     
    if ComputeCost
        residual = y - A*x;
        cost(i) = sum(abs(residual(:)).^2) + sum(abs(lambda * x(:))); 
    end
end
