function [IM4_LSQR,I_Lsqr] = bpd_lsqr(A_b,sdn2_v)

%%%%%%%%%%%%%%Code for Model Resolution based BPD using LSQR based image reconstruction using
%%%%%%%%%%%%%% Optimal Choice of Regularization as and k=25 regularization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  for i=2:50
tic;
l_step =40;


[U1,V1,B1] =lsqr_b_hybrid(A_b,sdn2_v,l_step,1);
T = eye(l_step);
%  lambda=5e-1;
lambda = 0.3;
%   lambda = fminbnd(@(lambda) opt_lambda_cw(B1(1:l_step,1:l_step+1),sdn2_v, lambda, V1, A_b*0.15, T, U1), 0, lambda, optimset( 'MaxIter',1000, 'TolX', 1e-16));

reg = lambda;
T = eye(l_step);
Hess = B1(1:l_step,1:l_step+1)'*B1(1:l_step,1:l_step+1);
yk = (Hess + reg.*eye(size(B1(1:l_step,1:l_step+1),2)))\(norm(sdn2_v,2).*B1(1:l_step,1:l_step+1)'*T(:,1));


%% changed
I_Lsqr = V1(:,1:length(yk))*yk;
toc
I_Lsqr = reshape(I_Lsqr,201,201);


%%

K = B1'*B1;
N3 = (K + (lambda.*eye(size((K),2))))\(K);


lambda1 = 0.00001;                                      % lambda : regularization parameter
Nit = 10000;                                           % Nit : number of iterations
mu1 = 0.01;                                           % mu : ADMM parameter
tic
% Run BPD algorithm
[x_BPD, cost] = bpd_salsa_sparsemtx(yk, N3, lambda1, mu1, Nit);


foo = V1(:,1:length(yk))*x_BPD;


IM4_LSQR = reshape(foo,2*100+1,2*100+1);
toc;
