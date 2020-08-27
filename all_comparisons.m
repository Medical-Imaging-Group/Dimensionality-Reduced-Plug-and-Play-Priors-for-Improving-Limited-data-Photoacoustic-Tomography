
%% code for generating all the results for a single phantom for the paper " 
%  Dimensionality Reduced Plug and Play Priors for Improving Limited data Photoacoustic Tomography "
%  Author - Navchetan Awasthi
%  Date Written - 30/August/2019
%  Date Modified - 21/November/2019

%%
function [I_Lsqr, IM4_LSQR, x, x_k ] = all_comparisons(A_b, b, tol, U, S, V, inverse_term)
%% BPD+LSQR
% This will give the results using base function of lpd_lsqr.
% Used for Comparisons

%%
[IM4_LSQR,I_Lsqr] = bpd_lsqr(A_b,b);

%% Total Variation
% This will give the results using base function of tv.
% Used for Comparisons

%%
tic
[x,~,~] = tv(b(:),A_b,inverse_term,tol);
toc

%% IDBP
% This will give the results using plug and play priors
% Load first the SVD matrix for initialization of U, V and S and the
% measurement vector b;

%%

close all
Number_of_Singular_Values = 40401;
k = 200;
sigma = zeros(size(S));
sigma(1:Number_of_Singular_Values,:)=S(1:Number_of_Singular_Values,:);

% Inverse SVD Solution
G = U'*b(:);
H = sigma'*G;
I = V*H;
imshow(reshape(I,201,201),[]);

params.lambda = 1 ;%Lagrange multiplier
params.beta = 0.5; %The regularization parameter 
kparams.sigma=sqrt(params.beta/params.lambda);
kparams.verb = 0;
kparams.display = 0;
kparams.niter = 50;  % number of iterations
kparams.c_TV = .018;
kparams.lambda = 2*kparams.sigma^2*kparams.c_TV;% initial regularization
    
Y_tilde1 = I;
Y_tilde = I;
tic
for i=1:k
    i
    [x_k,err,my_tv,lalist] = perform_tv_denoising(reshape(Y_tilde1,201,201),kparams);
    int_k = A_b*reshape(x_k,40401,1);
    int_k_n = U'*int_k;
    Y_tilde1 = Y_tilde + reshape(x_k,40401,1)-V*((sigma'*int_k_n));
end
toc
    figure
    imshow(imcomplement(reshape(x_k,201,201)),[]);
    figure
    imshow((reshape(x_k,201,201)),[]);



