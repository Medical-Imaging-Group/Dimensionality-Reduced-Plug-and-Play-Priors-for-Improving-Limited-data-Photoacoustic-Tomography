
%% code for generating all the results for the phantom data for the paper " Dimensionality Reduced Plug and Play Priors for Improving Limited data Photoacoustic Tomography "
%   Author - Navchetan Awasthi
%  Date Written - 30/August/2019
%  Date Modified - 21/November/2019

%% PAT phantom
display('PAT_phantom')

load('P1.mat')
b = Rref(:);
tol = 0.1;
[I_Lsqr_PAT_20db, IM4_LSQR_PAT_20db, x_PAT_20db,...
    x_k_PAT_20db ] = all_comparisons(A_b, b, tol, U, S, V, inverse_term);

load('P2.mat')
b = Rref(:);
tol = 0.01;
[I_Lsqr_PAT_40db, IM4_LSQR_PAT_40db, x_PAT_40db,...
    x_k_PAT_40db ] = all_comparisons(A_b, b, tol, U, S, V, inverse_term);

load('P3.mat')
b = Rref(:);
tol = 0.001;
[I_Lsqr_PAT_60db, IM4_LSQR_PAT_60db, x_PAT_60db,...
    x_k_PAT_60db ] = all_comparisons(A_b, b, tol, U, S, V, inverse_term);

%% Blood _Vessel
display('Blood_Vessel_phantom')

load('P4.mat')
b = Rref(:);
tol = 0.1;
[I_Lsqr_blood_vessel_20db, IM4_LSQR_blood_vessel_20db, x_blood_vessel_20db,...
    x_k_blood_vessel_20db ] = all_comparisons(A_b, b, tol, U, S, V, inverse_term);

load('P5.mat')
b = Rref(:);
tol = 0.01;
[I_Lsqr_blood_vessel_40db, IM4_LSQR_blood_vessel_40db, x_blood_vessel_40db,...
    x_k_blood_vessel_40db ] = all_comparisons(A_b, b, tol, U, S, V, inverse_term);

load('P6.mat')
b = Rref(:);
tol = 0.001;
[I_Lsqr_blood_vessel_60db, IM4_LSQR_blood_vessel_60db, x_blood_vessel_60db,...
    x_k_blood_vessel_60db ] = all_comparisons(A_b, b, tol, U, S, V, inverse_term);

%% Derenzo Phantom
display('Derenzo_phantom')

load('P7.mat')
b = Rref(:);
tol = 0.1;
[I_Lsqr_Derenzo_20db, IM4_LSQR_Derenzo_20db, x_Derenzo_20db,...
    x_k_Derenzo_20db ] = all_comparisons(A_b, b, tol, U, S, V, inverse_term);

load('P8.mat')
b = Rref(:);
tol = 0.01;
[I_Lsqr_Derenzo_40db, IM4_LSQR_Derenzo_40db, x_Derenzo_40db,...
    x_k_Derenzo_40db ] = all_comparisons(A_b, b, tol, U, S, V, inverse_term);

load('P9.mat')
b = Rref(:);
tol = 0.001;
[I_Lsqr_Derenzo_60db, IM4_LSQR_Derenzo_60db, x_Derenzo_60db,...
    x_k_Derenzo_60db ] = all_comparisons(A_b, b, tol, U, S, V, inverse_term);
