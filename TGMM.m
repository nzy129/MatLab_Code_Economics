function [beta,VCV,error,W_new]=TGMM(Y,X,Z,W0)
%%function for two-step linear GMM, input is data Y, X,and IV Z, initial
%%weighting matrix W0

%%output is second-step GMM coefficient beta,second-step VCV matrix,
%%their predicted errors and new weighting matrix in second stage.

[beta1,VCV1,error1]=LinGMM(Y,X,Z,W0);
M1=error1.*Z;
N=size(Y,1);
Omega_new=1/N*(M1'*M1);%get the new efficient weighting matrix.
W_new=inv(Omega_new);
[beta,VCV,error]=LinGMM(Y,X,Z,W_new);






end
