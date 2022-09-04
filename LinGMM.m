function [beta, VCV, error]=LinGMM(Y,X,Z,W)
%%function for linear GMM, input is data Y, X,and IV Z, weighting matrix W.
%%Output is GMM coefficient beta, VCV matrix and predicted errors.
N=size(Y,1);
beta=(X'*Z*W*Z'*X)\(X'*Z*W*Z'*Y);
error=Y-X*beta;
M=error.*Z;
Omega=1/N*(M'*M);
G=1/N*Z'*X;

VCV=(G'*W*G)\(G'*W*Omega*W*G)/(G'*W*G);








end
