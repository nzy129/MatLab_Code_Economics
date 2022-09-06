% function for linear GMM, input is data Y, X,and IV Z, weighting matrix W.
% Written by Zheyu Ni at the Ohio State University
% Output is GMM coefficient beta, VCV matrix and predicted errors.

function [beta, VCV, error]=LinGMM(Y,X,Z,W)

N=size(Y,1);
beta=(X'*Z*W*Z'*X)\(X'*Z*W*Z'*Y);
error=Y-X*beta;
M=error.*Z;
Omega=1/N*(M'*M);
G=1/N*Z'*X;

VCV=(G'*W*G)\(G'*W*Omega*W*G)/(G'*W*G);


end
