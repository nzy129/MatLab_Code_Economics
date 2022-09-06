function [beta,VCV,error,W_new,J]=TGMM_cluster(Y,X,Z,W0,g)
% function for two-step linear GMM, input is data Y, X,and IV Z, initial
% weighting matrix W0, group information g

% output is second-step GMM coefficient beta,second-step VCV matrix,
% their predicted errors and new weighting matrix in second stage.
% Clustered Standard Error
[beta1,VCV1,error1]=LinGMM(Y,X,Z,W0);

B=zeros(k,k);
for i = 1:nofg
    X_g=X(g==i,:);
    error_g=error1(g==i);
    %B=B+X_g'*error_g*error_g'*X_g;
    
    B=B+(error_g'*X_g)'*(error_g'*X_g);
   
end

N=size(Y,1);
%get the new efficient weighting matrix.
W_new=inv(1/N*B);

[beta,VCV,error]=LinGMM(Y,X,Z,W_new);

%J=1/N*error'*Z*W_new*(error'*Z);
J=error'*Z\B/(error'*Z);%N is cancelled out. 
end
