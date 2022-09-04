% cluster OLS
%based on Cameron_Miller 2015
%Written by Zheyu Ni in the Ohio State Univeristy
%Y is the dependent variable X is dependent variable, g contains group
%information

function [beta, se_c]=OLS_Cluster(Y,X,g)


beta=(X'*X)\X'*Y;
error=Y-X*beta;
k=size(X,2);
n=size(X,1);
nofg=size(unique(g),1);
B=[];
for i = 1:nofg
    X_g=X(g==i,:);
    error_g=error(g==i);
    B=B+X_g'*(error_g*error_g')*X_g;
   
end

var_c=(X'*X)\B/(X'*X);
se_c=sqrt(diag(var_c));
