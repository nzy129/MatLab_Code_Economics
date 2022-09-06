function [beta,se_r,se,RSS]=OLS_r(Y,X)

beta=(X'*X)\X'*Y;
error=Y-X*beta;
k=size(X,2);
n=size(X,1);
var_r=(X'*X)\(X'.*repmat((error.^2)',k,1)*X)/(X'*X);
se_r=sqrt(diag(var_r));
var=sum(error.^2)/(n-k)*inv(X'*X);
se=sqrt(diag(var));
RSS=sum(error.^2);

end
