function [A,G]=fa(X,err)
%
if (nargin==1), err=1e-4;end

[n,m]=size(X);

Sigma=X'*X/(n-1);
S=inv(Sigma);
G=inv(diag(diag(S)));

GSigma=G.^(-1/2)*Sigma*G.^(-1/2);
[U,L]=eig(GSigma);
A=G.^(1/2)*U*(L-)