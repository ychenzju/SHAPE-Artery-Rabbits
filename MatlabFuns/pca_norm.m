function  [M,S,Y]=pca_norm(X)
% Normalization of matrix X
% Mean(Y(:,k)) = 0; Var(Y(:,k)) = 1

[n,m]=size(X);
M=mean(X);
S=var(X);
for (j=1:m)
   Y(:,j)=(X(:,j)-M(j))/sqrt(S(j));
end
