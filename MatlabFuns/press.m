function R=press(X)
% Sum of square of a matrix elements
% A=(aij); A is a matrix of n*m
% R=Sigma(aij^2)

[n,m]=size(X);
sum=0;
for i=1:n
   for j=1:m
      sum=sum+X(i,j)^2;
   end
end
R=sum;