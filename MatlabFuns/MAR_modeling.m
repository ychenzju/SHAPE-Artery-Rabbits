function [Beta,F]=MAR_modeling(X,p)

% Set up a p-dimensional multivariable AR model to fit the relation
% between X and Y.
% 

% calculate the covariance of X (with -p+1 to p-1 grade)

[n,m]=size(X);

YY=X(p+1:n,:);
XX=[];
for i=1:p
   XX=[XX,X(p-i+1:n-i,:)];
end
   
TAOp=[];
for h=1:p+1
   temp=zeros(m,m);
   for t=1:n-(h-1)
      temp=temp+X(t+(h-1),:)'*X(t,:);
   end
   TAOp(:,:,h)=temp/n;
end

GAMA=[];
TAO=[];
for i=1:p
   temp_T=[];
   for j=1:p
      if (i<=j)
         temp_T=[temp_T,TAOp(:,:,j-i+1)];
      else
         temp_T=[temp_T,TAOp(:,:,i-j+1)];
      end
   end
   TAO=[TAO;temp_T];
   GAMA=[GAMA,TAOp(:,:,i+1)];
end

Beta=inv(TAO)*GAMA';
F=YY-XX*Beta;

rate=press(F)/press(YY)

