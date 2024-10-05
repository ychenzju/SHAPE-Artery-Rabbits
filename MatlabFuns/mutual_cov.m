function [cov,coeff]=mutual_cov(v1,v2,grade)

% covs=E{V1(t)'V2(t+grade)}
% v1 and v2 should be a column vector
% and the length should be same

if (nargin==1), v2=v1;end
if (nargin<3), grade=0;end

[n1,m1]=size(v1);
if (m1<1)|(m1>1),
   disp('The format of the first data incorrect!');
   return;
end
[n2,m2]=size(v2);
if (m2<1)|(m2>1),
   disp('The format of the second data incorrect!');
   return;
end
if (n1<n2)|(n1>n2),
   disp('The length of two vector should be same');
   return;
end

f_temp=0;
if (grade==0)
   for f_i=1:n1
      f_temp=f_temp+v1(f_i)*v2(f_i);
   end
   cov=f_temp/n1;
elseif (grade>0),
   for f_i=1:(n1-grade)
      f_temp=f_temp+v1(f_i)*v2(f_i+grade);
   end
   cov=f_temp/(n1-grade);
else
   for f_i=1:(n1+grade)
      f_temp=f_temp+v1(f_i-grade)*v2(f_i);
   end
   cov=f_temp/(n1+grade);
end

if nargout>1,coeff=cov/sqrt(var(v1)*var(v2));end


