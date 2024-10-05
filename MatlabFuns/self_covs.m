function [covs,corrs]=self_covs(DAT,grade)

%autocorrelation calculation
%
%DAT: should be a vector
%grade: the grade for covariances

[n,m]=size(DAT);
%m must be 1
if ((m<1)|(m>1))&((n<1)|(n>1)),
   disp('The format of input data is incorrect! It sould be a vector!');
   return;
end

if (m==1),
   f_DAT=DAT;
else
   f_DAT=DAT';
end;
l=length(f_DAT);

if (nargin==1),grade=0;end

covs=[];
for i=0:grade
   r_sum=0;
   for k=1:(l-i)
      r_sum=r_sum+f_DAT(k)*f_DAT(k+i);
   end
   covs=[covs,r_sum/l];
end

if (nargout==2)
   corrs=covs/covs(1);
end

   