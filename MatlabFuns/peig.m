function [eigVector,eigValue,validity]=peig(A,tol,tolnum,U)
%
[n,m]=size(A);
if (n>m)|(n<m)
   disp('The matrix A should be a Square Matrix!');
   valid=0;
   return;
end

if (nargin<=3),   U=ones(n);end
if (nargin<=2),   tolnum=10000;end
if (nargin<=1),   tol=10e-6;end

valid=1;

miu0=0.0;
miu1=0.0;
miu2=0.0;

miu=max(U);
U=U/miu;
num=1;
err=10e6;
while (num<tolnum)&(err>tol)
   V=A*U;
   miu=max(V);
   
   lamda=miu0-(miu1-miu0).^2/(miu-2*miu1+miu0);
   if (miu==0),
      disp('The principle eigenvalue is 0!');
      disp('You should choose a new initial value for U');
      valid=0;
      return;
   end
   
   U=V/miu;
   err=abs(miu-miu2);
   
   miu0=miu1;
   miu1=miu2;
   miu2=miu;
   num=num+1;
end

if (num==tolnum)
   valid=0;
   disp('The maximum number of iteration exceeded!');
   disp('The result is invalid!');
end

if (nargout>=3),
   validity=valid;
end

if (nargout>=2),
   eigValue=lamda;
end

if (nargout>=1),
   eigVector=U;
end

