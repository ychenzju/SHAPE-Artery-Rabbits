function [V]=xTy(X,Y);

% V=X'Y
% X and Y are variables with zero mean

[nX,mX]=size(X);
if (nargin>1),
   [nY,mY]=size(Y);
   if (nX<nY)|(nX>nY),
      disp('WRONG: The length of X should be equal to Y!');
      return;
   end
end


if (nargin==1),
   V=X'*X;
elseif (nargin==2),
   V=X'*Y;
end

   

