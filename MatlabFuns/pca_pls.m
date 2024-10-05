function [P,T,R]=pca_pls(X,d,err)
% @ function [P,T,R]=pca_pls(X,d,err)
% Principal Component Analysis via a nonlinear iterative method (often used in PLS)
%
% X is a (n * m) matrix; 
%       n -- the number of samples; m -- the number of variables in one sample
% d is the dimension of implicit subspace
%       on which projection will be conducted
% P is the projection matrix, (m * d)
% T is the score matrix, (n * d); T = XP
% R is the residuals, R = X - XPP'

%   Y. Chen 05-09-00
%   Revised Y. Chen 02-11-03
%   Copyright General Electric, Inc. 
%   $Revision: 1.1 $  $Date: 2003/02/12 15:16:04 $

if (nargin==2),err=1e-6;end
[n,m]=size(X);

P=[];
T_h=[];
E=X;
h=0;
rand_col=random('unid',m);
th=E(:,rand_col);

while (h<d)
   h=h+1;
   delta=2*err;
   th_PRE=0;
   while (delta>err)
      ph=(th'*E)'/(th'*th);
      ph=ph/norm(ph);
      th=E*ph;
      delta=norm(th-th_PRE);
      th_PRE=th;
   end
   E=E-th*ph';
   T_h=[T_h,th];
   P=[P,ph];
end

if (nargout==2),
   T=T_h;
elseif (nargout==3),
   T=T_h;
   R=E;
end