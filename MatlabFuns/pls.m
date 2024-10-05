function [W,P,Q,B,RY,RX]=pls(X,Y,d,err)
% Partial Least Square Method
% As a biased regression method, PLS can be used to
% estimate the parameters of multi-dimensional
% regression equation.
% Y=X*Beta+F
% where Y is (n,m1) matrix, X is (n,m2) and Beta is 
% (m2,m1). n is the number of samples.
% Beta=W*inv(P'W)*Q'
% The residues F=RY

%   Y. Chen 05-09-00
%   Revised Y. Chen 02-11-03
%   Copyright General Electric, Inc. 
%   $Revision: 1.1 $  $Date: 2003/02/12 15:16:04 $

if (nargin==3), err=1e-6;end
[n_X,m_X]=size(X);
[n_Y,m_Y]=size(Y);
if (n_X<n_Y)
   Y=Y(1:n_X,:);
   n=n_X;
elseif (n_X>n_Y)
   X=X(1:n_Y,:);
   n=n_Y;
end

E=X;F=Y;
W=[];P=[];Q=[];B=[];C=[];modi_Q=[];
h=0;
rand_col_Y=random('unid',m_Y);
uh=F(:,rand_col_Y);
while (h<d)
   h=h+1;
   delta=2*err;
   th_PRE=0;
   uh=rand(n_Y,1);
	while (delta>err)
      wh=(E'*uh)/(uh'*uh);%[m_X,1]
      wh=wh/norm(wh);
      th=(E*wh)/(wh'*wh);%[n,1]
      qh=(F'*th)/(th'*th);%[m_Y,1]
      qh=qh/norm(qh);
      uh=(F*qh)/(qh'*qh);%[n,1]
      delta=norm(th-th_PRE);
      th_PRE=th;
   end
   ph=(E'*th)/(th'*th);%[m_X,1]
   %wh=wh*norm(ph);
   wh=wh/norm(wh);
   ph=ph/norm(ph);
   modi_qh=(F'*uh)/(uh'*uh);
   bh=(uh'*th)'/(th'*th);
   
   E=E-th*ph';
   F=F-bh*th*qh';
   
   W=[W,wh];
   P=[P,ph];
   Q=[Q,qh];
   B=[B,bh];
   modi_Q=[modi_Q,modi_qh];
end
B=diag(B);

if (nargout==5),
   RY=F;
elseif (nargout==6),
   RX=E;
   RY=F;
end




   