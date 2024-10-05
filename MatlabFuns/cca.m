function [A,B,V]=cca(X,Y,d);

%Canonical Correlation Analysis
%X and Y should be with mean zero

[nX,mX]=size(X);
[nY,mY]=size(Y);
if ((mX<mY)|(nX<nY)|(nX>nY)),
   disp('WRONG: The varaible number in X should be bigger than Y!');
   disp('---OR--- ');
   disp('WRONG: The length of X doesnot equal to Y!');
   return;
end
if (nargin<3),
   d=mY;
end

%estimation of auto- and cross- correlation of X and Y
Lxx=xTy(X);
Lyy=xTy(Y);
Lxy=xTy(X,Y);
Lyx=Lxy';

pinv_Lxx=pinv(Lxx);
pinv_Lyy=pinv(Lyy);

L1=pinv_Lxx*Lxy;
L2=pinv_Lyy*Lyx;

M1=L1*L2;
M2=L2*L1;

V=zeros(d);
A=[];
B=[];
for k=1:d
   [vecA,valA]=peig(M1);
   [vecB,valB]=peig(M2);
   
   A=[A,vecA]
   B=[B,vecB]
   V(k,k)=valA
   M1=M1-valA*vecA*vecA';
   M2=M2-valB*vecB*vecB';
end





