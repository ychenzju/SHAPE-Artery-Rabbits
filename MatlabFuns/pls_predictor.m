function [Y_prediction]=pls_predictor(X,W,P,B,Q);

%
[nX,mX]=size(X);
[nP,mP]=size(P);
[nB,mB]=size(B);
[nQ,mQ]=size(Q);
if (mX<nP)|(mX>nP)|(mP<nB)|(mP>nB)|(nB<mB)|(nB>mB)|(nB<mQ)|(nB>mQ),
   disp('The Input Data Structure Does not Match with Each Other!');
   return;
end

dim=mP;
X_R=X;
Y_prediction=zeros(nX,nQ);
for k=1:dim
   T_p(:,k)=(X_R*W(:,k))/(W(:,k)'*W(:,k));
   X_R=X_R-T_p(:,k)*P(:,k)';
   Y_prediction=Y_prediction+B(k,k)*T_p(:,k)*Q(:,k)';
end
