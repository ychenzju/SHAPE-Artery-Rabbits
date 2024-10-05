function [P,T,R]=pca_cv(X,part,err)
% @function [P,T,R]=pca_cv(X,part,err)
% Principal Component Analysis with Cross-Validation
% The reference samples are divided into 'part' blocks. Use 'part-1'
% blocks of the sample data to construct PCA model and use the remaining
% block for validation. The dimension of implicit subspace is determined
% through the PRESS (compare the errors among different PCA 
% model with different implicit space dimension)
%
% X is a (n * m) matrix; 
%       n -- the number of samples; m -- the number of variables in one sample
% d is the dimension of implicit subspace
%       on which projection will be conducted
% P is the projection matrix, (m * d)
% T is the score matrix, (n * d); T = XP
% R is the residuals, R = X - XPP'
%
%   See also pca_pls

%   Y. Chen 5-10-00
%   Revised Y. Chen 2-11-03
%   Copyright General Electric, Inc. 
%   $Revision: 1.1 $  $Date: 2003/02/12 15:16:04 $


if (nargin==1),err=1e-6;part=5;end
if (nargin==2),err=1e-6;end

doPlot = 1;

[n,m]=size(X);
partition(1)=1;
for i=2:part
   partition(i)=round((i-1)*n/part);
end
partition(part+1)=n;

press_total=10000000;
press_total_PRE=inf;
h=0;
while ((press_total<press_total_PRE)&(h<m))
   h=h+1;
   sumOfPress=0;
   for j=1:part
      if (j==1),
         XX=X(partition(2):partition(part+1),:);
      elseif (j<part)
         XX=[X(partition(1):partition(j),:);X(partition(j+1):partition(part+1),:)];
      else
         XX=X(1:partition(part),:);
      end
      [PP,TT]=pca_pls(XX,h);
      YY=X(partition(j):partition(j+1),:);
      sumOfPress=sumOfPress+press(YY-YY*PP*PP');
   end
   press_total_PRE=press_total;
   press_total=sumOfPress;
   
   if doPlot==1
       press_total_PRE
       press_total
   end
end

if (h<2)|(press_total<press_total_PRE)
   [P,T,R]=pca_pls(X,h,err);
else
   [P,T,R]=pca_pls(X,h-1,err);
end

if doPlot==1
    h
    P
    press(R)
end
